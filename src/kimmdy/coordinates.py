"""coordinate, topology and plumed modification functions"""

from dataclasses import dataclass, field
from enum import Enum
import logging
import multiprocessing as mp  # to start helper processes
import os  # to read the number of helpers from the environment
from concurrent.futures import ProcessPoolExecutor  # runs the work on several helpers at once
from copy import deepcopy
from pathlib import Path
from typing import Optional, Union
from math import sqrt

import MDAnalysis as mda
import numpy as np

from kimmdy.constants import DEFAULT_EDISSOC, FFFUNC
from kimmdy.parsing import read_plumed, write_plumed
from kimmdy.recipe import Place
from kimmdy.tasks import TaskFiles
from kimmdy.topology.atomic import (
    Angle,
    AngleId,
    AngleType,
    Bond,
    BondId,
    BondType,
    Dihedral,
    DihedralId,
    DihedralType,
    Pair,
    Exclusion,
    Interaction,
    AtomicType,
    InteractionTypes,
    MultipleDihedrals,
)
from kimmdy.topology.ff import FF
from kimmdy.topology.topology import MoleculeType, Topology
from kimmdy.topology.utils import match_atomic_item_to_atomic_type
from enum import Enum, auto

logger = logging.getLogger(__name__)


class PairTransition(Enum):
    Morph = auto()
    Create = auto()
    Vanish = auto()


class NonbondedState(Enum):
    Excluded = auto()
    OneFour = auto()
    Full = auto()


def to_pairkey(key: BondId | AngleId | DihedralId) -> tuple[str, str]:
    """Turn key into a pair key
    Takes the key for a higher ordre interaction (bond, angle, dihedral)
    and from this takes the outermost atoms and returns them as a key for the corresponding
    pair interaction by sorting them in ascending order.
    """
    pairkey: tuple[str, str] = tuple(sorted([key[0], key[-1]], key=int))  # type: ignore
    return pairkey


# coordinates
def place_atom(
    files: TaskFiles, step: Place, ttime: Optional[float] = None
) -> TaskFiles:
    """Place an atom to new coords at the last time point of the trajectory"""
    logger = files.logger
    logger.info("Starting place_atom task")
    logger.debug(step)
    logger.debug(f"ttime {ttime}")
    trr = files.input["trr"]
    tpr = files.input["tpr"]

    # if place_atom is called multiple times in one _apply_recipe task
    if "trr" in files.output.keys():
        trr = files.output["trr"]

    u = mda.Universe(str(tpr), str(trr), topology_format="tpr", format="trr")

    for ts in u.trajectory[::-1]:
        if ttime is not None:
            if abs(ts.time - ttime) > 1e-5:  # 0.01 fs
                continue
        atm_move = u.select_atoms(f"index {step.ix_to_place}")
        atm_move[0].position = step.new_coords

        break
    else:
        raise LookupError(
            f"Did not find time {ttime} in trajectory "
            f"with length {u.trajectory[-1].time}"
        )

    trr_out = files.outputdir / "coord_mod.trr"
    gro_out = files.outputdir / "coord_mod.gro"

    u.atoms.write(trr_out)  # type: ignore
    u.atoms.write(gro_out)  # type: ignore

    files.output["trr"] = trr_out
    files.output["gro"] = gro_out

    logger.debug(
        "Exit place_atom, final coordinates written to "
        f"{'/'.join(trr_out.parts[-2:])}"
    )
    return files


def is_parameterized(entry: Interaction):
    """Parameterized topology entries have c0 and c1 attributes != None"""
    return entry.c0 is not None and entry.c1 is not None


def get_explicit_or_type(
    key: tuple[str, ...],
    interaction: Optional[Interaction],
    interaction_types: InteractionTypes,
    mol: MoleculeType,
    periodicity: str = "",
) -> Union[Interaction, AtomicType, None]:
    """Takes an Interaction and associated key, InteractionTypes, Topology
    and Periodicity (for dihedrals) and returns an object with the parameters of this Interaction
    """
    if not interaction:
        return None

    if is_parameterized(interaction):
        return interaction

    type_key = [mol.atoms[x].type for x in key]

    atomic_type = match_atomic_item_to_atomic_type(
        type_key, interaction_types, periodicity
    )

    if atomic_type is None:
        m = f"Could not find explicit parameters for {[type_key,interaction]} in the FF. Consider setting the parameterizer to `grappa` to generate missing parameters in the fly."
        logger.error(m)
        raise ValueError(m)

    return atomic_type


@dataclass
class PairSets:
    added: set[tuple[str, str]] = field(default_factory=set)
    removed: set[tuple[str, str]] = field(default_factory=set)
    morphed: set[tuple[str, str]] = field(default_factory=set)


class BondsPairSets(PairSets):
    pass


class AnglesPairSets(PairSets):
    pass


class DihedralsPairSets(PairSets):
    pass


@dataclass
class AffectedInteractions:
    """keeping track of affected interactions during the merge process to add the correct helper pairs"""

    atoms: set[str] = field(default_factory=set)
    bonds: BondsPairSets = field(default_factory=BondsPairSets)
    angles: AnglesPairSets = field(default_factory=AnglesPairSets)
    dihedrals: DihedralsPairSets = field(default_factory=DihedralsPairSets)


# the main process puts its merger here so the helpers can read it
_PARALLEL_DIHEDRAL_MERGER = None


def _merge_dihedral_worker(payload):
    """Work out the relax parameters of one dihedral (runs in a helper process).

    Input: the dihedral's four atoms, and whether it is an improper.
    It looks up the dihedral before and after the reaction and mixes them
    into one dihedral that changes from "before" to "after" during the relax.
    It only returns the result; the main process puts it into the topology.

    Returns: the atoms, the mixed dihedral (None if nothing changes), whether
    the dihedral is new, removed or changed, and whether it is an improper.
    """
    key, use_improper = payload  # the four atoms of one dihedral, and whether it is an improper
    merger = _PARALLEL_DIHEDRAL_MERGER  # the main process's merger (only read here)
    if merger is None:
        raise RuntimeError("parallel dihedral worker was not initialized")
    # look up the dihedral before (A) and after (B) the reaction
    dihedral_a = merger._get_explicit_MultipleDihedrals(
        key=key, use_state_b=False, use_improper=use_improper
    )
    dihedral_b = merger._get_explicit_MultipleDihedrals(
        key=key, use_state_b=True, use_improper=use_improper
    )

    # note whether the dihedral is new, removed or changed
    category = None
    if not use_improper:
        if dihedral_a is None and dihedral_b is None:
            raise ValueError(
                f"Can't find parameters for dihedral with key {key}, got None for both dihedrals."
            )
        pairkey = to_pairkey(key)
        if dihedral_a is not None and dihedral_b is not None:
            if (
                pairkey[0] in merger.affected_interactions.atoms
                or pairkey[1] in merger.affected_interactions.atoms
            ):
                category = "morphed"
        elif dihedral_a is not None:
            category = "removed"
        else:
            category = "added"
        if dihedral_a == dihedral_b:
            return key, None, category, use_improper  # unchanged dihedral: nothing to do

    # mix every periodicity of A and B into one before-to-after dihedral
    keys_a = set(dihedral_a.dihedrals.keys()) if dihedral_a else set()
    keys_b = set(dihedral_b.dihedrals.keys()) if dihedral_b else set()
    merged = MultipleDihedrals(
        *key,
        funct=(
            FFFUNC["mult_improper_dihedral"]
            if use_improper
            else FFFUNC["mult_proper_dihedral"]
        ),
        dihedrals={},
    )
    for periodicity in sorted(
        keys_a | keys_b, key=lambda value: int(value) if value else -1
    ):
        interaction_a = dihedral_a.dihedrals.get(periodicity) if dihedral_a else None
        interaction_b = dihedral_b.dihedrals.get(periodicity) if dihedral_b else None
        merged.dihedrals[periodicity] = merger._interpolate_two_dihedrals(
            key=key,
            dihedral_a=interaction_a,
            dihedral_b=interaction_b,
            funct=merged.funct,
        )
    return key, merged, category, use_improper  # the main process writes the result


class MoleculeTypeMerger:
    """Takes two MoleculeTypes and joins them for a smooth free-energy like parameter transition simulation.

    mol_a is modified in place.
    """

    def __init__(
        self,
        mol_a: MoleculeType,
        mol_b: MoleculeType,
        ff: FF,
        morse_only: bool = False,
        retype_transferred_h: bool = True,
        parallel_dihedral_merge: bool = False,
        parallel_dihedral_workers: int | None = None,
        harmonic_acceptor_h_relax_term: bool = False,
        harmonic_acceptor_h_relax_k: float = 50000.0,
    ) -> None:
        self.mol_a = mol_a
        self.mol_b = mol_b
        self.ff = ff
        self.morse_only = morse_only
        self.retype_transferred_h = bool(retype_transferred_h)
        self.harmonic_acceptor_h_relax_term = bool(harmonic_acceptor_h_relax_term)
        self.harmonic_acceptor_h_relax_k = float(harmonic_acceptor_h_relax_k)
        # created bond key -> its state-B equilibrium length, filled by merge_bonds
        self._created_bond_lengths: dict[tuple[str, str], str] = {}
        self.parallel_dihedral_merge = bool(parallel_dihedral_merge)
        if parallel_dihedral_workers is None:  # no number in kimmdy.yml: take it from the environment
            try:
                parallel_dihedral_workers = int(
                    os.environ.get("KIMMDY_DIHEDRAL_WORKERS", "1")
                )
            except ValueError:
                logger.warning("Ignoring invalid KIMMDY_DIHEDRAL_WORKERS; using one worker")
                parallel_dihedral_workers = 1
        self.parallel_dihedral_workers = max(1, int(parallel_dihedral_workers))  # at least one helper
        self.transferred_h_type = self._find_lj_hydrogen_type()
        self.default_morse_well_depth = 300  # [kJ mol^-1]
        self.default_beta_for_lj = 19  # matches LJ steepness
        self.default_morse_dist_factor = 1.12  # sigmaij* 1.12 = LJ minimum
        logger.info(
            f"Using default morse well depth of {self.default_morse_well_depth} and beta of {self.default_beta_for_lj}"
        )
        logger.info(
            f"Using default morse distance factor of {self.default_morse_dist_factor}"
        )

        self.affected_interactions = AffectedInteractions()
        # Keep the physical endpoint graphs before mol_a is mutated into the
        # union topology used for slow growth.
        self._bond_graph_a = self._build_bond_graph(mol_a)
        self._bond_graph_b = self._build_bond_graph(mol_b)
        # hydroxyl or water H atoms (type HO or HW) whose bonded partners change
        # in this reaction, i.e. the moving H; used by _get_LJ_parameters
        self.moved_lj_free_hydrogens = {
            atom_nr
            for atom_nr, atom in mol_a.atoms.items()
            if atom.type in ("HO", "HW")
            and self._bond_graph_a[atom_nr] != self._bond_graph_b.get(atom_nr, set())
        }
        self._pairs_a = {to_pairkey(key) for key in mol_a.pairs}
        self._pairs_b = {to_pairkey(key) for key in mol_b.pairs}
        self._explicit_exclusions_a = self._expand_exclusions(mol_a)
        self._explicit_exclusions_b = self._expand_exclusions(mol_b)
        self.union_exclusion_helper_pairs: set[tuple[str, str]] = set()

    def merge(self):
        """modiefies mol_a, the reactive moleculetype of of top_a, in place"""
        logger.info(f"Merging topologies for slow growth")
        self.merge_atoms()
        self.merge_bonds()
        self.merge_angles()
        self.merge_dihedrals()
        self.merge_exclusions()

        if not self.morse_only:
            self.merge_pairs()
            self.add_helper_pairs()
            self.add_union_exclusion_helper_pairs()
            # sort pairs because their order can get messed up
            # by adding the helper pairs after the already existing pairs
            self.mol_a.pairs = dict(
                sorted(
                    self.mol_a.pairs.items(),
                    key=lambda item: (int(item[0][0]), int(item[0][1])),
                )
            )
            # sort exclusions after adding helper pairs because they
            # come with their own exclusions
            self.mol_a.exclusions = dict(
                sorted(
                    self.mol_a.exclusions.items(),
                    key=lambda item: [int(i) for i in item[0]],
                )
            )
        else:
            # same as in (keysA - keysB) in v1
            for key in self.affected_interactions.bonds.removed:
                self.mol_a.pairs.pop(key, None)
                self.mol_a.exclusions[key] = Exclusion(*key)

        # after the bonded terms are merged: every interaction that touches a
        # transferred H is written with explicit parameters by then, so the new
        # type needs no <type>-H1 entry in the force field
        self._retype_transferred_hydrogens()
        self._add_harmonic_acceptor_h_terms()

        logger.info(f"Finished merging topologies")

    @staticmethod
    def _build_bond_graph(mol: MoleculeType) -> dict[str, set[str]]:
        graph = {atom_nr: set() for atom_nr in mol.atoms}
        for atom_1, atom_2 in mol.bonds:
            graph[atom_1].add(atom_2)
            graph[atom_2].add(atom_1)
        return graph

    @staticmethod
    def _expand_exclusions(mol: MoleculeType) -> set[tuple[str, str]]:
        pairs: set[tuple[str, str]] = set()
        for exclusion in mol.exclusions.values():
            atom_nrs = exclusion.key()
            if len(atom_nrs) < 2:
                continue
            for excluded_atom in atom_nrs[1:]:
                pairs.add(to_pairkey((atom_nrs[0], excluded_atom)))
        return pairs

    @staticmethod
    def _nrexcl(mol: MoleculeType) -> int:
        try:
            return int(mol.nrexcl)
        except (TypeError, ValueError) as error:
            raise ValueError(
                f"Invalid nrexcl value {mol.nrexcl!r} for moleculetype "
                f"{mol.name!r}."
            ) from error

    @staticmethod
    def _atoms_within_bond_distance(
        graph: dict[str, set[str]],
        initial_atoms: set[str],
        max_distance: int,
    ) -> set[str]:
        atoms = set(initial_atoms)
        frontier = set(initial_atoms)
        for _ in range(max_distance):
            next_frontier = {
                neighbor
                for atom_nr in frontier
                for neighbor in graph.get(atom_nr, set())
                if neighbor not in atoms
            }
            if not next_frontier:
                break
            atoms.update(next_frontier)
            frontier = next_frontier
        return atoms

    @staticmethod
    def _pairs_within_bond_distance(
        graph: dict[str, set[str]],
        initial_atoms: set[str],
        max_distance: int,
    ) -> set[tuple[str, str]]:
        pairs: set[tuple[str, str]] = set()
        for atom_1 in initial_atoms:
            visited = {atom_1}
            frontier = {atom_1}
            for _ in range(max_distance):
                next_frontier = {
                    neighbor
                    for atom_nr in frontier
                    for neighbor in graph.get(atom_nr, set())
                    if neighbor not in visited
                }
                if not next_frontier:
                    break
                pairs.update(
                    to_pairkey((atom_1, atom_2)) for atom_2 in next_frontier
                )
                visited.update(next_frontier)
                frontier = next_frontier
        return pairs

    @staticmethod
    def _is_within_bond_distance(
        graph: dict[str, set[str]],
        atom_1: str,
        atom_2: str,
        max_distance: int,
    ) -> bool:
        if max_distance < 1:
            return False
        visited = {atom_1}
        frontier = {atom_1}
        for _ in range(max_distance):
            next_frontier = {
                neighbor
                for atom_nr in frontier
                for neighbor in graph.get(atom_nr, set())
                if neighbor not in visited
            }
            if atom_2 in next_frontier:
                return True
            if not next_frontier:
                return False
            visited.update(next_frontier)
            frontier = next_frontier
        return False

    def _nonbonded_state(
        self,
        pair_key: tuple[str, str],
        graph: dict[str, set[str]],
        explicit_pairs: set[tuple[str, str]],
        explicit_exclusions: set[tuple[str, str]],
        nrexcl: int,
    ) -> NonbondedState:
        if pair_key in explicit_pairs:
            return NonbondedState.OneFour
        if pair_key in explicit_exclusions or self._is_within_bond_distance(
            graph, *pair_key, max_distance=nrexcl
        ):
            return NonbondedState.Excluded
        return NonbondedState.Full

    def _get_explicit_MultipleDihedrals(
        self,
        key: tuple[str, str, str, str],
        use_state_b: bool,
        use_improper: bool = False,
        periodicity_max: int = 6,
    ) -> Optional[MultipleDihedrals]:
        """Takes a valid dihedral key and returns explicit
        dihedral parameters for a given topology
        """
        if use_improper:
            funct = FFFUNC["mult_improper_dihedral"]
            dihedraltypes = self.ff.improper_dihedraltypes
        else:
            dihedraltypes = self.ff.proper_dihedraltypes
            funct = FFFUNC["mult_proper_dihedral"]

        if use_state_b:
            mol = self.mol_b
        else:
            mol = self.mol_a

        if use_improper:
            dihedrals_in = mol.improper_dihedrals.get(key)
        else:
            dihedrals_in = mol.proper_dihedrals.get(key)

        if not dihedrals_in:
            return None

        if "" not in dihedrals_in.dihedrals.keys():
            # empty string means implicit parameters
            # not having that means the dihedral is already
            # explicitly parameterized
            return dihedrals_in

        type_key = []
        for id in key:
            t = mol.atoms[id].type
            type_key.append(t)

        multiple_dihedrals = MultipleDihedrals(*key, funct=funct, dihedrals={})
        for periodicity in range(1, periodicity_max + 1):
            p = str(periodicity)
            match_obj = match_atomic_item_to_atomic_type(
                id=type_key, types=dihedraltypes, periodicity=p
            )
            if match_obj:
                assert isinstance(match_obj, DihedralType)
                multiple_dihedrals.dihedrals[p] = Dihedral(
                    *key,
                    funct=funct,
                    c0=match_obj.c0,
                    c1=match_obj.c1,
                    periodicity=match_obj.periodicity,
                )

        if not multiple_dihedrals.dihedrals:
            return None

        return multiple_dihedrals

    def _get_LJ_parameters(
        self,
        id1: str,
        id2: str,
        use_state_b: bool,
        is_1_4: bool = False,
        lj_free_h_in_a: bool = False,
    ):
        """Calculate LJ terms sigma and epsilon from atom types.

        A hydroxyl or water H (type HO or HW) has no LJ size, so it is normally
        looked up as an H type with a size (H1). With lj_free_h_in_a and
        retype_transferred_h on, the moving H keeps its real, size-free type in
        state A instead, so a pair term with it starts at zero and grows in
        with lambda. Returns sigma and epsilon for the atom pair.
        """

        # defaults top line:
        # non-bonded function type; combination rule; generate pairs (no/yes); fudge LJ (); fudge QQ ()
        comb_rule = self.ff.defaults[0][1]
        fudgeLJ = float(self.ff.defaults[0][3])
        fudgeQQ = float(self.ff.defaults[0][4])  # not used for now

        if use_state_b:
            mol = self.mol_b
        else:
            mol = self.mol_a

        t1 = mol.atoms[id1].type
        t2 = mol.atoms[id2].type

        # Protect breaking/binding hydrogen types without LJ potential. Keep
        # the historical H1 lookup when the opt-in relaxation retyping is off.
        replacement_h_type = (
            self.transferred_h_type
            if self.retype_transferred_h and self.transferred_h_type is not None
            else "H1"
        )
        # with retype_transferred_h, the moving H keeps its size-free type in
        # state A of a pair term, like its atom type in the relax topology
        keep_real_a_type = (
            lj_free_h_in_a and self.retype_transferred_h and not use_state_b
        )
        if t1 in ["HW", "HO"] and not (
            keep_real_a_type and id1 in self.moved_lj_free_hydrogens
        ):
            t1 = replacement_h_type
        if t2 in ["HW", "HO"] and not (
            keep_real_a_type and id2 in self.moved_lj_free_hydrogens
        ):
            t2 = replacement_h_type

        type1 = self.ff.atomtypes[t1]
        type2 = self.ff.atomtypes[t2]

        # see <https://manual.gromacs.org/current/reference-manual/topologies/parameter-files.html#nbpar>
        if comb_rule == "1":
            sigma = np.sqrt(float(type1.sigma) * float(type2.sigma))
            epsilon = np.sqrt(float(type1.epsilon) * float(type2.epsilon))
        elif comb_rule == "2":
            sigma = 0.5 * (float(type1.sigma) + float(type2.sigma))
            epsilon = np.sqrt(float(type1.epsilon) * float(type2.epsilon))
        elif comb_rule == "3":
            sigma = np.sqrt(float(type1.sigma) * float(type2.sigma))
            epsilon = np.sqrt(float(type1.epsilon) * float(type2.epsilon))
        else:
            raise ValueError("Unknown combination rule of forcefield")

        if is_1_4:
            # scale interaction for 1-4 interactions
            epsilon = epsilon * fudgeLJ

        return sigma, epsilon

    def _get_morse_parameters(
        self, atomtypes: tuple[str, str], bond: Bond | None = None
    ):
        # fall back to defaults if no parameters are found
        well_depth = self.default_morse_well_depth
        beta = self.default_beta_for_lj
        type_key: tuple[str, str] = tuple(sorted(atomtypes))  # type: ignore
        d = DEFAULT_EDISSOC.get(type_key)
        if d is not None:
            well_depth = d
        k = bond.c1 if bond else None
        if k is not None:
            # like in gmx /src/gromacs/gmxpreprocess/tomorse.cpp
            # https://github.com/gromacs/gromacs/blob/beca834f6269e0f6028cc7b62a0cd72e9ff446cf/src/gromacs/gmxpreprocess/tomorse.cpp#L195
            beta = sqrt(float(k) / (2 * well_depth))
            # k = beta^2 * (2 * well_depth)
        return well_depth, beta

    # 400.0 19.0

    def _make_pair(
        self,
        id1: str,
        id2: str,
        transition: PairTransition,
        from_1_4: bool = False,
        to_1_4: bool = False,
        force_explicit: bool = False,
    ) -> Pair:
        """Generates morphing pair interaction

        If it is for a binding event, the pair is vanishing, as it will be an
        exclusion once bound. If a bond is breaking, the pair interaction
        is slowly turned on, as it was excluded previously.

        Parameters
        ----------
        id1
            Atom one
        id2
            Atom two
        transition
            Is the pair turning on or off or morphing
        from_1_4
            Is the pair a 1-4 interaction in the A state
        to_1_4
            Will the pair be a 1-4 interaction in the B state
        force_explicit
            Always write A/B LJ parameters, even when they are identical.

        Returns
        -------
        Pair
            Morphing pair
        """
        # state A: the moving H is size-free here (only pairs, not Morse bond lengths)
        sigmaij_a, epsilonij_a = self._get_LJ_parameters(
            id1=id1, id2=id2, use_state_b=False, is_1_4=from_1_4, lj_free_h_in_a=True
        )
        sigmaij_b, epsilonij_b = self._get_LJ_parameters(
            id1=id1, id2=id2, use_state_b=True, is_1_4=to_1_4
        )
        comment = ""

        # if all parameters are the same, the pair doesn't change
        # so we can keep the default from the forcefield
        if (
            transition == PairTransition.Morph
            and sigmaij_a == sigmaij_b
            and epsilonij_a == epsilonij_b
            and not force_explicit
        ):
            return Pair(id1, id2, funct=self.ff.defaults[0][0])

        # defaults are 0.0 for all c0-c3
        c_kwargs = dict(
            zip(
                [f"c{i}" for i in range(4)],
                [f"{0.0:.5f}" for _ in range(4)],
            )
        )
        if transition == PairTransition.Morph:
            comment += " morphing pair"
        elif transition == PairTransition.Vanish:
            comment += " vanishing pair"
        elif transition == PairTransition.Create:
            comment += " creating pair"

        if from_1_4:
            comment += " from 1-4 interaction in A"
        if to_1_4:
            comment += " to 1-4 interaction in B"

        if transition == PairTransition.Vanish or transition == PairTransition.Morph:
            # pair interaction turning off (bind or morph)
            c_kwargs["c0"] = f"{sigmaij_a:.5f}"
            c_kwargs["c1"] = f"{epsilonij_a:.5f}"
        if transition == PairTransition.Create or transition == PairTransition.Morph:
            # pair interaction turning on (break or morph)
            c_kwargs["c2"] = f"{sigmaij_b:.5f}"
            c_kwargs["c3"] = f"{epsilonij_b:.5f}"

        return Pair(
            id1, id2, funct=self.ff.defaults[0][0], **c_kwargs, comment=f"; {comment}"
        )

    def _interpolate_two_dihedrals(
        self,
        key: tuple[str, str, str, str],
        dihedral_a: Optional[Dihedral],
        dihedral_b: Optional[Dihedral],
        funct: str,
    ) -> Dihedral:
        """Merge one to two Dihedrals into a Dihedral in free-energy syntax.

        Only one of the dihedrals can be None at any given time.
        """
        if funct not in [
            FFFUNC["mult_proper_dihedral"],
            FFFUNC["mult_improper_dihedral"],
        ]:
            m = f"Can only interpolate between proper (type {FFFUNC['mult_proper_dihedral']}) or improper (type {FFFUNC['mult_proper_dihedral']}) dihedrals"
            logger.error(m)
            raise ValueError(m)

        # construct parameterized Dihedral
        if dihedral_a is not None and dihedral_b is not None:
            # both parameters are given
            # transition between parameters
            # this can be the case when both dihedrals exist
            # but the atomtypes of participating atoms change from A to B state
            # or the parameters are specified explicitly
            dihedralmerge = Dihedral(
                *key,
                funct=funct,
                c0=dihedral_a.c0,
                c1=dihedral_a.c1,
                periodicity=dihedral_a.periodicity,
                c3=dihedral_b.c0,
                c4=dihedral_b.c1,
                c5=dihedral_b.periodicity,
                comment=f"; morphing dihedral",
            )
        elif dihedral_a is not None:
            # dihedral only in A
            # vanish it by transitiononing to 0
            dihedralmerge = Dihedral(
                *key,
                funct=funct,
                c0=dihedral_a.c0,
                c1=dihedral_a.c1,
                periodicity=dihedral_a.periodicity,
                c3=dihedral_a.c0,
                c4="0.00",
                c5=dihedral_a.periodicity,
                comment=f"; vanishing dihedral",
            )
        elif dihedral_b is not None:
            # dihedral only in B
            # create it by transitioning from 0
            dihedralmerge = Dihedral(
                *key,
                funct=funct,
                c0=dihedral_b.c0,
                c1="0.00",
                periodicity=dihedral_b.periodicity,
                c3=dihedral_b.c0,
                c4=dihedral_b.c1,
                c5=dihedral_b.periodicity,
                comment=f"; creating dihedral",
            )
        else:
            m = f"Tried to merge two dihedrals of {key} but no parameterized dihedrals found!"
            logger.error(m)
            raise ValueError(m)
        return dihedralmerge

    def _find_lj_hydrogen_type(self) -> str | None:
        """Pick a hydrogen atom type of the force field that has a Lennard-Jones size.

        A hydroxyl or water H (type HO or HW) has no size of its own. When such
        an H moves onto another atom it needs a type with a size; this returns
        the first one found (H1 in Amber force fields), or None if there is none.
        """
        # try the usual names first, then every other type whose name starts with H
        preferred = ("H1", "H", "HC", "HGA")
        candidates = list(preferred) + sorted(
            name for name in self.ff.atomtypes if name not in preferred and name.upper().startswith("H")
        )

        # return the first candidate with a non-zero size (epsilon)
        for name in candidates:
            atomtype = self.ff.atomtypes.get(name)
            if atomtype is None:
                continue
            try:
                if float(atomtype.epsilon) != 0.0:
                    return name
            except (TypeError, ValueError):
                continue
        return None

    def _add_harmonic_acceptor_h_terms(self):
        """Add a spring between the moving H and its new atom, for the relax only.

        The moving H is the hydrogen that loses one bond and gains a new one in
        the reaction. The spring has the normal length of the new bond and its
        strength grows from 0 at the start of the relax to
        harmonic_acceptor_h_relax_k at the end, so it pulls the H onto its new
        atom while the old bond fades. It is written as bond type 6, a plain
        spring that does not switch off any other forces between the atoms.
        Does nothing when the option is off.
        """
        if not self.harmonic_acceptor_h_relax_term:
            return

        # atoms that lost a bond in this reaction (the moving H is one of them)
        broken_atoms = {nr for key in self.affected_interactions.bonds.removed for nr in key}

        # for every new bond whose H also lost a bond, add the spring
        added = 0
        for bond_key, b0 in self._created_bond_lengths.items():
            moving_h = [
                nr
                for nr in bond_key
                if nr in broken_atoms
                and nr in self.mol_a.atoms
                and self._is_hydrogen(self.mol_a.atoms[nr])
            ]
            if len(moving_h) != 1:
                continue
            # third key element keeps the helper next to the Morse bond of the
            # same atom pair; only the values are written to the topology
            self.mol_a.bonds[(*bond_key, "harmonic_acceptor_h")] = Bond(
                *bond_key,
                funct="6",
                c0=f"{float(b0):.5f}",  # length at the start of the relax
                c1=f"{0.0:.5f}",  # strength at the start: zero
                c2=f"{float(b0):.5f}",  # length at the end (same)
                c3=f"{self.harmonic_acceptor_h_relax_k:.5f}",  # strength at the end
                comment="; harmonic acceptor-H helper growing with lambda (relax only)",
            )
            added += 1
        logger.info("Added %d harmonic acceptor-H relax helper bond(s)", added)

    @staticmethod
    def _is_hydrogen(atom) -> bool:
        """Tell whether an atom is a hydrogen: by its mass, or by its type name if the mass is missing."""
        try:
            return float(atom.mass) < 1.5
        except (TypeError, ValueError):
            return str(atom.type).upper().startswith("H")

    def _retype_transferred_hydrogens(self):
        """Give a moved hydroxyl or water H a Lennard-Jones size at the end of the relax.

        A hydroxyl or water H (type HO or HW) has no size of its own. This
        keeps its old type at the start of the relax and sets the type with a
        size (from _find_lj_hydrogen_type) only for the end state, so the H's
        size grows in gradually. Only the relax topology is changed. Does
        nothing when the option is off.
        """
        if not self.retype_transferred_h:
            return

        # atoms whose bonds changed in this reaction
        candidate_atoms = set(self.affected_interactions.atoms)
        for pair in (
            self.affected_interactions.bonds.added
            | self.affected_interactions.bonds.removed
        ):
            candidate_atoms.update(pair)
        for bond_key in set(self.mol_a.bonds) ^ set(self.mol_b.bonds):
            candidate_atoms.update(bond_key)
        if self.transferred_h_type is None:
            logger.warning("No LJ-bearing hydrogen type found in the force field; skipping transferred-H retyping")
            return

        # retype every hydroxyl or water H among them, in the end state only
        retyped = 0
        for atom_nr in candidate_atoms:
            atom = self.mol_a.atoms.get(atom_nr)
            if atom is None or atom.type not in {"HW", "HO"}:
                continue
            implicit = self._implicit_interactions(atom_nr)
            if implicit:
                logger.warning(
                    "Not retyping transferred H %s: %s still take their parameters "
                    "from the force field by atom type, which has no entry for %s",
                    atom_nr,
                    ", ".join(implicit),
                    self.transferred_h_type,
                )
                continue
            # state B only: state A keeps the original (LJ-free) type, so an H
            # that sits in a tight hydrogen bond at the reaction frame does not
            # get its full Lennard-Jones radius at lambda 0 and blow the relax
            # up; its LJ parameters grow in linearly with lambda instead
            atom.typeB = self.transferred_h_type
            atom.comment = (
                f"; transferred H retyped {atom.type} -> "
                f"{self.transferred_h_type} in state B only"
            )
            retyped += 1
        if retyped:
            logger.info(
                "Retyped %d transferred H atom(s) in relaxation topology (state B only)",
                retyped,
            )

    def _implicit_interactions(self, atom_nr: str) -> list[str]:
        """List the terms of one atom that still get their numbers from the force field by type.

        Retyping an H is only safe when every term it is part of already has
        explicit numbers: the force field has no bond, angle or dihedral entry
        for the new H type on a hydroxyl oxygen. Returns names such as
        "bond 12-13" for every term that still has no numbers.
        """
        # look through bonds, angles, dihedrals and pairs that contain the atom
        implicit: list[str] = []
        sections = (
            ("bond", self.mol_a.bonds),
            ("angle", self.mol_a.angles),
            ("dihedral", self.mol_a.proper_dihedrals),
            ("improper dihedral", self.mol_a.improper_dihedrals),
            ("pair", self.mol_a.pairs),
        )
        for name, interactions in sections:
            for key, interaction in interactions.items():
                if atom_nr not in key:
                    continue
                # proper and improper dihedrals are grouped by periodicity
                terms = list(getattr(interaction, "dihedrals", {}).values()) or [
                    interaction
                ]
                if any(getattr(term, "c0", None) is None for term in terms):
                    implicit.append(f"{name} {'-'.join(key)}")
        return implicit

    def merge_atoms(self):
        # only iterating over the keys of mol_a
        # because KIMMDY can't add or remove atoms,
        # the the keys should be the same in a and b
        logger.info("Merging atoms")
        for nr in self.mol_a.atoms.keys():
            atomA = self.mol_a.atoms[nr]
            atomB = self.mol_b.atoms[nr]
            if atomA != atomB and not (
                self.morse_only and atomA.charge == atomB.charge
            ):
                # the simplified version doesn't change if the charges are the same
                # set B parameters first
                # to not overwrite them
                # keep the name of the atom
                # but transition the type, charge and mass
                atomA.typeB = deepcopy(atomB.type)
                atomA.chargeB = deepcopy(atomB.charge)
                atomA.massB = deepcopy(atomB.mass)
                atomA.comment = "; morphing atomtypes"
                self.affected_interactions.atoms.add(nr)

    def merge_bonds(self):
        logger.info("Merging bonds")
        bond_keys_a = set(self.mol_a.bonds.keys())
        bond_keys_b = set(self.mol_b.bonds.keys())
        bond_keys = bond_keys_a | bond_keys_b

        for bond_key in bond_keys:
            bond_a = self.mol_a.bonds.get(bond_key)
            bond_b = self.mol_b.bonds.get(bond_key)
            # the bond can either be in A, B or both
            # in the first two cases the respective other will be None

            if bond_a is None and bond_b is None:
                m = f"Can't find parameters for bond with key {bond_key}, got None for both bonds. This should be an impossible state."
                logger.error(m)
                raise ValueError(m)

            bond_a = get_explicit_or_type(
                bond_key, bond_a, self.ff.bondtypes, self.mol_a
            )
            bond_b = get_explicit_or_type(
                bond_key, bond_b, self.ff.bondtypes, self.mol_b
            )

            if isinstance(bond_a, BondType):
                bond_a = Bond(
                    bond_a.i, bond_a.j, funct=bond_a.funct, c0=bond_a.c0, c1=bond_a.c1
                )
            if isinstance(bond_b, BondType):
                bond_b = Bond(
                    bond_b.i, bond_b.j, funct=bond_b.funct, c0=bond_b.c0, c1=bond_b.c1
                )
            assert isinstance(bond_a, Bond) or bond_a is None
            assert isinstance(bond_b, Bond) or bond_b is None
            if bond_a == bond_b:
                # bonds are equal because none of their atomtypes changed,
                # and the bond is still there
                # nothing to be done
                continue

            atomtypes_a: tuple[str, str] = tuple([self.mol_a.atoms[atom_id].type for atom_id in bond_key])  # type: ignore
            atomtypes_b: tuple[str, str] = tuple([self.mol_b.atoms[atom_id].type for atom_id in bond_key])  # type: ignore

            sigmaij_a, epsilonij_a = self._get_LJ_parameters(
                *bond_key, use_state_b=False
            )
            sigmaij_b, epsilonij_b = self._get_LJ_parameters(
                *bond_key, use_state_b=True
            )

            depth_a, beta_a = self._get_morse_parameters(
                atomtypes=atomtypes_a, bond=bond_a
            )
            depth_b, beta_b = self._get_morse_parameters(
                atomtypes=atomtypes_b, bond=bond_b
            )

            # both bonds exist-> transition between parameters
            if bond_a is not None and bond_b is not None:
                if not (
                    isinstance(bond_a, (Bond, BondType))
                    and isinstance(bond_b, (Bond, BondType))
                ):
                    m = f"Can't find parameters for bond with key {bond_key}, got explicits/types: {bond_a} and {bond_b} for a bonds: {bond_a}, {bond_b}"
                    logger.error(m)
                    raise ValueError(m)
                if not (bond_a.funct == bond_b.funct == FFFUNC["harmonic_bond"]):
                    m = f"In slow-growth, bond functionals need to be harmonic, but {bond_key} is not. It is in A: {bond_a} and B: {bond_b}"
                    logger.error(m)
                    raise ValueError(m)
                self.mol_a.bonds[bond_key] = Bond(
                    *bond_key,
                    funct=bond_b.funct,
                    c0=bond_a.c0,
                    c1=bond_a.c1,
                    c2=bond_b.c0,
                    c3=bond_b.c1,
                    comment=f"; morphing harmonic bond",
                )
            elif isinstance(bond_a, (Bond, BondType)):
                # bond only exists in A -> vanish bond
                pairkey = to_pairkey(bond_key)
                self.affected_interactions.bonds.removed.add(pairkey)
                if self.morse_only:
                    # use a bond that pretends to be LJ via morse
                    # because `simplified` will not touch the pairs
                    self.mol_a.bonds[bond_key] = Bond(
                        *bond_key,
                        funct=FFFUNC["morse_bond"],
                        c0=bond_a.c0,  # original bond distance in A
                        c1=f"{depth_a:.5f}",  # original well depth in A
                        c2=f"{beta_a:.5f}",  # original steepness in A
                        c3=f"{sigmaij_b*self.default_morse_dist_factor:.5f}",  # r0 goes to eq distance for LJ
                        c4=f"{epsilonij_b:.5f}",  # well depth is epsilon
                        c5=f"{self.default_beta_for_lj:.5f}",  # default beta for LJ
                        comment=f"; vanishing bond from morse to LJ",
                    )
                else:
                    # just remove the bond, LJ interactions will be done by the Pair
                    self.mol_a.bonds[bond_key] = Bond(
                        *bond_key,
                        funct=FFFUNC["morse_bond"],
                        c0=bond_a.c0,  # original bond distance in A
                        c1=f"{depth_a:.5f}",  # original well depth in A
                        c2=f"{beta_a:.5f}",  # original steepness in A
                        c3=f"{bond_a.c0}",  # new r0 in B (=stays the same, because bond vanishes anyways)
                        c4=f"{0.0:.5f}",  # well depth goes to 0
                        c5=f"{self.default_beta_for_lj:.5f}",  # default beta for LJ
                        comment=f"; vanishing bond from morse to nothing",
                    )

            elif isinstance(bond_b, (Bond, BondType)):
                # bond only exists in B -> create bond
                if epsilonij_a == 0:
                    m = f"epsilonij_a is 0 for {bond_key}"
                    logger.error(m)
                    raise ValueError(m)
                pairkey = to_pairkey(bond_key)
                self.affected_interactions.bonds.added.add(pairkey)
                self._created_bond_lengths[bond_key] = bond_b.c0
                if self.morse_only:
                    # use a bond that pretends to be LJ via morse
                    self.mol_a.bonds[bond_key] = Bond(
                        *bond_key,
                        funct=FFFUNC["morse_bond"],
                        # starts further away, so it can pull in atoms of the bond
                        c0=f"{sigmaij_a*self.default_morse_dist_factor:.5f}",  # sigmaij* 1.12 = LJ minimum at lambda 0, where the atom types are still from A
                        c1=f"{epsilonij_a:.5f}",  # epsilon is well depth of LJ
                        c2=f"{self.default_beta_for_lj:.5f}",  # beta from default beta for LJ
                        c3=bond_b.c0,  # final bond distance
                        c4=f"{depth_b:.5f}",  # final well depth
                        c5=f"{beta_b:.5f}",  # final steepness (beta)
                        comment=f"; creating bond from LJ to morse",
                    )
                else:
                    # just add the bond, LJ interactions will be done by the Pair
                    self.mol_a.bonds[bond_key] = Bond(
                        *bond_key,
                        funct=FFFUNC["morse_bond"],
                        # starts further away, so it can pull in atoms of the bond
                        c0=f"{sigmaij_a*self.default_morse_dist_factor:.5f}",  # sigmaij* 1.12 = LJ minimum at lambda 0, where the atom types are still from A
                        c1=f"{0.0:.5f}",  # well depth from 0
                        c2=f"{0.0:.5f}",  # beta from 0 = flat
                        c3=bond_b.c0,  # final bond distance
                        c4=f"{depth_b:.5f}",  # final morse well depth
                        c5=f"{beta_b:.5f}",  # final morse steepness (beta)
                        comment=f"; creating bond from nothing to morse",
                    )

    @staticmethod
    def _urey_bradley_state(
        angle: Angle | AngleType,
        fallback_r13: str,
        disabled: bool = False,
    ) -> tuple[str, str, str, str]:
        if angle.c0 is None or angle.c1 is None:
            raise ValueError(f"Angle is missing explicit parameters: {angle}")
        if angle.funct == FFFUNC["harmonic_angle"]:
            return (
                angle.c0,
                "0.00" if disabled else angle.c1,
                fallback_r13,
                "0.00",
            )
        if angle.funct == FFFUNC["urey_bradley_angle"]:
            if angle.c2 is None or angle.c3 is None:
                raise ValueError(
                    f"Urey-Bradley angle is missing r13 or kUB: {angle}"
                )
            if disabled:
                return angle.c0, "0.00", angle.c2, "0.00"
            return angle.c0, angle.c1, angle.c2, angle.c3
        raise ValueError(
            "Slow growth supports only harmonic (1) and Urey-Bradley (5) "
            f"angles, got function {angle.funct}: {angle}"
        )

    def _merge_urey_bradley_angles(
        self,
        key: tuple[str, str, str],
        angle_a: Angle | AngleType | None,
        angle_b: Angle | AngleType | None,
    ) -> Angle:
        urey_angle = (
            angle_b
            if angle_b is not None
            and angle_b.funct == FFFUNC["urey_bradley_angle"]
            else angle_a
        )
        if urey_angle is None or urey_angle.c2 is None:
            raise ValueError(
                f"Cannot establish Urey-Bradley r13 for angle {key}: "
                f"A={angle_a}, B={angle_b}"
            )
        fallback_r13 = urey_angle.c2

        if angle_a is None:
            assert angle_b is not None
            state_a = self._urey_bradley_state(
                angle_b, fallback_r13, disabled=True
            )
        else:
            state_a = self._urey_bradley_state(angle_a, fallback_r13)

        if angle_b is None:
            assert angle_a is not None
            state_b = self._urey_bradley_state(
                angle_a, fallback_r13, disabled=True
            )
        else:
            state_b = self._urey_bradley_state(angle_b, fallback_r13)

        return Angle(
            *key,
            funct=FFFUNC["urey_bradley_angle"],
            c0=state_a[0],
            c1=state_a[1],
            c2=state_a[2],
            c3=state_a[3],
            c4=state_b[0],
            c5=state_b[1],
            c6=state_b[2],
            c7=state_b[3],
            comment="; morphing Urey-Bradley angle",
        )

    def merge_angles(self):
        logger.info("Merging angles")
        keys = set(self.mol_a.angles.keys()) | set(self.mol_b.angles.keys())
        for key in keys:
            angle_a = self.mol_a.angles.get(key)
            angle_b = self.mol_b.angles.get(key)

            if angle_a == angle_b:
                self.affected_interactions.angles.morphed.add(to_pairkey(key))
                continue

            parameterized_a = get_explicit_or_type(
                key, angle_a, self.ff.angletypes, self.mol_a
            )
            parameterized_b = get_explicit_or_type(
                key, angle_b, self.ff.angletypes, self.mol_b
            )
            if parameterized_a is not None and not isinstance(
                parameterized_a, (Angle, AngleType)
            ):
                raise ValueError(
                    f"Cannot resolve A-state angle parameters for {key}: "
                    f"{parameterized_a}"
                )
            if parameterized_b is not None and not isinstance(
                parameterized_b, (Angle, AngleType)
            ):
                raise ValueError(
                    f"Cannot resolve B-state angle parameters for {key}: "
                    f"{parameterized_b}"
                )

            functions = {
                angle.funct
                for angle in (parameterized_a, parameterized_b)
                if angle is not None
            }
            if FFFUNC["urey_bradley_angle"] in functions:
                self.mol_a.angles[key] = self._merge_urey_bradley_angles(
                    key, parameterized_a, parameterized_b
                )
                pairkey = to_pairkey(key)
                if parameterized_a is None:
                    self.affected_interactions.angles.added.add(pairkey)
                elif parameterized_b is None:
                    self.affected_interactions.angles.removed.add(pairkey)
                else:
                    self.affected_interactions.angles.morphed.add(pairkey)
                continue

            unsupported = functions - {FFFUNC["harmonic_angle"]}
            if unsupported:
                raise ValueError(
                    f"Unsupported slow-growth angle functions {unsupported} "
                    f"for angle {key}."
                )

            if parameterized_a is not None and parameterized_b is not None:
                self.mol_a.angles[key] = Angle(
                    *key,
                    funct=FFFUNC["harmonic_angle"],
                    c0=parameterized_a.c0,
                    c1=parameterized_a.c1,
                    c2=parameterized_b.c0,
                    c3=parameterized_b.c1,
                    comment="; morphing angle",
                )
                self.affected_interactions.angles.morphed.add(to_pairkey(key))
            elif parameterized_a is not None:
                self.mol_a.angles[key] = Angle(
                    *key,
                    funct=FFFUNC["harmonic_angle"],
                    c0=parameterized_a.c0,
                    c1=parameterized_a.c1,
                    c2=parameterized_a.c0,
                    c3="0.00",
                    comment="; vanishing angle",
                )
                self.affected_interactions.angles.removed.add(to_pairkey(key))
            elif parameterized_b is not None:
                self.mol_a.angles[key] = Angle(
                    *key,
                    funct=FFFUNC["harmonic_angle"],
                    c0=parameterized_b.c0,
                    c1="0.00",
                    c2=parameterized_b.c0,
                    c3=parameterized_b.c1,
                    comment="; creating angle",
                )
                self.affected_interactions.angles.added.add(to_pairkey(key))
            else:
                logger.warning("Could not parameterize angle %s.", key)

    def merge_dihedrals(self):
        # parallel code when switched on; otherwise the one-CPU code below runs
        if self.parallel_dihedral_merge and self.parallel_dihedral_workers > 1:
            self._merge_dihedrals_parallel()
            return
        # proper dihedrals
        # proper dihedrals have a nested structure and need a different treatment from bonds, angles and improper dihedrals
        logger.info("Merging dihedrals")
        keys = set(self.mol_a.proper_dihedrals.keys()) | set(
            self.mol_b.proper_dihedrals.keys()
        )
        for key in keys:
            multiple_dihedralsA = self._get_explicit_MultipleDihedrals(
                key=key, use_state_b=False, use_improper=False
            )
            multiple_dihedralsB = self._get_explicit_MultipleDihedrals(
                key=key, use_state_b=True, use_improper=False
            )

            if multiple_dihedralsA is None and multiple_dihedralsB is None:
                m = f"Can't find parameters for dihedral with key {key}, got None for both dihedrals. This should be an impossible state."
                logger.error(m)
                raise ValueError(m)
            elif multiple_dihedralsA is not None and multiple_dihedralsB is not None:
                # dihedrals exist in both A and B -> transition between parameters for the 1-4 interactions
                # in case their atomtypes changed
                pairkey = to_pairkey(key)
                if (
                    pairkey[0] in self.affected_interactions.atoms
                    or pairkey[1] in self.affected_interactions.atoms
                ):
                    self.affected_interactions.dihedrals.morphed.add(pairkey)
            elif isinstance(multiple_dihedralsA, MultipleDihedrals):
                # dihedral only exists in A -> vanishing dihedral
                # pair transitions from half strength 1-4 interaction to full strength
                pairkey = to_pairkey(key)
                self.affected_interactions.dihedrals.removed.add(pairkey)
            elif isinstance(multiple_dihedralsB, MultipleDihedrals):
                # dihedral only exists in B -> creating dihedral
                # corresponding pair transitions from full strength to half strength 1-4 interaction
                pairkey = to_pairkey(key)
                self.affected_interactions.dihedrals.added.add(pairkey)

            if multiple_dihedralsA == multiple_dihedralsB:
                # dihedrals are exactly the same, nothing further to be done
                continue

            keysA = (
                set(multiple_dihedralsA.dihedrals.keys())
                if multiple_dihedralsA
                else set()
            )
            keysB = (
                set(multiple_dihedralsB.dihedrals.keys())
                if multiple_dihedralsB
                else set()
            )

            self.mol_a.proper_dihedrals[key] = MultipleDihedrals(
                *key, funct=FFFUNC["mult_proper_dihedral"], dihedrals={}
            )

            periodicity_keys = keysA | keysB
            for periodicity_key in periodicity_keys:
                interactionA = (
                    multiple_dihedralsA.dihedrals.get(periodicity_key)
                    if multiple_dihedralsA
                    else None
                )
                interactionB = (
                    multiple_dihedralsB.dihedrals.get(periodicity_key)
                    if multiple_dihedralsB
                    else None
                )
                self.mol_a.proper_dihedrals[key].dihedrals[periodicity_key] = (
                    self._interpolate_two_dihedrals(
                        key=key,
                        dihedral_a=interactionA,
                        dihedral_b=interactionB,
                        funct=FFFUNC["mult_proper_dihedral"],
                    )
                )

        # improper dihedrals
        keys = set(self.mol_a.improper_dihedrals.keys()) | set(
            self.mol_b.improper_dihedrals.keys()
        )
        for key in keys:
            multiple_dihedralsA = self._get_explicit_MultipleDihedrals(
                key=key, use_state_b=False, use_improper=True
            )
            multiple_dihedralsB = self._get_explicit_MultipleDihedrals(
                key=key, use_state_b=True, use_improper=True
            )
            keysA = (
                set(multiple_dihedralsA.dihedrals.keys())
                if multiple_dihedralsA
                else set()
            )
            keysB = (
                set(multiple_dihedralsB.dihedrals.keys())
                if multiple_dihedralsB
                else set()
            )

            self.mol_a.improper_dihedrals[key] = MultipleDihedrals(
                *key, funct=FFFUNC["mult_improper_dihedral"], dihedrals={}
            )

            periodicity_keys = keysA | keysB
            for periodicity_key in periodicity_keys:
                assert isinstance(periodicity_key, str)
                interactionA = (
                    multiple_dihedralsA.dihedrals.get(periodicity_key)
                    if multiple_dihedralsA
                    else None
                )
                interactionB = (
                    multiple_dihedralsB.dihedrals.get(periodicity_key)
                    if multiple_dihedralsB
                    else None
                )

                self.mol_a.improper_dihedrals[key].dihedrals[periodicity_key] = (
                    self._interpolate_two_dihedrals(
                        key=key,
                        dihedral_a=interactionA,
                        dihedral_b=interactionB,
                        funct=FFFUNC["mult_improper_dihedral"],
                    )
                )

    def _merge_dihedrals_parallel(self):
        """Do the work of merge_dihedrals on several CPUs at once.

        Every dihedral is sent to _merge_dihedral_worker in a helper process.
        They are sent in a fixed order (sorted by atom number) and the results
        come back in the same order, so every run gives the same topology.
        The main process then writes the results into the relax topology.
        If helper processes cannot be started, the normal one-CPU code is used.
        """
        global _PARALLEL_DIHEDRAL_MERGER
        logger.info(
            "Merging dihedrals in parallel with %d workers",
            self.parallel_dihedral_workers,
        )
        # list all dihedrals, sorted by atom number so the order never changes
        proper_keys = sorted(
            set(self.mol_a.proper_dihedrals) | set(self.mol_b.proper_dihedrals),
            key=lambda key: tuple(int(atom) for atom in key),
        )
        improper_keys = sorted(
            set(self.mol_a.improper_dihedrals) | set(self.mol_b.improper_dihedrals),
            key=lambda key: tuple(int(atom) for atom in key),
        )
        payload = [(key, False) for key in proper_keys] + [
            (key, True) for key in improper_keys
        ]
        if not payload:
            return
        _PARALLEL_DIHEDRAL_MERGER = self  # let the helpers read this merger
        try:
            context = mp.get_context("fork")
        except ValueError:
            logger.warning("Fork workers are unavailable; using serial dihedral merging")
            _PARALLEL_DIHEDRAL_MERGER = None
            self.parallel_dihedral_merge = False
            self.merge_dihedrals()  # helpers not possible: use the normal one-CPU code
            return
        try:
            with ProcessPoolExecutor(
                max_workers=self.parallel_dihedral_workers, mp_context=context
            ) as executor:
                # results come back in the order they were sent
                results = list(executor.map(_merge_dihedral_worker, payload))
        finally:
            _PARALLEL_DIHEDRAL_MERGER = None  # clean up when the helpers are done

        # write each result into the topology, in the same sorted order
        for key, merged, category, use_improper in results:
            if category == "morphed":
                self.affected_interactions.dihedrals.morphed.add(to_pairkey(key))
            elif category == "added":
                self.affected_interactions.dihedrals.added.add(to_pairkey(key))
            elif category == "removed":
                self.affected_interactions.dihedrals.removed.add(to_pairkey(key))
            if merged is not None:
                target = (
                    self.mol_a.improper_dihedrals
                    if use_improper
                    else self.mol_a.proper_dihedrals
                )
                target[key] = merged

    def merge_pairs(self):
        """Merge pairs that are from the respective pairs sections of the topologies"""
        logger.info("Merging pairs")
        new_pairs = {}
        keys = set(self.mol_a.pairs.keys()) | set(self.mol_b.pairs.keys())
        for key in keys:
            pair_a = self.mol_a.pairs.get(key)
            pair_b = self.mol_b.pairs.get(key)

            if pair_a is not None and pair_b is not None:
                if is_parameterized(pair_a) or is_parameterized(pair_b):
                    new_pairs[key] = self._make_pair(
                        key[0], key[1], PairTransition.Morph, from_1_4=True, to_1_4=True
                    )
                else:
                    new_pairs[key] = Pair(
                        key[0],
                        key[1],
                        funct=self.ff.defaults[0][0],
                    )
            elif pair_a is not None:
                new_pairs[key] = self._make_pair(
                    key[0], key[1], PairTransition.Morph, from_1_4=True
                )
            elif pair_b is not None:
                new_pairs[key] = self._make_pair(
                    key[0], key[1], PairTransition.Morph, to_1_4=True
                )

        # NOTE: this reqires re-linking top.pairs to top.reactive_molecule.pairs
        # after wards, because the pairs are not the same objects anymore
        self.mol_a.pairs = {}
        for key, pair in new_pairs.items():
            self.mol_a.pairs[key] = pair

    def merge_exclusions(self):
        """Merge exclusions by adding the exclusions from the B state to the A state"""
        logger.info("Merging exclusions")
        for key in self.mol_b.exclusions.keys():
            self.mol_a.exclusions[key] = Exclusion(*key)

    def add_helper_pairs(self):
        """
        growing/shrinking pairs are used to turn on/off
        LJ / non-bonded interactions.
        If it is for a binding event, the pair is vanishing, as it will be an
        exclusion once bound. If a bond is breaking, the pair interaction
        is slowly turned on, as it was excluded previously.
        """
        logger.info(f"Adding helper pairs")
        # pairs that stem from bonds breaking or forming
        # collected by looking at the changing bonds, angles and dihedrals (=higher order interactions)
        self.helper_pairs: dict[tuple[str, str], Pair] = {}

        # a bond is either added or removed, but can't be both.
        # but for higher order interactions the end-atoms (of the angle or dihedral) can be involved
        # in multiple changing higher order interactions
        # e.g. two atoms may loose a dihedral between them and gain a different one
        added_bonds = (
            self.affected_interactions.bonds.added
            - self.affected_interactions.bonds.removed
        )
        removed_bonds = (
            self.affected_interactions.bonds.removed
            - self.affected_interactions.bonds.added
        )
        added_angles = (
            self.affected_interactions.angles.added
            - self.affected_interactions.angles.removed
        )
        removed_angles = (
            self.affected_interactions.angles.removed
            - self.affected_interactions.angles.added
        )
        swapping_angles = (
            self.affected_interactions.angles.added
            & self.affected_interactions.angles.removed
        )
        morphing_angles = self.affected_interactions.angles.morphed
        all_angles = {
            to_pairkey(key)
            for key in set(self.mol_a.angles.keys()) | set(self.mol_b.angles.keys())
        }
        added_dihedrals = (
            self.affected_interactions.dihedrals.added
            - self.affected_interactions.dihedrals.removed
        )
        removed_dihedrals = (
            self.affected_interactions.dihedrals.removed
            - self.affected_interactions.dihedrals.added
        )
        swapping_dihedrals = (
            self.affected_interactions.dihedrals.added
            & self.affected_interactions.dihedrals.removed
        )
        morphing_dihedrals = self.affected_interactions.dihedrals.morphed

        for pair_key in added_bonds:
            if pair_key in removed_dihedrals:
                # the bond is added, but the atoms where in a dihedral in A
                # so they had 1-4 interactions in A and will be excluded in B
                self.helper_pairs[pair_key] = self._make_pair(
                    *pair_key, transition=PairTransition.Vanish, from_1_4=True
                )
            if (
                pair_key
                not in removed_angles
                | removed_dihedrals
                | morphing_angles
                | morphing_dihedrals
            ):
                # the bond is added and the atoms didn't see each other in A
                # so the full strength LJ interaction is turned off to 0
                self.helper_pairs[pair_key] = self._make_pair(
                    *pair_key, transition=PairTransition.Vanish
                )

        for pair_key in removed_bonds:
            if (
                pair_key
                not in added_angles
                | added_dihedrals
                | morphing_angles
                | morphing_angles
            ):
                # the bond is removed and the atoms will not see each other in B
                # so the pair interaction is turned on (at full strength)
                self.helper_pairs[pair_key] = self._make_pair(
                    *pair_key, transition=PairTransition.Create
                )

        # pairs of atoms that end up with an angle between them where there was none before
        for pairkey in added_angles:
            if pairkey in removed_dihedrals:
                # these atoms shared a dihedral in A, so they start as a 1-4 interaction
                self.helper_pairs[pairkey] = self._make_pair(
                    *pairkey, transition=PairTransition.Vanish, from_1_4=True
                )
            if (
                pairkey
                not in removed_bonds
                | removed_dihedrals
                | morphing_dihedrals
                | morphing_angles
            ):
                # the angle is added and the atoms saw each other in A, LJ is turned off
                self.helper_pairs[pairkey] = self._make_pair(
                    *pairkey, transition=PairTransition.Vanish
                )

        # pairs that end up with no angle between them where there was one before
        for pairkey in removed_angles:
            if (
                not pairkey
                in added_bonds | added_dihedrals | morphing_dihedrals | morphing_angles
            ):
                # the atoms where in an angle A and will not see each other in B
                # we morph from 0 to full strength
                self.helper_pairs[pairkey] = self._make_pair(
                    *pairkey, transition=PairTransition.Create
                )

        # atoms that end up with a dihedral between them where there was none before
        for pairkey in added_dihedrals:
            if pairkey in removed_bonds | removed_angles:
                # if the atoms had some lower order interaction in A, they where excluded before
                # and should transition from 0 to 1-4
                self.helper_pairs[pairkey] = self._make_pair(
                    *pairkey, transition=PairTransition.Create, to_1_4=True
                )
            if pairkey not in removed_bonds | removed_angles | morphing_angles:
                # if the atoms didn't see each other in A, we morph from full strength LJ to 1-4
                # dihedral only exists in B -> creating dihedral
                # corresponding pair transitions from full strength to half strength 1-4 interaction
                self.helper_pairs[pairkey] = self._make_pair(
                    *pairkey, transition=PairTransition.Morph, to_1_4=True
                )

        # atoms that end up with no dihedral between them where there was one before
        for pairkey in removed_dihedrals:
            if not pairkey in added_bonds | added_angles:
                # dihedral only exists in A -> vanishing dihedral
                # pair transitions from half strength 1-4 interaction to full strength
                # but only if it's not already handled by a bond or angle
                self.helper_pairs[pairkey] = self._make_pair(
                    *pairkey, transition=PairTransition.Morph, from_1_4=True
                )

        for pairkey in morphing_dihedrals:
            if pairkey not in all_angles:
                # dihedrals exist in both A and B -> transition between parameters for the 1-4 interactions
                # in case their atomtypes changed
                self.helper_pairs[pairkey] = self._make_pair(
                    *pairkey,
                    transition=PairTransition.Morph,
                    from_1_4=True,
                    to_1_4=True,
                )

        # these pairs of atoms are first involved on one dihedral and end up still involved together in a dihedral
        # but it's a different dihedral
        for pairkey in swapping_dihedrals:
            if pairkey not in swapping_angles:
                # dihedrals exist in both A and B -> transition between parameters for the 1-4 interactions
                self.helper_pairs[pairkey] = self._make_pair(
                    *pairkey,
                    transition=PairTransition.Morph,
                    from_1_4=True,
                    to_1_4=True,
                )

        for key, pair in self.helper_pairs.items():
            if pair is not None:
                self.mol_a.pairs[key] = pair

            # add general exclusions for each pair
            self.mol_a.exclusions[key] = Exclusion(*key)

    def add_union_exclusion_helper_pairs(self):
        """Restore endpoint interactions hidden by the merged bond graph.

        Vanishing and creating bonds coexist in the slow-growth topology.
        Together they can form short paths that exist in neither physical
        endpoint. GROMACS applies nrexcl to that union graph, so these paths
        can accidentally suppress otherwise full nonbonded interactions.
        """
        changed_bonds = (
            self.affected_interactions.bonds.added
            | self.affected_interactions.bonds.removed
        )
        if not changed_bonds:
            return

        merged_graph = self._build_bond_graph(self.mol_a)
        merged_nrexcl = self._nrexcl(self.mol_a)
        changed_atoms = {atom_nr for bond in changed_bonds for atom_nr in bond}
        local_atoms = self._atoms_within_bond_distance(
            merged_graph,
            changed_atoms,
            merged_nrexcl,
        )
        merged_excluded_pairs = self._pairs_within_bond_distance(
            merged_graph,
            local_atoms,
            merged_nrexcl,
        )

        nrexcl_a = self._nrexcl(self.mol_a)
        nrexcl_b = self._nrexcl(self.mol_b)
        for pair_key in sorted(
            merged_excluded_pairs,
            key=lambda pair: (int(pair[0]), int(pair[1])),
        ):
            # Explicit endpoint pairs and the established helper-pair logic
            # have already been merged and remain authoritative.
            if pair_key in self.mol_a.pairs:
                continue

            state_a = self._nonbonded_state(
                pair_key,
                self._bond_graph_a,
                self._pairs_a,
                self._explicit_exclusions_a,
                nrexcl_a,
            )
            state_b = self._nonbonded_state(
                pair_key,
                self._bond_graph_b,
                self._pairs_b,
                self._explicit_exclusions_b,
                nrexcl_b,
            )
            if (
                state_a == NonbondedState.Excluded
                and state_b == NonbondedState.Excluded
            ):
                continue

            if state_a == NonbondedState.Excluded:
                transition = PairTransition.Create
            elif state_b == NonbondedState.Excluded:
                transition = PairTransition.Vanish
            else:
                transition = PairTransition.Morph

            pair = self._make_pair(
                *pair_key,
                transition=transition,
                from_1_4=state_a == NonbondedState.OneFour,
                to_1_4=state_b == NonbondedState.OneFour,
                force_explicit=(
                    state_a == NonbondedState.Full
                    or state_b == NonbondedState.Full
                ),
            )
            pair.comment = (
                f"{pair.comment or ''} "
                "; restoring interaction excluded only by merged bond graph"
            ).strip()
            self.helper_pairs[pair_key] = pair
            self.union_exclusion_helper_pairs.add(pair_key)
            self.mol_a.pairs[pair_key] = pair
            self.mol_a.exclusions[pair_key] = Exclusion(*pair_key)

        if self.union_exclusion_helper_pairs:
            logger.info(
                "Added %d helper pairs for interactions excluded only by the "
                "merged slow-growth bond graph.",
                len(self.union_exclusion_helper_pairs),
            )


def merge_top_slow_growth(
    top_a: Topology,
    top_b: Topology,
    morse_only: bool = False,
    retype_transferred_h: bool = True,
    parallel_dihedral_merge: bool = False,
    parallel_dihedral_workers: int | None = None,
    harmonic_acceptor_h_relax_term: bool = False,
    harmonic_acceptor_h_relax_k: float = 50000.0,
):
    """Takes two Topologies and joins them for a smooth free-energy like parameter transition simulation.
    Modifies topB in place.
    All changes are contained to the `Reactive` moleculetype.
    """

    MoleculeTypeMerger(
        mol_a=top_a.reactive_molecule,
        mol_b=top_b.reactive_molecule,
        ff=top_a.ff,
        morse_only=morse_only,
        retype_transferred_h=retype_transferred_h,
        parallel_dihedral_merge=parallel_dihedral_merge,
        parallel_dihedral_workers=parallel_dihedral_workers,
        harmonic_acceptor_h_relax_term=harmonic_acceptor_h_relax_term,
        harmonic_acceptor_h_relax_k=harmonic_acceptor_h_relax_k,
    ).merge()

    # because we clear some sections of the moleculetype,
    # like self.mol_a.pairs = {},
    # breaking the link between
    # `top.pairs` and `top.reactive_molecule.pairs`,
    # we need to re-link them for convenient access
    top_a._link_atomics()

    return


def break_bond_plumed(
    files: TaskFiles, breakpair: tuple[str, str], newplumed: Path
) -> None:
    """Break bond in plumed configuration file.

    Parameters
    ----------
    files:

    breakpair:
    """

    plumed_path = files.input["plumed"]
    # if break_bond_plumed is called multiple times in one _apply_recipe task
    if "plumed" in files.output.keys():
        plumed_path = files.output["plumed"]
    assert plumed_path is not None, "plumed_path should not be None"
    plumed_dict = read_plumed(plumed_path)

    files.output["plumed"] = newplumed

    logger.debug(
        f"Reading: {files.input['plumed']} and writing modified plumed input to "
        f"{files.output['plumed']}."
    )

    new_distances = {}
    broken_distances = []
    for k, v in plumed_dict["labeled_action"].items():
        if all(x in v["atoms"] for x in breakpair):
            broken_distances.append(k)
        else:
            new_distances[k] = v

    plumed_dict["labeled_action"] = new_distances

    for line in plumed_dict["prints"]:
        line["ARG"] = [id for id in line["ARG"] if id not in broken_distances]
        line["FILE"] = files.input["plumed_out"]

    write_plumed(plumed_dict, newplumed)
