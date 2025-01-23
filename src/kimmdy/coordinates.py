"""coordinate, topology and plumed modification functions"""

from enum import Enum
import logging
from copy import deepcopy
from pathlib import Path
from typing import Optional, Union
from math import sqrt

import MDAnalysis as mda
import numpy as np

from kimmdy.constants import DEFAULT_EDISSOC, REACTIVE_MOLECULEYPE, FFFUNC
from kimmdy.parsing import read_plumed, write_plumed
from kimmdy.recipe import Place
from kimmdy.tasks import TaskFiles
from kimmdy.topology.atomic import (
    Angle,
    AngleType,
    Bond,
    BondType,
    Dihedral,
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
    use_state_b: bool = False,
) -> Union[Interaction, AtomicType, None]:
    """Takes an Interaction and associated key, InteractionTypes, Topology
    and Periodicity (for dihedrals) and returns an object with the parameters of this Interaction
    """
    if not interaction:
        return None

    if is_parameterized(interaction):
        return interaction

    if use_state_b:
        type_key = []
        for id in key:
            t = getattr(mol.atoms[id], "typeB", None)
            t = t if t else mol.atoms[id].type
            type_key.append(t)
    else:
        type_key = [mol.atoms[x].type for x in key]

    atomic_type = match_atomic_item_to_atomic_type(
        type_key, interaction_types, periodicity
    )

    if atomic_type is None:
        m = f"Could not find explicit parameters for {[type_key,interaction]} in the FF. Consider setting the parameterizer to `grappa` to generate missing parameters in the fly."
        logger.error(m)
        raise ValueError(m)

    return atomic_type


class MoleculeTypeMerger:
    """Takes two Topologies and joins them for a smooth free-energy like parameter transition simulation"""

    def __init__(
        self,
        mol_a: MoleculeType,
        mol_b: MoleculeType,
        ff: FF,
        enable_morph_pairs: bool,
    ) -> None:
        self.mol_a = mol_a
        self.mol_b = mol_b
        self.ff = ff
        self.enable_morph_pairs = enable_morph_pairs
        self.default_morse_well_depth = 300  # [kJ mol^-1]
        self.default_beta_for_lj = 19  # matches LJ steepness
        self.stateA_morse_dist_factor = 1.12  # sigmaij* 1.12 = LJ minimum
        self.stateB_morse_dist_factor = 1.12  # sigmaij* 1.12 = LJ minimum
        logger.info(
            f"Using default morse well depth of {self.default_morse_well_depth} and beta of {self.default_beta_for_lj}"
        )
        logger.info(
            f"Using default morse distance factor of {self.stateA_morse_dist_factor}"
        )

        self.affected_pairs = {}

    def merge(self) -> MoleculeType:
        """modiefies mol_b, the reactive moleculetype of of top_b, in place"""
        logger.info(f"Merging topologies for slow growth")
        self.merge_atoms()
        self.merge_bonds()
        self.merge_angles()
        self.merge_dihedrals()
        self.merge_pairs()
        self.mol_b.find_radicals()
        return self.mol_b

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
            if use_improper:
                dihedrals_in = self.mol_b.improper_dihedrals.get(key)
            else:
                dihedrals_in = self.mol_b.proper_dihedrals.get(key)
        else:
            if use_improper:
                dihedrals_in = self.mol_a.improper_dihedrals.get(key)
            else:
                dihedrals_in = self.mol_a.proper_dihedrals.get(key)

        if not dihedrals_in:
            return None

        # not using implicit parameters
        # if "" not in dihedrals_in.dihedrals.keys():
        #     # empty string means implicit parameters
        #     return dihedrals_in

        type_key = []
        for id in key:
            if use_state_b:
                t = getattr(self.mol_b.atoms[id], "typeB", None)
                t = t if t else self.mol_b.atoms[id].type
            else:
                t = self.mol_a.atoms[id].type

            type_key.append(t)

        multiple_dihedrals = MultipleDihedrals(*key, funct=funct, dihedrals={})
        for periodicity in range(1, periodicity_max + 1):
            match_obj = match_atomic_item_to_atomic_type(
                id=type_key, types=dihedraltypes, periodicity=str(periodicity)
            )
            if match_obj:
                assert isinstance(match_obj, DihedralType)
                multiple_dihedrals.dihedrals[str(periodicity)] = Dihedral(
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
        self, id1: str, id2: str, use_state_b: bool, is_1_4: bool = False
    ):
        """Calculate LJ terms sigma and epsilon from atom types"""

        # defaults top line:
        # non-bonded function type; combination rule; generate pairs (no/yes); fudge LJ (); fudge QQ ()
        comb_rule = self.ff.defaults[0][1]
        fudgeLJ = float(self.ff.defaults[0][3])
        fudgeQQ = float(self.ff.defaults[0][4])  # not used for now

        if use_state_b:
            t1 = getattr(self.mol_b.atoms[id1], "typeB", None)
            t2 = getattr(self.mol_b.atoms[id2], "typeB", None)
            t1 = t1 if t1 else self.mol_b.atoms[id1].type
            t2 = t2 if t2 else self.mol_b.atoms[id2].type
            type1 = self.ff.atomtypes[t1]
            type2 = self.ff.atomtypes[t2]
        else:
            t1 = self.mol_a.atoms[id1].type
            t2 = self.mol_a.atoms[id2].type

        # amber fix for breaking/binding atom types without LJ potential
        if t1 in ["HW", "HO"]:
            logger.debug(f"amber fix for {id1} of type {t1} typeA set to H1")
            t1 = "H1"
        if t2 in ["HW", "HO"]:
            logger.debug(f"amber fix for {id2} of type {t2} typeA set to H1")
            t2 = "H1"

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

        logger.debug(
            f"Got LJ parameters for {id1}-{id2} of types {type1.type}-{type2.type}: sigma {sigma} and epsilon {epsilon} from mol_b ({use_state_b})"
        )

        if is_1_4:
            # apply fudgeLJ
            epsilon = epsilon * fudgeLJ

        return sigma, epsilon

    def _get_morse_parameters(
        self, atomtypes: tuple[str, str], bond: Bond | None = None
    ):
        well_depth = self.default_morse_well_depth
        type_key: tuple[str, str] = tuple(sorted(atomtypes))  # type: ignore
        d = DEFAULT_EDISSOC.get(type_key)
        if d is not None:
            well_depth = d
        k = bond.c1 if bond else None
        if k is not None:
            # like in gmx /src/gromacs/gmxpreprocess/tomorse.cpp
            beta = sqrt(float(k) / (2 * well_depth))
        else:
            beta = self.default_beta_for_lj
        return well_depth, beta

    def _make_pair(
        self,
        id1: str,
        id2: str,
        transition: PairTransition,
        from_1_4: bool = False,
        to_1_4: bool = False,
    ) -> Pair:
        """Generates morphing pair interaction

        If it is for a binding event, the pair is vanishing, as it will be an
        exclusion once bound. If a bond is breaking, the pair interaction
        is slowly turned on, as it was excluded previously.

        Parameters
        ----------
        id1 : str
            Atom one
        id2 : str
            Atom two

        Returns
        -------
        Pair
            Morphing pair
        """
        sigmaij_a, epsilonij_a = self._get_LJ_parameters(
            id1=id1, id2=id2, use_state_b=False, is_1_4=from_1_4
        )
        sigmaij_b, epsilonij_b = self._get_LJ_parameters(
            id1=id1, id2=id2, use_state_b=True, is_1_4=to_1_4
        )

        # if all parameters are the same, the pair doesn't change
        # so we can keep the default from the forcefield
        if (
            transition == PairTransition.Morph
            and sigmaij_a == sigmaij_b
            and epsilonij_a == epsilonij_b
        ):
            logger.debug(f"Pair {id1}-{id2} doesn't change, keeping default")
            return Pair(id1, id2, funct=self.ff.defaults[0][0])

        # defaults are 0.0 for all c0-c3
        c_kwargs = dict(
            zip(
                [f"c{i}" for i in range(4)],
                [f"{0.0:.5f}" for _ in range(4)],
            )
        )
        if transition == PairTransition.Vanish or transition == PairTransition.Morph:
            # pair interaction turning off (bind or morph)
            c_kwargs["c0"] = f"{sigmaij_a:.5f}"
            c_kwargs["c1"] = f"{epsilonij_a:.5f}"
        if transition == PairTransition.Create or transition == PairTransition.Morph:
            # pair interaction turning on (break or morph)
            c_kwargs["c2"] = f"{sigmaij_b:.5f}"
            c_kwargs["c3"] = f"{epsilonij_b:.5f}"

        return Pair(id1, id2, funct=self.ff.defaults[0][0], **c_kwargs)

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
            dihedralmerge = Dihedral(
                *key,
                funct=funct,
                c0=dihedral_a.c0,
                c1=dihedral_a.c1,
                periodicity=dihedral_a.periodicity,
                c3=dihedral_b.c0,
                c4=dihedral_b.c1,
                c5=dihedral_b.periodicity,
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
            )
        else:
            m = f"Tried to merge two dihedrals of {key} but no parameterized dihedrals found!"
            logger.error(m)
            raise ValueError(m)
        return dihedralmerge

    def merge_atoms(self):
        # only iterating over the keys of mol_a
        # because KIMMDY can't add or remove atoms,
        # the the keys should be the same in a and b
        for nr in self.mol_a.atoms.keys():
            atomA = self.mol_a.atoms[nr]
            atomB = self.mol_b.atoms[nr]
            if not atomA.eq_by_params(atomB):
                # set B parameters first
                # to not overwrite them
                atomB.typeB = deepcopy(atomB.type)
                atomB.type = deepcopy(atomA.type)
                atomB.chargeB = deepcopy(atomB.charge)
                atomB.charge = deepcopy(atomA.charge)
                atomB.massB = deepcopy(atomB.mass)
                atomB.mass = deepcopy(atomA.mass)
                if atomA.charge == atomB.charge:
                    logger.debug(
                        f"Atom {nr} changed but not the charges!\n\tA:{atomA}\n\tB:{atomB} "
                    )

    def merge_bonds(self):
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
                bond_key, bond_b, self.ff.bondtypes, self.mol_b, use_state_b=True
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
            atomtypes_b: tuple[str, str] = tuple(
                [
                    (
                        t
                        if (t := getattr(self.mol_b.atoms[id], "typeB", None))
                        else self.mol_b.atoms[id].type
                    )
                    for id in bond_key
                ]
            )  # type: ignore

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
                self.mol_b.bonds[bond_key] = Bond(
                    *bond_key,
                    funct=bond_b.funct,
                    c0=bond_a.c0,
                    c1=bond_a.c1,
                    c2=bond_b.c0,
                    c3=bond_b.c1,
                )
                logger.debug(
                    f"Created harmonic bond to bond transition for {bond_key} (transition): {self.mol_b.bonds[bond_key]}"
                )

                # make sure the actual bond is excluded just in case
                # None is a dummy value for the pair
                # the just generates no pair but the exclusion anyways
                # in self.merge_pairs later
                self.affected_pairs[bond_key] = None

            elif isinstance(bond_a, (Bond, BondType)):
                # bond only exists in A -> vanish bond
                self.mol_b.bonds[bond_key] = Bond(
                    *bond_key,
                    funct=FFFUNC["morse_bond"],
                    c0=bond_a.c0,  # bond distance
                    c1=f"{depth_a:.5f}",  # D
                    c2=f"{beta_a:.5f}",  # beta
                    c3=f"{bond_a.c0}",  # new r0 (=stays the same, because bond vanishes anyways)
                    c4=f"{0.0:.5f}",  # well depth goes to 0
                    c5=f"{self.default_beta_for_lj:.5f}",  # default beta for LJ
                )
                logger.debug(
                    f"Created Morse to LJ transition bond for {bond_key} (vanish): {self.mol_b.bonds[bond_key]}"
                )
                # the bond is vanishing, so the pair interaction is turned on (at full strength)
                # this can be overwritten if e.g. a pair for the same atoms is created later
                # due to affected angles or dihedrals
                # this should win over the bonds because if the atoms are still affected
                # by dihedrals after the bond is gone, the pair should only be turned on
                # to half strength
                logger.debug(
                    f"Adding pair for bond {bond_key} breaking to affected_pairs"
                )
                self.affected_pairs[bond_key] = self._make_pair(
                    *bond_key, PairTransition.Create
                )

            elif isinstance(bond_b, (Bond, BondType)):
                # bond only exists in B -> create bond
                if epsilonij_a == 0:
                    m = f"epsilonij_a is 0 for {bond_key}"
                    logger.error(m)
                    raise ValueError(m)
                self.mol_b.bonds[bond_key] = Bond(
                    *bond_key,
                    funct=FFFUNC["morse_bond"],
                    # starts further away, so it can pull in atoms of the bond
                    c0=f"{sigmaij_a*self.stateA_morse_dist_factor:.5f}",  # sigmaij* 1.12 = LJ minimum at lambda 0, where the atom types are still from A
                    c1=f"{0.0:.5f}",  # well depth from 0
                    c2=f"{0.0:.5f}",  # beta from 0 = flat
                    c3=bond_b.c0,  # final bond distance
                    c4=f"{depth_b:.5f}",  # final well depth
                    c5=f"{beta_b:.5f}",  # final steepness (beta)
                )
                logger.debug(
                    f"Created LJ to Morse transition bond for {bond_key} (create): {self.mol_b.bonds[bond_key]}"
                )
                # the bond is created, so the pair interaction is turned off
                logger.debug(
                    f"Adding pair for bond {bond_key} being created to affected_pairs"
                )
                self.affected_pairs[bond_key] = self._make_pair(
                    *bond_key, PairTransition.Vanish
                )

    def merge_angles(self):
        keys = set(self.mol_a.angles.keys()) | set(self.mol_b.angles.keys())
        for key in keys:
            angle_a = self.mol_a.angles.get(key)
            angle_b = self.mol_b.angles.get(key)

            if angle_a == angle_b:
                # angles are the same, nothing to be done.
                # atoms will also continue to be automatically excluded
                # so nothing to be done with pairs
                continue

            parameterizedA = get_explicit_or_type(
                key, angle_a, self.ff.angletypes, self.mol_a
            )
            parameterizedB = get_explicit_or_type(
                key, angle_b, self.ff.angletypes, self.mol_b
            )

            if parameterizedA is not None and parameterizedB is not None:
                if not (
                    isinstance(parameterizedA, (Angle, AngleType))
                    and isinstance(parameterizedB, (Angle, AngleType))
                ):
                    m = f"Can't find parameters for bond with key {key}, got explicits/types: {parameterizedA} and {parameterizedB} for a bonds: {angle_a}, {angle_b}"
                    logger.error(m)
                    raise ValueError(m)
                # both angles exist -> transition between parameters
                self.mol_b.angles[key] = Angle(
                    *key,
                    funct=parameterizedB.funct,
                    c0=parameterizedA.c0,
                    c1=parameterizedA.c1,
                    c2=parameterizedB.c0,
                    c3=parameterizedB.c1,
                )
                # nothing to be done with pairs
                # atoms will also continue to be automatically excluded
            elif isinstance(parameterizedA, (Angle, AngleType)):
                # angle only exists in A -> vanish angle
                self.mol_b.angles[key] = Angle(
                    *key,
                    funct=FFFUNC["harmonic_angle"],
                    c0=parameterizedA.c0,
                    c1=parameterizedA.c1,
                    c2=parameterizedA.c0,
                    c3="0.00",
                )
                # the angle (1 - 2 - 3) is vanished, so the pair interaction
                # between the outer atoms (1 -- 3) is turned on
                # (1 -- 2 and 2 -- 3 are already handled during the bond handling)
                logger.debug(f"Adding pair for angle {key} breaking to affected_pairs")
                self.affected_pairs[(key[0], key[2])] = self._make_pair(
                    id1=key[0], id2=key[2], transition=PairTransition.Create
                )
            elif isinstance(parameterizedB, (Angle, AngleType)):
                # angle only exists in B -> create angle
                self.mol_b.angles[key] = Angle(
                    *key,
                    funct=FFFUNC["harmonic_angle"],
                    c0=parameterizedB.c0,
                    c1="0.00",
                    c2=parameterizedB.c0,
                    c3=parameterizedB.c1,
                )
                # the angle is created, so the non-bonded interaction
                # between the outer atoms is turned off via a pair
                logger.debug(
                    f"Adding pair for angle {key} being created to affected_pairs"
                )
                self.affected_pairs[(key[0], key[2])] = self._make_pair(
                    id1=key[0], id2=key[2], transition=PairTransition.Vanish
                )
            else:
                logger.warning(f"Could not parameterize angle {key}.")

    def merge_dihedrals(self):
        # proper dihedrals
        # proper dihedrals have a nested structure and need a different treatment from bonds, angles and improper dihedrals
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
                logger.debug(
                    f"Adding pair for dihedral {key} being morphed to affected_pairs"
                )
                self.affected_pairs[(key[0], key[3])] = self._make_pair(
                    id1=key[0],
                    id2=key[3],
                    transition=PairTransition.Morph,
                    from_1_4=True,
                    to_1_4=True,
                )
            elif isinstance(multiple_dihedralsA, MultipleDihedrals):
                # dihedral only exists in A -> vanishing dihedral
                # pair transitions from half strength 1-4 interaction to full strength
                logger.debug(
                    f"Adding pair for dihedral {key} breaking to affected_pairs"
                )
                self.affected_pairs[(key[0], key[3])] = self._make_pair(
                    id1=key[0],
                    id2=key[3],
                    transition=PairTransition.Morph,
                    from_1_4=True,
                )
            elif isinstance(multiple_dihedralsB, MultipleDihedrals):
                # dihedral only exists in B -> creating dihedral
                # pair transitions from full strength to half strength 1-4 interaction
                logger.debug(
                    f"Adding pair for dihedral {key} being created to affected_pairs"
                )
                self.affected_pairs[(key[0], key[3])] = self._make_pair(
                    id1=key[0], id2=key[3], transition=PairTransition.Morph, to_1_4=True
                )

            if key == ("15", "17", "19", "24"):
                print(key, multiple_dihedralsA, multiple_dihedralsB)
            if multiple_dihedralsA == multiple_dihedralsB:
                # dihedrals are the same, nothing further to be done
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

            self.mol_b.proper_dihedrals[key] = MultipleDihedrals(
                *key, funct=FFFUNC["mult_proper_dihedral"], dihedrals={}
            )
            periodicities = keysA | keysB
            for periodicity in periodicities:
                interactionA = (
                    multiple_dihedralsA.dihedrals.get(periodicity)
                    if multiple_dihedralsA
                    else None
                )
                interactionB = (
                    multiple_dihedralsB.dihedrals.get(periodicity)
                    if multiple_dihedralsB
                    else None
                )

                self.mol_b.proper_dihedrals[key].dihedrals[periodicity] = (
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

            self.mol_b.improper_dihedrals[key] = MultipleDihedrals(
                *key, funct=FFFUNC["mult_improper_dihedral"], dihedrals={}
            )
            periodicities = keysA | keysB
            for periodicity in periodicities:
                assert isinstance(periodicity, str)
                interactionA = (
                    multiple_dihedralsA.dihedrals.get(periodicity)
                    if multiple_dihedralsA
                    else None
                )
                interactionB = (
                    multiple_dihedralsB.dihedrals.get(periodicity)
                    if multiple_dihedralsB
                    else None
                )

                self.mol_b.improper_dihedrals[key].dihedrals[periodicity] = (
                    self._interpolate_two_dihedrals(
                        key=key,
                        dihedral_a=interactionA,
                        dihedral_b=interactionB,
                        funct=FFFUNC["mult_improper_dihedral"],
                    )
                )

    def merge_pairs(self):
        """
        If it is for a binding event, the pair is vanishing, as it will be an
        exclusion once bound. If a bond is breaking, the pair interaction
        is slowly turned on, as it was excluded previously.

        growing/shrinking pairs are used to turn on/off
        LJ / non-bonded interactions
        """
        if not self.enable_morph_pairs:
            return

        keys = set(self.mol_a.pairs.keys()) | set(self.mol_b.pairs.keys())
        for key in keys:
            pair_a = self.mol_a.pairs.get(key)
            pair_b = self.mol_b.pairs.get(key)

            if pair_a is None and pair_b is None:
                m = f"Can't find parameters for bond with key {key}, got None for both bonds. This should be an impossible state."
                logger.error(m)
                raise ValueError(m)

            # transitioning between actual pairs
            if pair_a is not None and pair_b is not None:
                logger.debug(
                    f"Pair {key} is in both topologies, transitioning between them"
                )
                # TODO: handle explicitly parameterized pairs
                self.affected_pairs[key] = self._make_pair(
                    key[0], key[1], PairTransition.Morph, from_1_4=True, to_1_4=True
                )
            elif pair_a is not None:
                if self.affected_pairs.get(key) is None:
                    logger.debug(
                        f"Pair {key} is only in A, going from 1-4 to full strength LJ"
                    )
                    self.affected_pairs[key] = self._make_pair(
                        key[0], key[1], PairTransition.Morph, from_1_4=True
                    )
                else:
                    # pair is already in there based on dihedrals or angles
                    pass
            elif pair_b is not None:
                if self.affected_pairs.get(key) is None:
                    logger.debug(
                        f"Pair {key} is only in B, going from full strength LJ to 1-4"
                    )
                    self.affected_pairs[key] = self._make_pair(
                        key[0], key[1], PairTransition.Morph, to_1_4=True
                    )
                else:
                    # pair is already in there based on dihedrals or angles
                    pass

            # Add general exclusions for each pair
            # such that automatic non-bonded interactions
            # don't interfer with our efforts
            self.mol_b.exclusions[key] = Exclusion(*key)

        # pairs that stem from bonds breaking or forming
        # collected by looking at the changing bonds, angles (=like second order bonds)
        # and dihedrals
        for key, pair in self.affected_pairs.items():
            if pair is not None:
                self.mol_b.pairs[key] = pair

            # gromacs does not automatically create exclusions for our our new morse bonds (type 3).
            # see https://manual.gromacs.org/current/reference-manual/topologies/molecule-definition.html#excl
            # "Particles are considered bonded when they are connected by “chemical” bonds ([ bonds ] types 1 to 5, 7 or 8)"
            self.mol_b.exclusions[key] = Exclusion(*key)


def merge_top_slow_growth(
    top_a: Topology, top_b: Topology, enable_morph_pairs: bool
) -> Topology:
    """Takes two Topologies and joins them for a smooth free-energy like parameter transition simulation.
    Modifies topB in place.
    All changes are contained to the `Reactive` moleculetype.
    """

    MoleculeTypeMerger(
        mol_a=top_a.moleculetypes[REACTIVE_MOLECULEYPE],
        mol_b=top_b.moleculetypes[REACTIVE_MOLECULEYPE],
        ff=top_b.ff,
        enable_morph_pairs=enable_morph_pairs,
    ).merge()
    return top_b


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
