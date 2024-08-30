"""coordinate, topology and plumed modification functions"""

import logging
from copy import deepcopy
from pathlib import Path
from typing import Optional, Union

import MDAnalysis as mda
import numpy as np

from kimmdy.constants import REACTIVE_MOLECULEYPE, FFFUNC
from kimmdy.parsing import read_plumed, write_plumed
from kimmdy.recipe import Place
from kimmdy.tasks import TaskFiles
from kimmdy.topology.atomic import (
    Angle,
    Bond,
    Dihedral,
    DihedralType,
    InteractionType,
    Pair,
    Exclusion,
    ImproperDihedralId,
    Interaction,
    AtomicType,
    InteractionTypes,
    MultipleDihedrals,
    ProperDihedralId,
)
from kimmdy.topology.ff import FF
from kimmdy.topology.topology import MoleculeType, Topology
from kimmdy.topology.utils import match_atomic_item_to_atomic_type

logger = logging.getLogger(__name__)


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


def get_explicit_MultipleDihedrals(
    dihedral_key: tuple[str, str, str, str],
    mol: MoleculeType,
    dihedrals_in: Optional[MultipleDihedrals],
    ff: FF,
    periodicity_max: int = 6,
) -> Optional[MultipleDihedrals]:
    """Takes a valid dihedral key and returns explicit
    dihedral parameters for a given topology
    """
    if not dihedrals_in:
        return None

    if "" not in dihedrals_in.dihedrals.keys():
        # empty string means implicit parameters
        return dihedrals_in

    type_key = [mol.atoms[id].type for id in dihedral_key]

    multiple_dihedrals = MultipleDihedrals(
        *dihedral_key, FFFUNC["mult_proper_dihedral"], {}
    )
    for periodicity in range(1, periodicity_max + 1):
        match_obj = match_atomic_item_to_atomic_type(
            type_key, ff.proper_dihedraltypes, str(periodicity)
        )
        if match_obj:
            assert isinstance(match_obj, DihedralType)
            multiple_dihedrals.dihedrals[str(periodicity)] = Dihedral(
                *dihedral_key,
                funct=FFFUNC["mult_proper_dihedral"],
                c0=match_obj.c0,
                c1=match_obj.c1,
                periodicity=match_obj.periodicity,
            )

    if not multiple_dihedrals.dihedrals:
        return None

    return multiple_dihedrals


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
    match_obj = match_atomic_item_to_atomic_type(
        type_key, interaction_types, periodicity
    )

    if match_obj:
        assert isinstance(match_obj, InteractionType)
        # FIXME: This is a property of `match_atomic_item_to_atomic_type` and should
        # be tested there.
        return match_obj
    else:
        raise ValueError(
            f"Could not find explicit parameters for {[type_key,interaction]} in line or in {interaction_types}."
        )


def merge_dihedrals(
    dihedral_key: tuple[str, str, str, str],
    dihedral_a: Optional[Dihedral],
    dihedral_b: Optional[Dihedral],
    dihedral_types_a: Union[
        dict[ProperDihedralId, DihedralType], dict[ImproperDihedralId, DihedralType]
    ],
    dihedral_types_b: Union[
        dict[ProperDihedralId, DihedralType], dict[ImproperDihedralId, DihedralType]
    ],
    molA: MoleculeType,
    molB: MoleculeType,
    funct: str,
    periodicity: str,
) -> Dihedral:
    """Merge one to two Dihedrals or -Types into a Dihedral in free-energy syntax"""
    # convert implicit standard ff parameters to explicit, if necessary
    if dihedral_a:
        parameterizedA = get_explicit_or_type(
            dihedral_key,
            dihedral_a,
            dihedral_types_a,
            molA,
            periodicity,
        )
    else:
        parameterizedA = None

    if dihedral_b:
        parameterizedB = get_explicit_or_type(
            dihedral_key,
            dihedral_b,
            dihedral_types_b,
            molB,
            periodicity,
        )
    else:
        parameterizedB = None

    # construct parameterized Dihedral
    if parameterizedA is not None and parameterizedB is not None:
        # same

        assert type(parameterizedA) == Dihedral or type(parameterizedA) == DihedralType
        assert type(parameterizedB) == Dihedral or type(parameterizedB) == DihedralType

        dihedralmerge = Dihedral(
            *dihedral_key,
            funct=funct,
            c0=parameterizedA.c0,
            c1=parameterizedA.c1,
            periodicity=parameterizedA.periodicity,
            c3=parameterizedB.c0,
            c4=parameterizedB.c1,
            c5=parameterizedB.periodicity,
        )
    elif parameterizedA is not None:
        # breaking
        if not (
            type(parameterizedA) == Dihedral or type(parameterizedA) == DihedralType
        ):
            m = f"parameterizedA {parameterizedA} is not a Dihedral or DihedralType"
            logger.warning(m)
            raise ValueError(m)
        dihedralmerge = Dihedral(
            *dihedral_key,
            funct=funct,
            c0=parameterizedA.c0,
            c1=parameterizedA.c1,
            periodicity=parameterizedA.periodicity,
            c3=parameterizedA.c0,
            c4="0.00",
            c5=parameterizedA.periodicity,
        )
    elif parameterizedB is not None:
        # binding
        assert type(parameterizedB) == Dihedral or type(parameterizedB) == DihedralType
        dihedralmerge = Dihedral(
            *dihedral_key,
            funct=funct,
            c0=parameterizedB.c0,
            c1="0.00",
            periodicity=parameterizedB.periodicity,
            c3=parameterizedB.c0,
            c4=parameterizedB.c1,
            c5=parameterizedB.periodicity,
        )
    else:
        m = f"Tried to merge two dihedrals of {dihedral_key} but no parameterized dihedrals found!"
        logger.error(m)
        raise ValueError(m)
    return dihedralmerge


def merge_top_moleculetypes_slow_growth(
    molA: MoleculeType,
    molB: MoleculeType,
    ff: FF,
    morph_pairs: bool,
) -> MoleculeType:
    """Takes two Topologies and joins them for a smooth free-energy like parameter transition simulation"""

    def get_LJ_parameters(idx1: str, idx2: str):
        """Calculate LJ terms sigma and epsilon from atom types"""

        comb_rule = ff.defaults[0][1]
        # fudgeLJ = ff.defaults[0][3] # could be used to compensate fudge in pairs

        type1 = ff.atomtypes[molB.atoms[idx1].type]
        type2 = ff.atomtypes[molB.atoms[idx2].type]

        if comb_rule == "1":
            v = np.sqrt(float(type1.sigma) * float(type2.sigma))
            w = np.sqrt(float(type1.epsilon) * float(type2.epsilon))
        elif comb_rule == "2":
            v = 0.5 * (float(type1.sigma) + float(type2.sigma))
            w = np.sqrt(float(type1.epsilon) * float(type2.epsilon))
        elif comb_rule == "3":
            v = np.sqrt(float(type1.sigma) * float(type2.sigma))
            w = np.sqrt(float(type1.epsilon) * float(type2.epsilon))
        else:
            raise ValueError("Unknown combination rule of forcefield")

        return v, w

    def make_pair(idx1: str, idx2: str, bind: bool) -> Pair:
        """Generates morphing pair interaction

        If it is for a binding event, the pair is vanishing, as it will be an
        exclusion once bound. If a bond is breaking, the pair interaction
        is slowly turned on, as it was excluded previously.

        Parameters
        ----------
        idx1 : str
            Atom one
        idx2 : str
            Atom two
        bind : bool
            Binding or breaking event. Determines morphing direction.

        Returns
        -------
        Pair
            Morphing pair
        """
        v, w = get_LJ_parameters(idx1, idx2)
        c_kwargs = dict(
            zip(
                [f"c{i}" for i in range(4)],
                [f"{0.0:.5f}" for _ in range(4)],
            )
        )
        if morph_pairs:
            # Bind: pair interaction turning off
            if bind:
                c_kwargs["c0"] = f"{v:.5f}"
                c_kwargs["c1"] = f"{w:.5f}"
            # Break: pair interaction turning on
            else:
                c_kwargs["c2"] = f"{v:.5f}"
                c_kwargs["c3"] = f"{w:.5f}"

        return Pair(idx1, idx2, funct=ff.defaults[0][0], **c_kwargs)

    def merge_pairs(break_pair: Optional[Pair], bind_pair: Optional[Pair]):
        """Merges two morphing pairs into a morphing pair for slow-growth.
        Takes starting state of break_pair and end state of bind_pair.

        Parameters
        ----------
        break_pair : Optional[Pair]
            Pair for breaking event
        bind_pair : Optional[Pair]
            Pair for binding event

        Returns
        -------
        Pair
            Merged pair containing four parameters.
        """

        assert len(ps := list(filter(lambda x: x, (break_pair, bind_pair)))) > 0
        ai = ps[0].ai  # type: ignore
        aj = ps[0].aj  # type: ignore
        funct = ps[0].funct  # type: ignore

        for p in ps:
            assert p.c0 is not None, "Pair must contain c0"  # type: ignore
            assert p.c1 is not None, "Pair must contain c1"  # type: ignore
            assert p.c2 is not None, "Pair must contain c2"  # type: ignore
            assert p.c3 is not None, "Pair must contain c3"  # type: ignore

        if break_pair and bind_pair:
            assert (
                break_pair.funct == bind_pair.funct
            ), "Pair functional must be the same to interpolate"
            assert (break_pair.ai == bind_pair.ai) or (
                break_pair.ai == bind_pair.aj
            ), "Atoms must be the same in pairs to interpolate between"
            assert (break_pair.aj == bind_pair.ai) or (
                break_pair.aj == bind_pair.aj
            ), "Atoms must be the same in pairs to interpolate between"

        return Pair(
            ai,
            aj,
            funct=funct,
            c0=break_pair.c0 if break_pair else bind_pair.c0,  # type: ignore
            c1=break_pair.c1 if break_pair else bind_pair.c1,  # type: ignore
            c2=bind_pair.c2 if bind_pair else break_pair.c2,  # type: ignore
            c3=bind_pair.c3 if bind_pair else break_pair.c3,  # type: ignore
        )

    hyperparameters = {
        "morse_well_depth": 300,  # [kJ mol-1]
        "beta": 19,  # matches LJ steepness
    }

    # atoms
    for nr in molA.atoms.keys():
        atomA = molA.atoms[nr]
        atomB = molB.atoms[nr]
        if atomA != atomB:
            if atomA.charge != atomB.charge:
                atomB.typeB = deepcopy(atomB.type)
                atomB.type = deepcopy(atomA.type)
                atomB.chargeB = deepcopy(atomB.charge)
                atomB.charge = deepcopy(atomA.charge)
                atomB.massB = deepcopy(atomB.mass)
                atomB.mass = deepcopy(atomA.mass)
            else:
                logger.debug(
                    f"Atom {nr} changed but not the charges!\n\tA:{atomA}\n\tB:{atomB} "
                )

    break_bind_atoms = {}

    # bonds
    keysA = set(molA.bonds.keys())
    keysB = set(molB.bonds.keys())
    keys = keysA | keysB

    new_pairs_break = {}
    new_pairs_bind = {}

    for key in keys:
        interactionA = molA.bonds.get(key)
        interactionB = molB.bonds.get(key)

        if interactionA != interactionB:
            parameterizedA = get_explicit_or_type(key, interactionA, ff.bondtypes, molA)
            parameterizedB = get_explicit_or_type(key, interactionB, ff.bondtypes, molB)

            # both bonds exist
            if parameterizedA and parameterizedB:
                if not (
                    parameterizedA.funct
                    == parameterizedB.funct
                    == FFFUNC["harmonic_bond"]
                ):
                    m = f"In slow-growth, bond functionals need to be harmonic, but {key} is not. It is in A: {parameterizedA} and B: {parameterizedB}"
                    logger.error(m)
                    raise ValueError(m)
                molB.bonds[key] = Bond(
                    *key,
                    funct=parameterizedB.funct,
                    c0=parameterizedA.c0,
                    c1=parameterizedA.c1,
                    c2=parameterizedB.c0,
                    c3=parameterizedB.c1,
                )

            # bond only exists in A -> vanish bond
            elif parameterizedA:
                atomtypes = [molA.atoms[atom_id].type for atom_id in key]
                break_bind_atoms[key[0]] = atomtypes[0]
                break_bind_atoms[key[1]] = atomtypes[1]
                sigmaij, epsilonij = get_LJ_parameters(*key)

                molB.bonds[key] = Bond(
                    *key,
                    funct=FFFUNC["morse_bond"],
                    c0=parameterizedA.c0,  # b
                    c1=f"{hyperparameters['morse_well_depth']:.5f}",  # D
                    c2=f"{hyperparameters['beta']:.5f}",  # beta
                    c3=f"{sigmaij*1.12:.5f}",  # sigmaij* 1.12 = LJ minimum
                    c4=f"{0.0:.5f}",  # well depth -> zero
                    c5=f"{0.0:.5f}",
                )

                # update bound_to
                atompair = [molB.atoms[key[0]], molB.atoms[key[1]]]
                atompair[0].bound_to_nrs.append(atompair[1].nr)
                atompair[1].bound_to_nrs.append(atompair[0].nr)

                # Find all neighbors of the bond
                neighbor_sides = []
                for idx in key:
                    bonds = list(filter(lambda b: idx in b, keysB))
                    neighbor_set = set()
                    for b in bonds:
                        neighbor_set.add(b[0])
                        neighbor_set.add(b[1])
                    neighbor_sides.append(neighbor_set)
                neighbor_set.discard(idx)  # avoid double counting central bond

                # handle vanishing exclusions neighbors1 - key[1]
                for a1 in neighbor_sides[0]:
                    pk = tuple(sorted((a1, key[1]), key=int))
                    new_pairs_break[pk] = make_pair(*pk, bind=False)

                # handle vanishing exclusions neighbors2 - key[0]
                for a2 in neighbor_sides[1]:
                    pk = tuple(sorted((a2, key[0]), key=int))
                    new_pairs_break[pk] = make_pair(*pk, bind=False)

            # bond only exists in B -> create bond
            elif parameterizedB:
                atomtypes = [molB.atoms[atom_id].type for atom_id in key]
                break_bind_atoms[key[0]] = atomtypes[0]
                break_bind_atoms[key[1]] = atomtypes[1]
                sigmaij, epsilonij = get_LJ_parameters(*key)

                molB.bonds[key] = Bond(
                    *key,
                    funct=FFFUNC["morse_bond"],
                    c0=f"{sigmaij*1.12:.5f}",  # sigmaij* 1.12 = LJ minimum
                    c1=f"{0.0:.5f}",  # well depth is epsilonij
                    c2=f"{0.0:.5f}",
                    c3=parameterizedB.c0,  # b
                    c4=f"{hyperparameters['morse_well_depth']:.5f}",  # D
                    c5=f"{hyperparameters['beta']:.5f}",  # beta
                )

                # Find all neighbors of the bond
                neighbor_sides = []
                for idx in key:
                    bonds = list(filter(lambda b: idx in b, keysA))
                    neighbor_set = set()
                    for b in bonds:
                        neighbor_set.add(b[0])
                        neighbor_set.add(b[1])
                    neighbor_sides.append(neighbor_set)
                neighbor_set.discard(idx)  # type: ignore  # avoid double counting central bond

                # handle growing exclusions neighbors1 - key[1]
                for a1 in neighbor_sides[0]:
                    pk = tuple(sorted((a1, key[1]), key=int))
                    new_pairs_bind[pk] = make_pair(*pk, bind=True)

                # handle growing exclusions neighbors2 - key[0]
                for a2 in neighbor_sides[1]:
                    pk = tuple(sorted((a2, key[0]), key=int))
                    new_pairs_bind[pk] = make_pair(*pk, bind=True)

    for key in set(new_pairs_bind.keys()) | set(new_pairs_break.keys()):
        break_pair = new_pairs_break.get(key)
        bind_pair = new_pairs_bind.get(key)

        morph_pair = merge_pairs(break_pair, bind_pair)

        molB.pairs[key] = morph_pair
        molB.exclusions[key] = Exclusion(*key)

    # angles
    keys = set(molA.angles.keys()) | set(molB.angles.keys())
    for key in keys:
        interactionA = molA.angles.get(key)
        interactionB = molB.angles.get(key)

        if interactionA != interactionB:
            parameterizedA = get_explicit_or_type(
                key, interactionA, ff.angletypes, molA
            )
            parameterizedB = get_explicit_or_type(
                key, interactionB, ff.angletypes, molB
            )
            if parameterizedA and parameterizedB:
                molB.angles[key] = Angle(
                    *key,
                    funct=parameterizedB.funct,
                    c0=parameterizedA.c0,
                    c1=parameterizedA.c1,
                    c2=parameterizedB.c0,
                    c3=parameterizedB.c1,
                )
            elif parameterizedA:
                molB.angles[key] = Angle(
                    *key,
                    funct=FFFUNC["harmonic_angle"],
                    c0=parameterizedA.c0,
                    c1=parameterizedA.c1,
                    c2=parameterizedA.c0,
                    c3="0.00",
                )
            elif parameterizedB:
                molB.angles[key] = Angle(
                    *key,
                    funct=FFFUNC["harmonic_angle"],
                    c0=parameterizedB.c0,
                    c1="0.00",
                    c2=parameterizedB.c0,
                    c3=parameterizedB.c1,
                )
            else:
                logger.warning(f"Could not parameterize angle {key}.")

    # dihedrals
    # proper dihedrals
    # proper dihedrals have a nested structure and need a different treatment from bonds, angles and improper dihedrals
    # if indices change atomtypes and parameters change because of that, it will ignore these parameter change

    keys = set(molA.proper_dihedrals.keys()) | set(molB.proper_dihedrals.keys())
    for key in keys:
        multiple_dihedralsA = molA.proper_dihedrals.get(key)
        multiple_dihedralsB = molB.proper_dihedrals.get(key)

        if multiple_dihedralsA != multiple_dihedralsB:
            multiple_dihedralsA = get_explicit_MultipleDihedrals(
                key, molA, multiple_dihedralsA, ff
            )
            multiple_dihedralsB = get_explicit_MultipleDihedrals(
                key, molB, multiple_dihedralsB, ff
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

            molB.proper_dihedrals[key] = MultipleDihedrals(
                *key, FFFUNC["mult_proper_dihedral"], {}
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

                molB.proper_dihedrals[key].dihedrals[periodicity] = merge_dihedrals(
                    key,
                    interactionA,
                    interactionB,
                    ff.proper_dihedraltypes,
                    ff.proper_dihedraltypes,
                    molA,
                    molB,
                    FFFUNC["mult_proper_dihedral"],
                    periodicity,
                )

    # improper dihedrals
    keys = set(molA.improper_dihedrals.keys()) | set(molB.improper_dihedrals.keys())
    for key in keys:
        multiple_dihedralsA = molA.improper_dihedrals.get(key)
        multiple_dihedralsB = molB.improper_dihedrals.get(key)

        if multiple_dihedralsA != multiple_dihedralsB:
            multiple_dihedralsA = get_explicit_MultipleDihedrals(
                key, molA, multiple_dihedralsA, ff
            )
            multiple_dihedralsB = get_explicit_MultipleDihedrals(
                key, molB, multiple_dihedralsB, ff
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

            molB.improper_dihedrals[key] = MultipleDihedrals(
                *key, FFFUNC["mult_improper_dihedral"], {}
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

                molB.improper_dihedrals[key].dihedrals[periodicity] = merge_dihedrals(
                    key,
                    interactionA,
                    interactionB,
                    ff.improper_dihedraltypes,
                    ff.improper_dihedraltypes,
                    molA,
                    molB,
                    FFFUNC["mult_improper_dihedral"],
                    periodicity,
                )

    # amber fix for breaking/binding atom types without LJ potential
    for k, v in break_bind_atoms.items():
        if v in ["HW", "HO"]:
            molB.atoms[k].type = "H1"
            if molB.atoms[k].typeB is not None:
                molB.atoms[k].typeB = "H1"

    # update is_radical attribute of Atom objects in topology
    molB.find_radicals()

    return molB


def merge_top_slow_growth(
    topA: Topology, topB: Topology, morph_pairs: bool
) -> Topology:
    """Takes two Topologies and joins them for a smooth free-energy like parameter transition simulation.
    For now this assumes that only one moleculeype is of interest.
    """

    molA = topA.moleculetypes[REACTIVE_MOLECULEYPE]
    molB = topB.moleculetypes[REACTIVE_MOLECULEYPE]
    molB = merge_top_moleculetypes_slow_growth(
        molA=molA, molB=molB, ff=topB.ff, morph_pairs=morph_pairs
    )
    # not necessary, will be updated automatically on `to_dict()`
    # topB._update_dict()

    return topB


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
