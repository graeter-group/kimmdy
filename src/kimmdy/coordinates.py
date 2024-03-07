"""coordinate, topology and plumed modification functions"""

import logging
from copy import deepcopy
from pathlib import Path
from typing import Optional, Union

import MDAnalysis as mda
import numpy as np

from kimmdy.constants import REACTIVE_MOLECULEYPE
from kimmdy.parsing import read_plumed, write_plumed
from kimmdy.recipe import Place
from kimmdy.tasks import TaskFiles
from kimmdy.topology.atomic import (
    Angle,
    Bond,
    Dihedral,
    DihedralType,
    Exclusion,
    ImproperDihedralId,
    Interaction,
    InteractionType,
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

    u.atoms.write(trr_out)
    u.atoms.write(gro_out)

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

    multiple_dihedrals = MultipleDihedrals(*dihedral_key, "9", {})
    for periodicity in range(1, periodicity_max + 1):
        match_obj = match_atomic_item_to_atomic_type(
            type_key, ff.proper_dihedraltypes, str(periodicity)
        )
        if match_obj:
            assert isinstance(match_obj, DihedralType)
            multiple_dihedrals.dihedrals[str(periodicity)] = Dihedral(
                *dihedral_key,
                funct="9",
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
) -> Union[Interaction, InteractionType, None]:
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
        assert type(parameterizedA) == Dihedral or type(parameterizedA) == DihedralType
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
            funct="9",
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
) -> MoleculeType:
    """Takes two Topologies and joins them for a smooth free-energy like parameter transition simulation"""
    hyperparameters = {"morse_well_depth": 300}  # [kJ mol-1]

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

    for key in keys:
        interactionA = molA.bonds.get(key)
        interactionB = molB.bonds.get(key)

        if interactionA != interactionB:
            parameterizedA = get_explicit_or_type(key, interactionA, ff.bondtypes, molA)
            parameterizedB = get_explicit_or_type(key, interactionB, ff.bondtypes, molB)
            if parameterizedA and parameterizedB:
                molB.bonds[key] = Bond(
                    *key,
                    funct=parameterizedB.funct,
                    c0=parameterizedA.c0,
                    c1=parameterizedA.c1,
                    c2=parameterizedB.c0,
                    c3=parameterizedB.c1,
                )
            elif parameterizedA:
                atomtypes = [molA.atoms[atom_id].type for atom_id in key]
                break_bind_atoms[key[0]] = atomtypes[0]
                break_bind_atoms[key[1]] = atomtypes[1]
                # use combination rule type 2 (typically used by amber force fields)
                sigmas = [ff.atomtypes[at].sigma for at in atomtypes]
                sigmaij = 0.5 * (float(sigmas[0]) + float(sigmas[1]))
                epsilons = [ff.atomtypes[at].epsilon for at in atomtypes]
                epsilonij = np.sqrt(float(epsilons[0]) * float(epsilons[1]))

                # morse well steepness
                assert parameterizedA.c1 is not None, f"parameterizedA.c1 is not set"
                beta = np.sqrt(
                    float(parameterizedA.c1) / (2 * hyperparameters["morse_well_depth"])
                )

                molB.bonds[key] = Bond(
                    *key,
                    funct="3",
                    c0=parameterizedA.c0,
                    c1=f"{hyperparameters['morse_well_depth']:5.3f}",
                    c2=f"{beta:5.3f}",
                    c3=f"{sigmaij*1.12:7.5f}",  # sigmaij* 1.12 = LJ minimum
                    c4=f"{epsilonij:7.5f}",  # well depth is epsilonij
                    c5=f"{beta:5.3f}",
                )

                # update bound_to
                atompair = [molB.atoms[key[0]], molB.atoms[key[1]]]
                atompair[0].bound_to_nrs.append(atompair[1].nr)
                atompair[1].bound_to_nrs.append(atompair[0].nr)

            elif parameterizedB:
                atomtypes = [molB.atoms[atom_id].type for atom_id in key]
                break_bind_atoms[key[0]] = atomtypes[0]
                break_bind_atoms[key[1]] = atomtypes[1]
                # use combination rule type 2 (typically used by amber force fields)
                sigmas = [ff.atomtypes[at].sigma for at in atomtypes]
                sigmaij = 0.5 * (float(sigmas[0]) + float(sigmas[1]))
                epsilons = [ff.atomtypes[at].epsilon for at in atomtypes]
                epsilonij = np.sqrt(float(epsilons[0]) * float(epsilons[1]))
                # morse well steepness
                assert parameterizedB.c1 is not None, f"parameterizedB.c1 is not set"
                beta = np.sqrt(
                    float(parameterizedB.c1) / (2 * hyperparameters["morse_well_depth"])
                )
                molB.bonds[key] = Bond(
                    *key,
                    funct="3",
                    c0=f"{sigmaij*1.12:7.5f}",  # sigmaij* 1.12 = LJ minimum
                    c1=f"{epsilonij:7.5f}",  # well depth is epsilonij
                    c2=f"{beta:5.3f}",
                    c3=parameterizedB.c0,
                    c4=f"{hyperparameters['morse_well_depth']:5.3f}",
                    c5=f"{beta:5.3f}",
                )

    # pairs and exclusions
    exclusions = molB.exclusions
    for key in keysA - keysB:
        molB.pairs.pop(key, None)
        exclusion = Exclusion(*key)
        exclusions[key] = exclusion

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
                    funct="1",
                    c0=parameterizedA.c0,
                    c1=parameterizedA.c1,
                    c2=parameterizedA.c0,
                    c3="0.00",
                )
            elif parameterizedB:
                molB.angles[key] = Angle(
                    *key,
                    funct="1",
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

            molB.proper_dihedrals[key] = MultipleDihedrals(*key, "9", {})
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
                    "9",
                    periodicity,
                )

    # improper dihedrals
    # TODO: duplicate of proper dihedrals, could refactor

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

            molB.improper_dihedrals[key] = MultipleDihedrals(*key, "9", {})
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
                    "9",
                    periodicity,
                )

    # amber fix for breaking/binding atom types without LJ potential
    # breakpoint()

    for k, v in break_bind_atoms.items():
        if v in ["HW", "HO"]:
            molB.atoms[k].type = "H1"
            if molB.atoms[k].typeB is not None:
                molB.atoms[k].typeB = "H1"

    # update is_radical attribute of Atom objects in topology
    molB.find_radicals()

    return molB


def merge_top_slow_growth(topA: Topology, topB: Topology) -> Topology:
    """Takes two Topologies and joins them for a smooth free-energy like parameter transition simulation.


    TODO: for now this assumes that only one moleculeype (the first, index 0) is of interest.
    """

    molA = topA.moleculetypes[REACTIVE_MOLECULEYPE]
    molB = topB.moleculetypes[REACTIVE_MOLECULEYPE]
    molB = merge_top_moleculetypes_slow_growth(molA, molB, topB.ff)

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
