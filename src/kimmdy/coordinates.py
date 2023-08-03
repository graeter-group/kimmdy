"""
TODO: WIP
"""
import logging
from typing import Union
from copy import deepcopy

from kimmdy.parsing import read_top
from kimmdy.tasks import TaskFiles
from kimmdy.topology.topology import Topology
from kimmdy.topology.atomic import (
    Bond,
    Angle,
    Dihedral,
    DihedralType,
    MultipleDihedrals,
    Pair,
    Interaction,
    Interactions,
    InteractionType,
    InteractionTypes,
    InteractionIds,
    BondId,
)
from kimmdy.topology.utils import (
    match_atomic_item_to_atomic_type,
    get_protein_section,
    set_protein_section,
)


def is_parameterized(entry: Interaction):
    """Parameterized topology entries have c0 and c1 attributes != None"""
    return entry.c0 != None and entry.c1 != None


def get_keys(
    interactionA: Interactions, interactionB: Interactions
) -> tuple[list[InteractionIds], list[InteractionIds], list[InteractionIds]]:
    """For a pair of Interactions, returns same keys, keys only in
    interactionA == breaking and keys only in interactionB == binding
    """
    # also works for MultipleDihedrals.dihedrals
    keysA = set(interactionA.keys())
    keysB = set(interactionB.keys())

    same: set[InteractionIds] = set.intersection(keysA, keysB)
    breaking = keysA - keysB
    binding = keysB - keysA

    return list(same), list(breaking), list(binding)


def get_explicit_MultipleDihedrals(
    dihedral_key: tuple[str, str, str, str], top: Topology, periodicity_max: int = 6
):
    """Takes a valid dihedral key and returns explicit
    dihedral parameters for a given topology
    """
    type_key = [top.atoms[x].type for x in dihedral_key]

    multiple_dihedrals = MultipleDihedrals(*dihedral_key, "9", {})
    for periodicity in range(1, periodicity_max + 1):
        match_obj = match_atomic_item_to_atomic_type(
            type_key, top.ff.proper_dihedraltypes, str(periodicity)
        )
        if match_obj:
            assert isinstance(match_obj, DihedralType)
            l = [*dihedral_key, "9", match_obj.c0, match_obj.c1, match_obj.periodicity]
            multiple_dihedrals.dihedrals[str(periodicity)] = Dihedral.from_top_line(l)

    return multiple_dihedrals


def get_explicit_or_type(
    key: tuple[str, ...],
    interaction: Interaction,
    interaction_types: InteractionTypes,
    top: Topology,
    periodicity: str = "",
) -> Union[Interaction, InteractionType]:
    """Takes an Interaction and associated key, InteractionTypes, Topology
    and Periodicity (for dihedrals) and returns an object with the parameters of this Interaction
    """
    if is_parameterized(interaction):
        return interaction

    type_key = [top.atoms[x].type for x in key]
    match_obj = match_atomic_item_to_atomic_type(
        type_key, interaction_types, periodicity
    )

    if match_obj:
        assert isinstance(match_obj, InteractionType)
        return match_obj
    else:
        raise ValueError(
            f"Could not find explicit parameters for {[type_key,interaction]} in line or in force field attached to {top}."
        )


def merge_parameterized_dihedrals(
    dihedral_key: tuple[str, str, str, str],
    dihedralA: Union[Dihedral, DihedralType, None],
    dihedralB: Union[Dihedral, DihedralType, None],
    funct: str,
) -> Dihedral:
    """Merge one to two Dihedrals or -Types into a Dihedral in free-energy syntax"""
    if dihedralA and dihedralB:
        # same
        dihedralmerge = Dihedral(
            *dihedral_key,
            funct=funct,
            c0=dihedralA.c0,
            c1=dihedralA.c1,
            periodicity=dihedralA.periodicity,
            c3=dihedralB.c0,
            c4=dihedralB.c1,
            c5=dihedralB.periodicity,
        )
    elif dihedralA:
        # breaking
        dihedralmerge = Dihedral(
            *dihedral_key,
            funct=funct,
            c0=dihedralA.c0,
            c1=dihedralA.c1,
            periodicity=dihedralA.periodicity,
            c3=dihedralA.c0,
            c4="0.00",
            c5=dihedralA.periodicity,
        )
    elif dihedralB:
        # binding
        dihedralmerge = Dihedral(
            *dihedral_key,
            funct="9",
            c0=dihedralB.c0,
            c1="0.00",
            periodicity=dihedralB.periodicity,
            c3=dihedralB.c0,
            c4=dihedralB.c1,
            c5=dihedralB.periodicity,
        )
    else:
        raise ValueError(
            f"Tried to merge two dihedrals of {dihedral_key} but no parameterized dihedrals found!"
        )
    return dihedralmerge


def merge_top_parameter_growth(
    topA: Topology, topB: Topology, focus_nr: Union[list[str], None] = None
) -> Topology:
    """Takes two Topologies and joins them for a smooth free-energy like parameter transition simulation"""
    hyperparameters = {
        "morse_well_depth": "400.0",
        "morse_steepness": "10.0",
        "morse_dist_factor": 5,
    }  # well_depth D [kJ/mol], steepness [nm-1], dist_factor for bond length
    logging.info(f"Merging topologies {topA} and {topB}")

    # TODO:
    # think about how to bring focus_nr into this

    ## atoms
    for nr in topA.atoms.keys():
        atomA = topA.atoms[nr]
        atomB = topB.atoms[nr]
        if atomA != atomB:
            if atomA.charge != atomB.charge:
                atomB.typeB = deepcopy(atomB.type)
                atomB.type = deepcopy(atomA.type)
                atomB.chargeB = deepcopy(atomB.charge)
                atomB.charge = deepcopy(atomA.charge)
                atomB.massB = deepcopy(atomB.mass)
                atomB.mass = deepcopy(atomA.mass)
            else:
                logging.debug(
                    f"Atom {nr} with A:{atomA} and B:{atomB} changed during changemanager step but not the charges!"
                )

    ## bonds
    same, breaking, binding = get_keys(topA.bonds, topB.bonds)

    for bond_key in same:
        interactionA = topA.bonds[bond_key]
        interactionB = topB.bonds[bond_key]

        if interactionA != interactionB:
            parameterizedA = get_explicit_or_type(
                bond_key, interactionA, topA.ff.bondtypes, topA
            )
            parameterizedB = get_explicit_or_type(
                bond_key, interactionB, topB.ff.bondtypes, topB
            )

            topB.bonds[bond_key] = Bond(
                *bond_key,
                funct=parameterizedB.funct,
                c0=parameterizedA.c0,
                c1=parameterizedA.c1,
                c2=parameterizedB.c0,
                c3=parameterizedB.c1,
            )

    for bond_key in breaking:
        interactionA = topA.bonds[bond_key]
        parameterizedA = get_explicit_or_type(
            bond_key, interactionA, topA.ff.bondtypes, topA
        )

        topB.bonds[bond_key] = Bond(
            *bond_key,
            funct="3",
            c0=parameterizedA.c0,
            c1=hyperparameters["morse_well_depth"],
            c2=hyperparameters["morse_steepness"],
            c3=f"{float(parameterizedA.c0)*hyperparameters['morse_dist_factor']:7.5f}",
            c4="0.00",
            c5=hyperparameters["morse_steepness"],
        )

        # update bound_to
        atompair = [topB.atoms[bond_key[0]], topB.atoms[bond_key[1]]]
        atompair[0].bound_to_nrs.append(atompair[1].nr)
        atompair[1].bound_to_nrs.append(atompair[0].nr)

    for bond_key in binding:
        interactionB = topB.bonds[bond_key]
        parameterizedB = get_explicit_or_type(
            bond_key, interactionB, topB.ff.bondtypes, topB
        )

        topB.bonds[bond_key] = Bond(
            *bond_key,
            funct="3",
            c0=f"{float(parameterizedB.c0)*hyperparameters['morse_dist_factor']:7.5f}",
            c1="0.00",
            c2=hyperparameters["morse_steepness"],
            c3=parameterizedB.c0,
            c4=hyperparameters["morse_well_depth"],
            c5=hyperparameters["morse_steepness"],
        )

    ## pairs and exclusions
    exclusions_content = get_protein_section(topB.top, "exclusions")
    if exclusions_content is None:
        # maybe hook this up to empty_sections if it gets accessible
        exclusions = {"content": [], "else_content": [], "extra": [], "condition": None}
        topB.top["moleculetype_0"]["subsections"]["exclusions"] = exclusions
        exclusions_content = exclusions["content"]

    for bond_key in breaking:
        topB.pairs.pop(bond_key, None)
        exclusions_content.append(list(bond_key))

    set_protein_section(topB.top, "exclusions", exclusions_content)

    ## angles
    same, breaking, binding = get_keys(topA.angles, topB.angles)

    for angle_key in same:
        interactionA = topA.angles[angle_key]
        interactionB = topB.angles[angle_key]

        if interactionA != interactionB:
            parameterizedA = get_explicit_or_type(
                angle_key, interactionA, topA.ff.angletypes, topA
            )
            parameterizedB = get_explicit_or_type(
                angle_key, interactionB, topB.ff.angletypes, topB
            )
            topB.angles[angle_key] = Angle(
                *angle_key,
                funct=parameterizedB.funct,
                c0=parameterizedA.c0,
                c1=parameterizedA.c1,
                c2=parameterizedB.c0,
                c3=parameterizedB.c1,
            )

    for angle_key in breaking:
        interactionA = topA.angles[angle_key]
        parameterizedA = get_explicit_or_type(
            angle_key, interactionA, topA.ff.angletypes, topA
        )
        topB.angles[angle_key] = Angle(
            *angle_key,
            funct="1",
            c0=parameterizedA.c0,
            c1=parameterizedA.c1,
            c2=parameterizedA.c0,
            c3="0.00",
        )

    for angle_key in binding:
        interactionB = topB.angles[angle_key]
        parameterizedB = get_explicit_or_type(
            angle_key, interactionB, topB.ff.angletypes, topB
        )
        topB.angles[angle_key] = Angle(
            *angle_key,
            funct="1",
            c0=parameterizedB.c0,
            c1="0.00",
            c2=parameterizedB.c0,
            c3=parameterizedB.c1,
        )

    ## dihedrals
    ## proper dihedrals
    # if indices change atomtypes and parameters change because of that, it will ignore these parameter change

    same, breaking, binding = get_keys(topA.proper_dihedrals, topB.proper_dihedrals)

    for dihedral_key in same:
        # dihedrals have a nested structure and need a different treatment from bonds and angles
        multiple_dihedralsA = topA.proper_dihedrals[dihedral_key]
        multiple_dihedralsB = topB.proper_dihedrals[dihedral_key]
        if multiple_dihedralsA != multiple_dihedralsB:
            # convert implicit standard ff parameters to explicit
            multiple_dihedralsA = (
                get_explicit_MultipleDihedrals(dihedral_key, topA)
                if "" in multiple_dihedralsA.dihedrals.keys()
                else multiple_dihedralsA
            )
            multiple_dihedralsB = (
                get_explicit_MultipleDihedrals(dihedral_key, topB)
                if "" in multiple_dihedralsB.dihedrals.keys()
                else multiple_dihedralsB
            )

            # periodicity-wise construction of merged parameters
            periodicity_same, periodicity_breaking, periodicity_binding = get_keys(
                multiple_dihedralsA.dihedrals, multiple_dihedralsB.dihedrals
            )

            for periodicity in periodicity_same:
                interactionA = multiple_dihedralsA.dihedrals[periodicity]
                interactionB = multiple_dihedralsB.dihedrals[periodicity]
                if interactionA != interactionB:
                    parameterizedA = get_explicit_or_type(
                        dihedral_key,
                        interactionA,
                        topA.ff.proper_dihedraltypes,
                        topA,
                        periodicity,
                    )
                    parameterizedB = get_explicit_or_type(
                        dihedral_key,
                        interactionB,
                        topB.ff.proper_dihedraltypes,
                        topB,
                        periodicity,
                    )
                    multiple_dihedralsB.dihedrals[
                        periodicity
                    ] = merge_parameterized_dihedrals(
                        dihedral_key, parameterizedA, parameterizedB, "9"
                    )

            for periodicity in periodicity_breaking:
                interactionA = multiple_dihedralsA.dihedrals[periodicity]
                parameterizedA = get_explicit_or_type(
                    dihedral_key,
                    interactionA,
                    topA.ff.proper_dihedraltypes,
                    topA,
                    periodicity,
                )
                multiple_dihedralsB.dihedrals[
                    periodicity
                ] = merge_parameterized_dihedrals(
                    dihedral_key, parameterizedA, None, "9"
                )

            for periodicity in periodicity_binding:
                interactionB = multiple_dihedralsB.dihedrals[periodicity]
                parameterizedB = get_explicit_or_type(
                    dihedral_key,
                    interactionB,
                    topB.ff.proper_dihedraltypes,
                    topB,
                    periodicity,
                )
                multiple_dihedralsB.dihedrals[
                    periodicity
                ] = merge_parameterized_dihedrals(
                    dihedral_key, None, parameterizedB, "9"
                )
        topB.proper_dihedrals[dihedral_key] = multiple_dihedralsB

    for dihedral_key in breaking:
        multiple_dihedralsA = topA.proper_dihedrals.get(dihedral_key)
        multiple_dihedralsA = (
            get_explicit_MultipleDihedrals(dihedral_key, topA)
            if "" in multiple_dihedralsA.dihedrals.keys()
            else multiple_dihedralsA
        )
        multiple_dihedralsB = MultipleDihedrals(*dihedral_key, "9", {})
        for periodicity, parameterizedA in multiple_dihedralsA.dihedrals.items():
            multiple_dihedralsB.dihedrals[periodicity] = merge_parameterized_dihedrals(
                dihedral_key, parameterizedA, None, "9"
            )
        topB.proper_dihedrals[dihedral_key] = multiple_dihedralsB

    for dihedral_key in binding:
        multiple_dihedralsB = topB.proper_dihedrals.get(dihedral_key)
        multiple_dihedralsB = (
            get_explicit_MultipleDihedrals(dihedral_key, topB)
            if "" in multiple_dihedralsB.dihedrals.keys()
            else multiple_dihedralsB
        )
        for periodicity, parameterizedB in multiple_dihedralsB.dihedrals.items():
            multiple_dihedralsB.dihedrals[periodicity] = merge_parameterized_dihedrals(
                dihedral_key, None, parameterizedB, "9"
            )
        topB.proper_dihedrals[dihedral_key] = multiple_dihedralsB

    ## improper dihedrals
    # all impropers in amber99SB ffbonded.itp have a periodicity of 2
    # but not the ones defined in aminoacids.rtp. For now, I am assuming
    # a periodicity of 2 in this section
    same, breaking, binding = get_keys(topA.improper_dihedrals, topB.improper_dihedrals)

    for dihedral_key in same:
        interactionA = topA.improper_dihedrals[dihedral_key]
        interactionB = topB.improper_dihedrals[dihedral_key]
        if interactionA != interactionB:
            # convert implicit standard ff parameters to explicit, if necessary
            parameterizedA = get_explicit_or_type(
                dihedral_key, interactionA, topA.ff.improper_dihedraltypes, topA, "2"
            )
            parameterizedB = get_explicit_or_type(
                dihedral_key, interactionB, topB.ff.improper_dihedraltypes, topB, "2"
            )
            topB.improper_dihedrals[dihedral_key] = merge_parameterized_dihedrals(
                dihedral_key, parameterizedA, parameterizedB, "4"
            )

    for dihedral_key in breaking:
        interactionA = topA.improper_dihedrals[dihedral_key]
        parameterizedA = get_explicit_or_type(
            dihedral_key, interactionA, topA.ff.improper_dihedraltypes, topA, "2"
        )
        topB.improper_dihedrals[dihedral_key] = merge_parameterized_dihedrals(
            dihedral_key, parameterizedA, None, "4"
        )

    for dihedral_key in binding:
        interactionB = topB.improper_dihedrals[dihedral_key]
        parameterizedB = get_explicit_or_type(
            dihedral_key, interactionB, topB.ff.improper_dihedraltypes, topB, "2"
        )
        topB.improper_dihedrals[dihedral_key] = merge_parameterized_dihedrals(
            dihedral_key, None, parameterizedB, "4"
        )

    ## update is_radical attribute of Atom objects in topology
    topB._test_for_radicals()

    return topB
