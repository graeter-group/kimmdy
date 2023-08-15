import logging
from typing import Union
from copy import deepcopy
from kimmdy.topology.topology import Topology
from kimmdy.topology.atomic import (
    Bond,
    Angle,
    Dihedral,
    DihedralType,
    MultipleDihedrals,
    Interaction,
    InteractionType,
    InteractionTypes,
)
from kimmdy.topology.utils import (
    match_atomic_item_to_atomic_type,
    get_protein_section,
    set_protein_section,
)

logger = logging.getLogger(__name__)


def is_parameterized(entry: Interaction):
    """Parameterized topology entries have c0 and c1 attributes != None"""
    return entry.c0 != None and entry.c1 != None


def get_explicit_MultipleDihedrals(
    dihedral_key: tuple[str, str, str, str],
    top: Topology,
    dihedrals_in: Union[MultipleDihedrals, None],
    periodicity_max: int = 6,
) -> Union[MultipleDihedrals, None]:
    """Takes a valid dihedral key and returns explicit
    dihedral parameters for a given topology
    """
    if not dihedrals_in:
        return None

    if not "" in dihedrals_in.dihedrals.keys():
        # empty string means implicit parameters
        return dihedrals_in

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

    if not multiple_dihedrals.dihedrals:
        return None

    return multiple_dihedrals


def get_explicit_or_type(
    key: tuple[str, ...],
    interaction: Union[Interaction, None],
    interaction_types: InteractionTypes,
    top: Topology,
    periodicity: str = "",
) -> Union[Interaction, InteractionType, None]:
    """Takes an Interaction and associated key, InteractionTypes, Topology
    and Periodicity (for dihedrals) and returns an object with the parameters of this Interaction
    """
    if not interaction:
        return None

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


def merge_dihedrals(
    dihedral_key: tuple[str, str, str, str],
    interactionA: Union[Dihedral, None],
    interactionB: Union[Dihedral, None],
    interaction_typesA: dict[tuple[str, ...], DihedralType],
    interaction_typesB: dict[tuple[str, ...], DihedralType],
    topA: Topology,
    topB: Topology,
    funct: str,
    periodicity: str,
) -> Dihedral:
    """Merge one to two Dihedrals or -Types into a Dihedral in free-energy syntax"""
    # convert implicit standard ff parameters to explicit, if necessary
    if interactionA:
        parameterizedA = get_explicit_or_type(
            dihedral_key,
            interactionA,
            interaction_typesA,
            topA,
            periodicity,
        )
    else:
        parameterizedA = None

    if interactionB:
        parameterizedB = get_explicit_or_type(
            dihedral_key,
            interactionB,
            interaction_typesB,
            topB,
            periodicity,
        )
    else:
        parameterizedB = None

    # construct parameterized Dihedral
    if parameterizedA and parameterizedB:
        # same
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
    elif parameterizedA:
        # breaking
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
    elif parameterizedB:
        # binding
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
                logger.debug(
                    f"Atom {nr} with A:{atomA} and B:{atomB} changed during changemanager step but not the charges!"
                )

    ## bonds
    keysA = set(topA.bonds.keys())
    keysB = set(topB.bonds.keys())
    keys = keysA | keysB

    for key in keys:
        interactionA = topA.bonds.get(key)
        interactionB = topB.bonds.get(key)

        if interactionA != interactionB:
            parameterizedA = get_explicit_or_type(
                key, interactionA, topA.ff.bondtypes, topA
            )
            parameterizedB = get_explicit_or_type(
                key, interactionB, topB.ff.bondtypes, topB
            )
            if parameterizedA and parameterizedB:
                topB.bonds[key] = Bond(
                    *key,
                    funct=parameterizedB.funct,
                    c0=parameterizedA.c0,
                    c1=parameterizedA.c1,
                    c2=parameterizedB.c0,
                    c3=parameterizedB.c1,
                )
            elif parameterizedA:
                topB.bonds[key] = Bond(
                    *key,
                    funct="3",
                    c0=parameterizedA.c0,
                    c1=hyperparameters["morse_well_depth"],
                    c2=hyperparameters["morse_steepness"],
                    c3=f"{float(parameterizedA.c0)*hyperparameters['morse_dist_factor']:7.5f}",
                    c4="0.00",
                    c5=hyperparameters["morse_steepness"],
                )

                # update bound_to
                atompair = [topB.atoms[key[0]], topB.atoms[key[1]]]
                atompair[0].bound_to_nrs.append(atompair[1].nr)
                atompair[1].bound_to_nrs.append(atompair[0].nr)

            elif parameterizedB:
                topB.bonds[key] = Bond(
                    *key,
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
    for key in keysA - keysB:
        topB.pairs.pop(key, None)
        exclusions_content.append(list(key))

    set_protein_section(topB.top, "exclusions", exclusions_content)

    ## angles
    keys = set(topA.angles.keys()) | set(topB.angles.keys())
    for key in keys:
        interactionA = topA.angles.get(key)
        interactionB = topB.angles.get(key)

        if interactionA != interactionB:
            parameterizedA = get_explicit_or_type(
                key, interactionA, topA.ff.angletypes, topA
            )
            parameterizedB = get_explicit_or_type(
                key, interactionB, topB.ff.angletypes, topB
            )
            if parameterizedA and parameterizedB:
                topB.angles[key] = Angle(
                    *key,
                    funct=parameterizedB.funct,
                    c0=parameterizedA.c0,
                    c1=parameterizedA.c1,
                    c2=parameterizedB.c0,
                    c3=parameterizedB.c1,
                )
            elif parameterizedA:
                topB.angles[key] = Angle(
                    *key,
                    funct="1",
                    c0=parameterizedA.c0,
                    c1=parameterizedA.c1,
                    c2=parameterizedA.c0,
                    c3="0.00",
                )
            elif parameterizedB:
                topB.angles[key] = Angle(
                    *key,
                    funct="1",
                    c0=parameterizedB.c0,
                    c1="0.00",
                    c2=parameterizedB.c0,
                    c3=parameterizedB.c1,
                )
            else:
                logger.warning(f"Could not parameterize angle {key}.")

    ## dihedrals
    ## proper dihedrals
    # proper dihedrals have a nested structure and need a different treatment from bonds, angles and improper dihedrals
    # if indices change atomtypes and parameters change because of that, it will ignore these parameter change

    keys = set(topA.proper_dihedrals.keys()) | set(topB.proper_dihedrals.keys())
    for key in keys:
        multiple_dihedralsA = topA.proper_dihedrals.get(key)
        multiple_dihedralsB = topB.proper_dihedrals.get(key)

        if multiple_dihedralsA != multiple_dihedralsB:
            multiple_dihedralsA = get_explicit_MultipleDihedrals(
                key, topA, multiple_dihedralsA
            )
            multiple_dihedralsB = get_explicit_MultipleDihedrals(
                key, topB, multiple_dihedralsB
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

            topB.proper_dihedrals[key] = MultipleDihedrals(*key, "9", {})
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

                topB.proper_dihedrals[key].dihedrals[periodicity] = merge_dihedrals(
                    key,
                    interactionA,
                    interactionB,
                    topA.ff.proper_dihedraltypes,
                    topB.ff.proper_dihedraltypes,
                    topA,
                    topB,
                    "9",
                    periodicity,
                )

    ## improper dihedrals
    # all impropers in amber99SB ffbonded.itp have a periodicity of 2
    # but not the ones defined in aminoacids.rtp. For now, I am assuming
    # a periodicity of 2 in this section

    keys = set(topA.improper_dihedrals.keys()) | set(topB.improper_dihedrals.keys())
    for key in keys:
        interactionA = topA.improper_dihedrals.get(key)
        interactionB = topB.improper_dihedrals.get(key)

        if interactionA != interactionB:
            topB.improper_dihedrals[key] = merge_dihedrals(
                key,
                interactionA,
                interactionB,
                topA.ff.improper_dihedraltypes,
                topB.ff.improper_dihedraltypes,
                topA,
                topB,
                "4",
                "2",
            )

    ## update is_radical attribute of Atom objects in topology
    topB._test_for_radicals()

    return topB
