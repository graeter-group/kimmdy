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
    MultipleDihedrals,
    Pair,
    Atomic,
    AtomicType,
)
from kimmdy.topology.utils import (
    match_atomic_item_to_atomic_type,
    get_protein_section,
    set_protein_section,
)


def is_parameterized(entry: Atomic):
    """Parameterized topology entries have c0 and c1 attributes != None"""
    # TODO: this will fail for Atom (no c0 and c1)
    # should we rename mass and charge of Atom to c0 and c1
    # instead of patching our way aroiund this here?
    return entry.c0 != None and entry.c1 != None


def get_explicit_MultipleDihedrals(dihedral_key: list[str], top: Topology):
    type_key = [top.atoms[x].type for x in dihedral_key]
    multiple_dihedrals = MultipleDihedrals(*dihedral_key, "9", {})
    for periodicity in range(0, 7):
        match_obj = match_atomic_item_to_atomic_type(
            type_key, top.ff.proper_dihedraltypes, str(periodicity)
        )
        if match_obj:
            multiple_dihedrals.dihedrals[str(periodicity)] = Dihedral(
                *dihedral_key, "9", match_obj.c0, match_obj.c1, match_obj.periodicity
            )

    return multiple_dihedrals


def get_atomicobj(key: list[str], type: Atomic, top: Topology, periodicity: str = ""):
    # TODO: ugly
    key = tuple(key)
    type_key = [top.atoms[x].type for x in key]
    if type == Bond:
        instance = top.bonds.get(key)
        match_obj = match_atomic_item_to_atomic_type(type_key, top.ff.bondtypes)
    elif type == Angle:
        instance = top.angles.get(key)
        match_obj = match_atomic_item_to_atomic_type(type_key, top.ff.angletypes)
    elif type == Dihedral:
        if periodicity:
            multiple_instance = top.proper_dihedrals.get(key)
            instance = multiple_instance.dihedrals.get(periodicity)
            match_obj = match_atomic_item_to_atomic_type(
                type_key, top.ff.proper_dihedraltypes, periodicity
            )
        else:
            instance = top.improper_dihedrals.get(key)
            match_obj = match_atomic_item_to_atomic_type(
                type_key, top.ff.improper_dihedraltypes, "2"
            )
    else:
        raise ValueError(f"Could not match type {type} of atomic object.")
    if instance:
        if is_parameterized(instance):
            return instance

    if match_obj:
        return match_obj
    else:
        raise ValueError(
            f"Match object is {match_obj} for {[type_key,instance]} of type {type}."
        )


def get_keys(atomicA: dict, atomicB: dict) -> tuple[list, list, list]:
    keysA = set(atomicA.keys())
    keysB = set(atomicB.keys())

    same = set.intersection(keysA, keysB)
    breaking = keysA - keysB
    binding = keysB - keysA

    return list(same), list(breaking), list(binding)


def merge_same(same: list, topA: Topology, topB: Topology, type: Atomic):
    # TODO:
    # also ugly
    # should also be able to deal with dihedral
    for key in same:
        if type is Bond:
            instanceA = topA.bonds.get(key)
            instanceB = topB.bonds.get(key)
        elif type is Angle:
            instanceA = topA.angles.get(key)
            instanceB = topB.angles.get(key)
        if instanceA != instanceB:
            # assuming no instance has explicit standard ff parameters
            instance_objA = get_atomicobj(key, type, topA)
            instance_objB = get_atomicobj(key, type, topB)
            instanceB.c2 = deepcopy(instance_objB.c0)
            instanceB.c3 = deepcopy(instance_objB.c1)
            instanceB.c0 = deepcopy(instance_objA.c0)
            instanceB.c1 = deepcopy(instance_objA.c1)


def merge_top_parameter_growth(
    topA: Topology, topB: Topology, focus_nr: Union[list[str], None] = None
) -> Topology:
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
                    f"Atom {nr} changed during changemanager step but not the charges!"
                )

    ## bonds
    same, breaking, binding = get_keys(topA.bonds, topB.bonds)

    merge_same(same, topA, topB, Bond)

    for bond_key in breaking:
        topB.bonds[bond_key] = Bond(*bond_key, "1")
        bond_objA = get_atomicobj(bond_key, Bond, topA)

        topB.bonds[(bond_key)] = Bond(
            *bond_key,
            funct="3",
            c0=deepcopy(bond_objA.c0),
            c1=deepcopy(hyperparameters["morse_well_depth"]),
            c2=deepcopy(hyperparameters["morse_steepness"]),
            c3=(
                f"{float(deepcopy(bond_objA.c0))*hyperparameters['morse_dist_factor']:7.5f}"
            ),
            c4="0.00",
            c5=deepcopy(hyperparameters["morse_steepness"]),
        )

        # update bound_to
        atompair = [topB.atoms[bond_key[0]], topB.atoms[bond_key[1]]]
        atompair[0].bound_to_nrs.append(atompair[1].nr)
        atompair[1].bound_to_nrs.append(atompair[0].nr)

    for bond_key in binding:
        bond_objB = get_atomicobj(bond_key, Bond, topB)

        topB.bonds[(bond_key)] = Bond(
            *bond_key,
            funct="3",
            c0=(
                f"{float(deepcopy(bond_objB.c0))*hyperparameters['morse_dist_factor']:7.5f}"
            ),
            c1="0.00",
            c2=deepcopy(hyperparameters["morse_steepness"]),
            c3=deepcopy(bond_objB.c0),
            c4=deepcopy(hyperparameters["morse_well_depth"]),
            c5=deepcopy(hyperparameters["morse_steepness"]),
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

    merge_same(same, topA, topB, Angle)

    for angle_key in breaking:
        topB.angles[angle_key] = Angle(*angle_key, "1")
        angle_objA = get_atomicobj(angle_key, Angle, topA)

        topB.angles[(angle_key)] = Angle(
            *angle_key,
            funct="1",
            c0=deepcopy(angle_objA.c0),
            c1=deepcopy(angle_objA.c1),
            c2=deepcopy(angle_objA.c0),
            c3="0.00",
        )

    for angle_key in binding:
        angle_objB = get_atomicobj(angle_key, Angle, topB)

        topB.angles[(angle_key)] = Angle(
            *angle_key,
            funct="1",
            c0=deepcopy(angle_objB.c0),
            c1="0.00",
            c2=deepcopy(angle_objB.c0),
            c3=deepcopy(angle_objB.c1),
        )

    ## dihedrals
    ## porper dihedrals
    # if indices change atomtypes and parameters change because of that, it will ignore these parameter change

    same, breaking, binding = get_keys(topA.proper_dihedrals, topB.proper_dihedrals)

    for dihedral_key in same:
        multiple_dihedralsA = topA.proper_dihedrals.get(dihedral_key)
        multiple_dihedralsB = topB.proper_dihedrals.get(dihedral_key)
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

            # assuming no dihedral has explicit standard ff parameters

            periodicity_same, periodicity_breaking, periodicity_binding = get_keys(
                multiple_dihedralsA.dihedrals, multiple_dihedralsB.dihedrals
            )

            for periodicity in periodicity_same:
                dihedralA = multiple_dihedralsA.dihedrals[periodicity]
                dihedralB = multiple_dihedralsB.dihedrals[periodicity]
                if dihedralA != dihedralB:
                    dihedral_objA = get_atomicobj(
                        dihedral_key, Dihedral, topA, periodicity
                    )
                    dihedral_objB = get_atomicobj(
                        dihedral_key, Dihedral, topB, periodicity
                    )

                    dihedralB.c3 = deepcopy(dihedral_objB.c0)
                    dihedralB.c4 = deepcopy(dihedral_objB.c1)
                    dihedralB.c5 = deepcopy(dihedral_objB.periodicity)
                    dihedralB.c0 = deepcopy(dihedral_objA.c0)
                    dihedralB.c1 = deepcopy(dihedral_objA.c1)
                    dihedralB.periodicity = deepcopy(dihedral_objA.periodicity)

            for periodicity in periodicity_breaking:
                dihedral_objA = get_atomicobj(dihedral_key, Dihedral, topA, periodicity)

                multiple_dihedralsB.dihedrals[periodicity] = Dihedral(
                    *dihedral_key,
                    funct="9",
                    c0=deepcopy(dihedral_objA.c0),
                    c1=deepcopy(dihedral_objA.c1),
                    periodicity=periodicity,
                    c3=deepcopy(dihedral_objA.c0),
                    c4="0.00",
                    c5=periodicity,
                )

            for periodicity in periodicity_binding:
                dihedral_objB = get_atomicobj(dihedral_key, Dihedral, topB, periodicity)

                multiple_dihedralsB.dihedrals[periodicity] = Dihedral(
                    *dihedral_key,
                    funct="9",
                    c0=deepcopy(dihedral_objB.c0),
                    c1="0.00",
                    periodicity=periodicity,
                    c3=deepcopy(dihedral_objB.c0),
                    c4=deepcopy(dihedral_objB.c1),
                    c5=periodicity,
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
        for periodicity, dihedral_objA in multiple_dihedralsA.dihedrals.items():
            multiple_dihedralsB.dihedrals[periodicity] = Dihedral(
                *dihedral_key,
                funct="9",
                c0=deepcopy(dihedral_objA.c0),
                c1=deepcopy(dihedral_objA.c1),
                periodicity=periodicity,
                c3=deepcopy(dihedral_objA.c0),
                c4="0.00",
                c5=periodicity,
            )
        topB.proper_dihedrals[dihedral_key] = multiple_dihedralsB

    for dihedral_key in binding:
        multiple_dihedralsB = topB.proper_dihedrals.get(dihedral_key)
        multiple_dihedralsB = (
            get_explicit_MultipleDihedrals(dihedral_key, topB)
            if "" in multiple_dihedralsB.dihedrals.keys()
            else multiple_dihedralsB
        )
        for periodicity, dihedral_objB in multiple_dihedralsB.dihedrals.items():
            multiple_dihedralsB.dihedrals[periodicity] = Dihedral(
                *dihedral_key,
                funct="9",
                c0=deepcopy(dihedral_objB.c0),
                c1="0.00",
                periodicity=periodicity,
                c3=deepcopy(dihedral_objB.c0),
                c4=deepcopy(dihedral_objB.c1),
                c5=periodicity,
            )
        topB.proper_dihedrals[dihedral_key] = multiple_dihedralsB

    ## update is_radical attribute of Atom objects in topology
    topB._test_for_radicals()

    ## improper dihedrals
    same, breaking, binding = get_keys(topA.improper_dihedrals, topB.improper_dihedrals)

    for dihedral_key in same:
        dihedralA = topA.improper_dihedrals.get(dihedral_key)
        dihedralB = topB.improper_dihedrals.get(dihedral_key)
        if dihedralA != dihedralB:
            # convert implicit standard ff parameters to explicit, if necessary
            dihedral_objA = get_atomicobj(dihedral_key, Dihedral, topA)
            dihedral_objB = get_atomicobj(dihedral_key, Dihedral, topB)

            dihedralB.c3 = deepcopy(dihedral_objB.c0)
            dihedralB.c4 = deepcopy(dihedral_objB.c1)
            dihedralB.c5 = deepcopy(dihedral_objB.periodicity)
            dihedralB.c0 = deepcopy(dihedral_objA.c0)
            dihedralB.c1 = deepcopy(dihedral_objA.c1)
            dihedralB.periodicity = deepcopy(dihedral_objB.periodicity)

    for dihedral_key in breaking:
        dihedral_objA = get_atomicobj(dihedral_key, Dihedral, topA)

        topB.improper_dihedrals[(dihedral_key)] = Dihedral(
            *dihedral_key,
            funct="4",
            c0=deepcopy(dihedral_objA.c0),
            c1=deepcopy(dihedral_objA.c1),
            periodicity=deepcopy(dihedral_objA.periodicity),
            c3=deepcopy(dihedral_objA.c0),
            c4="0.00",
            c5=deepcopy(dihedral_objA.periodicity),
        )

    for dihedral_key in binding:
        dihedral_objB = get_atomicobj(dihedral_key, Dihedral, topB)

        topB.improper_dihedrals[(dihedral_key)] = Dihedral(
            *dihedral_key,
            funct="4",
            c0=deepcopy(dihedral_objB.c0),
            c1="0.00",
            periodicity=deepcopy(dihedral_objB.periodicity),
            c3=deepcopy(dihedral_objB.c0),
            c4=deepcopy(dihedral_objB.c1),
            c5=deepcopy(dihedral_objB.periodicity),
        )
    return topB
