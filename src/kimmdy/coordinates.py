import logging
from typing import Union
from copy import deepcopy

from kimmdy.parsing import read_top
from kimmdy.tasks import TaskFiles
from kimmdy.topology.topology import Topology
from kimmdy.topology.atomic import Bond, Angle, Dihedral, Pair, Atomic, AtomicType
from kimmdy.topology.utils import (
    match_atomic_item_to_atomic_type,
    get_protein_section,
    set_protein_section,
)


def is_parameterized(entry: Atomic):
    """Parameterized topology entries have c0 and c1 attributes != None"""
    return entry.c0 != None and entry.c1 != None


def get_atomicobj(key: list[str], type: AtomicType, top: Topology):
    # ugly
    key = tuple(key)
    type_key = [top.atoms[x].type for x in key]
    if type == Bond:
        instance = top.bonds.get(key)
        match_obj = match_atomic_item_to_atomic_type(type_key, top.ff.bondtypes)
    elif type == Angle:
        instance = top.angles.get(key)
        match_obj = match_atomic_item_to_atomic_type(type_key, top.ff.angletypes)
    elif type == Dihedral:
        instance = top.dihedrals.get(key)
        if instance.funct == "4":
            match_obj = match_atomic_item_to_atomic_type(
                type_key, top.ff.improper_dihedraltypes
            )
        elif instance.funct == "9":
            match_obj = match_atomic_item_to_atomic_type(
                type_key, top.ff.proper_dihedraltypes
            )
    else:
        raise ValueError(f"Could not match type {type} of atomic object.")
    if is_parameterized(instance):
        return instance
    else:
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


def merge_top_prmgrowth(
    files: TaskFiles, focus_nr: Union[list[str], None] = None
) -> Topology:
    # not the most robust way to get topA and topB
    hyperprms = {
        "morse_well_depth": "400.0",
        "morse_steepness": "10.0",
        "morse_dist_factor": 5,
    }  # well_depth D [kJ/mol], steepness [nm-1], dist_factor for bond length
    logging.info(f"Merging topologies {files.input['top']} and {files.output['top']}")
    topADict = read_top(files.input["top"])
    topBDict = read_top(files.output["top"])
    topA = Topology(topADict)
    topB = Topology(topBDict)

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
            c1=deepcopy(hyperprms["morse_well_depth"]),
            c2=deepcopy(hyperprms["morse_steepness"]),
            c3=(f"{float(deepcopy(bond_objA.c0))*hyperprms['morse_dist_factor']:7.5f}"),
            c4="0.00",
            c5=deepcopy(hyperprms["morse_steepness"]),
        )

    for bond_key in binding:
        bond_objB = get_atomicobj(bond_key, Bond, topB)

        topB.bonds[(bond_key)] = Bond(
            *bond_key,
            funct="3",
            c0=(f"{float(deepcopy(bond_objB.c0))*hyperprms['morse_dist_factor']:7.5f}"),
            c1="0.00",
            c2=deepcopy(hyperprms["morse_steepness"]),
            c3=deepcopy(bond_objB.c0),
            c4=deepcopy(hyperprms["morse_well_depth"]),
            c5=deepcopy(hyperprms["morse_steepness"]),
        )

    ## pairs and exclusions
    try:
        exclusions_content = get_protein_section(topB.top, "exclusions")

    except ValueError:
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
    # dihedraltypes does not work at the moment
    # proper_dihedraltypes also has the periodicity as key
    # maybe match_atomic_item_to_atomic_type could give us dihedraltypes with all periodicities
    # compare entry-wise

    # same, breaking, binding = get_keys(topA.dihedrals, topB.dihedrals)

    # for dihedral_key in same:
    #     dihedralA = topA.proper_dihedrals.get(dihedral_key)
    #     dihedralB = topB.proper_dihedrals.get(dihedral_key)
    #     if dihedralA != dihedralB:
    #         # assuming no dihedral has explicit standard ff parameters
    #         dihedral_objA = get_atomicobj(dihedral_key, dihedralA, topA)
    #     dihedral_objB = get_atomicobj(dihedral_key, Dihedral, topB)
    #     dihedralB.c0 = deepcopy(dihedral_objB.c0)
    #     dihedralB.c1 = deepcopy(dihedral_objB.c1)
    #     dihedralB.c2 = deepcopy(dihedral_objB.c2)
    #     dihedralB.c3 = deepcopy(dihedral_objB.c3)
    #     dihedralB.c4 = deepcopy(dihedral_objB.c4)
    #     dihedralB.c5 = deepcopy(dihedral_objB.c5)

    return topB
