import MDAnalysis as MDA
from pathlib import Path
import numpy as np
import logging
from typing import Union
from copy import deepcopy

from kimmdy.parsing import read_top
from kimmdy.tasks import TaskFiles
from kimmdy.topology.topology import Topology
from kimmdy.topology.atomic import Bond, Angle, Dihedral, Pair, Atomic
from kimmdy.topology.utils import match_atomic_item_to_atomic_type, get_protein_section, set_protein_section


## copied from changemanager. should be put into utils
def get_ff_sections(ffdir):
    return read_top(ffdir / "ffbonded.itp")


def parameterize_bonded_terms(ffprm, terms_atomtype, prop, terms):
    """
    takes a term (bond or angle) and adds the parameters from
    the force field ffbonded.itp file that match the atomtypes
    """
    print(f"Looking for atom parameters in atoms section of topology for {terms}")
    # prop must be 'bondtypes' or 'angletypes'
    terms_prm = []
    for i, term in enumerate(terms_atomtype):
        for entry in ffprm[prop]:
            if term == entry[: len(term)] or term[::-1] == entry[: len(term)]:
                stop = entry.index(";")
                terms_prm.append([*terms[i], *entry[len(term) : stop]])
                break
        else:
            print(f"No parameters found for {term}/{terms[i]}!")
    print(f"Parameterized these terms: {terms_prm}")
    return terms_prm


# def is_parameterized(term):
#     """does not work for dihedrals
#     basically checks whether there are floats in the term
#     """
#     for idx in term:
#         if not idx.isdigit():
#             if idx.replace(".", "", 1).isdigit():
#                 return True
#     return False
def is_parameterized(entry: Atomic):
    """Parameterized topology entries have c0 and c1 attributes != None"""
    return entry.c0 != None and entry.c1 != None


def get_bondobj(bond_key: list[str], bond: Bond, top: Topology):
    if is_parameterized(bond):
        return bond
    else:
        bondtype_key = [top.atoms[bond_key[0]].type, top.atoms[bond_key[1]].type]
        return match_atomic_item_to_atomic_type(bondtype_key, top.ff.bondtypes)

def get_angleobj(angle_key: list[str], angle: Angle, top: Topology):
    if is_parameterized(angle):
        return angle
    else:
        angletype_key = [top.atoms[x].type for x in angle_key]
        return match_atomic_item_to_atomic_type(angletype_key, top.ff.angletypes)

def get_dihedralobj(dihedral_key: list[str], dihedral: Dihedral, top: Topology):
    if is_parameterized(dihedral):
        return dihedral
    else:
        dihedraltype_key = [top.atoms[x].type for x in dihedral_key]
        return match_atomic_item_to_atomic_type(dihedraltype_key, top.ff.proper_dihedraltypes)

##


def merge_top_prmgrowth(
    files: TaskFiles, focus_nr: Union[list[str], None] = None
) -> Topology:
    # ffdir = (files.input["ff"],)
    # not the most robust way to get topA and topB
    hyperprms = {
        "morse_well_depth": "400.0",
        "morse_steepness": "10.0",
        'morse_dist_factor': 3
    }  # well_depth D [kJ/mol], steepness [nm-1]
    topADict = read_top(files.input["top"])
    topBDict = read_top(files.output["top"])
    topA = Topology(topADict)
    topB = Topology(topBDict)

    # ToDo: what about implicit parameters?? especially dihedrals
    # think about how to bring focus_nr into this
    # for nr in topB.atoms.keys():

    #     # atoms
    #     atomA = topA.atoms[nr]
    #     atomB = topB.atoms[nr]
    #     if atomA != atomB:
    #         if atomA.charge != atomB.charge:
    #             atomB.typeB = deepcopy(atomB.type)
    #             atomB.type = deepcopy(atomA.type)
    #             atomB.chargeB = deepcopy(atomB.charge)
    #             atomB.charge = deepcopy(atomA.charge)
    #             atomB.massB = deepcopy(atomB.mass)
    #             atomB.mass = deepcopy(atomA.mass)
    #         else:
    #             logging.debug(
    #                 f"Atom {nr} changed during changemanager step but not the charges!"
    #             )

    #ToDo: Generalize

    # # bonds
    bondA_keys = set(topA.bonds.keys())
    bondB_keys = set(topB.bonds.keys())

    same = set.intersection(bondA_keys, bondB_keys)
    breaking = bondA_keys - bondB_keys
    binding = bondB_keys - bondA_keys

    for bond_key in same:
        bondA = topA.bonds.get(bond_key)
        bondB = topB.bonds.get(bond_key)
        if bondA != bondB:
            # assuming no bond has explicit standard ff parameters
            bond_objA = get_bondobj(bond_key, bondA, topA)
            bond_objB = get_bondobj(bond_key, bondB, topB)
            bondB.c2 = deepcopy(bond_objB.c0)
            bondB.c3 = deepcopy(bond_objB.c1)
            bondB.c0 = deepcopy(bond_objA.c0)
            bondB.c1 = deepcopy(bond_objA.c1)

    for bond_key in breaking:
        topB.bonds[bond_key] = Bond(*bond_key, "1")
        # topB.bind_bond(bond_key)      # this doesn't work for this use case
        bondB = topB.bonds.get(bond_key)
        bond_objA = get_bondobj(bond_key, bondA, topA)

        bondB.funct = "3"  # Morse potential
        bondB.c0 = deepcopy(bond_objA.c0)
        bondB.c1 = deepcopy(hyperprms["morse_well_depth"])
        bondB.c2 = deepcopy(hyperprms["morse_steepness"])
        bondB.c3 = f"{float(deepcopy(bond_objA.c0))*hyperprms['morse_dist_factor']:7.5f}"
        bondB.c4 = "0.00"
        bondB.c5 = deepcopy(hyperprms["morse_steepness"])

    # deal with pairs and exclusions
    try:
        exclusions_content = get_protein_section(topB.top,"exclusions")

    except ValueError:
        # maybe hook this up to empty_sections if it gets accessible
        exclusions = {"content": [], "else_content": [], "extra": [], "condition": None}
        topB.top["moleculetype_0"]["subsections"]["exclusions"] = exclusions
        exclusions_content = exclusions["content"]

    for bond_key in breaking:
        topB.pairs.pop(bond_key,None)
        exclusions_content.append(list(bond_key))

    set_protein_section(topB.top,"exclusions",exclusions_content)


    for bond_key in binding:
        bondB = topB.bonds.get(bond_key)
        bond_objB = get_bondobj(bond_key, bondB, topB)

        bondB.funct = "3"  # Morse potential
        bondB.c0 = deepcopy(bond_objB.c0)
        bondB.c1 = "0.00"
        bondB.c2 = deepcopy(hyperprms['morse_steepness'])
        bondB.c3 = deepcopy(bond_objB.c0)
        bondB.c4 = deepcopy(hyperprms['morse_well_depth'])
        bondB.c5 = deepcopy(hyperprms['morse_steepness'])

    # # angles
    angleA_keys = set(topA.angles.keys())
    angleB_keys = set(topB.angles.keys())

    same = set.intersection(angleA_keys, angleB_keys)
    breaking = angleA_keys - angleB_keys
    binding = angleB_keys - angleA_keys

    for angle_key in same:
        angleA = topA.angles.get(angle_key)
        angleB = topB.angles.get(angle_key)
        if angleA != angleB:
            # assuming no angle has explicit standard ff parameters
            angle_objA = get_angleobj(angle_key, angleA, topA)
            angle_objB = get_angleobj(angle_key, angleB, topB)
            angleB.c2 = deepcopy(angle_objB.c0)
            angleB.c3 = deepcopy(angle_objB.c1)
            angleB.c0 = deepcopy(angle_objA.c0)
            angleB.c1 = deepcopy(angle_objA.c1)

    for angle_key in breaking:
        topB.angles[angle_key] = Angle(*angle_key, "1")
        # topB.bind_angle(angle_key)      # this doesn't work for this use case
        angleB = topB.angles.get(angle_key)
        angle_objA = get_angleobj(angle_key, angleA, topA)

        angleB.c0 = deepcopy(angle_objA.c0)
        angleB.c1 = deepcopy(angle_objA.c1)
        angleB.c2 = deepcopy(angle_objA.c0)
        angleB.c3 = "0.00"

    for angle_key in binding:
        angleB = topB.angles.get(angle_key)
        angle_objB = get_angleobj(angle_key, angleB, topB)

        angleB.c0 = deepcopy(angle_objB.c0)
        angleB.c1 = "0.00"
        angleB.c2 = deepcopy(angle_objB.c0)
        angleB.c3 = deepcopy(angle_objB.c1)

    # # dihedrals
    # dihedraltypes does not work at the moment
    # proper_dihedraltypes also has the periodicity as key
    # maybe match_atomic_item_to_atomic_type could give us dihedraltypes with all periodicities
    # compare entry-wise

    # dihedralA_keys = set(topA.proper_dihedrals.keys())
    # dihedralB_keys = set(topB.proper_dihedrals.keys())

    # same = set.intersection(dihedralA_keys, dihedralB_keys)
    # breaking = dihedralA_keys - dihedralB_keys
    # binding = dihedralB_keys - dihedralA_keys
    # for dihedral_key in same:
    #     dihedralA = topA.proper_dihedrals.get(dihedral_key)
    #     dihedralB = topB.proper_dihedrals.get(dihedral_key)
    #     # if dihedralA != dihedralB:
    #     #     # assuming no dihedral has explicit standard ff parameters
    #     #     dihedral_objA = get_dihedralobj(dihedral_key, dihedralA, topA)
    #     dihedral_objB = get_dihedralobj(dihedral_key, dihedralB, topB)
    #     dihedralB.c0 = deepcopy(dihedral_objB.c0)
    #     dihedralB.c1 = deepcopy(dihedral_objB.c1)
    #     dihedralB.c2 = deepcopy(dihedral_objB.c2)
    #     dihedralB.c3 = deepcopy(dihedral_objB.c3)
    #     dihedralB.c4 = deepcopy(dihedral_objB.c4)
    #     dihedralB.c5 = deepcopy(dihedral_objB.c5)


    return topB



# # pairs
#             if name == "pairs":
#                 if not csplit[:2] == list(CR[0].atom_idx):
#                     clist.append(csplit)
#         else:
#             clist.append(csplit)
#     content = clist
#     return content
