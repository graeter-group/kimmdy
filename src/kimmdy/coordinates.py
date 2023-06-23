import MDAnalysis as MDA
from pathlib import Path
import numpy as np
import logging
from typing import Union

from kimmdy.parsing import read_topol
from kimmdy.tasks import TaskFiles
from kimmdy.topology.topology import Topology
from kimmdy.topology.atomic import Bond, Angle, Dihedral, Pair, Atomic
from kimmdy.topology.utils import match_atomic_item_to_atomic_type


## copied from changemanager. should be put into utils
def get_ff_sections(ffdir):
    return read_topol(ffdir / "ffbonded.itp")


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


##


def merge_top_prmgrowth(
    files: TaskFiles, focus_nr: Union[list[str], None] = None
) -> Topology:
    # ffdir = (files.input["ff"],)
    # not the most robust way to get topA and topB
    hyperprms = {
        "morse_well_depth": "400",
        "morse_steepness": "10",
    }  # well_depth D [kJ/mol], steepness [nm-1]
    topADict = read_topol(files.input["top"])
    topBDict = read_topol(files.output["top"])
    topA = Topology(topADict)
    topB = Topology(topBDict)

    # ToDo: what about implicit parameters?? especially dihedrals
    # think about how to bring focus_nr into this
    for nr in topB.atoms.keys():

        # atoms
        atomA = topA.atoms[nr]
        atomB = topB.atoms[nr]
        if atomA != atomB:
            if atomA.charge != atomB.charge:
                atomB.chargeB = atomB.charge
                atomB.charge = atomA.charge
            else:
                logging.debug(
                    f"Atom {nr} changed during changemanager step but not the charges!"
                )

        # # bonds
        bondA_keys = set(topA.bonds.keys())
        bondB_keys = set(topB.bonds.keys())

        same = set.intersection(bondA_keys, bondB_keys)
        broken = bondA_keys - bondB_keys
        bound = bondB_keys - bondA_keys

        for bond_key in same:
            bondA = topA.bonds.get(bond_key)
            bondB = topB.bonds.get(bond_key)
            if bondA != bondB:
                # assuming no bond has explicit standard ff parameters
                bond_objA = get_bondobj(bond_key, bondA, topA)
                bond_objB = get_bondobj(bond_key, bondB, topB)
                bondB.c2 = bond_objB.c0
                bondB.c3 = bond_objB.c3
                bondB.c0 = bond_objA.c0
                bondB.c1 = bond_objA.c1

        for bond_key in broken:
            topB.bind_bond(bond_key)
            bondB = topB.bonds.get(bond_key)
            bond_objA = get_bondobj(bond_key, bondA, topA)

            bondB.funct = "3"  # Morse potential
            bondB.c0 = bond_objA.c0
            bondB.c1 = hyperprms["morse_well_depth"]
            bondB.c2 = hyperprms["morse_steepness"]
            bondB.c3 = "0.00"
            # bondB.c4 = '0.00'
            # bondB.c5 = '0.00'

        for bond_key in bound:
            bondB = topB.bonds.get(bond_key)
            bond_objB = get_bondobj(bond_key, bondB, topB)

            bondB.funct = "3"  # Morse potential
            bondB.c0 = "0.00"
            bondB.c1 = "0.00"
            bondB.c2 = "0.00"
            bondB.c3 = bond_objA.c0
            # bondB.c4 = hyperprms['morse_well_depth']
            # bondB.c5 = hyperprms['morse_steepness']
    return topB
    # # angles

    # angle_keys = topB._get_atom_angles(focus_nr[0]) + topB._get_atom_angles(
    #     focus_nr[1]
    # )
    # for key in angle_keys:
    #     angle = topB.angles.get(key)
    #     if (
    #         angle is None
    #         or topB.ffpatches is None
    #         or topB.ffpatches.anglepatches is None
    #     ):
    #         continue
    #     id = [topB.atoms[i].radical_type() for i in key]
    #     patch = match_id_to_patch(id, topB.ffpatches.anglepatches)
    #     if patch is None:
    #         continue
    #     id_base = [topB.atoms[i].radical_type() for i in key]
    #     topB._apply_param_patch(angle, id_base, patch, topB.ff.angletypes)

    # # proper dihedrals and pairs
    # dihedral_keys = self._get_atom_proper_dihedrals(
    #     focus_nr[0]
    # ) + self._get_atom_proper_dihedrals(focus_nr[1])
    # for key in dihedral_keys:
    #     dihedral = self.proper_dihedrals.get(key)
    #     if (
    #         dihedral is None
    #         or self.ffpatches is None
    #         or self.ffpatches.anglepatches is None
    #     ):
    #         continue
    #     id = [self.atoms[i].radical_type() for i in key]
    #     patch = match_id_to_patch(id, self.ffpatches.dihedralpatches)
    #     if patch is None:
    #         continue
    #     id_base = [self.atoms[i].radical_type() for i in key]
    #     self._apply_param_patch(
    #         dihedral, id_base, patch, self.ff.proper_dihedraltypes


# get toppath_A from runmanager iterating through self.filehist['n'] to find self.filehist['x']['_run_recipe']['output']['top']
# def merge_section_slowgrowth(
#     name: str,
#     content: list[str],
#     CR: ConversionRecipe,
#     state_A_reduced: Topology,
#     ffpath: Path,
# ):
#     # name is 'bonds' or 'angles'
#     logging.debug("CR", CR)
#     clist = []
#     recipeatoms = [x for y in CR for x in y.atom_idx]
#     recipeatoms = list(set(recipeatoms))

#     # holds true for stateB as well except for the HAT atom
#     atom_idx_at = {
#         key: value
#         for (key, value) in map(lambda a: [a[0], a[1]], state_A_reduced["atoms"])
#     }

#     inserted_angles = False

#     ffprm = get_ff_sections(ffpath)
#     for c in content:
#         # TODO: make merge work
#         csplit = c.split()
#         if any([idx in csplit for idx in recipeatoms]):
#             # move check up
# # bonds
#             if name == "bonds":
#                 atoms = csplit[0:2]
#                 logging.debug(atoms)
#                 # atom_idx of bind (exists only in stateB)
#                 if atoms == list(CR[1].atom_idx):
#                     CR_bonds = {
#                         "break": [
#                             x
#                             for x in state_A_reduced["bonds"]
#                             if CR[0].atom_idx == tuple(x[:2])
#                         ][0],
#                         "bind": csplit,
#                     }
#                     for key, bond in CR_bonds.items():
#                         if not is_parameterized(bond):
#                             bond_at = [atom_idx_at[x] for x in bond[:2]]
#                             CR_bonds[key] = parameterize_bonded_terms(
#                                 ffprm, [bond_at], "bondtypes", [bond[:2]]
#                             )
#                             # unpack
#                             CR_bonds[key] = CR_bonds[key][0]
#                     CR_bonds["break"] = [
#                         *CR_bonds["break"][:2],
#                         "3",
#                         CR_bonds["break"][3],
#                         "400.0",
#                         "10.0",
#                         "0.00",
#                         "0.00",
#                         "0.00",
#                     ]
#                     CR_bonds["bind"] = [
#                         *CR_bonds["bind"][:2],
#                         "3",
#                         "0.00",
#                         "0.00",
#                         "0.00",
#                         CR_bonds["bind"][3],
#                         "400.0",
#                         "10.0",
#                     ]
#                     clist.append(CR_bonds["bind"])
#                     clist.append(CR_bonds["break"])
#                 # bonded terms that are both in stateA and stateB
#                 else:
#                     [bond_A] = [x for x in state_A_reduced["bonds"] if atoms == x[:2]]
#                     bond_B = csplit
#                     bonds = [bond_A, bond_B]
#                     logging.debug(bonds)

#                     for i, bond in enumerate(bonds):
#                         if not is_parameterized(bond):
#                             bond_at = [atom_idx_at[x] for x in bond[:2]]
#                             [bonds[i]] = parameterize_bonded_terms(
#                                 ffprm, [bond_at], "bondtypes", [bond[:2]]
#                             )
#                             logging.debug("!!", bond)
#                     logging.debug(bonds)
#                     if not bonds[0] == bonds[1]:
#                         bond_slowgrowth = [*bonds[0][:5], *bonds[1][3:5]]
#                     else:
#                         bond_slowgrowth = bonds[1]
#                     logging.debug(bond_slowgrowth)
#                     clist.append(bond_slowgrowth)
#                 # atom_A = [x if idxs[2:4] in x else None for x in state_A_reduced["bonds"]]
#             # logging.debug(atoms,atom_A)

# #angles
#             if name == "angles":
#                 atoms = csplit[:3]
#                 if atoms[1] == CR[0].atom_idx[0]:
#                     # potentially broken angle -> not in state B
#                     if not inserted_angles:
#                         # definitively broken
#                         addterms = [
#                             x
#                             for x in state_A_reduced["angles"]
#                             if all([y in x[:4] for y in CR[0].atom_idx])
#                         ]
#                         logging.debug(CR[0].atom_idx, addterms)
#                         for i, term in enumerate(addterms):
#                             if not is_parameterized(term):
#                                 term_at = [atom_idx_at[x] for x in term[:3]]
#                                 [term] = parameterize_bonded_terms(
#                                     ffprm, [term_at], "angletypes", [term[:3]]
#                                 )
#                                 logging.debug(term)
#                                 addterms[i] = [*term[:6], term[4], "0.00"]
#                             else:
#                                 addterms[i] = [*term[:6], term[4], "0.00"]
#                         logging.debug(addterms)
#                         clist.extend(addterms)
#                         inserted_angles = True

#                     # angle should exist in A
#                     [angle_A] = [x for x in state_A_reduced["angles"] if atoms == x[:3]]
#                     angle_B = csplit
#                     if angle_A == angle_B:
#                         clist.append(csplit)
#                     elif is_parameterized(angle_A) != is_parameterized(angle_B):
#                         if is_parameterized(angle_A):
#                             angle_B_at = [atom_idx_at[x] for x in angle_B[:3]]
#                             [angle_B] = parameterize_bonded_terms(
#                                 ffprm, [angle_B_at], "angletypes", [angle_B[:3]]
#                             )
#                         elif is_parameterized(angle_B):
#                             angle_A_at = [atom_idx_at[x] for x in angle_A[:3]]
#                             [angle_A] = parameterize_bonded_terms(
#                                 ffprm, [angle_A_at], "angletypes", [angle_A[:3]]
#                             )
#                         if angle_A == angle_B:
#                             clist.append(csplit)
#                         else:
#                             angle_slowgrowth = [*angle_A[:6], *angle_B[4:6]]
#                             clist.append(angle_slowgrowth)
#                     else:
#                         angle_slowgrowth = [*angle_A[:6], *angle_B[4:6]]
#                         clist.append(angle_slowgrowth)

#                 elif atoms[1] == CR[1].atom_idx[0]:
#                     # this line, among others assumes a difference between atom_idx[0] and atom_idx[1] which is the case for HAT but not other reactions
#                     angle_B = csplit
#                     if CR[1].atom_idx[1] in atoms:
#                         # doesn't exist in state_A
#                         if not is_parameterized(angle_B):
#                             angle_B_at = [atom_idx_at[x] for x in angle_B[:3]]
#                             [angle_B] = parameterize_bonded_terms(
#                                 ffprm, [angle_B_at], "angletypes", [angle_B[:3]]
#                             )
#                         angle_slowgrowth = [*angle_B[:5], "0.00", *angle_B[4:6]]
#                         clist.append(angle_slowgrowth)
#                     else:
#                         # copied
#                         # angle should exist in A
#                         [angle_A] = [
#                             x for x in state_A_reduced["angles"] if atoms == x[:3]
#                         ]
#                         angle_B = csplit
#                         if angle_A == angle_B:
#                             clist.append(csplit)
#                         elif is_parameterized(angle_A) != is_parameterized(angle_B):
#                             if is_parameterized(angle_A):
#                                 angle_B_at = [atom_idx_at[x] for x in angle_B[:3]]
#                                 [angle_B] = parameterize_bonded_terms(
#                                     ffprm, [angle_B_at], "angletypes", [angle_B[:3]]
#                                 )
#                             elif is_parameterized(angle_B):
#                                 angle_A_at = [atom_idx_at[x] for x in angle_A[:3]]
#                                 [angle_A] = parameterize_bonded_terms(
#                                     ffprm, [angle_A_at], "angletypes", [angle_A[:3]]
#                                 )
#                             if angle_A == angle_B:
#                                 clist.append(csplit)
#                             else:
#                                 angle_slowgrowth = [*angle_A[:6], *angle_B[4:6]]
#                                 clist.append(angle_slowgrowth)
#                         else:
#                             angle_slowgrowth = [*angle_A[:6], *angle_B[4:6]]
#                             clist.append(angle_slowgrowth)
#                 else:
#                     clist.append(csplit)
# # pairs
#             if name == "pairs":
#                 if not csplit[:2] == list(CR[0].atom_idx):
#                     clist.append(csplit)
#         else:
#             clist.append(csplit)
#     content = clist
#     return content
