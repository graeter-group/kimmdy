import MDAnalysis as MDA
from pathlib import Path
import numpy as np
import logging

from kimmdy.utils import find_radical_pos
from kimmdy.parsing import read_topol, Topology
from kimmdy.reaction import ReactionResult, ConversionRecipe, ConversionType


def place_hydrogen(tpr: Path, trr: Path, recipe: ConversionRecipe):
    """use find_radical_pos to change the coordinates of the HAT hydrogen according to the reaction recipe"""
    # current recipe format: ConversionRecipe(type=[ConversionType.MOVE], atom_idx=[[from_H_nr, rad_nr]])                   )
    u = MDA.Universe(str(tpr), str(trr), topology_format="tpr", format="trr")
    # set timestep to last one in trajectory
    u.trajectory[-1]
    rad = u.select_atoms(f"bynum {recipe['atom_idx'][1]}")[0]
    bonded_rad = rad.bonded_atoms
    H_options = find_radical_pos(rad, bonded_rad)

    H_initial = u.select_atoms(f"bynum {recipe['atom_idx'][0]}")[0]

    # dealing with cases where several positions are returned as H_options
    # continues with the positions thats closest to the initial H_initial position
    distances = [np.sum(np.absolute(x - H_initial.position)) for x in H_options]
    H_final = H_options[np.argmin(distances)]

    system = u.select_atoms("all")
    H_initial.position = H_final
    system.write(trr.with_name(f"{trr.stem}_mod.trr"))
    return


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


def is_parameterized(term):
    """does not work for dihedrals
    basically checks whether there are floats in the term
    """
    for idx in term:
        if not idx.isdigit():
            if idx.replace(".", "", 1).isdigit():
                return True
    return False


##

# get toppath_A from runmanager iterating through self.filehist['n'] to find self.filehist['x']['_run_recipe']['output']['top']
def merge_section_slowgrowth(
    name: str,
    content: list[str],
    CR: ConversionRecipe,
    state_A_reduced: Topology,
    ffpath: Path,
):
    # name is 'bonds' or 'angles'
    logging.debug("CR", CR)
    clist = []
    recipeatoms = [x for y in CR for x in y.atom_idx]
    recipeatoms = list(set(recipeatoms))

    # holds true for stateB as well except for the HAT atom
    atom_idx_at = {
        key: value
        for (key, value) in map(lambda a: [a[0], a[1]], state_A_reduced["atoms"])
    }

    inserted_angles = False

    ffprm = get_ff_sections(ffpath)
    for c in content:
        # TODO: make merge work
        csplit = c.split()
        if any([idx in csplit for idx in recipeatoms]):
            # move check up
            if name == "bonds":
                atoms = csplit[0:2]
                logging.debug(atoms)
                # atom_idx of bind (exists only in stateB)
                if atoms == list(CR[1].atom_idx):
                    CR_bonds = {
                        "break": [
                            x
                            for x in state_A_reduced["bonds"]
                            if CR[0].atom_idx == tuple(x[:2])
                        ][0],
                        "bind": csplit,
                    }
                    for key, bond in CR_bonds.items():
                        if not is_parameterized(bond):
                            bond_at = [atom_idx_at[x] for x in bond[:2]]
                            CR_bonds[key] = parameterize_bonded_terms(
                                ffprm, [bond_at], "bondtypes", [bond[:2]]
                            )
                            # unpack
                            CR_bonds[key] = CR_bonds[key][0]
                    CR_bonds["break"] = [
                        *CR_bonds["break"][:2],
                        "3",
                        CR_bonds["break"][3],
                        "400.0",
                        "10.0",
                        "0.00",
                        "0.00",
                        "0.00",
                    ]
                    CR_bonds["bind"] = [
                        *CR_bonds["bind"][:2],
                        "3",
                        "0.00",
                        "0.00",
                        "0.00",
                        CR_bonds["bind"][3],
                        "400.0",
                        "10.0",
                    ]
                    clist.append(CR_bonds["bind"])
                    clist.append(CR_bonds["break"])
                # bonded terms that are both in stateA and stateB
                else:
                    [bond_A] = [x for x in state_A_reduced["bonds"] if atoms == x[:2]]
                    bond_B = csplit
                    bonds = [bond_A, bond_B]
                    logging.debug(bonds)

                    for i, bond in enumerate(bonds):
                        if not is_parameterized(bond):
                            bond_at = [atom_idx_at[x] for x in bond[:2]]
                            [bonds[i]] = parameterize_bonded_terms(
                                ffprm, [bond_at], "bondtypes", [bond[:2]]
                            )
                            logging.debug("!!", bond)
                    logging.debug(bonds)
                    if not bonds[0] == bonds[1]:
                        bond_slowgrowth = [*bonds[0][:5], *bonds[1][3:5]]
                    else:
                        bond_slowgrowth = bonds[1]
                    logging.debug(bond_slowgrowth)
                    clist.append(bond_slowgrowth)
                # atom_A = [x if idxs[2:4] in x else None for x in state_A_reduced["bonds"]]
            # logging.debug(atoms,atom_A)
            if name == "angles":
                atoms = csplit[:3]
                if atoms[1] == CR[0].atom_idx[0]:
                    # potentially broken angle -> not in state B
                    if not inserted_angles:
                        # definitively broken
                        addterms = [
                            x
                            for x in state_A_reduced["angles"]
                            if all([y in x[:4] for y in CR[0].atom_idx])
                        ]
                        logging.debug(CR[0].atom_idx, addterms)
                        for i, term in enumerate(addterms):
                            if not is_parameterized(term):
                                term_at = [atom_idx_at[x] for x in term[:3]]
                                [term] = parameterize_bonded_terms(
                                    ffprm, [term_at], "angletypes", [term[:3]]
                                )
                                logging.debug(term)
                                addterms[i] = [*term[:6], term[4], "0.00"]
                            else:
                                addterms[i] = [*term[:6], term[4], "0.00"]
                        logging.debug(addterms)
                        clist.extend(addterms)
                        inserted_angles = True

                    # angle should exist in A
                    [angle_A] = [x for x in state_A_reduced["angles"] if atoms == x[:3]]
                    angle_B = csplit
                    if angle_A == angle_B:
                        clist.append(csplit)
                    elif is_parameterized(angle_A) != is_parameterized(angle_B):
                        if is_parameterized(angle_A):
                            angle_B_at = [atom_idx_at[x] for x in angle_B[:3]]
                            [angle_B] = parameterize_bonded_terms(
                                ffprm, [angle_B_at], "angletypes", [angle_B[:3]]
                            )
                        elif is_parameterized(angle_B):
                            angle_A_at = [atom_idx_at[x] for x in angle_A[:3]]
                            [angle_A] = parameterize_bonded_terms(
                                ffprm, [angle_A_at], "angletypes", [angle_A[:3]]
                            )
                        if angle_A == angle_B:
                            clist.append(csplit)
                        else:
                            angle_slowgrowth = [*angle_A[:6], *angle_B[4:6]]
                            clist.append(angle_slowgrowth)
                    else:
                        angle_slowgrowth = [*angle_A[:6], *angle_B[4:6]]
                        clist.append(angle_slowgrowth)

                elif atoms[1] == CR[1].atom_idx[0]:
                    # this line, among others assumes a difference between atom_idx[0] and atom_idx[1] which is the case for HAT but not other reactions
                    angle_B = csplit
                    if CR[1].atom_idx[1] in atoms:
                        # doesn't exist in state_A
                        if not is_parameterized(angle_B):
                            angle_B_at = [atom_idx_at[x] for x in angle_B[:3]]
                            [angle_B] = parameterize_bonded_terms(
                                ffprm, [angle_B_at], "angletypes", [angle_B[:3]]
                            )
                        angle_slowgrowth = [*angle_B[:5], "0.00", *angle_B[4:6]]
                        clist.append(angle_slowgrowth)
                    else:
                        # copied
                        # angle should exist in A
                        [angle_A] = [
                            x for x in state_A_reduced["angles"] if atoms == x[:3]
                        ]
                        angle_B = csplit
                        if angle_A == angle_B:
                            clist.append(csplit)
                        elif is_parameterized(angle_A) != is_parameterized(angle_B):
                            if is_parameterized(angle_A):
                                angle_B_at = [atom_idx_at[x] for x in angle_B[:3]]
                                [angle_B] = parameterize_bonded_terms(
                                    ffprm, [angle_B_at], "angletypes", [angle_B[:3]]
                                )
                            elif is_parameterized(angle_B):
                                angle_A_at = [atom_idx_at[x] for x in angle_A[:3]]
                                [angle_A] = parameterize_bonded_terms(
                                    ffprm, [angle_A_at], "angletypes", [angle_A[:3]]
                                )
                            if angle_A == angle_B:
                                clist.append(csplit)
                            else:
                                angle_slowgrowth = [*angle_A[:6], *angle_B[4:6]]
                                clist.append(angle_slowgrowth)
                        else:
                            angle_slowgrowth = [*angle_A[:6], *angle_B[4:6]]
                            clist.append(angle_slowgrowth)
                else:
                    clist.append(csplit)

            if name == "pairs":
                if not csplit[:2] == list(CR[0].atom_idx):
                    clist.append(csplit)
        else:
            clist.append(csplit)
    content = clist
    return content
