import logging
from kimmdy.reaction import ConversionRecipe, ConversionType
from kimmdy.parsing import (
    read_plumed,
    write_plumed,
    read_topol,
    write_topol,
    topol_split_dihedrals,
    topol_merge_propers_impropers,
    Topology,
)
from kimmdy.utils import (
    str_to_int,check_idx,sort_bond,sort_angle,sort_dihedral,sort_improper
)
from pathlib import Path

from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from enum import Enum, auto
from collections.abc import Iterable
from typing import Generator
from copy import deepcopy
import itertools
from operator import itemgetter
import re


def modify_top(recipe: ConversionRecipe, oldtop: Path, newtop: Path, ffdir: Path):
    logging.info(f"Reading: {oldtop} and writing modified topology to {newtop}.")
    topology = read_topol(oldtop)

    for type, pair in zip(recipe.type, recipe.atom_idx):
        if type == ConversionType.BREAK:
            topology = break_bond_top(topology, pair)
        elif type == ConversionType.MOVE:
            topology = move_bond_top(topology, pair, ffdir)

    write_topol(topology, newtop)


def break_bond_top(topology: Topology, breakpair: tuple[int, int]) -> Topology:
    """Break bonds in topology.
    removes bond, angles and dihedrals where breakpair was involved
    """
    breakpair = (str(breakpair[0]),str(breakpair[1]))
    topology["bonds"] = [
        bond
        for bond in topology["bonds"]
        if not all(x in [bond[0], bond[1]] for x in breakpair)
    ]
    topology["angles"] = [
        angle
        for angle in topology["angles"]
        if not all(x in [angle[0], angle[1], angle[2]] for x in breakpair)
    ]
    topology["dihedrals"] = [
        dihedral
        for dihedral in topology["dihedrals"]
        if not all(
            x in [dihedral[0], dihedral[1], dihedral[2], dihedral[3]] for x in breakpair
        )
    ]

    header = topology["pairs"][0]
    dihpairs = [[d[0], d[3]] for d in topology["dihedrals"]]
    for pair in dihpairs:  # keeping pairs in proper order (which one goes first in the definition)
        if str_to_int(pair[0]) != 0 and str_to_int(pair[1]) != 0:
            pair_int = [int(pair[0]), int(pair[1])]
            if pair_int[0] > pair_int[1]:
                pair.reverse()
    dihpairs = sorted(dihpairs, key=sort_bond)  
    
    topology["pairs"] = [pair for pair in topology["pairs"] if pair[:2] in dihpairs]
    topology["pairs"].insert(0,header)

    return topology


def move_bond_top(
    topology: Topology, movepair: tuple[int, int], ffdir: Path
) -> Topology:

    topology = topol_split_dihedrals(topology)
    heavy_idx = find_heavy(topology["bonds"], movepair[0])
    logging.debug(f"Heavy atom bound to HAT hydrogen has idx {heavy_idx}")

    logging.info(f"---- Dealing with the 'from' part ----")

    # build localGraph and fill it 
    fromGraph = localGraph(topology, heavy_idx, ffdir)
    fromGraph.construct_graph()
    fromGraph.order_lists()
    fromGraph.update_atoms_list()
    fromGraph.update_bound_to()
    fromGraph.build_PADs()
    fromGraph.order_lists()

    # remove terms with from_H 
    termdict = fromGraph.get_terms_with_atom(movepair[0])
    fromGraph.remove_terms(termdict)
    topology = topol_remove_terms(topology, termdict)  

    # parameterize around the new radical from_heavy
    atom_terms = fromGraph.parameterize_around_atom(heavy_idx)  
    topology = topol_add_terms(topology, atom_terms)  


    logging.info(f"---- Dealing with the 'to' part ----\n")

    # build localGraph and fill it 
    toGraph = localGraph(topology, movepair[1], ffdir)
    toGraph.construct_graph()
    toGraph.order_lists()
    toGraph.update_atoms_list()
    toGraph.update_bound_to()
    toGraph.add_Atom(Atom(movepair[0]))
    toGraph.order_lists()
    toGraph.update_atoms_list()
    toGraph.add_Bond(deepcopy(movepair))
    toGraph.update_bound_to()
    toGraph.order_lists()
    toGraph.build_PADs()
    toGraph.order_lists()

    # set the correct atomtype,resname for the HAT hydrogen based on the ff definition
    atoms_idxs_from = fromGraph.atoms_idx
    heavy_resname = fromGraph.AtomList[atoms_idxs_from.index(heavy_idx)].resname  
    atmdef = toGraph.correct_atomprops(movepair,heavy_resname)
    topol_change_at_an(topology,atmdef) 

    # this is to stop overwriting bonds,angles for the newly parameterized "from" part if they are next to each other
    atom_terms = toGraph.parameterize_around_atom(movepair[1])
    atom_terms = terms_keep_only(movepair)
    topology = topol_add_terms(topology, atom_terms)

    # add pairs of the from_H at the new position
    atom_terms_H = toGraph.get_terms_with_atom(movepair[0], add_function=True)
    for section in ["bonds", "angles", "propers", "impropers"]:
        atom_terms_H[section].clear()
    topology = topol_add_terms(topology, atom_terms_H)

    # have the right impropers at the to_heavy atom
    atom_terms_add, atom_terms_remove = toGraph.compare_ff_impropers(heavy_idx, movepair[1])
    topology = topol_add_terms(topology, atom_terms_add)  # no need, yet
    topology = topol_remove_terms(topology, atom_terms_remove)

    topology = topol_merge_propers_impropers(topology)
    return topology


class Atom:
    """
    A class containing atom information as in the atoms section of the topology
    and a bound_to variable
    """

    def __init__(self, idx: str, atomtype=None, atomname=None, resname=None):
        self.idx = idx
        self.atomtype = atomtype
        self.atomname = atomname
        self.resname = resname
        self.bound_to = []


class localGraph:
    def __init__(self, topology, heavy_idx, ffdir):
        self.topology = topology
        self.heavy_idx = heavy_idx
        self.radicals = []
        self.ffdir = ffdir
        self.ff = {}

        self.AtomList = [Atom(self.heavy_idx)]
        self.BondList = []
        self.PairList = []
        self.AngleList = []
        self.ProperList = []
        self.ImproperList = []
        self.atoms_idx = [self.heavy_idx]
        self.atoms_atomtype = []
        self.atoms_atomname = []
        self.atoms_resname = []

## Propagate changes to all attributes 
    def update_atoms_list(self):
        self.atoms_idx = self.get_AtomList_property("idx")
        self.get_atomprop("atomtype")
        self.get_atomprop("atomname")
        self.get_atomprop("resname")
        self.atoms_atomtype = self.get_AtomList_property("atomtype")
        self.atoms_atomname = self.get_AtomList_property("atomname")
        self.atoms_resname = self.get_AtomList_property("resname")

    def update_bound_to(self):
        bound_dict = {}
        for atom_idx in self.atoms_idx:
            bound_dict[atom_idx] = []
        for bond in self.BondList:
            if not bond[1] in bound_dict[bond[0]]:
                bound_dict[bond[0]].append(bond[1])
            if not bond[0] in bound_dict[bond[1]]:
                bound_dict[bond[1]].append(bond[0])
        for atom_idx in self.atoms_idx:
            self.AtomList[self.atoms_idx.index(atom_idx)].bound_to = sorted(
                bound_dict[atom_idx], key=str_to_int
            )

    def order_lists(self):
        self.AtomList = sorted(self.AtomList, key=check_idx)
        self.atoms_idx = sorted(self.atoms_idx, key=str_to_int)
        self.atoms_atomtype = [
            x
            for _, x in sorted(zip(self.atoms_idx, self.atoms_atomtype), key=sort_bond)
        ]  # not so great solution
        self.atoms_atomname = [
            x
            for _, x in sorted(zip(self.atoms_idx, self.atoms_atomname), key=sort_bond)
        ]
        self.atoms_resname = [
            x for _, x in sorted(zip(self.atoms_idx, self.atoms_resname), key=sort_bond)
        ]
        self.BondList = sorted(self.BondList, key=sort_bond)
        self.PairList = sorted(self.PairList, key=sort_bond)
        self.AngleList = sorted(self.AngleList, key=sort_angle)
        self.ProperList = sorted(self.ProperList, key=sort_dihedral)
        self.ImproperList = sorted(self.ImproperList, key=sort_improper)

        for (
            bond
        ) in (
            self.BondList
        ):  # keeping bonds in proper order (which one goes first in the definition)
            bond_int = [int(bond[0]), int(bond[1])]
            if bond_int[0] > bond_int[1]:
                bond.reverse()

## Retrieve properties
    def get_atomprop(self, property):
        for entry in self.topology["atoms"]:
            if entry[0] in self.atoms_idx:
                if property == "atomtype":
                    self.AtomList[self.atoms_idx.index(entry[0])].atomtype = entry[1]
                elif property == "atomname":
                    self.AtomList[self.atoms_idx.index(entry[0])].atomname = entry[4]
                elif property == "resname":
                    self.AtomList[self.atoms_idx.index(entry[0])].resname = entry[3]

    def get_AtomList_property(self, property):
        if property == "idx":
            return [x.idx for x in self.AtomList]
        elif property == "bound_to":
            return [x.bound_to for x in self.AtomList]
        elif property == "atomtype":
            return [x.atomtype for x in self.AtomList]
        elif property == "atomname":
            return [x.atomname for x in self.AtomList]
        elif property == "resname":
            return [x.resname for x in self.AtomList]
        else:
            raise ValueError("Can't find property {} in object".format(property))

## Methods related to filling the LocalGraph attributes from the topology dict
    def construct_graph(self, depth=3):
        """
        searches the bonds section of a topology up to depth bonds deep
        to create a local graph i.e. filling the AtomList and BondList
        """
        curratoms = [self.heavy_idx]
        for i in range(depth):
            addlist = []
            for bond in self.topology["bonds"]:
                if any(x in bond[:2] for x in curratoms):
                    if not bond[:2] in self.BondList:
                        self.BondList.append(bond[:2])
                    for j in [0, 1]:
                        if not bond[j] in [x.idx for x in self.AtomList]:
                            addlist.append(bond[j])
                            self.AtomList.append(Atom(bond[j]))
                            self.atoms_idx.append(bond[j])
            curratoms = deepcopy(addlist)

    def build_PADs(self):
        """
        in:
        - localGraph with AtomList that contains atom objects with bound_to attribute that contains the bonded atom idxs
        out:
        - filled AngleList, ProperList 
        - filled PairList as ProperList outer atoms
        - ImproperList as defined in the topology
        """
        # defining angles
        for atom in self.AtomList:
            if len(atom.bound_to) >= 2:
                for comb in itertools.combinations(atom.bound_to, 2):
                    self.AngleList.append([comb[0], atom.idx, comb[1]])

        # defining propers
        for atom1 in self.AtomList:
            if len(atom1.bound_to) >= 2:
                for partner in atom1.bound_to:
                    atom2 = self.AtomList[self.atoms_idx.index(partner)]
                    if len(atom2.bound_to) >= 2 and int(atom2.idx) > int(atom1.idx):
                        atom1_periphery = [x for x in atom1.bound_to if x != atom2.idx]
                        atom2_periphery = [x for x in atom2.bound_to if x != atom1.idx]
                        for prod in itertools.product(atom1_periphery, atom2_periphery):
                            self.ProperList.append(
                                [prod[0], atom1.idx, atom2.idx, prod[1]]
                            )

        # taking impropers from topology
        for term in self.topology["impropers"]:
            if all([x in self.atoms_idx for x in term[:4]]):
                self.ImproperList.append(term)

        # all pairs have a dihedral between them
        for dihedral in self.ProperList:  # this looks stupid, needs to check for existing pairs in the case of five membered rings           
            if int(dihedral[0]) < int(dihedral[3]):                
                self.PairList.append([dihedral[0], dihedral[3]])
            else:
                self.PairList.append([dihedral[3], dihedral[0]])

    def correct_atomprops(self, movepair, heavy_resname):
        atoms_idxs = self.atoms_idx
        H_atomtype, H_atomname = self.get_H_ff_at_an(movepair[1])
        self.AtomList[atoms_idxs.index(movepair[0])].atomtype = H_atomtype  
        self.AtomList[atoms_idxs.index(movepair[0])].atomname = H_atomname
        self.AtomList[atoms_idxs.index(movepair[0])].resname = heavy_resname
        self.update_atoms_list()
        return [movepair[0],H_atomtype,H_atomname,self.AtomList[atoms_idxs.index(movepair[0])].resname]

## Methods related to retrieving ff terms from a LocalGraph that contain a certain atom (search for index)
    def search_terms(self, TypeList, idx: str, center):
        if TypeList == []:
            return []
        CopyList = deepcopy(TypeList)
        hits = []
        if center and len(TypeList[0]) > 2:
            start = 1
            end = -1
        else:
            start = 0
            end = None
        for entry in CopyList:
            if idx in entry[start:end]:
                hits.append(entry)
        return hits

    def add_function_to_terms(self, termdict):
        termdict_funct = deepcopy(termdict)
        funct_dict = {
            "bonds": "1",
            "pairs": "1",
            "angles": "1",
            "propers": "9",
            "impropers": "4",
        }
        for section in ["bonds", "pairs", "angles", "propers", "impropers"]:
            termdict_funct[section] = [
                [*x, funct_dict[section]] for x in termdict_funct[section]
            ]
        return termdict_funct

    def get_terms_with_atom(self, idx: str, center=False, add_function=False):
        """
        in: 
        - idx:str           = atom index 
        - center:bool       = take angle/dihedral terms with idx in center positions only
        - add_function:bool = add GROMACS function # to the terms
        out:
        - termdict:dict     = bonds,pairs,angles,propers,impropers containing atom idx
        """
        if not idx in self.atoms_idx:
            return

        termdict = {}
        termdict["bonds"] = self.search_terms(self.BondList, idx, center)
        termdict["pairs"] = self.search_terms(self.PairList, idx, center)
        termdict["angles"] = self.search_terms(self.AngleList, idx, center)
        termdict["propers"] = self.search_terms(self.ProperList, idx, center)
        termdict["impropers"] = self.search_terms(self.ImproperList, idx, center)

        if add_function:
            termdict_funct = self.add_function_to_terms(termdict)
            return termdict_funct
        return termdict

## Manipulate attributes
    def add_Atom(self, atom: Atom):
        if not atom in self.AtomList:
            self.AtomList.append(atom)
            self.atoms_idx.append(atom.idx)
        else:
            logging.warning(f"Atom with idx {Atom.idx} already exists!")

    def add_Bond(self, bond: list):
        if all(x in self.atoms_idx for x in bond):
            self.BondList.append(bond)
        else:
            logging.warning("Could not establish bond between idxs {bond}!")

    def remove_terms(self, termdict):
        logging.debug(f"Removing terms {termdict} from local graph.")
        graphdict = {
            "bonds": self.BondList,
            "pairs": self.PairList,
            "angles": self.AngleList,
            "propers": self.ProperList,
            "impropers": self.ImproperList,
        }
        for section in ["pairs", "bonds", "angles", "propers", "impropers"]:
            try:
                for term in termdict[section]:
                    graphdict[section].remove(term)
                    logging.debug(f"removed term {term} from graph {section}")
            except ValueError:
                logging.warning(f"Couldn't find term {term} in TypeList")

        self.update_bound_to()

#obsolete?
    def compare_section(self, section, TypeList):
        """
        in:
        - topology section
        - TypeList (e.g, BondList, PairList)
        out:
        - ff terms that are in the TypeList but not in the topology section
        """
        typelen = len(TypeList[0])
        missing_terms = deepcopy(TypeList)
        for entry in section:
            if entry[:typelen] in missing_terms:
                missing_terms.remove(entry[:typelen])
        return missing_terms
## obsolete?
    def find_missing_terms(self):
        """
        returns ff terms that are not part of the topology dict but in the TypeLists 
        """
        termdict = {}
        termdict["bonds"] = self.compare_section(self.topology["bonds"], self.BondList)
        termdict["pairs"] = self.compare_section(self.topology["pairs"], self.PairList)
        termdict["angles"] = self.compare_section(self.topology["angles"], self.AngleList)
        termdict["propers"] = self.compare_section(self.topology["propers"], self.ProperList)
        return termdict


## Method that interact with the force field files 
    def is_radical(self, atom_idx):
        """
        in:
        - atom_idx              = index of the atom that is checked for being a radical
        - self.atoms_atomtype   = atomtypes in the LocalGraph (comes from the topology dict)
        - nbonds_dict           = # of bonds of amber FF99SB-ILDNP* atomtypes (as written in the paper)
        out:
        - bool whether atom_idx is a radical (-> has fewer bonds than would be expected for the atomtype)
        """

        nbonds_dict = {
            ("MG", "NA", "CO"): 0,
            (
                "H",
                "HW",
                "HO",
                "HS",
                "HA",
                "HC",
                "H1",
                "H2",
                "H3",
                "HP",
                "H4",
                "H5",
                "HO",
                "H0",
                "HP",
                "O",
                "O2",
                "Cl",
                "Na",
                "I",
                "F",
                "Br",
            ): 1,
            ("NB", "NC", "OW", "OH", "OS", "SH", "S"): 2,
            (
                "C",
                "CN",
                "CB",
                "CR",
                "CK",
                "CC",
                "CW",
                "CV",
                "C*",
                "CQ",
                "CM",
                "CA",
                "CD",
                "CZ",
                "N",
                "NA",
                "N*",
                "N2",
            ): 3,
            ("CT", "N3", "P", "SO"): 4,
        }  # compare to atom type perception paper (2006) same as in HAT_utils.py
        atom_type = self.atoms_atomtype[self.atoms_idx.index(atom_idx)]
        try:
            nbonds = [v for k, v in nbonds_dict.items() if atom_type in k][0]
        except IndexError:
            raise IndexError(
                "{} not in atomtype dictionary nbonds_dict".format(atom_type)
            )

        logging.debug(
            f"Potential radical with index {atom_idx} is bound to {self.AtomList[self.atoms_idx.index(atom_idx)].bound_to} other atoms"
        )
        if len(self.AtomList[self.atoms_idx.index(atom_idx)].bound_to) < nbonds:
            logging.info(f"{atom_idx} is a radical")
            return True
        else:
            logging.info(f"{atom_idx} is not a radical")
            return False

    def ff_AA_todict(self):
        with open(self.ffdir / "aminoacids.rtp", "r") as ff:
            foo = ff.readlines()
        ffaminoacids = {}
        key = "header"
        ffaminoacids[key] = {}
        subkey = "header"
        ffaminoacids[key][subkey] = []
        for line in foo:
            if len(line.split()) > 0:
                if line.startswith("["):
                    key = line.split()[1].strip()
                    ffaminoacids[key] = {}
                    subkey = "header"
                    ffaminoacids[key][subkey] = []
                elif line.startswith(" ["):
                    subkey = line.split()[1].strip()
                    ffaminoacids[key][subkey] = []
                else:
                    line = re.sub("[+-]", "", line)  # useful?
                    ffaminoacids[key][subkey].append(line.split())
        return ffaminoacids

    def compare_ff_impropers(self, heavy_idx, to_idx):
        """
        looks for obsolete patched impropers
        and impropers in the force field that should be in the topology
        """
        # improper_bbdict = {('C','CA','N','H'):['180.00','4.60240','2'],('CA','N','C','O'):['180.00','43.93200','2'],('N','CA','C','N'):['105.4','0.75','1']}
        ffaminoacids = self.ff_AA_todict()
        logging.debug(ffaminoacids["ALA"]["impropers"])

        rmvdict = {"impropers": []}
        adddict = {"impropers": []}
        # heavy_res = self.atoms_resname[self.atoms_idx.index(heavy_idx)]
        # heavy_atomname = self.atoms_atomname[self.atoms_idx.index(heavy_idx)]
        to_res = self.atoms_resname[self.atoms_idx.index(to_idx)]
        to_atomname = self.atoms_atomname[self.atoms_idx.index(to_idx)]

        # from_impropers = []
        to_impropers = []
        for improper in self.ImproperList:
            # if heavy_idx == improper[2]:
            #     from_impropers.append(improper)
            if to_idx == improper[2]:
                to_impropers.append(improper)

        impropers_atomnames = [
            [
                *[self.atoms_atomname[self.atoms_idx.index(x)] for x in improper[:4]],
                *improper[4:],
            ]
            for improper in to_impropers
        ]
        for i, improper_atomname in enumerate(impropers_atomnames):
            if not improper_atomname[:4] in ffaminoacids[to_res]["impropers"]:
                rmvdict["impropers"].append(to_impropers[i])

        for improper_res in ffaminoacids[to_res][
            "impropers"
        ]:  # probably need to work on this
            if to_atomname == improper_res[2]:
                improper_res_idx = [
                    self.atoms_idx[self.atoms_atomname.index(x)] for x in improper_res
                ]

                logging.debug(
                    f"Added improper {improper_res}/{improper_res_idx} from FF"
                )
                adddict["impropers"].append([*improper_res_idx, "4"])

        return adddict, rmvdict

    def get_H_ff_at_an(
        self, heavy_idx
    ):  # might have to be refined to find the exact atom name
        """
        Takes a heavy atom index and searches in
        the force field aminoacids.rtp for the residue definition of that heavy atom.
        Then, a bond of that heavy atom to a hydrogen is searched in the residue definition and its
        atomtype and atomname returned
        """
        heavy_atomname = self.atoms_atomname[self.atoms_idx.index(heavy_idx)]
        heavy_res = self.atoms_resname[self.atoms_idx.index(heavy_idx)]
        logging.debug([heavy_atomname, heavy_res])
        ffaminoacids = self.ff_AA_todict()

        for bond in ffaminoacids[heavy_res]["bonds"]:
            if heavy_atomname in bond:
                if any([x.startswith("H") for x in bond]):
                    H_atomname = bond[0] if bond[0] != heavy_atomname else bond[1]
                    for atom in ffaminoacids[heavy_res]["atoms"]:
                        if atom[0] == H_atomname:
                            logging.info(
                                f"HAT hydrogen atomtype has been found to be {atom[1]} at its new position."
                            )
                            return (atom[1], H_atomname)
        logging.warn(f"Found no new atomtype for HAT hydrogen!")

    def get_ff_sections(self):
        ffbonded = read_topol(self.ffdir / "ffbonded.itp")  # not intended use of read_topol but should work fine
        self.ff["bondtypes"] = ffbonded["bondtypes"]
        self.ff["angletypes"] = ffbonded["angletypes"]

    def parameterize_bonded_terms(self, prop, terms):
        """
        takes a term (bond or angle) and adds the parameters from
        the force field ffbonded.itp file that match the atomtypes
        """
        logging.info(
            f"Looking for atom parameters in atoms section of topology for {terms}"
        )
        terms = deepcopy(terms)
        terms_prm = []

        terms_mapped = [
            [self.atoms_atomtype[self.atoms_idx.index(x)] for x in term]
            for term in terms
        ]

        for i, term in enumerate(terms_mapped):
            for entry in self.ff[prop]:
                if term == entry[: len(term)] or term[::-1] == entry[: len(term)]:
                    stop = entry.index(";")
                    terms_prm.append([*terms[i], *entry[len(term) : stop]])
                    break
            else:
                logging.warning(f"No parameters found for {term}/{terms[i]}!")
        logging.debug(f"Parameterized these terms: {terms_prm}")
        return terms_prm


## Parameterization
    def patch_bond(self, bonds, atom_idx, newfrac = 0.98, SOfrac = 0.955, NCaR_offset = 0.006, CNR_offset = 0.008):
        logging.debug(f"Patching bonds {bonds}")
        # factor by which all non-hydrogen bonds of the atom_idx are corrected

        for bond in bonds:
            partnerpos = 0 if str(atom_idx) == bond[1] else 1
            radicalpos = bond[1 - partnerpos]
            partner_idx = bond[partnerpos]
            partner_atomtype = self.atoms_atomtype[self.atoms_idx.index(partner_idx)]
            radical_atomtype = self.atoms_atomtype[self.atoms_idx.index(atom_idx)]
            radicalname = self.atoms_atomname[self.atoms_idx.index(atom_idx)]

            if partner_atomtype.startswith("H"):
                # logging.debug(f"{partner_idx} is a hydrogen, no change of parameters necessary!")
                continue

            req = float(bond[3]) * newfrac  # change to r_eq happens here
            k = float(bond[4])

            # special fixes to bond length that are added on top of the general *0.98 factor
            if (
                partner_atomtype == "N"
                and radical_atomtype == "CT"
                and radicalname == "CA"
            ):  # N-CA of C-alpha atom_idx
                logging.debug("!!Found CaR!!")
                req -= NCaR_offset

            if (
                partner_atomtype == "C"
                and radical_atomtype == "N"
                and radicalname == "N"
            ):  # C-N of N atom_idx
                logging.debug("!!Found bb NR!!")
                req += CNR_offset

            if any(at in ["S", "O"] for at in [radical_atomtype, partner_atomtype]):
                req = req * (SOfrac / newfrac)  # correcting by a different value

            bond[3] = "{:7.5f}".format(req)  # will carry back to atom_terms
            bond[4] = "{:13.6f}".format(k)
            bond.append(" ; patched parameter")

    def patch_angle(self, angles, atom_idx, newtheteq = 117, aromatic_offset = 10  ):
        logging.debug(f"Patching angles {angles}")        

        for angle in angles:
            radical_atomtype = self.atoms_atomtype[self.atoms_idx.index(atom_idx)]

            theteq = float(angle[4])
            k = float(angle[5])

            if radical_atomtype in ["CR", "NA", "CW", "CC", "CA", "CN", "C*"]:
                theteq += aromatic_offset
            else:
                theteq = newtheteq

            angle[4] = "{:11.7f}".format(theteq)  # will carry back to atom_terms
            angle[5] = "{:10.6f}".format(k)
            angle.append(" ; patched parameter")

    def patch_dihedral_CA(self, propers, atom_idx, phivals = ["1.6279944", "21.068532", "1.447664"], psivals = ["6.556746", "20.284450", "0.297901"]):
         # phivals and psivals from own MCSA
        logging.debug(f"Attempting to patch propers {propers}")

        radical_resname = self.atoms_resname[self.atoms_idx.index(atom_idx)]
        radical_atomname = self.atoms_atomname[self.atoms_idx.index(atom_idx)]

        if radical_resname in ["Pro", "Hyp"] or not radical_atomname == "CA":
            logging.debug(
                f"Did not add a dihedral term because the atom_idx is in Pro or Hyp (here {radical_resname}) or not a C-alpha (here {radical_atomname})"
            )
            return propers

        phi = []
        psi = []
        propers_mapped = [
            [self.atoms_atomname[self.atoms_idx.index(x)] for x in dihedral[:4]]
            for dihedral in propers
        ]
        logging.debug(f"dihedrals by atomname {propers_mapped}")
        for i, dihedral in enumerate(propers_mapped):
            if dihedral == ["C", "N", "CA", "C"]:  # Phi
                for ii in range(3):
                    phi.append(
                        [
                            *propers[i],
                            "180.000000",
                            phivals[ii],
                            str(ii + 1),
                            " ; patched parameter",
                        ]
                    )

            elif dihedral == ["N", "CA", "C", "N"]:  # Psi
                for ii in range(3):
                    psi.append(
                        [
                            *propers[i],
                            "180.000000",
                            psivals[ii],
                            str(ii + 1),
                            " ; patched parameter",
                        ]
                    )

        return [*phi, *psi]

    def patch_improper(self, atom_idx,newphik = "43.93200"):
        # same as backbone N improper or 4.6024000
        logging.debug(f"Attempting to patch improper")
          
        radical_Atom = self.AtomList[self.atoms_idx.index(atom_idx)]

        logging.debug(
            f"Potential improper center {atom_idx} is bound to {len(radical_Atom.bound_to)} Atoms."
        )
        # this would mean atom_idx went from 4 to 3 partners -> tetrahedral to planar
        if (len(radical_Atom.bound_to) == 3):  
            improper = [
                radical_Atom.bound_to[0],
                radical_Atom.bound_to[1],
                str(atom_idx),
                radical_Atom.bound_to[2],
                "4",
                "180.0000000",
                newphik,
                "2",
                " ; patched parameter",
            ]  # improper entry
            return [improper]
        else:
            return []

    def parameterize_around_atom(self, atom_idx):
        logging.info(f" --- Parameterizing around atom {atom_idx} ---\n")
        logging.debug(self.atoms_idx)
        logging.debug(self.atoms_atomtype)
        logging.debug(self.atoms_atomname)
        logging.debug(self.atoms_resname)

        atom_terms = self.get_terms_with_atom(atom_idx, center=True)
        atom_terms["pairs"].clear()

        # deal with ff
        self.get_ff_sections()
        atom_terms["bonds"] = self.parameterize_bonded_terms(
            "bondtypes", atom_terms["bonds"]
        )
        atom_terms["angles"] = self.parameterize_bonded_terms(
            "angletypes", atom_terms["angles"]
        )
         # adding function to propers, impropers
        atom_terms["propers"] = [[*x, "9"] for x in atom_terms["propers"]] 
        atom_terms["impropers"] = [[*x, "4"] if len(x) == 4 else x for x in atom_terms["impropers"]]

        # apply patches
        if self.is_radical(atom_idx):
            self.patch_bond(atom_terms["bonds"], atom_idx)
            self.patch_angle(atom_terms["angles"], atom_idx)
            atom_terms["propers"] = self.patch_dihedral_CA(
                atom_terms["propers"], atom_idx
            )
            atom_terms["impropers"] = self.patch_improper(atom_idx)
        return atom_terms


def find_heavy(bonds, H_idx: str):
    """
    takes the bonds section of a topology object and an atom index
    to return the index of the first found bond partner.
    """
    for bond in bonds:
        if H_idx in bond[:2]:
            if H_idx == bond[0]:
                return bond[1]
            elif H_idx == bond[1]:
                return bond[0]
    raise ValueError(
        'Atom from ConversionRecipe not found in "bonds" section of the topology'
    )

def terms_keep_only(movepair,heavy_idx,atom_terms):
    """
    this is to stop overwriting bonds,angles for the newly parameterized "from" part if they are next to each other
    """
    mapsection = {"bonds": slice(0, 2), "angles": slice(1, 2), "propers": slice(1, 3)}
    rmvdict = {"bonds": [], "angles": [], "propers": []}
    for section in rmvdict.keys():  
        for term in atom_terms[section]:
            if movepair[0] in term:
                continue
            if heavy_idx in term[mapsection[section]]:
                rmvdict[section].append(term)
        for entry in rmvdict[section][::-1]:
            atom_terms[section].remove(entry)
    return atom_terms

## operations on the topology dict
def topol_add_terms(topology, adddict):
    """
    adds terms from adddict to the topologydict and returns the modified dictionary
    this could mean replacing or adding terms
    """
    logging.info(f"Adding/replacing following terms: {adddict}")
    copydict = deepcopy(topology)
    mapsection = {"pairs": 2, "bonds": 2, "angles": 3, "propers": 4, "impropers": 4}
    sortkeys = {
        "pairs": sort_bond,
        "bonds": sort_bond,
        "angles": sort_angle,
        "propers": sort_dihedral,
        "impropers": sort_improper,
    }
    for section in ["pairs", "bonds", "angles", "propers", "impropers"]:
        mapval = mapsection[section]
        if section in adddict.keys():
            for term in adddict[section]:
                try:
                    searchlist = list(
                        map(itemgetter(*list(range(mapval))), copydict[section])
                    )  # creating list of
                    if (
                        section == "propers" and len(term) > 5 and term[7] in ["2", "3"]
                    ):  # unhappy that I have to add this line to get all three dihedral terms
                        copydict[section].append(term)
                        logging.debug(
                            f"appended term {term} to section {section} because it is thought to be part of a multi-dihedral term"
                        )
                    else:
                        copydict[section][searchlist.index(tuple(term[:mapval]))] = term
                        logging.debug(f"replaced previous term with {term}")
                except ValueError:
                    logging.debug(
                        f"term {term} not in topology section {section}, appending it at the end"
                    )
                    copydict[section].append(term)
            copydict[section] = sorted(copydict[section], key=sortkeys[section])
    logging.debug(f"sorted section: {copydict[section]}")

    rmv_dih = []
    prevdih = []
    for i, entry in enumerate(copydict["propers"]):
        if entry[:4] == prevdih[:4]:
            if len(prevdih) == 5 and len(entry) > 5:
                rmv_dih.append(i)
                continue
        else:
            prevdih = entry

    logging.debug(f"Removing propers at these positions: {rmv_dih}")
    for val in rmv_dih[::-1]:
        del copydict["propers"][val]

    return copydict


def topol_remove_terms(topology, rmvdict):
    """
    removes terms in rmvdict that are part of the topologydict and returns the modified dictionary
    """
    logging.info(f"Deleting following terms: {rmvdict}")
    copydict = deepcopy(topology)
    mapsection = {"pairs": 2, "bonds": 2, "angles": 3, "propers": 4, "impropers": 4}
    for section in ["pairs", "bonds", "angles", "propers", "impropers"]:
        mapval = mapsection[section]
        if section in rmvdict.keys():
            for term in rmvdict[section]:
                try:
                    searchlist = list(
                        map(itemgetter(*list(range(mapval))), copydict[section])
                    )  # creating list of
                    del copydict[section][searchlist.index(tuple(term[:mapval]))]
                    logging.debug(f"removed term {term}")
                except ValueError:
                    logging.debug(
                        f"Couldn't remove term {term} because it is not part of the topology section {section}"
                    )

    return copydict


def topol_change_at_an(topology, atom):
    """
    changing the atoms entry of topology to give the HAT hydrogen the right atomtype and atomname
    """
    for entry in topology["atoms"]:
        if atom[2] == entry[4] and atom[3] == entry[3] and str_to_int(atom[2][-1]) != 0:
            atom[2] = atom[2][:-1] + str(int(atom[2][-1]) + 1)
    for entry in topology["atoms"]:
        if atom[0] == entry[0]:
            logging.debug(f"changed at/an of {entry} to {atom}")
            entry[1] = atom[1]
            entry[4] = atom[2]
    return

def modify_plumed(
    recipe: ConversionRecipe,
    oldplumeddat: Path,
    newplumeddat: Path,
    plumeddist: str,
):

    logging.info(
        f"Reading: {oldplumeddat} and writing modified plumed input to {newplumeddat}."
    )
    plumeddat = read_plumed(oldplumeddat)

    for type, pair in zip(recipe.type, recipe.atom_idx):
        if type == ConversionType.BREAK:
            plumeddat = break_bond_plumed(plumeddat, pair, plumeddist)
        elif type == ConversionType.MOVE:
            plumeddat = move_bond_plumed(plumeddat, pair, plumeddist)

    write_plumed(plumeddat, newplumeddat)

def break_bond_plumed(plumeddat, breakpair, plumeddist):
    new_distances = []
    broken_distances = []
    for line in plumeddat["distances"]:
        if all(x in line["atoms"] for x in breakpair):
            broken_distances.append(line["id"])
        else:
            new_distances.append(line)

    plumeddat["distances"] = new_distances

    for line in plumeddat["prints"]:
        line["ARG"] = [id for id in line["ARG"] if not id in broken_distances]
        line["FILE"] = plumeddist

    return plumeddat


def move_bond_plumed(plumeddat, movepair, plumeddist):
    raise NotImplementedError(
        "Plumeddat Changer for moving Atoms is not implemented yet."
    )
