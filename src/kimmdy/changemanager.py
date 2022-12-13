from __future__ import annotations
from dataclasses import dataclass, field
import logging
from typing import Any, Optional
from kimmdy.reaction import ConversionRecipe, ConversionType
from kimmdy.constants import ATOMTYPE_BONDORDER
from kimmdy.parsing import (
    read_plumed,
    write_plumed,
    read_topol,
    write_topol,
    split_dihedrals,
    merge_propers_impropers,
    TopologyDict,
)
from kimmdy.topology import Topology
from kimmdy.utils import (
    str_to_int_or_0,
    check_idx,
    sort_bond,
    sort_angle,
    sort_dihedral,
    sort_improper,
)
from pathlib import Path

from copy import deepcopy
import itertools
from operator import itemgetter
import re


def modify_top(
    recipe: ConversionRecipe, oldtop: Path, newtop: Path, ffdir: Path, ffpatch: Path
):
    logging.info(f"Reading: {oldtop} and writing modified topology to {newtop}.")
    topologyDict = read_topol(oldtop)
    topology = Topology(topologyDict, ffdir, ffpatch)

    for conversion in recipe:
        if conversion.type == ConversionType.BREAK:
            topology.break_bond(conversion.atom_idx)
        elif conversion.type == ConversionType.BIND:
            topology.bind_bond(conversion.atom_idx)

    write_topol(topology.top, newtop)


def bind_bond(topology: TopologyDict, atompair: tuple[str, str], ffdir: Path):
    """Add a bond in topology.

    Move an atom (typically H for Hydrogen Atom Transfer) to a new location.
    Modifies the topology dictionary in place.
    It keeps track of affected terms in the topology via a graph representation of the topology
    and applies the necessary changes to bonds, angles and dihedrals (proper and improper).
    Furthermore, it modifies to function types in the topology to account for radicals.

    Parameters
    ----------
    topology: dict
        dictionary representation of the topology
    movepair: tuple[int, int]
        a tuple of integers with the atoms indices
        `from`, the atom being moved and
        `to`, the atom to which the `from` atom will be bound
    ffdir: Path
        path to the forcefield directory. Needs to contain aminoacids.rtp.
    """
    atompair_str = [str(x) for x in atompair]
    split_dihedrals(topology)

    # build localGraph and fill it
    to_graph = LocalGraph(topology, atompair_str[1], ffdir, add_bond=atompair_str)

    # set the correct atomtype,resname for the HAT hydrogen based on the ff definition
    heavy_resname = atompair_str[0]
    atmdef = to_graph.correct_atomprops(atompair_str, heavy_resname)
    topol_change_at_an(topology, atmdef)

    # this is to prevent overwriting bonds,angles for the newly parameterized "from" part,
    # if they are next to each other
    atom_terms = to_graph.parameterize_around_atom(atompair_str[1])
    atom_terms = terms_keep_only(atompair_str, heavy_idx, atom_terms)
    topology = topol_add_terms(topology, atom_terms)

    # add pairs of the from_atom at the new position
    atom_terms_from = to_graph.get_terms_with_atom(atompair_str[0], add_function=True)
    for section in ["bonds", "angles", "propers", "impropers"]:
        atom_terms_from[section].clear()
    topology = topol_add_terms(topology, atom_terms_from)

    # have the right impropers at the to_heavy atom
    atom_terms_add, atom_terms_remove = to_graph.compare_ff_impropers(
        heavy_idx, atompair_str[1]
    )
    topology = topol_add_terms(topology, atom_terms_add)  # no need, yet
    topology = topol_remove_terms(topology, atom_terms_remove)

    topology = merge_propers_impropers(topology)


@dataclass
class Atom:
    """Information about one atom

    A class containing atom information as in the atoms section of the topology.
    An atom keeps a list of which atoms it is bound to.

    From gromacs topology:
    ; nr type resnr residue atom cgnr charge mass typeB chargeB massB
    Here: nr = idx, type = atomtype, atom = atomname, residue = resname
    """

    idx: str
    atomtype: Optional[str] = None
    atomname: Optional[str] = None
    resname: Optional[str] = None
    # bound_to: list[str] = field(default_factory=list)
    bound_to: list[Atom] = field(default_factory=list)


@dataclass
class Bond:
    """Information about one bond

    From gromacs topology:
    ; ai aj funct c0 c1 c2 c3
    """

    i: str
    j: str
    funct: Optional[str] = None


def bond_section_to_bond_dict(ls: list[list[str]]) -> dict:
    d = {}
    for l in ls:
        if l[0] == ";":
            continue
        i, j, f = l
        i, j = sorted([i, j])
        d[(i, j)] = Bond(i, j, f)
    return d


def atom_section_to_atom_dict(ls: list[list[str]]) -> dict:
    d = {}
    for l in ls:
        if l[0] == ";":
            continue
        d[l[0]] = Atom(l[0], l[1], l[4], l[3])
    return d


def reciprocal_bonds(d: dict) -> dict:
    reciproce_dict = {}
    for k in d.keys():
        reciproce_dict[k[::-1]] = d[k]
    return reciproce_dict


class LocalGraph:
    def __init__(
        self,
        topology: Topology,
        # central_atom_idx: str,
        central_pair: tuple[str, str],
        ffdir: Path,
        depth: int = 3,
    ):
        self.topology = topology
        self.central_pair = central_pair
        # self.central_atom_idx = central_atom_idx
        self.ffdir = ffdir
        self.ff = {}

        # self.atoms = [Atom(self.central_atom_idx)]
        self.atoms: list[Atom] = []
        self.bonds: list[Bond] = []
        self.pairs = []
        self.angles = []
        self.proper_dihedrals = []
        self.improper_dihedrals = []
        self.topology_bonds = bond_section_to_bond_dict(topology["bonds"])
        self.topology_bonds_reciprocal = reciprocal_bonds(self.topology_bonds)
        self.topology_atoms = atom_section_to_atom_dict(topology["atoms"])

        # add initial two atoms to the graph
        one, two = [self.topology_atoms[idx] for idx in central_pair]
        one.bound_to.append(two)
        two.bound_to.append(one)
        self.atoms.append(one)
        self.atoms.append(two)

        # add bond if it exists
        if bond := self.topology_bonds[tuple(sorted([one.idx, two.idx]))]:
            self.bonds.append(bond)

        for atom in self.atoms:
            for key, bond in self.topology_bonds.items():
                if atom.idx in key:
                    self.bonds.append(bond)

        # self.construct_graph(depth)
        # self.order_lists()
        # self.update_atoms_list()
        # self.update_bound_to()
        # self.build_PADs()
        # self.order_lists()

    def __repr__(self) -> str:
        s = "\n".join(
            [
                f"atoms: {self.atoms}",
                f"bonds: {self.bonds}",
                f"pairs: {self.pairs}",
                f"angles: {self.angles}",
                f"proper_dihedrals: {self.proper_dihedrals}",
                f"improper_dihedrals: {self.improper_dihedrals}",
                # f'atomdict: {self.topology_atoms}',
                # f'bondsdict: {self.topology_bonds}',
            ]
        )
        return s

    def construct_graph(self, depth=3):
        """
        searches the bonds section of a topology up to `depth` bonds deep
        to create a local graph i.e. filling the self.atoms list and self.bonds list
        """
        curratoms = [self.central_atom_idx]
        for _ in range(depth):
            addlist = []
            for bond in self.topology["bonds"]:
                if any(x in bond[:2] for x in curratoms):
                    if not bond[:2] in self.bonds:
                        self.bonds.append(bond[:2])
                    for j in [0, 1]:
                        if not bond[j] in [x.idx for x in self.atoms]:
                            addlist.append(bond[j])
                            self.atoms.append(Atom(bond[j]))
                            self.atoms_idx.append(bond[j])
            curratoms = deepcopy(addlist)

    def order_lists(self):
        # return early on empty topology
        if not self.atoms:
            return
        self.atoms = sorted(self.atoms, key=check_idx)
        self.atoms_idx = sorted(self.atoms_idx, key=str_to_int_or_0)
        # TODO: find more elegant solution for ordering multiple properties
        self.atoms_atomtype = [
            x
            for _, x in sorted(zip(self.atoms_idx, self.atoms_atomtype), key=sort_bond)
        ]
        self.atoms_atomname = [
            x
            for _, x in sorted(zip(self.atoms_idx, self.atoms_atomname), key=sort_bond)
        ]
        self.atoms_resname = [
            x for _, x in sorted(zip(self.atoms_idx, self.atoms_resname), key=sort_bond)
        ]
        self.bonds = sorted(self.bonds, key=sort_bond)
        self.pairs = sorted(self.pairs, key=sort_bond)
        self.angles = sorted(self.angles, key=sort_angle)
        self.proper_dihedrals = sorted(self.proper_dihedrals, key=sort_dihedral)
        self.improper_dihedrals = sorted(self.improper_dihedrals, key=sort_improper)

        # keep bonds in proper order (which one goes first in the definition)
        for bond in self.bonds:
            bond_int = [int(bond[0]), int(bond[1])]
            if bond_int[0] > bond_int[1]:
                bond.reverse()

    # Propagate changes to all attributes
    def update_atoms_list(self):
        self.atoms_idx = self.get_atoms_property("idx")
        self.get_atomprop("atomtype")
        self.get_atomprop("atomname")
        self.get_atomprop("resname")
        self.atoms_atomtype = self.get_atoms_property("atomtype")
        self.atoms_atomname = self.get_atoms_property("atomname")
        self.atoms_resname = self.get_atoms_property("resname")

    def update_bound_to(self):
        bound_dict = {}
        for atom_idx in self.atoms_idx:
            bound_dict[atom_idx] = []
        for bond in self.bonds:
            if not bond[1] in bound_dict[bond[0]]:
                bound_dict[bond[0]].append(bond[1])
            if not bond[0] in bound_dict[bond[1]]:
                bound_dict[bond[1]].append(bond[0])
        for atom_idx in self.atoms_idx:
            self.atoms[self.atoms_idx.index(atom_idx)].bound_to = sorted(
                bound_dict[atom_idx], key=str_to_int_or_0
            )

    # Retrieve properties
    def get_atomprop(self, property: str):
        for entry in self.topology["atoms"]:
            if entry[0] in self.atoms_idx:
                if property == "atomtype":
                    self.atoms[self.atoms_idx.index(entry[0])].atomtype = entry[1]
                elif property == "atomname":
                    self.atoms[self.atoms_idx.index(entry[0])].atomname = entry[4]
                elif property == "resname":
                    self.atoms[self.atoms_idx.index(entry[0])].resname = entry[3]

    def get_atoms_property(self, property) -> list[Any]:
        if property == "idx":
            return [x.idx for x in self.atoms]
        elif property == "bound_to":
            return [x.bound_to for x in self.atoms]
        elif property == "atomtype":
            return [x.atomtype for x in self.atoms]
        elif property == "atomname":
            return [x.atomname for x in self.atoms]
        elif property == "resname":
            return [x.resname for x in self.atoms]
        else:
            raise ValueError("Can't find property {} in object".format(property))

    def build_PADs(self):
        """
        in:
        - localGraph with AtomList that contains atom objects with bound_to attribute that contains the bonded atom idxs
        out:
        - filled AngleList, ProperList
        - filled PairList as ProperList outer atoms
        - ImproperList as defined in the topology
        """
        # define angles
        for atom in self.atoms:
            if len(atom.bound_to) >= 2:
                for comb in itertools.combinations(atom.bound_to, 2):
                    self.angles.append([comb[0], atom.idx, comb[1]])

        # define propers
        for atom1 in self.atoms:
            if len(atom1.bound_to) >= 2:
                for partner in atom1.bound_to:
                    atom2 = self.atoms[self.atoms_idx.index(partner)]
                    if len(atom2.bound_to) >= 2 and int(atom2.idx) > int(atom1.idx):
                        atom1_periphery = [x for x in atom1.bound_to if x != atom2.idx]
                        atom2_periphery = [x for x in atom2.bound_to if x != atom1.idx]
                        for prod in itertools.product(atom1_periphery, atom2_periphery):
                            self.proper_dihedrals.append(
                                [prod[0], atom1.idx, atom2.idx, prod[1]]
                            )

        # take impropers from topology
        for term in self.topology["impropers"]:
            if all([x in self.atoms_idx for x in term[:4]]):
                self.improper_dihedrals.append(term)

        # all pairs have a dihedral between them
        for dihedral in self.proper_dihedrals:
            # check for existing pairs in the case of five membered rings
            if int(dihedral[0]) < int(dihedral[3]):
                self.pairs.append([dihedral[0], dihedral[3]])
            else:
                self.pairs.append([dihedral[3], dihedral[0]])

    def correct_atomprops(self, movepair, heavy_resname):
        atoms_idxs = self.atoms_idx
        H_atomtype, H_atomname = self.get_H_ff_at_an(movepair[1])
        self.atoms[atoms_idxs.index(movepair[0])].atomtype = H_atomtype
        self.atoms[atoms_idxs.index(movepair[0])].atomname = H_atomname
        self.atoms[atoms_idxs.index(movepair[0])].resname = heavy_resname
        self.update_atoms_list()
        return (
            movepair[0],
            H_atomtype,
            H_atomname,
            self.atoms[atoms_idxs.index(movepair[0])].resname,
        )

    # Methods related to retrieving ff terms from a LocalGraph that contain a certain atom (search for index)
    def search_terms(self, entries, idx: str, center: bool) -> list[Any]:
        if entries == []:
            return []
        hits = []
        if center and len(entries[0]) > 2:
            start = 1
            end = -1
        else:
            start = 0
            end = None
        for entry in entries:
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

    def get_terms_with_atom(
        self, idx: str, center: bool = False, add_function: bool = False
    ) -> dict[str, Any]:
        """
        in:
        - idx:str           = atom index
        - center:bool       = take angle/dihedral terms with idx in center positions only
        - add_function:bool = add GROMACS function # to the terms
        out:
        - termdict:dict     = bonds,pairs,angles,propers,impropers containing atom idx
        """
        if idx not in self.atoms_idx:
            msg = f"Atom with index {idx} is not in the atoms of the local_graph {self}"
            logging.error(msg)
            raise IndexError(msg)

        termdict = {}
        termdict["bonds"] = self.search_terms(self.bonds, idx, center)
        termdict["pairs"] = self.search_terms(self.pairs, idx, center)
        termdict["angles"] = self.search_terms(self.angles, idx, center)
        termdict["propers"] = self.search_terms(self.proper_dihedrals, idx, center)
        termdict["impropers"] = self.search_terms(self.improper_dihedrals, idx, center)

        if add_function:
            termdict_funct = self.add_function_to_terms(termdict)
            return termdict_funct
        return termdict

    # Manipulate attributes
    def add_atom(self, atom: Atom):
        if atom.idx not in self.atoms_idx:
            self.atoms.append(atom)
            self.atoms_idx.append(atom.idx)
        else:
            logging.warning(f"Atom with idx {atom.idx} already exists!")

    def add_bond(self, bond: list[str]):
        if all(x in self.atoms_idx for x in bond):
            self.bonds.append(bond)
        else:
            logging.warning(f"Could not establish bond between idxs {bond}!")

    def remove_terms(self, termdict):
        logging.debug(f"Removing terms {termdict} from local graph.")
        graphdict = {
            "bonds": self.bonds,
            "pairs": self.pairs,
            "angles": self.angles,
            "propers": self.proper_dihedrals,
            "impropers": self.improper_dihedrals,
        }
        for section in ["pairs", "bonds", "angles", "propers", "impropers"]:
            for term in termdict[section]:
                try:
                    graphdict[section].remove(term)
                    logging.debug(f"removed term {term} from graph {section}")
                except ValueError:
                    logging.warning(f"Couldn't find term {term} in TypeList")

        self.update_bound_to()

    # Method that interact with the force field files
    def is_radical(self, atom_idx):
        """
        in:
        - atom_idx              = index of the atom that is checked for being a radical
        - self.atoms_atomtype   = atomtypes in the LocalGraph (comes from the topology dict)
        - nbonds_dict           = # of bonds of amber FF99SB-ILDNP* atomtypes (as written in the paper)
        out:
        - bool whether atom_idx is a radical (-> has fewer bonds than would be expected for the atomtype)
        """

        # compare to atom type perception paper (2006) same as in HAT_utils.py

        atom_type = self.atoms_atomtype[self.atoms_idx.index(atom_idx)]
        try:
            nbonds = [v for k, v in ATOMTYPE_BONDORDER.items() if atom_type in k][0]
        except IndexError:
            raise IndexError(
                "{} not in atomtype dictionary nbonds_dict".format(atom_type)
            )

        logging.debug(
            f"Potential radical with index {atom_idx} is bound to {self.atoms[self.atoms_idx.index(atom_idx)].bound_to} other atoms"
        )
        if len(self.atoms[self.atoms_idx.index(atom_idx)].bound_to) < nbonds:
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
        """look for obsolete patched improper dihedrals
        and improper dihedrals in the force field that should be in the topology
        """
        # improper_bbdict = {('C','CA','N','H'):['180.00','4.60240','2'],('CA','N','C','O'):['180.00','43.93200','2'],('N','CA','C','N'):['105.4','0.75','1']}
        ffaminoacids = self.ff_AA_todict()
        logging.debug(ffaminoacids["ALA"]["impropers"])

        rmvdict = {"impropers": []}
        adddict = {"impropers": []}
        to_res = self.atoms_resname[self.atoms_idx.index(to_idx)]
        to_atomname = self.atoms_atomname[self.atoms_idx.index(to_idx)]

        to_impropers = []
        for improper in self.improper_dihedrals:
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

        # TODO: probably need to work on this
        for improper_res in ffaminoacids[to_res]["impropers"]:
            if to_atomname == improper_res[2]:
                improper_res_idx = [
                    self.atoms_idx[self.atoms_atomname.index(x)] for x in improper_res
                ]

                logging.debug(
                    f"Added improper {improper_res}/{improper_res_idx} from FF"
                )
                adddict["impropers"].append([*improper_res_idx, "4"])

        return adddict, rmvdict

    # TODO: might have to be refined to find the exact atom name
    def get_H_ff_at_an(self, heavy_idx: str):
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
        logging.warn("Found no new atomtype for HAT hydrogen!")

    def get_ff_sections(self):
        ffbonded = read_topol(
            self.ffdir / "ffbonded.itp"
        )  # not intended use of read_topol but should work fine
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

    # Parameterization
    def patch_bond(
        self,
        bonds,
        atom_idx,
        newfrac=0.98,
        SOfrac=0.955,
        NCaR_offset=0.006,
        CNR_offset=0.008,
    ):
        """Correct all non-hydrogen bonds of the atom_idx by a factor.

        Parameters
        ----------
        newfrac:    factor for r_eq of C,N
        SOfrac:     factor for r_eq of S,O
        NCaR_offset:offset for r_eq of the N-Ca backbone bond for a radical Ca
        CNR_offset: offset for r_eq of the C-N backbone for for a radical N
        """
        logging.debug(f"Patching bonds {bonds}")

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

            # change to r_eq happens here
            req = float(bond[3]) * newfrac
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
            bond.extend([";", "patched", "parameter"])

    def patch_angle(self, angles, atom_idx, newtheteq=117, aromatic_offset=10):
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
            angle.extend([";", "patched", "parameter"])

    def patch_dihedral_CA(
        self,
        propers,
        atom_idx,
        phivals=["1.6279944", "21.068532", "1.447664"],
        psivals=["6.556746", "20.284450", "0.297901"],
    ):
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
                            ";",
                            "patched",
                            "parameter",
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
                            ";",
                            "patched",
                            "parameter",
                        ]
                    )

        return [*phi, *psi]

    def patch_improper(self, atom_idx, newphik="43.93200"):
        # same as backbone N improper or 4.6024000
        logging.debug(f"Attempting to patch improper")

        radical_Atom = self.atoms[self.atoms_idx.index(atom_idx)]

        logging.debug(
            f"Potential improper center {atom_idx} is bound to {len(radical_Atom.bound_to)} Atoms."
        )
        # this would mean atom_idx went from 4 to 3 partners -> tetrahedral to planar
        if len(radical_Atom.bound_to) == 3:
            improper = [
                radical_Atom.bound_to[0],
                radical_Atom.bound_to[1],
                str(atom_idx),
                radical_Atom.bound_to[2],
                "4",
                "180.0000000",
                newphik,
                "2",
                ";",
                "patched",
                "parameter",
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
        atom_terms["impropers"] = [
            [*x, "4"] if len(x) == 4 else x for x in atom_terms["impropers"]
        ]

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


def terms_keep_only(movepair, heavy_idx, atom_terms):
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
        if (
            atom[2] == entry[4]
            and atom[3] == entry[3]
            and str_to_int_or_0(atom[2][-1]) != 0
        ):
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
    breakpair = [str(x) for x in breakpair]
    for line in plumeddat["distances"]:
        if all(x in line["atoms"] for x in breakpair):
            broken_distances.append(line["id"])
        else:
            new_distances.append(line)

    plumeddat["distances"] = new_distances

    for line in plumeddat["prints"]:
        line["ARG"] = [id for id in line["ARG"] if id not in broken_distances]
        line["FILE"] = plumeddist

    return plumeddat


def move_bond_plumed(plumeddat, movepair, plumeddist):
    raise NotImplementedError(
        "Plumeddat Changer for moving Atoms is not implemented yet."
    )
