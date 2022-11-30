from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import Optional, Tuple
from xml.etree.ElementTree import Element
from kimmdy.parsing import TopologyDict, read_xml_ff
from itertools import takewhile
import re
import textwrap
import logging

from kimmdy.utils import sort_bond, str_to_int_or_0


@dataclass(order=True)
class Atom:
    """Information about one atom

    A class containing atom information as in the atoms section of the topology.
    An atom keeps a list of which atoms it is bound to.

    From gromacs topology:
    ; nr type resnr residue atom cgnr charge mass typeB chargeB massB
    """

    nr: str
    type: str
    resnr: str
    residue: str
    atom: str
    cgnr: str
    charge: str
    mass: str
    typeB: Optional[str] = None
    chargeB: Optional[str] = None
    massB: Optional[str] = None
    # TODO: use this with a local graph representation
    bound_to_nrs: list[str] = field(default_factory=list)

    @classmethod
    def from_top_line(cls, l: list[str]):
        l = list(takewhile(is_not_comment, l))

        return cls(
            nr=l[0],
            type=l[1],
            resnr=l[2],
            residue=l[3],
            atom=l[4],
            cgnr=l[5],
            charge=l[6],
            mass=l[7],
            typeB=field_or_none(l, 8),
            chargeB=field_or_none(l, 9),
            massB=field_or_none(l, 9),
        )


@dataclass(order=True)
class Bond:
    """Information about one bond

    A class containing bond information as in the bonds section of the topology.

    From gromacs topology:
    'ai', 'aj', 'funct', 'c0', 'c1', 'c2', 'c3
    """

    ai: str
    aj: str
    funct: str
    c0: Optional[str] = None
    c1: Optional[str] = None
    c2: Optional[str] = None
    c3: Optional[str] = None

    @classmethod
    def from_top_line(cls, l: list[str]):
        l = list(takewhile(is_not_comment, l))

        return cls(
            ai=l[0],
            aj=l[1],
            funct=l[2],
            c0=field_or_none(l, 3),
            c1=field_or_none(l, 4),
            c2=field_or_none(l, 5),
            c3=field_or_none(l, 6),
        )


@dataclass(order=True)
class Dihedral:
    """Information about one dihedral

    A class containing bond information as in the dihedrals section of the topology.
    Proper dihedrals have funct 9.
    Improper dihedrals have funct 4.

    From gromacs topology:
    ';', 'ai', 'aj', 'ak', 'al', 'funct', 'c0', 'c1', 'c2', 'c3', 'c4', 'c5'
    """

    ai: str
    aj: str
    ak: str
    al: str
    funct: str
    c0: Optional[str] = None
    c1: Optional[str] = None
    c2: Optional[str] = None
    c3: Optional[str] = None
    c4: Optional[str] = None
    c5: Optional[str] = None

    @classmethod
    def from_top_line(cls, l: list[str]):
        l = list(takewhile(is_not_comment, l))

        return cls(
            ai=l[0],
            aj=l[1],
            ak=l[2],
            al=l[3],
            funct=l[4],
            c0=field_or_none(l, 5),
            c1=field_or_none(l, 6),
            c2=field_or_none(l, 7),
            c3=field_or_none(l, 8),
            c4=field_or_none(l, 9),
            c5=field_or_none(l, 10),
        )


@dataclass(order=True)
class Angle:
    """Information about one angle

    A class containing angle information as in the angles section of the topology.

    From gromacs topology:
    ';', 'ai', 'aj', 'ak', 'funct', 'c0', 'c1', 'c2', 'c3'
    """

    ai: str
    aj: str
    ak: str
    funct: str
    c0: Optional[str] = None
    c1: Optional[str] = None
    c2: Optional[str] = None
    c3: Optional[str] = None

    @classmethod
    def from_top_line(cls, l: list[str]):
        l = list(takewhile(is_not_comment, l))

        return cls(
            ai=l[0],
            aj=l[1],
            ak=l[2],
            funct=l[3],
            c0=field_or_none(l, 4),
            c1=field_or_none(l, 5),
            c2=field_or_none(l, 6),
            c3=field_or_none(l, 7),
        )


@dataclass(order=True)
class Pair:
    """Information about one pair

    A class containing pair information as in the pair section of the topology.

    From gromacs topology:
    ai', 'aj', 'funct', 'c0', 'c1', 'c2', 'c3'
    """

    ai: str
    aj: str
    funct: str
    c0: Optional[str] = None
    c1: Optional[str] = None
    c2: Optional[str] = None
    c3: Optional[str] = None

    @classmethod
    def from_top_line(cls, l: list[str]):
        l = list(takewhile(is_not_comment, l))

        return cls(
            ai=l[0],
            aj=l[1],
            funct=l[2],
            c0=field_or_none(l, 3),
            c1=field_or_none(l, 4),
            c2=field_or_none(l, 5),
            c3=field_or_none(l, 6),
        )


class Topology:
    def __init__(
        self, top: TopologyDict, ffdir: Path, ffpatch: Optional[Path] = None
    ) -> None:
        self.top = top
        self.patch = None
        if ffpatch:
            self.patch = read_xml_ff(ffpatch)
        self.forcefield_directory = ffdir
        self.ff = {}

        self.atoms: dict[str, Atom] = {}
        self.bonds: list[Bond] = []
        self.dihedrals: list[Dihedral] = []
        self.pairs: list[Pair] = []
        self.angles: list[Angle] = []
        self.proper_dihedrals: list[Dihedral] = []
        self.improper_dihedrals: list[Dihedral] = []

        self._get_atoms()
        self._get_bonds()
        self._get_pairs()
        self._get_angles()
        self._get_dihedrals()

        self._update_dict()

        self._initialize_graph()

    def _update_dict(self):
        self.top["atoms"] = [attributes_to_list(x) for x in self.atoms.values()]
        self.top["bonds"] = [attributes_to_list(x) for x in self.bonds]
        self.top["pairs"] = [attributes_to_list(x) for x in self.pairs]
        self.top["angles"] = [attributes_to_list(x) for x in self.angles]
        self.top["dihedrals"] = [attributes_to_list(x) for x in self.dihedrals]

    def to_dict(self) -> TopologyDict:
        self._update_dict()
        return self.top

    def __repr__(self) -> str:
        return textwrap.dedent(
            f"""\
        Topology with
        {len(self.atoms)} atoms,
        {len(self.bonds)} bonds,
        {len(self.angles)} angles,
        {len(self.pairs)} pairs,
        {len(self.dihedrals)} dihedrals
        """
        )

    def __str__(self) -> str:
        return str(self.atoms)

    def _get_atoms(self):
        ls = self.top["atoms"]
        for l in ls:
            if l[0] != ";":
                atom = Atom.from_top_line(l)
                self.atoms[atom.nr] = atom

    def _get_bonds(self):
        ls = self.top["bonds"]
        for l in ls:
            if l[0] != ";":
                self.bonds.append(Bond.from_top_line(l))

    def _get_pairs(self):
        ls = self.top["pairs"]
        for l in ls:
            if l[0] != ";":
                self.pairs.append(Pair.from_top_line(l))

    def _get_angles(self):
        ls = self.top["angles"]
        for l in ls:
            if l[0] != ";":
                self.angles.append(Angle.from_top_line(l))

    def _get_dihedrals(self):
        ls = self.top["dihedrals"]
        for l in ls:
            if l[0] != ";":
                dihedral = Dihedral.from_top_line(l)
                self.dihedrals.append(dihedral)

    def get_proper_dihedrals(self):
        return [dihedral for dihedral in self.dihedrals if dihedral.funct == '9']

    def get_improper_dihedrals(self):
        return [dihedral for dihedral in self.dihedrals if dihedral.funct == '4']

    def _initialize_graph(self):
        for bond in self.bonds:
            i = bond.ai
            j = bond.aj
            self.atoms[i].bound_to_nrs.append(j)
            self.atoms[j].bound_to_nrs.append(i)

    def patch_parameters(self):
        pass

    def break_bond(self, atompair: tuple[str, str]):
        """Break bonds in topology.

        removes bond, angles and dihedrals where atompair was involved.
        Modifies the topology dictionary in place.
        It modifies to function types and parameters in the topology to account for radicals.
        """
        radical_pair = [atom for atom in self.atoms if atom.nr in atompair]
        # TODO: we can make this faster at some point by assuming
        # atoms are indexd by number.
        # Right now I don't think we can safely assume this from
        # just the topology file
        
        # atoms
        # let's make some radicals
        # maybe there are patches for atomtypes that become radicals
        if self.patch:
            if atompatches := self.patch.findall('Atoms/Atom[@class1]'):
                for atom in radical_pair:
                        logging.info(f"Adjust parameters for atom {atom.nr}.")

                        # don't turn a radical into a radical radical
                        if '_R' in atom.type: continue

                        atom.type = atom.type + '_R'
                        patch = match_attr(atompatches, 'class1', atom.type)
                        if patch is not None:
                            if mass_factor := patch.get('mass_factor'):
                                atom.mass = str(float(atom.mass) * float(mass_factor))
                            if charge_factor := patch.get('charge_factor'):
                                atom.charge = str(float(atom.charge) * float(charge_factor))

        # update bound_to
        # radical_pair[0].bound_to_nrs.remove(radical_pair[1].nr)
        # radical_pair[1].bound_to_nrs.remove(radical_pair[0].nr)

        # remove bonds
        self.bonds = [
            bond
            for bond in self.bonds
            if not all(x in [bond.ai, bond.aj] for x in atompair)
        ]

        # remove angles
        self.angles = [
            angle
            for angle in self.angles
            if not all(x in [angle.ai, angle.aj, angle.ak] for x in atompair)
        ]

        # remove proper and improper dihedrals
        self.dihedrals = [
            dihedral
            for dihedral in self.dihedrals
            if not all(x in [dihedral.ai, dihedral.aj, dihedral.ak, dihedral.al] for x in atompair)
        ]

        # remove pairs
        dihpairs = [[d.ai, d.al] if d.ai < d.al else [d.al, d.ai] for d in self.dihedrals]
        self.pairs = [pair for pair in self.pairs if [pair.ai, pair.aj] in dihpairs]

        ## TODO: matching ff patches to atoms might get ugly from here on out
        # get (unbroken) bonds that the now radicals in the atompair are still involved in
        for radical in radical_pair:
            for partner in radical.bound_to_nrs:
                print(radical)

        # update topology dictionary
        self._update_dict()


    def bind_bond(self, atompair: tuple[str, str]):
        """Add a bond in topology.

        Move an atom (typically H for Hydrogen Atom Transfer) to a new location.
        Modifies the topology dictionary in place.
        It keeps track of affected terms in the topology via a graph representation of the topology
        and applies the necessary changes to bonds, angles and dihedrals (proper and improper).
        Furthermore, it modifies to function types in the topology to account for radicals.

        Parameters
        ----------
        atompair: a tuple of integers with the atoms indices
            `from`, the atom being moved and
            `to`, the atom to which the `from` atom will be bound
        """
        raise NotImplemented("WIP")
        # TODO: remember to undo the ff patches when atoms are no longer radicals!

        # atompair_str = [str(x) for x in atompair]
        # split_dihedrals(topology)
        #
        # # build localGraph and fill it
        # to_graph = LocalGraph(topology, atompair_str[1], ffdir, add_bond=atompair_str)
        #
        # # set the correct atomtype,resname for the HAT hydrogen based on the ff definition
        # heavy_resname = atompair_str[0]
        # atmdef = to_graph.correct_atomprops(atompair_str, heavy_resname)
        # topol_change_at_an(topology, atmdef)
        #
        # # this is to prevent overwriting bonds,angles for the newly parameterized "from" part,
        # # if they are next to each other
        # atom_terms = to_graph.parameterize_around_atom(atompair_str[1])
        # atom_terms = terms_keep_only(atompair_str, heavy_idx, atom_terms)
        # topology = topol_add_terms(topology, atom_terms)
        #
        # # add pairs of the from_atom at the new position
        # atom_terms_from = to_graph.get_terms_with_atom(atompair_str[0], add_function=True)
        # for section in ["bonds", "angles", "propers", "impropers"]:
        #     atom_terms_from[section].clear()
        # topology = topol_add_terms(topology, atom_terms_from)
        #
        # # have the right impropers at the to_heavy atom
        # atom_terms_add, atom_terms_remove = to_graph.compare_ff_impropers(
        #     heavy_idx, atompair_str[1]
        # )
        # topology = topol_add_terms(topology, atom_terms_add)  # no need, yet
        # topology = topol_remove_terms(topology, atom_terms_remove)
        #
        # topology = merge_propers_impropers(topology)


def attributes_to_list(obj) -> list[str]:
    return list(takewhile(lambda x: x is not None, obj.__dict__.values()))


def field_or_none(l: list[str], i) -> Optional[str]:
    try:
        return l[i]
    except IndexError as _:
        return None


def is_not_comment(c: str):
    return c != ";"


def is_not_none(x) -> bool:
    return x is None

def match_attr(patches: list[Element], attr: str, m: str) -> Optional[Element]:
    matches = []
    for p in patches:
        if value := p.get(attr):
            if value == m: return p
            pattern = value.replace('*', r'.*').replace('+', r'\+')
            if re.match(pattern, m): matches.append(p)
    if matches:
        matches.sort(key=lambda x: x.get(attr))
        return matches[0]
    else:
        return None

