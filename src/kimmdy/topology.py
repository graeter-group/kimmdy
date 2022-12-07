from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import Any, Hashable, Optional, Tuple
from xml.etree.ElementTree import Element
from kimmdy.parsing import TopologyDict, read_topol, read_xml_ff, is_not_comment
from itertools import takewhile, permutations
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
class AtomType:
    """Information about one atom

    A class containing atom information as in the atoms section of the topology.
    An atom keeps a list of which atoms it is bound to.

    From gromacs version of the amber* ff:
    ; name      at.num  mass     charge ptype  sigma      epsilon
    """

    type: str
    at_num: str
    charge: str
    mass: str
    ptype: str
    sigma: str
    epsilon: str

    @classmethod
    def from_top_line(cls, l: list[str]):
        return cls(
            type=l[0],
            at_num=l[1],
            charge=l[2],
            mass=l[3],
            ptype=l[4],
            sigma=l[5],
            epsilon=l[6],
        )

@dataclass(order=True)
class Bond:
    """Information about one bond

    A class containing bond information as in the bonds section of the topology.
    From gromacs topology:
    'ai', 'aj', 'funct', 'c0', 'c1', 'c2', 'c3
    With ai < aj
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
    With ai < ak
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
        return cls(
            ai=l[0],
            aj=l[1],
            funct=l[2],
            c0=field_or_none(l, 3),
            c1=field_or_none(l, 4),
            c2=field_or_none(l, 5),
            c3=field_or_none(l, 6),
        )


class FF:
    """Conainer for parsed forcefield data."""

    def __init__(self, ffdir: Path):
        self.atomtypes: dict[str, AtomType] = {}
        self.bondtypes: dict[tuple[str, str], Bond] = {}
        self.angletypes: dict[tuple[str, str, str], Angle] = {}
        self.dihedraltypes: dict[tuple[str, str, str, str], Dihedral] = {}


        nonbonded_path = ffdir / "ffnonbonded.itp"
        nonbonded = read_topol(nonbonded_path)
        for l in nonbonded["atomtypes"]:
            if l[0][0] != ";":
                atomtype = AtomType.from_top_line(l)
                self.atomtypes[atomtype.type] = atomtype

        bonded_path = ffdir / "ffbonded.itp"
        bonded = read_topol(bonded_path)
        for l in bonded["bondtypes"]:
            bond = Bond.from_top_line(l)
            self.bondtypes[(bond.ai, bond.aj)] = bond
        for l in bonded["angletypes"]:
            angle = Angle.from_top_line(l)
            self.angletypes[(angle.ai, angle.aj, angle.ak)] = angle
        for l in bonded["dihedraltypes"]:
            dihedral = Dihedral.from_top_line(l)
            self.dihedraltypes[(dihedral.ai, dihedral.aj, dihedral.ak, dihedral.al)] = dihedral

    def __repr__(self) -> str:
        return textwrap.dedent(
            f"""\
        ForceField parameters with
        {len(self.atomtypes)} atomtypes,
        {len(self.bondtypes)} bondtypes,
        {len(self.angletypes)} angletypes,
        {len(self.dihedraltypes)} dihedraltypes
        """
        )

class FFPatches():
    """A container for forcefield patches
    """
    def __init__(self, path: Path) -> None:
        xml = read_xml_ff(path)
        self.atompatches = xml.findall("Atoms/Atom")
        self.bondpatches = xml.findall("Bonds/Bond")
        self.pairpatches = xml.findall("Pairs/Pair")
        self.anglepatches = xml.findall("Angles/Angle")

        # damn... <https://manual.gromacs.org/current/reference-manual/functions/bonded-interactions.html#proper-dihedrals-periodic-type>
        # dihedrals of type 9 can have multiple definitions
        # so the angle phi is both an identifier and parameter for the dihedral
        # e.g. in ffbonded
        # CT  CT  OS  CT    9       0.0      1.60247     3  ; 
        # CT  CT  OS  CT    9     180.0      0.41840     2  ; 
        # idea: first match up on atomtypes and phi
        # if phi does not match add a new dihedral with the new phi
        # and leave the old one in place
        self.dihedralpatches = xml.findall("Dihedrals/Dihedral")

    def __repr__(self) -> str:
        return textwrap.dedent(
            f"""\
        ForceField parameter patches with
        {len(self.atompatches)} atompatches,
        {len(self.bondpatches)} bondpatches,
        {len(self.pairpatches)} pairpatches,
        {len(self.anglepatches)} bondpatches,
        {len(self.dihedralpatches)} dihedralpatches,
        """
        )


class Topology:
    """Smart container for parsed topology data.

    A topology keeps track of connections and applies patches to parameters when bonds are broken or formed.
    """

    def __init__(
        self, top: TopologyDict, ffdir: Path, ffpatch: Optional[Path] = None
    ) -> None:
        self.top = top
        self.forcefield_directory = ffdir
        self.ff = FF(ffdir)
        if ffpatch:
            self.ffpatches = FFPatches(ffpatch)

        self.atoms: dict[str, Atom] = {}
        self.bonds: dict[tuple[str, str], Bond] = {}
        self.pairs: dict[tuple[str, str], Pair] = {}
        self.angles: dict[tuple[str, str, str], Angle] = {}
        self.dihedrals: dict[tuple[str, str, str, str], Dihedral] = {}
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
        self.top["bonds"] = [attributes_to_list(x) for x in self.bonds.values()]
        self.top["pairs"] = [attributes_to_list(x) for x in self.pairs.values()]
        self.top["angles"] = [attributes_to_list(x) for x in self.angles.values()]
        self.top["dihedrals"] = [attributes_to_list(x) for x in self.dihedrals.values()]

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
            atom = Atom.from_top_line(l)
            self.atoms[atom.nr] = atom

    def _get_bonds(self):
        ls = self.top["bonds"]
        for l in ls:
            bond = Bond.from_top_line(l)
            self.bonds[(bond.ai, bond.aj)] = bond

    def _get_pairs(self):
        ls = self.top["pairs"]
        for l in ls:
            pair = Pair.from_top_line(l)
            self.pairs[(pair.ai, pair.aj)] = pair

    def _get_angles(self):
        ls = self.top["angles"]
        for l in ls:
            angle = Angle.from_top_line(l)
            self.angles[(angle.ai, angle.aj, angle.ak)] = angle

    def _get_dihedrals(self):
        ls = self.top["dihedrals"]
        for l in ls:
            dihedral = Dihedral.from_top_line(l)
            self.dihedrals[(dihedral.ai, dihedral.aj, dihedral.ak, dihedral.al)] = dihedral

    def get_proper_dihedrals(self):
        return [dihedral for dihedral in self.dihedrals.values() if dihedral.funct == "9"]

    def get_improper_dihedrals(self):
        return [dihedral for dihedral in self.dihedrals.values() if dihedral.funct == "4"]

    def _initialize_graph(self):
        for bond in self.bonds.values():
            i = bond.ai
            j = bond.aj
            self.atoms[i].bound_to_nrs.append(j)
            self.atoms[j].bound_to_nrs.append(i)

    def break_bond(self, atompair_nrs: tuple[str, str]):
        """Break bonds in topology.

        removes bond, angles and dihedrals where atompair was involved.
        Modifies the topology dictionary in place.
        It modifies the function types and parameters in the topology to account for radicals.
        """
        atompair_nrs = tuple(sorted(atompair_nrs, key=str_to_int_or_0))
        atompair = [self.atoms[atompair_nrs[0]], self.atoms[atompair_nrs[1]]]

        # bonds
        # remove bonds
        removed = self.bonds.pop(atompair_nrs, None)
        logging.info(f"removed bond: {removed}")

        # update bound_to
        try:
            atompair[0].bound_to_nrs.remove(atompair[1].nr)
            atompair[1].bound_to_nrs.remove(atompair[0].nr)
        except ValueError as _:
            m = f"tried to remove bond between already disconnected atoms: {atompair}."
            logging.warning(m)

        # remove angles
        self.angles = {
                key:angle
                for key, angle in self.angles.items()
                if not all(x in [angle.ai, angle.aj, angle.ak] for x in atompair_nrs)
        }

        # remove proper and improper dihedrals
        self.dihedrals = {
                key:dihedral
                for key,dihedral in self.dihedrals.items()
                if not all( x in [dihedral.ai, dihedral.aj, dihedral.ak, dihedral.al] for x in atompair_nrs)
                }

        # remove pairs
        dihpairs = [
            tuple(sorted((d.ai, d.al), key=str_to_int_or_0)) for d in self.dihedrals.values()
        ]
        self.pairs = {
            key: value for key, value in self.pairs.items() if key in dihpairs
        }

        if self.ffpatches is not None:
            self._patch_parameters(atompair)

        self._update_dict()

    def bind_bond(self, atompair_nrs: tuple[str, str]):
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

        atompair_nrs = tuple(sorted(atompair_nrs, key=str_to_int_or_0))
        atompair = [self.atoms[atompair_nrs[0]], self.atoms[atompair_nrs[1]]]

        # bonds
        # add bond
        bond = Bond(atompair_nrs[0], atompair_nrs[1], '1')
        self.bonds[atompair_nrs] = bond
        logging.info(f"added bond: {bond}")

        # update bound_to
        atompair[0].bound_to_nrs.append(atompair[1].nr)
        atompair[1].bound_to_nrs.append(atompair[0].nr)

        # add angles
        all_angles = (self._get_atom_angles(atompair_nrs[0]) +
                      self._get_atom_angles(atompair_nrs[1]))
        for key in all_angles:
            if self.angles.get(key) is None:
                self.angles[key] = Angle(key[0], key[1], key[2], '1')

        # add proper and improper dihedrals
        self.dihedrals = {
                key:dihedral
                for key,dihedral in self.dihedrals.items()
                if not all( x in [dihedral.ai, dihedral.aj, dihedral.ak, dihedral.al] for x in atompair)
                }

        # add pairs
        dihpairs = [
            tuple(sorted((d.ai, d.al), key=str_to_int_or_0)) for d in self.dihedrals.values()
        ]
        self.pairs = {
            key: value for key, value in self.pairs.items() if key in dihpairs
        }

        # if there are no changed parameters for radicals, exit here
        if self.ffpatches is None:
            self._update_dict()
            return

        if self.ffpatches is not None:
            self._patch_parameters(atompair)

        self._update_dict()

    def _get_atom_angles(self, atom_nr):
        """
        each atom has a list of atoms it is bound to
        get a list of angles that one atom is involved in
        based in these lists.
        Angles between atoms ai, aj, ak satisfy ai < ak
        """
        return (self._get_center_atom_angles(atom_nr) +
                self._get_margin_atom_angles(atom_nr))

    def _get_center_atom_angles(self, atom_nr):
        # atom_nr in the middle of an angle
        angles = []
        aj = atom_nr
        partners = self.atoms[aj].bound_to_nrs
        for ai in partners:
            for ak in partners:
                if int(ai) < int(ak):
                    angles.append((ai, aj, ak))
        return angles

    def _get_margin_atom_angles(self, atom_nr):
        # atom_nr in the middle of an angle
        angles = []
        ai = atom_nr
        for aj in self.atoms[ai].bound_to_nrs:
                for ak in self.atoms[aj].bound_to_nrs:
                    if ai == ak: continue
                    if int(ai) < int(ak):
                        angles.append((ai, aj, ak))
                    else:
                        angles.append((ak, aj, ai))
        return angles

    def _patch_parameters(self, atompair):
        # Adjust parameters based on patch
        # atoms
        if atompatches := self.ffpatches.atompatches:
            for atom in atompair:
                logging.info(f"Adjust parameters for atom {atom}.")

                # don't turn a radical into a radical radical
                if "_R" in atom.type:
                    continue

                atom.type = atom.type + "_R"
                print(atom.type)
                patch = match_attr(atompatches, "class1", atom.type)
                if patch is not None:
                    if mass_factor := patch.get("mass_factor"):
                        atom.mass = str(float(atom.mass) * float(mass_factor))
                    if charge_factor := patch.get("charge_factor"):
                        atom.charge = str(float(atom.charge) * float(charge_factor))

        # get (unbroken) bonds that the now radicals in the atompair are still involved in
        if bondpatches := self.ffpatches.bondpatches:
            for radical in atompair:
                for partner in radical.bound_to_nrs:
                    bond_key = tuple(sorted([radical.nr, partner], key=str_to_int_or_0))
                    bond = self.bonds[bond_key]
                    atom_i = self.atoms[bond.ai]
                    atom_j = self.atoms[bond.aj]
                    patch = match_attr(bondpatches, "class1", atom_i.type)
                    if patch is not None:
                        logging.info(f"Adjust parameters for bond {bond}.")
                        c0 = patch.get("c0")
                        if c0 is not None:
                            bond.c0 = c0

                        c0_factor = patch.get("c0_factor")
                        if c0_factor is not None:
                            original_atomtypes = (
                                atom_i.type.removesuffix("_R"),
                                atom_j.type.removesuffix("_R"),
                            )
                            bondtype = get_by_permutations(
                                self.ff.bondtypes, original_atomtypes
                            )
                            if bondtype is not None and bondtype.c0 is not None:
                                bond.c0 = str(float(bondtype.c0) * float(c0_factor))

        # get (unbroken) angles that the now radicals in the atompair are still involved in
        if anglepatches := self.ffpatches.anglepatches:
            for radical in atompair:
                angles = []
                for partner in radical.bound_to_nrs:
                    for partner_partner in self.atoms[partner].bound_to_nrs:
                        key = (radical.nr, partner, partner_partner)
                        angle = get_by_permutations(self.angles, key)
                        if angle is not None:
                            logging.info(f"Adjust parameters for angle {angle}.")
                            angles.append(angle)

        # get (unbroken) angles that the now radicals in the atompair are still involved in
        if dihedralpatches := self.ffpatches.dihedralpatches:
            for radical in atompair:
                dihedrals = []
                for partner in radical.bound_to_nrs:
                    for partner_partner in self.atoms[partner].bound_to_nrs:
                        for partner_partner_partner in self.atoms[partner_partner].bound_to_nrs:
                            key = (radical.nr, partner, partner_partner, partner_partner_partner)
                            dihedral = get_by_permutations(self.dihedrals, key)
                            if dihedral is not None:
                                logging.info(f"Adjust parameters for dihedral {dihedral}.")
                                dihedrals.append(dihedral)



def attributes_to_list(obj) -> list[str]:
    return list(takewhile(lambda x: x is not None, obj.__dict__.values()))


def field_or_none(l: list[str], i) -> Optional[str]:
    try:
        return l[i]
    except IndexError as _:
        return None


def is_not_none(x) -> bool:
    return x is None


def match_attr(patches: list[Element], attr: str, m: str) -> Optional[Element]:
    matches = []
    for p in patches:
        if value := p.get(attr):
            if value == m:
                return p
            pattern = value.replace("*", r".*").replace("+", r"\+")
            if re.match(pattern, m):
                matches.append(p)
    if matches:
        matches.sort(key=lambda x: x.get(attr))
        return matches[0]
    else:
        return None


def match_multi_attr(
    patches: list[Element], attrs: list[str], m: list[str]
) -> Optional[Element]:
    raise NotImplementedError('WIP')



def get_by_permutations(d: dict, key) -> Optional[Any]:
    # TODO: oh, we might not even need this in all cases
    # not just bonds, but also
    # angles and dihedrals have a well defined order
    for k in permutations(key):
        value = d.get(k, None)
        if value is not None:
            return value
    return None

def get_element_id(e: Element) -> Optional[str]:
    id = None
    if e.tag == "Atom":
        id = ""
    elif e.tag == "Bond":
        id = ""
    return id

