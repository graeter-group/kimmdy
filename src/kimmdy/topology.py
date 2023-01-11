from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Optional
from xml.etree.ElementTree import Element
from kimmdy.parsing import TopologyDict, read_topol, read_xml_ff, read_rtp
from itertools import takewhile, permutations, combinations
from xml.etree.ElementTree import Element
import re
import textwrap
import logging


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
    bound_to_nrs: list[str] = field(default_factory=list)
    is_radical: bool = False

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
class ResidueAtomSpec:
    """Information about one atom in a residue
    ; name  type  charge  chargegroup
    """

    name: str
    type: str
    charge: str
    cgrp: str

    @classmethod
    def from_top_line(cls, l: list[str]):
        return cls(name=l[0], type=l[1], charge=l[2], cgrp=l[3])


@dataclass(order=True)
class ResidueBondSpec:
    """Information about one bond in a residue
    ; atom1 atom2      b0      kb
    """

    atom1: str
    atom2: str
    b0: Optional[str]
    kb: Optional[str]

    @classmethod
    def from_top_line(cls, l: list[str]):
        return cls(
            atom1=l[0], atom2=l[1], b0=field_or_none(l, 2), kb=field_or_none(l, 3)
        )


@dataclass(order=True)
class ResidueImroperSpec:
    """Information about one imroper dihedral in a residue
    ;atom1 atom2 atom3 atom4     q0     cq
    """

    atom1: str
    atom2: str
    atom3: str
    atom4: str
    q0: Optional[str]
    cq: Optional[str]

    @classmethod
    def from_top_line(cls, l: list[str]):
        return cls(
            atom1=l[0],
            atom2=l[1],
            atom3=l[2],
            atom4=l[3],
            q0=field_or_none(l, 4),
            cq=field_or_none(l, 5),
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
class BondType:
    """Information about one bondtype

    A class containing bond information as in the bonds section of the topology.
    From gromacs topology:
    'i', 'j', 'funct', 'c0', 'c1', 'c2', 'c3
    Where i and j are atomtypes
    """

    i: str
    j: str
    funct: str
    c0: Optional[str] = None
    c1: Optional[str] = None
    c2: Optional[str] = None
    c3: Optional[str] = None

    @classmethod
    def from_top_line(cls, l: list[str]):
        return cls(
            i=l[0],
            j=l[1],
            funct=l[2],
            c0=field_or_none(l, 3),
            c1=field_or_none(l, 4),
            c2=field_or_none(l, 5),
            c3=field_or_none(l, 6),
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


@dataclass(order=True)
class Angle:
    """Information about one angle

    A class containing angle information as in the angles section of the topology.

    From gromacs topology:
    ';', 'ai', 'aj', 'ak', 'funct', 'c0', 'c1', 'c2', 'c3'
    With aj < ai < ak
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
class AngleType:
    """Information about one angle

    A class containing angle information as in the angles section of the topology.

    From gromacs topology:
    ';', 'i', 'j', 'k', 'funct', 'c0', 'c1', 'c2', 'c3'
    where i,j,k are atomtypes
    """

    i: str
    j: str
    k: str
    funct: str
    c0: Optional[str] = None
    c1: Optional[str] = None
    c2: Optional[str] = None
    c3: Optional[str] = None

    @classmethod
    def from_top_line(cls, l: list[str]):
        return cls(
            i=l[0],
            j=l[1],
            k=l[2],
            funct=l[3],
            c0=field_or_none(l, 4),
            c1=field_or_none(l, 5),
            c2=field_or_none(l, 6),
            c3=field_or_none(l, 7),
        )


@dataclass(order=True)
class Dihedral:
    """Information about one dihedral

    A class containing bond information as in the dihedrals section of the topology.
    Proper dihedrals have funct 9.
    Improper dihedrals have funct 4.

    From gromacs topology:
    ';', 'ai', 'aj', 'ak', 'al', 'funct', 'c0', 'c1', 'c2', 'c3', 'c4', 'c5'
    For proper dihedrals (funct 9): aj < ak
    For improper dihedrals (funct 4): no guaranteed order
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
class DihedralType:
    """Information about one dihedral

    A class containing bond information as in the dihedrals section of the topology.
    Proper dihedrals have funct 9.
    Improper dihedrals have funct 4.

    From gromacs topology:
    ';', 'i', 'j', 'k', 'l', 'funct', 'c0', 'c1', 'c2', 'c3', 'c4', 'c5'
    Where i,j,k,l are atomtypes
    """

    i: str
    j: str
    k: str
    l: str
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
            i=l[0],
            j=l[1],
            k=l[2],
            l=l[3],
            funct=l[4],
            c0=field_or_none(l, 5),
            c1=field_or_none(l, 6),
            c2=field_or_none(l, 7),
            c3=field_or_none(l, 8),
            c4=field_or_none(l, 9),
            c5=field_or_none(l, 10),
        )


@dataclass(order=True)
class ResidueType:
    """Information about one residuetype"""

    residue: str
    atoms: dict[str, ResidueAtomSpec]
    bonds: dict[tuple[str, str], ResidueBondSpec]
    improper_dihedrals: dict[tuple[str, str, str, str], ResidueImroperSpec]

    @classmethod
    def from_section(cls, residue, d: dict[str, list[list[str]]]):
        atoms = {}
        bonds = {}
        impropers = {}
        if ls := d.get("atoms"):
            for l in ls:
                atom = ResidueAtomSpec.from_top_line(l)
                atoms[atom.name] = atom
        if ls := d.get("bonds"):
            for l in ls:
                bond = ResidueBondSpec.from_top_line(l)
                bonds[(bond.atom1, bond.atom2)] = bond
        if ls := d.get("impropers"):
            for l in ls:
                improper = ResidueImroperSpec.from_top_line(l)
                impropers[
                    (improper.atom1, improper.atom2, improper.atom3, improper.atom4)
                ] = improper

        return cls(residue, atoms, bonds, impropers)


class FF:
    """Conainer for parsed forcefield data."""

    def __init__(self, ffdir: Path):
        self.atomtypes: dict[str, AtomType] = {}
        self.bondtypes: dict[tuple[str, str], BondType] = {}
        self.angletypes: dict[tuple[str, str, str], AngleType] = {}
        self.proper_dihedraltypes: dict[
            tuple[str, str, str, str], list[DihedralType]
        ] = {}
        self.improper_dihedraltypes: dict[tuple[str, str, str, str], DihedralType] = {}
        self.residuetypes: dict[str, ResidueType]

        nonbonded_path = ffdir / "ffnonbonded.itp"
        nonbonded = read_topol(nonbonded_path)
        for l in nonbonded["atomtypes"]:
            if l[0][0] != ";":
                atomtype = AtomType.from_top_line(l)
                self.atomtypes[atomtype.type] = atomtype

        bonded_path = ffdir / "ffbonded.itp"
        bonded = read_topol(bonded_path)
        for l in bonded["bondtypes"]:
            bondtype = BondType.from_top_line(l)
            self.bondtypes[(bondtype.i, bondtype.j)] = bondtype
        for l in bonded["angletypes"]:
            angletype = AngleType.from_top_line(l)
            self.angletypes[(angletype.i, angletype.j, angletype.k)] = angletype
        for l in bonded["dihedraltypes"]:
            dihedraltype = DihedralType.from_top_line(l)
            # proper dihedrals can be defined multiple times
            # with a different phase
            if dihedraltype.funct == "4":
                self.improper_dihedraltypes[
                    (dihedraltype.i, dihedraltype.j, dihedraltype.k, dihedraltype.l)
                ] = dihedraltype
            elif dihedraltype.funct == "9":
                if (
                    self.proper_dihedraltypes.get(
                        (dihedraltype.i, dihedraltype.j, dihedraltype.k, dihedraltype.l)
                    )
                    is None
                ):
                    self.proper_dihedraltypes[
                        (dihedraltype.i, dihedraltype.j, dihedraltype.k, dihedraltype.l)
                    ] = [dihedraltype]
                else:
                    self.proper_dihedraltypes[
                        (dihedraltype.i, dihedraltype.j, dihedraltype.k, dihedraltype.l)
                    ].append(dihedraltype)

            # TODO
            self.residuetypes = {}
            aminoacids_path = ffdir / "aminoacids.rtp"
            aminoacids = read_rtp(aminoacids_path)
            for k, v in aminoacids.items():
                if k.startswith("BLOCK") or k == "bondedtypes":
                    continue
                self.residuetypes[k] = ResidueType.from_section(k, v)

    def __repr__(self) -> str:
        return textwrap.dedent(
            f"""\
        ForceField parameters with
        {len(self.atomtypes)} atomtypes,
        {len(self.bondtypes)} bondtypes,
        {len(self.angletypes)} angletypes,
        {len(self.proper_dihedraltypes)} dihedraltypes
        {len(self.residuetypes)} residuetypes
        """
        )


@dataclass
class ParamPatch:
    offset: Optional[float] = None
    factor: Optional[float] = None
    value: Optional[float] = None

    def update(self, new):
        self.__dict__.update(new)


def props_to_patches(props):
    patches = {}
    for k, v in props.items():
        if "_" not in k:
            continue
        key, suffix = k.split("_")
        if patches.get(key) == None:
            patches[key] = ParamPatch()
        if suffix == "factor":
            patches[key].update({"factor": float(v)})
        elif suffix == "offset":
            patches[key].update({"offset": float(v)})
        elif suffix == "value":
            patches[key].update({"value": float(v)})
    return patches


@dataclass(order=True)
class AtomPatch:
    """Instructions to patch one atom"""

    ai: str
    patches: dict[str, ParamPatch]

    @classmethod
    def from_element(cls, elem: Element):
        props = elem.attrib
        name = props.pop("ai", None)
        if name is None:
            raise ValueError("Atom patch must have an ai attribute")

        patches = props_to_patches(props)

        return cls(name, patches)


@dataclass(order=True)
class BondPatch:
    """Instructions to patch one bond"""

    ai: str
    aj: str
    patches: dict[str, ParamPatch]

    @classmethod
    def from_element(cls, elem: Element):
        props = elem.attrib
        ai = props.pop("ai", None)
        aj = props.pop("aj", None)
        if ai is None or aj is None:
            raise ValueError("Bond patch must have an ai and aj attribute")

        patches = props_to_patches(props)
        return cls(ai, aj, patches)


@dataclass(order=True)
class PairPatch:
    """Instructions to patch one pair"""

    ai: str
    aj: str
    patches: dict[str, ParamPatch]

    @classmethod
    def from_element(cls, elem: Element):
        props = elem.attrib
        ai = props.pop("ai", None)
        aj = props.pop("aj", None)
        if ai is None or aj is None:
            raise ValueError("Pair patch must have an ai and aj attribute")

        patches = props_to_patches(props)
        return cls(ai, aj, patches)


@dataclass(order=True)
class AnglePatch:
    """Instructions to patch one angle"""

    ai: str
    aj: str
    ak: str
    patches: dict[str, ParamPatch]

    @classmethod
    def from_element(cls, elem: Element):
        props = elem.attrib
        ai = props.pop("ai", None)
        aj = props.pop("aj", None)
        ak = props.pop("ak", None)
        if ai is None or aj is None or ak is None:
            raise ValueError("Angle patch must have an ai, aj and ak attribute")

        patches = props_to_patches(props)
        return cls(ai, aj, ak, patches)


@dataclass(order=True)
class DihedralPatch:
    """Instructions to patch one dihedral"""

    ai: str
    aj: str
    ak: str
    al: str
    periodicity: str
    patches: dict[str, ParamPatch]

    @classmethod
    def from_element(cls, elem: Element):
        props = elem.attrib
        ai = props.pop("ai", None)
        aj = props.pop("aj", None)
        ak = props.pop("ak", None)
        al = props.pop("al", None)
        periodicity = props.pop("periodicity", None)
        if ai is None or aj is None or ak is None or al is None or periodicity is None:
            raise ValueError(
                "Angle patch must have an ai, aj, ak, al and periodicity attribute"
            )

        patches = props_to_patches(props)
        return cls(ai, aj, ak, al, periodicity, patches)


class FFPatches:
    """A container for forcefield patches"""
    atompatches: list[AtomPatch]
    bondpatches: list[BondPatch]
    pairpatches: list[PairPatch]
    anglepatches: list[AnglePatch]
    dihedralpatches: list[DihedralPatch]

    def __init__(self, path: Path) -> None:
        xml = read_xml_ff(path)
        self.atompatches = []
        if elems := xml.findall("Atoms/Atom"):
            for elem in elems:
                self.atompatches.append(AtomPatch.from_element(elem))

        self.bondpatches = []
        if elems := xml.findall("Bonds/Bond"):
            for elem in elems:
                self.bondpatches.append(BondPatch.from_element(elem))

        self.pairpatches = []
        if elems := xml.findall("Pairs/Pair"):
            for elem in elems:
                self.pairpatches.append(PairPatch.from_element(elem))

        self.anglepatches = []
        if elems := xml.findall("Angles/Angle"):
            for elem in elems:
                self.anglepatches.append(AnglePatch.from_element(elem))

        # note... <https://manual.gromacs.org/current/reference-manual/functions/bonded-interactions.html#proper-dihedrals-periodic-type>
        # periodicity is also an identifier, not a parameter!
        self.dihedralpatches = []
        if elems := xml.findall("Dihedrals/Dihedral"):
            for elem in elems:
                self.dihedralpatches.append(DihedralPatch.from_element(elem))

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
        self,
        top: TopologyDict,
        ffdir: Optional[Path] = None,
        ffpatch: Optional[Path] = None,
    ) -> None:
        self.top = top
        self.forcefield_directory = ffdir
        self.atoms: dict[str, Atom] = {}
        self.bonds: dict[tuple[str, str], Bond] = {}
        self.pairs: dict[tuple[str, str], Pair] = {}
        self.angles: dict[tuple[str, str, str], Angle] = {}
        self.proper_dihedrals: dict[tuple[str, str, str, str], Dihedral] = {}
        self.improper_dihedrals: dict[tuple[str, str, str, str], Dihedral] = {}

        if ffdir:
            self.ff = FF(ffdir)
        self.ffpatches = None
        if ffpatch:
            self.ffpatches = FFPatches(ffpatch)

        # generate empty Topology if empty TopologyDict
        if self.top == {}:
            return

        self._parse_atoms()
        self._parse_bonds()
        self._parse_pairs()
        self._parse_angles()
        self._parse_dihedrals()

        # self._update_dict()

        self._initialize_graph()

    def _update_dict(self):
        self.top["atoms"] = [attributes_to_list(x) for x in self.atoms.values()]
        self.top["bonds"] = [attributes_to_list(x) for x in self.bonds.values()]
        self.top["pairs"] = [attributes_to_list(x) for x in self.pairs.values()]
        self.top["angles"] = [attributes_to_list(x) for x in self.angles.values()]
        self.top["dihedrals"] = [
            attributes_to_list(x) for x in self.proper_dihedrals.values()
        ]

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
        {len(self.proper_dihedrals)} proper dihedrals
        {len(self.improper_dihedrals)} improper dihedrals
        """
        )

    def __str__(self) -> str:
        return str(self.atoms)

    def _parse_atoms(self):
        ls = self.top["atoms"]
        for l in ls:
            atom = Atom.from_top_line(l)
            self.atoms[atom.nr] = atom

    def _parse_bonds(self):
        ls = self.top["bonds"]
        for l in ls:
            bond = Bond.from_top_line(l)
            self.bonds[(bond.ai, bond.aj)] = bond

    def _parse_pairs(self):
        ls = self.top["pairs"]
        for l in ls:
            pair = Pair.from_top_line(l)
            self.pairs[(pair.ai, pair.aj)] = pair

    def _parse_angles(self):
        ls = self.top["angles"]
        for l in ls:
            angle = Angle.from_top_line(l)
            self.angles[(angle.ai, angle.aj, angle.ak)] = angle

    def _parse_dihedrals(self):
        ls = self.top["dihedrals"]
        for l in ls:
            dihedral = Dihedral.from_top_line(l)
            if dihedral.funct == "4":
                self.improper_dihedrals[
                    (dihedral.ai, dihedral.aj, dihedral.ak, dihedral.al)
                ] = dihedral
            else:
                self.proper_dihedrals[
                    (dihedral.ai, dihedral.aj, dihedral.ak, dihedral.al)
                ] = dihedral

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
        atompair_nrs = tuple(sorted(atompair_nrs, key=int))
        atompair = [self.atoms[atompair_nrs[0]], self.atoms[atompair_nrs[1]]]

        def apply_param(old_value, new_params: ParamPatch):
            pass

        def patch_params(atom: Atom, patch: AtomPatch, spec: AtomType):
            for key, patch in patch.patches.items():
                # check if the key is already set in the topology

                # otherwise get the initial value from the FF

                # then patch it
                # new = apply_param(key, patch)
                pass
            

        # mark atoms as radicals
        for atom in atompair:
            atom.is_radical = True

            # patch parameters
            if self.ffpatches is None or self.ffpatches.atompatches is None: continue
            spec = self.ff.atomtypes.get(atom.type)
            for patch in self.ffpatches.atompatches:
                if patch.ai == atom.type + "_R":
                    patch_params(atom, patch, spec)

        # bonds
        # remove bonds
        removed = self.bonds.pop(atompair_nrs, None)
        logging.info(f"removed bond: {removed}")

        # remove angles
        angle_keys = self._get_atom_angles(atompair_nrs[0]) + self._get_atom_angles(
            atompair_nrs[1]
        )
        for key in angle_keys:
            if all([x in key for x in atompair_nrs]):
                self.angles.pop(key, None)

        # remove proper dihedrals
        # and pairs
        dihedral_keys = self._get_atom_proper_dihedrals(
            atompair_nrs[0]
        ) + self._get_atom_proper_dihedrals(atompair_nrs[1])
        for key in dihedral_keys:
            if all([x in key for x in atompair_nrs]):
                self.proper_dihedrals.pop(key, None)
                pairkey = tuple(sorted((key[0], key[3]), key=int))
                self.pairs.pop(pairkey, None)

        # and improper dihedrals
        dihedral_k_v = self._get_atom_improper_dihedrals(
            atompair_nrs[0]
        ) + self._get_atom_improper_dihedrals(atompair_nrs[1])
        for key, _ in dihedral_k_v:
            if all([x in key for x in atompair_nrs]):
                self.improper_dihedrals.pop(key, None)

        # update bound_to
        try:
            atompair[0].bound_to_nrs.remove(atompair[1].nr)
            atompair[1].bound_to_nrs.remove(atompair[0].nr)
        except ValueError as _:
            m = f"tried to remove bond between already disconnected atoms: {atompair}."
            logging.warning(m)


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

        atompair_nrs = tuple(sorted(atompair_nrs, key=int))
        atompair = [self.atoms[atompair_nrs[0]], self.atoms[atompair_nrs[1]]]

        # de-radialize if re-combining two radicals
        if all(map(lambda x: x.is_radical, atompair)):
            atompair[0].is_radical = False
            atompair[1].is_radical = False

        # update bound_to
        atompair[0].bound_to_nrs.append(atompair[1].nr)
        atompair[1].bound_to_nrs.append(atompair[0].nr)

        # bonds
        # add bond
        bond = Bond(atompair_nrs[0], atompair_nrs[1], "1")
        self.bonds[atompair_nrs] = bond
        logging.info(f"added bond: {bond}")

        # add angles
        angle_keys = self._get_atom_angles(atompair_nrs[0]) + self._get_atom_angles(
            atompair_nrs[1]
        )
        for key in angle_keys:
            if self.angles.get(key) is None:
                self.angles[key] = Angle(key[0], key[1], key[2], "1")

        # add proper and improper dihedrals
        # add proper dihedrals and pairs
        dihedral_keys = self._get_atom_proper_dihedrals(
            atompair_nrs[0]
        ) + self._get_atom_proper_dihedrals(atompair_nrs[1])
        for key in dihedral_keys:
            if self.proper_dihedrals.get(key) is None:
                self.proper_dihedrals[key] = Dihedral(
                    key[0], key[1], key[2], key[3], "9"
                )
            pairkey = tuple(str(x) for x in sorted([key[0], key[3]], key=int))
            if self.pairs.get(pairkey) is None:
                self.pairs[pairkey] = Pair(pairkey[0], pairkey[1], "1")

        # improper dihedral
        dihedral_k_v = self._get_atom_improper_dihedrals(
            atompair_nrs[0]
        ) + self._get_atom_improper_dihedrals(atompair_nrs[1])
        for key, value in dihedral_k_v:
            if self.improper_dihedrals.get(key) is None:
                # TODO: fix this after the demonstration
                c2 = None
                if value.q0 is not None:
                    c2 = "1"
                self.improper_dihedrals[key] = Dihedral(
                    key[0], key[1], key[2], key[3], "4", value.q0, value.cq, c2
                )

        # if there are no changed parameters for radicals, exit here
        if self.ffpatches is None:
            self._update_dict()
            return

    def _get_atom_bonds(self, atom_nr: str) -> list[tuple[str, str]]:
        ai = atom_nr
        bonds = []
        for aj in self.atoms[ai].bound_to_nrs:
            # TODO: maybe sort here and filter later instead
            if int(ai) < int(aj):
                bonds.append((ai, aj))
        return bonds

    def _get_atom_pairs(self, _: str) -> list[tuple[str, str]]:
        raise NotImplementedError(
            "get_atom_pairs is not implementes. Get the pairs as the endpoints of dihedrals instead."
        )

    def _get_atom_angles(self, atom_nr: str) -> list[tuple[str, str, str]]:
        """
        each atom has a list of atoms it is bound to
        get a list of angles that one atom is involved in
        based in these lists.
        Angles between atoms ai, aj, ak satisfy ai < ak
        """
        return self._get_center_atom_angles(atom_nr) + self._get_margin_atom_angles(
            atom_nr
        )

    def _get_center_atom_angles(self, atom_nr: str) -> list[tuple[str, str, str]]:
        # atom_nr in the middle of an angle
        angles = []
        aj = atom_nr
        partners = self.atoms[aj].bound_to_nrs
        for ai in partners:
            for ak in partners:
                if int(ai) < int(ak):
                    angles.append((ai, aj, ak))
        return angles

    def _get_margin_atom_angles(self, atom_nr: str) -> list[tuple[str, str, str]]:
        # atom_nr at the outer corner of angle
        angles = []
        ai = atom_nr
        for aj in self.atoms[ai].bound_to_nrs:
            for ak in self.atoms[aj].bound_to_nrs:
                if ai == ak:
                    continue
                if int(ai) < int(ak):
                    angles.append((ai, aj, ak))
                else:
                    angles.append((ak, aj, ai))
        return angles

    def _get_atom_proper_dihedrals(
        self, atom_nr: str
    ) -> list[tuple[str, str, str, str]]:
        """
        each atom has a list of atoms it is bound to.
        get a list of dihedrals that one atom is involved in
        based in these lists.
        """
        return self._get_center_atom_dihedrals(
            atom_nr
        ) + self._get_margin_atom_dihedrals(atom_nr)

    def _get_center_atom_dihedrals(
        self, atom_nr: str
    ) -> list[tuple[str, str, str, str]]:
        dihedrals = []
        aj = atom_nr
        partners = self.atoms[aj].bound_to_nrs
        for ai in partners:
            for ak in partners:
                if ai == ak:
                    continue
                for al in self.atoms[ak].bound_to_nrs:
                    if al == ak or aj == al:
                        continue
                    if int(aj) < int(ak):
                        dihedrals.append((ai, aj, ak, al))
                    else:
                        dihedrals.append((al, ak, aj, ai))
        return dihedrals

    def _get_margin_atom_dihedrals(
        self, atom_nr: str
    ) -> list[tuple[str, str, str, str]]:
        dihedrals = []
        ai = atom_nr
        for aj in self.atoms[ai].bound_to_nrs:
            for ak in self.atoms[aj].bound_to_nrs:
                if ai == ak:
                    continue
                for al in self.atoms[ak].bound_to_nrs:
                    if al == ak or aj == al:
                        continue
                    if int(aj) < int(ak):
                        dihedrals.append((ai, aj, ak, al))
                    else:
                        dihedrals.append((al, ak, aj, ai))
        return dihedrals

    def _get_atom_improper_dihedrals(self, atom_nr: str) -> list[tuple]:
        # TODO: cleanup and make more efficient
        # which improper dihedrals are used is defined for each residue
        # in aminoacids.rtp
        # get improper diheldrals from FF based on residue
        # TODO: handle impropers defined for the residue that
        # belongs to an adjacent atom, not just the the specied one
        atom = self.atoms[atom_nr]
        residue = self.ff.residuetypes.get(atom.residue)
        if residue is None:
            return []

        # <https://manual.gromacs.org/current/reference-manual/functions/bonded-interactions.html#improper-dihedrals>
        # atom in a line, like a regular dihedral:
        dihedrals = []
        dihedral_candidate_keys = self._get_margin_atom_dihedrals(
            atom_nr
        ) + self._get_center_atom_dihedrals(atom_nr)

        # atom in the center of a star/tetrahedron:
        ai = atom_nr
        partners = self.atoms[ai].bound_to_nrs
        if len(partners) >= 3:
            combs = combinations(partners, 3)
            for comb in combs:
                aj, ak, al = comb
                dihedral_candidate_keys.append((ai, aj, ak, al))

        # atom in corner of a star/tetrahedron:
        for a in self.atoms[atom_nr].bound_to_nrs:
            partners = self.atoms[a].bound_to_nrs
            if len(partners) >= 3:
                combs = permutations(partners + [a], 4)
                for comb in combs:
                    ai, aj, ak, al = comb
                    dihedral_candidate_keys.append((ai, aj, ak, al))

        # residues on aminoacids.rtp specify a dihedral to the next or previous
        # AA with -C and +N as the atomname
        for candidate in dihedral_candidate_keys:
            candidate_key = [self.atoms[atom_nr].atom for atom_nr in candidate]
            for i, nr in enumerate(candidate):
                if self.atoms[nr].resnr != atom.resnr:
                    if candidate_key[i] == "C":
                        candidate_key[i] = "-C"
                    elif candidate_key[i] == "N":
                        candidate_key[i] = "+N"

            candidate_key = tuple(candidate_key)
            dihedral = residue.improper_dihedrals.get(candidate_key)
            if dihedral:
                dihedrals.append((candidate, dihedral))

        return dihedrals

    def _patch_parameters(self, atompair: list[Atom]):
        if self.ffpatches is None:
            return
        # Adjust parameters based on patch
        # TODO: move to the section where the bonds are broken and made
        # to reduce iterations
        atom_nrs = [a.nr for a in atompair]

        # atoms
        if atompatches := self.ffpatches.atompatches:
            for atom in atompair:
                spec = self.ff.atomtypes.get(atom.type)
                print(spec)

        # get (unbroken) bonds that the now radicals in the atompair are still involved in
        if bondpatches := self.ffpatches.bondpatches:
            pass

        # get (unbroken) angles that the now radicals in the atompair are still involved in
        if anglepatches := self.ffpatches.anglepatches:
            pass

        if pairpatches := self.ffpatches.pairpatches:
            pass

        # get (unbroken) angles that the now radicals in the atompair are still involved in
        if dihedralpatches := self.ffpatches.dihedralpatches:
            pass


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
            pattern = value.replace("*", r".*").replace("+", r".+")
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
    raise NotImplementedError("WIP")


def get_by_permutations(d: dict, key) -> Optional[Any]:
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


def generate_topology_from_bound_to(
    atoms: list[Atom], ffdir: Path, ffpatch: Path
) -> Topology:
    top = Topology({}, ffdir, ffpatch)
    for atom in atoms:
        top.atoms[atom.nr] = atom

    # bonds
    keys = []
    for atom in top.atoms.values():
        keys = top._get_atom_bonds(atom.nr)
        for key in keys:
            top.bonds[key] = Bond(key[0], key[1], "1")

    # angles
    for atom in top.atoms.values():
        keys = top._get_atom_angles(atom.nr)
        for key in keys:
            top.angles[key] = Angle(key[0], key[1], key[2], "1")

    # dihedrals and pass
    for atom in top.atoms.values():
        keys = top._get_atom_proper_dihedrals(atom.nr)
        for key in keys:
            top.proper_dihedrals[key] = Dihedral(key[0], key[1], key[2], key[3], "9")
            pairkey = tuple(str(x) for x in sorted([key[0], key[3]], key=int))
            if top.pairs.get(pairkey) is None:
                top.pairs[pairkey] = Pair(pairkey[0], pairkey[1], "1")

    for atom in top.atoms.values():
        impropers = top._get_atom_improper_dihedrals(atom.nr)
        for key, improper in impropers:
            top.improper_dihedrals[key] = Dihedral(improper.atom1, improper.atom2, improper.atom3, improper.atom4, "4", improper.cq)

    return top


def match_typestring_to_patch(s: str, ps: list[Patch]):
    s = "CT_T"
    patch_ids = [p.id for p in ps]
    print(patch_ids)
    return ps[0]
