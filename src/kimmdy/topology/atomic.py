from dataclasses import dataclass, field
from typing import Optional, Union
from kimmdy.topology.utils import field_or_none


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

    def radical_type(self):
        if self.is_radical:
            return self.type + "_R"
        else:
            return self.type

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
    id_sym: str
    at_num: str
    charge: str
    mass: str
    ptype: str
    sigma: str
    epsilon: str
    id: str

    @classmethod
    def from_top_line(cls, l: list[str]):
        return cls(
            type=l[0],
            id_sym=l[0],
            at_num=l[1],
            charge=l[2],
            mass=l[3],
            ptype=l[4],
            sigma=l[5],
            epsilon=l[6],
            id=l[0],
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
    id: str
    id_sym: str
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
            id="---".join(l[:2]),
            id_sym="---".join(reversed(l[:2])),
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
    id: str
    id_sym: str
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
            id="---".join(l[:3]),
            id_sym="---".join(reversed(l[:3])),
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
    
    Note that proper dihedrals of type 9 can be defined multiple times, for different
    periodicities. This is why would-be parameter c2 is called periodicity and part of
    the `id`.

    From gromacs topology:
    ';', 'i', 'j', 'k', 'l', 'funct', 'c0', 'c1', 'periodicity', 'c3', 'c4', 'c5'
    Where i,j,k,l are atomtypes
    """

    i: str
    j: str
    k: str
    l: str
    id: str
    id_sym: str
    funct: str
    periodicity: str
    c0: Optional[str] = None
    c1: Optional[str] = None
    c3: Optional[str] = None
    c4: Optional[str] = None
    c5: Optional[str] = None

    @classmethod
    def from_top_line(cls, l: list[str]):
        periodicity = field_or_none(l, 7)
        if periodicity is None:
            periodicity = '2'
        return cls(
            i=l[0],
            j=l[1],
            k=l[2],
            l=l[3],
            id="---".join(l[:4]) + ':::' + periodicity,
            id_sym="---".join(reversed(l[:4])) + ':::' + periodicity,
            funct=l[4],
            periodicity=periodicity,
            c0=field_or_none(l, 5),
            c1=field_or_none(l, 6),
            c3=field_or_none(l, 8),
            c4=field_or_none(l, 9),
            c5=field_or_none(l, 10),
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
    b0: Optional[str] = None
    kb: Optional[str] = None

    @classmethod
    def from_top_line(cls, l: list[str]):
        return cls(
            atom1=l[0], atom2=l[1], b0=field_or_none(l, 2), kb=field_or_none(l, 3)
        )


@dataclass(order=True)
class ResidueImproperSpec:
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
class ResidueProperSpec:
    """Information about one imroper dihedral in a residue
    ;atom1 atom2 atom3 atom4     q0     cq
    """

    atom1: str
    atom2: str
    atom3: str
    atom4: str
    q0: Optional[str]

    @classmethod
    def from_top_line(cls, l: list[str]):
        return cls(
            atom1=l[0],
            atom2=l[1],
            atom3=l[2],
            atom4=l[3],
            q0=field_or_none(l, 4),
        )

@dataclass(order=True)
class ResidueType:
    """Information about one residuetype"""

    residue: str
    atoms: dict[str, ResidueAtomSpec]
    bonds: dict[tuple[str, str], ResidueBondSpec]
    proper_dihedrals: dict[tuple[str, str, str, str], ResidueProperSpec]
    improper_dihedrals: dict[tuple[str, str, str, str], ResidueImproperSpec]

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
                improper = ResidueImproperSpec.from_top_line(l)
                impropers[
                    (improper.atom1, improper.atom2, improper.atom3, improper.atom4)
                ] = improper

        return cls(residue, atoms, bonds, impropers)


AtomId = str
BondId = tuple[str, str]
AngleId = tuple[str, str, str]
DihedralId = tuple[str, str, str, str, str]
Atomic = Union[Atom, Bond, Pair, Angle, Dihedral]
AtomicType = Union[AtomType, BondType, AngleType, DihedralType]
AtomicTypes = Union[
    dict[AtomId, AtomType],
    dict[BondId, BondType],
    dict[AngleId, AngleType],
    dict[DihedralId, DihedralType],
]
