"""
Atomic datatypes for the topology such as Atom, Bond, Angle, Dihedral, etc.
The order of the fields comes from the gromacs topology file format.
See [gromacs manual](https://manual.gromacs.org/current/reference-manual/topologies/topology-file-formats.html#topology-file)
"""

from dataclasses import dataclass, field
from typing import Optional, Union

from kimmdy.constants import FFFUNC
from kimmdy.utils import field_or_none


@dataclass
class MoleculeTypeHeader:
    name: str
    nrexcl: str


@dataclass()
class Atom:
    """Information about one atom

    A class containing atom information as in the atoms section of the topology.
    An atom keeps a list of which atoms it is bound to and its radical state.

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
    mass: Optional[str] = None
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
            mass=field_or_none(l, 7),
            typeB=field_or_none(l, 8),
            chargeB=field_or_none(l, 9),
            massB=field_or_none(l, 10),
        )


@dataclass()
class PositionRestraint:
    """Information about one position restraint.

    A class containing information as in the position_restraints section of the topology.

    From gromacs topology:
    ; ai   funct    fc(x,y,z)
    """

    ai: str
    funct: str
    fc: tuple[str, str, str]
    condition: Optional[str] = None

    @classmethod
    def from_top_line(cls, l: list[str], condition=None):
        return cls(
            ai=l[0],
            funct=l[1],
            fc=(l[2], l[3], l[4]),
            condition=condition,
        )


# [ position_restraints ]
# ; you wouldn't normally use this for a molecule like Urea,
# ; but we include it here for didactic purposes
# ; ai   funct    fc
#    1     1     1000    1000    1000 ; Restrain to a point
#    2     1     1000       0    1000 ; Restrain to a line (Y-axis)
#    3     1     1000       0       0 ; Restrain to a plane (Y-Z-plane)
# [ dihedral_restraints ]
# ; ai   aj    ak    al  type  phi  dphi  fc
#     3    6     1    2     1  180     0  10
#     1    4     3    5     1  180     0  10


@dataclass()
class DihedralRestraint:
    """Information about one dihedral restraint.

    A class containing information as in the dihedral_restraints section of the topology.

    From gromacs topology:
    ; ai   aj    ak    al  type  phi  dphi  fc
    """

    ai: str
    aj: str
    ak: str
    al: str
    type: str
    phi: str
    dphi: str
    fc: str

    @classmethod
    def from_top_line(cls, l: list[str]):
        return cls(
            ai=l[0],
            aj=l[1],
            ak=l[2],
            al=l[3],
            type=l[4],
            phi=l[5],
            dphi=l[6],
            fc=l[7],
        )


@dataclass()
class Settle:
    """Information about one settles

    A class containing atom information as in the settle section of the topology.

    From gromacs topology:
    ; nr funct doh dhh
    """

    nr: str
    funct: str
    doh: str
    dhh: str

    @classmethod
    def from_top_line(cls, l: list[str]):
        return cls(
            nr=l[0],
            funct=l[1],
            doh=l[2],
            dhh=l[3],
        )


@dataclass()
class Exclusion:
    """Information about one exclusion

    A class containing atom information as in the exclusions section of the topology.

    It's unlikey we need this many atomnumbers in a single exclusion, but just in case.
    Because the gromacs manuals just says

    > Each line should start with one atom index, followed by one or more atom indices.
    > All non-bonded interactions between the first atom and the other atoms will be excluded.
    > -- https://manual.gromacs.org/current/reference-manual/topologies/molecule-definition.html#exclusions

    From gromacs topology:
    ; ai aj ak al am an ao ap
    """

    ai: str
    aj: Optional[str] = None
    ak: Optional[str] = None
    al: Optional[str] = None
    am: Optional[str] = None
    an: Optional[str] = None
    ao: Optional[str] = None
    ap: Optional[str] = None

    @classmethod
    def from_top_line(cls, l: list[str]):
        return cls(
            ai=l[0],
            aj=field_or_none(l, 1),
            ak=field_or_none(l, 2),
            al=field_or_none(l, 3),
            am=field_or_none(l, 4),
            an=field_or_none(l, 5),
            ao=field_or_none(l, 6),
            ap=field_or_none(l, 7),
        )

    def key(self):
        return tuple(
            field
            for field in [self.ai, self.aj, self.ak, self.al, self.am, self.an, self.ao]
            if field is not None
        )


@dataclass()
class AtomType:
    """Information about one atom type

    A class containing atom type information as in the atomtypes section of the forcefield.

    From gromacs version of the amber* ff:
    ; name at.num mass charge ptype sigma epsilon
    """

    type: str
    id: str
    id_sym: str
    at_num: str
    mass: str
    charge: str
    ptype: str
    sigma: str
    epsilon: str

    @classmethod
    def from_top_line(cls, l: list[str]):
        length = len(l)
        has_at_num = length >= 7
        if has_at_num:
            offset = 1
        else:
            offset = 0
        at_num = l[1] if has_at_num else ""
        return cls(
            type=l[0],
            id=l[0],
            id_sym=l[0],
            at_num=at_num,
            mass=l[1 + offset],
            charge=l[2 + offset],
            ptype=l[3 + offset],
            sigma=l[4 + offset],
            epsilon=l[5 + offset],
        )


@dataclass()
class Bond:
    """Information about one bond

    A class containing bond information as in the bonds section of the topology.

    From gromacs topology:
    ; ai aj funct c0 c1 c2 c3 c4 c5
    With ai < aj
    """

    ai: str
    aj: str
    funct: str
    c0: Optional[str] = None
    c1: Optional[str] = None
    c2: Optional[str] = None
    c3: Optional[str] = None
    c4: Optional[str] = None
    c5: Optional[str] = None

    @classmethod
    def from_top_line(cls, l: list[str]):
        funct = field_or_none(l, 2)
        if funct is None:
            funct = FFFUNC["harmonic_bond"]
        return cls(
            ai=l[0],
            aj=l[1],
            funct=funct,
            c0=field_or_none(l, 3),
            c1=field_or_none(l, 4),
            c2=field_or_none(l, 5),
            c3=field_or_none(l, 6),
            c4=field_or_none(l, 7),
            c5=field_or_none(l, 8),
        )


@dataclass()
class BondType:
    """Information about one bondtype

    A class containing bond information as in the bondtype section of the forcefield.

    From gromacs version of the amber* ff:
    ; i j func b0 kb
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
        funct = field_or_none(l, 2)
        if funct is None:
            funct = FFFUNC["harmonic_bond"]
        return cls(
            i=l[0],
            j=l[1],
            id="---".join(l[:2]),
            id_sym="---".join(reversed(l[:2])),
            funct=funct,
            c0=field_or_none(l, 3),
            c1=field_or_none(l, 4),
            c2=field_or_none(l, 5),
            c3=field_or_none(l, 6),
        )


@dataclass()
class NonbondParamType:
    """Information about one nonbonded parameterize

    typical in coarse grained models.
    A class containing nonbonded information as in the nonbond_params section of the forcefield.

    From gromacs:
    ; Lennard jones between beads
    ; i j funda sigma(nm) epsilon (kmol/mol)
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
        )


@dataclass()
class Pair:
    """Information about one pair

    A class containing pair information as in the pair section of the topology.

    From gromacs topology:
    ; ai aj funct c0 c1 c2 c3
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


@dataclass()
class Angle:
    """Information about one angle

    A class containing angle information as in the angles section of the topology.

    From gromacs topology:
    ; ai aj ak funct c0 c1 c2 c3
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


@dataclass()
class AngleType:
    """Information about one angletype

    A class containing angle type information as in the angletypes section of the forcefield.

    From gromacs version of the amber* ff:
    ; i j k func th0 cth
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


@dataclass()
class Dihedral:
    """Information about one proper or improper dihedral

    A class containing dihedral information as in the dihedrals section of the topology.
    Improper dihedrals have funct 4.
    Proper dihedrals have funct != 4, mostly funct 9.

    Note that proper dihedrals of type 9 can be defined multiple times, for different
    periodicities. This is why would-be parameter c2 is called periodicity.

    From gromacs topology:
    ; ai aj ak al funct c0 c1 c2 c3 c4 c5
    For proper dihedrals (funct 9): aj < ak
    """

    ai: str
    aj: str
    ak: str
    al: str
    funct: str
    c0: Optional[str] = None
    c1: Optional[str] = None
    periodicity: str = ""
    c3: Optional[str] = None
    c4: Optional[str] = None
    c5: Optional[str] = None

    @classmethod
    def from_top_line(cls, l: list[str]):
        periodicity = field_or_none(l, 7)
        if periodicity is None:
            periodicity = ""
        return cls(
            ai=l[0],
            aj=l[1],
            ak=l[2],
            al=l[3],
            funct=l[4],
            c0=field_or_none(l, 5),
            c1=field_or_none(l, 6),
            periodicity=periodicity,
            c3=field_or_none(l, 8),
            c4=field_or_none(l, 9),
            c5=field_or_none(l, 10),
        )


@dataclass
class MultipleDihedrals:
    """
    Multiple ``Dihedral``s with the same ai, aj, ak, al
    but different periodicities.
    funct should always be "9" when the length of dihedrals is > 1.
    The key of the dihedrals dict is the periodicity (c2).
    """

    ai: str
    aj: str
    ak: str
    al: str
    funct: str
    dihedrals: dict[str, Dihedral]


@dataclass()
class DihedralType:
    """Information about one dihedraltype

    A class containing dihedral type information as in the dihedraltypes
    section of the forcefield.
    Improper dihedrals have funct 4. Proper dihedrals have funct 9.

    Note that proper dihedrals of type 9 can be defined multiple times, for different
    periodicities. This is why would-be parameter c2 is called periodicity and part of
    the `id`.

    From gromacs version of the amber* ff:
    ; i j k l func phase kd pn
    """

    i: str
    j: str
    k: str
    l: str
    id: str
    id_sym: str
    funct: str
    c0: str
    c1: str
    periodicity: str
    c3: Optional[str] = None
    c4: Optional[str] = None
    c5: Optional[str] = None

    @classmethod
    def from_top_line(cls, l: list[str]):
        periodicity = field_or_none(l, 7)
        if periodicity is None:
            periodicity = "2"
        return cls(
            i=l[0],
            j=l[1],
            k=l[2],
            l=l[3],
            id="---".join(l[:4]) + ":::" + periodicity,
            id_sym="---".join(reversed(l[:4])) + ":::" + periodicity,
            funct=l[4],
            c0=l[5],
            c1=l[6],
            periodicity=periodicity,
            c3=field_or_none(l, 8),
            c4=field_or_none(l, 9),
            c5=field_or_none(l, 10),
        )


@dataclass()
class MultipleDihedralTypes:
    """
    Multiple ``DihedralTypes``s with the same ai, aj, ak, al
    but different periodicities.
    funct should always be "9" when the length of dihedrals is > 1.
    The key of the dihedral_types dict is the periodicity (c2).
    """

    ai: str
    aj: str
    ak: str
    al: str
    funct: str
    dihedral_types: dict[str, DihedralType]


@dataclass()
class ResidueAtomSpec:
    """Information about one atom in a residue

    ; name type charge chargegroup
    """

    name: str
    type: str
    charge: str
    cgrp: str

    @classmethod
    def from_top_line(cls, l: list[str]):
        return cls(name=l[0], type=l[1], charge=l[2], cgrp=l[3])


@dataclass()
class ResidueBondSpec:
    """Information about one bond in a residue

    ; atom1 atom2 b0 kb
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


@dataclass()
class ResidueImproperSpec:
    """Information about one improper dihedral in a residue

    ; atom1 atom2 atom3 atom4 c0(q0) c1(cp) c2(mult)
    """

    atom1: str
    atom2: str
    atom3: str
    atom4: str
    c0: Optional[str]
    c1: Optional[str]
    c2: Optional[str]

    @classmethod
    def from_top_line(cls, l: list[str]):
        return cls(
            atom1=l[0],
            atom2=l[1],
            atom3=l[2],
            atom4=l[3],
            c0=field_or_none(l, 4),
            c1=field_or_none(l, 5),
            c2=field_or_none(l, 6),
        )


@dataclass()
class ResidueProperSpec:
    """Information about one proper dihedral in a residue

    ; atom1 atom2 atom3 atom4 c0(q0) c1(cq) c2
    """

    atom1: str
    atom2: str
    atom3: str
    atom4: str
    c0: Optional[str]
    c1: Optional[str]
    c2: Optional[str]

    @classmethod
    def from_top_line(cls, l: list[str]):
        return cls(
            atom1=l[0],
            atom2=l[1],
            atom3=l[2],
            atom4=l[3],
            c0=field_or_none(l, 4),
            c1=field_or_none(l, 5),
            c2=field_or_none(l, 6),
        )


@dataclass()
class ResidueType:
    """Information about one residuetype from aminoacids.rtp"""

    residue: str
    atoms: dict[str, ResidueAtomSpec]
    bonds: dict[tuple[str, str], ResidueBondSpec]
    proper_dihedrals: dict[tuple[str, str, str, str], ResidueProperSpec]
    improper_dihedrals: dict[tuple[str, str, str, str], ResidueImproperSpec]

    @classmethod
    def from_section(cls, residue, d: dict[str, dict[str, list[list[str]]]]):
        atoms = {}
        bonds = {}
        propers = {}
        impropers = {}
        if ls := d.get("atoms"):
            for l in ls["content"]:
                atom = ResidueAtomSpec.from_top_line(l)
                atoms[atom.name] = atom
        if ls := d.get("bonds"):
            for l in ls["content"]:
                bond = ResidueBondSpec.from_top_line(l)
                bonds[(bond.atom1, bond.atom2)] = bond
        if ls := d.get("dihedrals"):
            for l in ls["content"]:
                proper = ResidueProperSpec.from_top_line(l)
                propers[(proper.atom1, proper.atom2, proper.atom3, proper.atom4)] = (
                    proper
                )
        if ls := d.get("impropers"):
            for l in ls["content"]:
                improper = ResidueImproperSpec.from_top_line(l)
                impropers[
                    (improper.atom1, improper.atom2, improper.atom3, improper.atom4)
                ] = improper

        return cls(residue, atoms, bonds, propers, impropers)


AtomId = str
BondId = tuple[str, str]
AngleId = tuple[str, str, str]
ProperDihedralId = tuple[str, str, str, str, str]
ImproperDihedralId = tuple[str, str, str, str]
Atomic = Union[Atom, Bond, Pair, Angle, Dihedral]
AtomicType = Union[AtomType, BondType, AngleType, DihedralType]
AtomicTypes = Union[
    dict[AtomId, AtomType],
    dict[BondId, BondType],
    dict[AngleId, AngleType],
    dict[ProperDihedralId, DihedralType],
    dict[ImproperDihedralId, DihedralType],
]
InteractionIds = Union[BondId, AngleId, ProperDihedralId, ImproperDihedralId]
Interaction = Union[Bond, Pair, Angle, Dihedral]
Interactions = Union[
    dict[BondId, Bond],
    dict[AngleId, Angle],
    dict[ProperDihedralId, Dihedral],
    dict[ImproperDihedralId, Dihedral],
]
InteractionType = Union[BondType, AngleType, DihedralType]
InteractionTypes = Union[
    dict[BondId, BondType],
    dict[AngleId, AngleType],
    dict[ProperDihedralId, DihedralType],
    dict[ImproperDihedralId, DihedralType],
]
