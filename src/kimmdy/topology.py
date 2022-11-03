from dataclasses import field, dataclass
from pathlib import Path
from typing import Optional
from kimmdy.parsing import TopologyDict
from itertools import takewhile


@dataclass
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


@dataclass
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



@dataclass
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


@dataclass
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


@dataclass
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
    def __init__(self, top: TopologyDict, ffdir: Path) -> None:
        self.top = top
        self.forcefield_directory = ffdir
        self.ff = {}

        self.atoms = []
        self.bonds = []
        self.dihedrals = []
        self.pairs = []
        self.angles = []
        self.proper_dihedrals = []
        self.improper_dihedrals = []

        self._get_atoms()
        self._get_bonds()
        self._get_pairs()
        self._get_angles()
        self._get_dihedrals()

        self._initialize_graph()

    def to_dict(self) -> TopologyDict:
        return self.top

    def __repr__(self) -> str:
        return self.top.__repr__()

    def __str__(self) -> str:
        return self.top.__str__()

    def _get_atoms(self):
        ls = self.top["atoms"]
        for l in ls:
            if l[0] != ';':
                self.atoms.append(Atom.from_top_line(l))

    def _get_bonds(self):
        ls = self.top["bonds"]
        for l in ls:
            if l[0] != ';':
                self.bonds.append(Bond.from_top_line(l))

    def _get_pairs(self):
        ls = self.top["pairs"]
        for l in ls:
            if l[0] != ';':
                self.pairs.append(Pair.from_top_line(l))
                
    def _get_angles(self):
        ls = self.top["angles"]
        for l in ls:
            if l[0] != ';':
                self.angles.append(Angle.from_top_line(l))

    def _get_dihedrals(self):
        ls = self.top["dihedrals"]
        for l in ls:
            if l[0] != ';':
                dihedral = Dihedral.from_top_line(l)
                self.dihedrals.append(dihedral)
                if dihedral.funct == '9': self.proper_dihedrals.append(dihedral)
                if dihedral.funct == '4': self.improper_dihedrals.append(dihedral)


    def _initialize_graph(self):
        pass


def field_or_none(l: list[str], i) -> Optional[str]:
    try:
        return l[i]
    except IndexError as _:
        return None


def is_not_comment(c: str):
    return c != ";"

