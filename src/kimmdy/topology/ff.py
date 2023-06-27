from __future__ import annotations
import textwrap
from pathlib import Path
from xml.etree.ElementTree import Element
from kimmdy.topology.atomic import *
from kimmdy.parsing import read_top, read_xml_ff, read_rtp
from typing import Union

from kimmdy.topology.utils import get_top_section


class FF:
    """Conainer for parsed forcefield data."""

    def __init__(self, top):
        self.atomtypes: dict[str, AtomType] = {}
        self.bondtypes: dict[tuple[str, str], BondType] = {}
        self.angletypes: dict[tuple[str, str, str], AngleType] = {}
        self.proper_dihedraltypes: dict[
            tuple[str, str, str, str, str], DihedralType
        ] = {}
        self.improper_dihedraltypes: dict[tuple[str, str, str, str], DihedralType] = {}
        self.residuetypes: dict[str, ResidueType]

        ffdir = top['ffdir']

        atomtypes = get_top_section(top, "atomtypes")
        if atomtypes is None:
            raise ValueError("atomtypes not found in top file")
        for l in atomtypes:
            atomtype = AtomType.from_top_line(l)
            self.atomtypes[atomtype.type] = atomtype

        bondtypes = get_top_section(top, "bondtypes")
        if bondtypes is None:
            raise ValueError("bondtypes not found in top file")
        for l in bondtypes:
            bondtype = BondType.from_top_line(l)
            self.bondtypes[(bondtype.i, bondtype.j)] = bondtype

        angletypes = get_top_section(top, "angletypes")
        if angletypes is None:
            raise ValueError("angletypes not found in top file")
        for l in angletypes:
            angletype = AngleType.from_top_line(l)
            self.angletypes[(angletype.i, angletype.j, angletype.k)] = angletype

        dihedraltypes = get_top_section(top, "dihedraltypes")
        if dihedraltypes is None:
            raise ValueError("dihedraltypes not found in top file")
        for l in dihedraltypes:
            dihedraltype = DihedralType.from_top_line(l)
            # proper dihedrals can be defined multiple times
            # with a different phase
            if dihedraltype.funct == "4":
                self.improper_dihedraltypes[
                    (dihedraltype.i, dihedraltype.j, dihedraltype.k, dihedraltype.l)
                ] = dihedraltype
            else:
                # e.g. proper dihedrals with dihedraltype.funct == "9":
                if (
                    self.proper_dihedraltypes.get(
                        (
                            dihedraltype.i,
                            dihedraltype.j,
                            dihedraltype.k,
                            dihedraltype.l,
                            dihedraltype.periodicity,
                        )
                    )
                    is None
                ):
                    self.proper_dihedraltypes[
                        (
                            dihedraltype.i,
                            dihedraltype.j,
                            dihedraltype.k,
                            dihedraltype.l,
                            dihedraltype.periodicity,
                        )
                    ] = dihedraltype

            self.residuetypes = {}
            if ffdir is None:
                return
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
        {len(self.proper_dihedraltypes)} proper dihedraltypes
        {len(self.improper_dihedraltypes)} improper dihedraltypes
        {len(self.residuetypes)} residuetypes
        """
        )


@dataclass
class ParamPatch:
    value: Optional[float] = None
    offset: Optional[float] = None
    factor: Optional[float] = None

    def update(self, new):
        self.__dict__.update(new)

    def apply(self, initial: float):
        result = initial
        if self.value is not None:
            result = self.value
        if self.offset is not None:
            result += self.offset
        if self.factor is not None:
            result *= self.factor
        return result


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
    id: str
    params: dict[str, ParamPatch]

    def __init__(self, elem: Element):
        props = elem.attrib
        ai = props.pop("ai", None)
        if ai is None:
            raise ValueError("Atom patch must have an ai attribute")

        self.ai = ai
        self.params = props_to_patches(props)
        self.id = self.ai
        self.id_sym = self.id


@dataclass(order=True)
class BondPatch:
    """Instructions to patch one bond"""

    ai: str
    aj: str
    id: str
    params: dict[str, ParamPatch]

    def __init__(self, elem: Element):
        props = elem.attrib
        ai = props.pop("ai", None)
        aj = props.pop("aj", None)
        if ai is None or aj is None:
            raise ValueError("Bond patch must have an ai and aj attribute")

        self.params = props_to_patches(props)
        self.ai = ai
        self.aj = aj
        self.id = "---".join([ai, aj])
        self.id_sym = "---".join(reversed([ai, aj]))


@dataclass(order=True)
class PairPatch:
    """Instructions to patch one pair"""

    ai: str
    aj: str
    id: str
    params: dict[str, ParamPatch]

    def __init__(self, elem: Element):
        props = elem.attrib
        ai = props.pop("ai", None)
        aj = props.pop("aj", None)
        if ai is None or aj is None:
            raise ValueError("Pair patch must have an ai and aj attribute")

        self.params = props_to_patches(props)
        self.ai = ai
        self.aj = aj
        self.id = "---".join([ai, aj])
        self.id_sym = "---".join(reversed([ai, aj]))


@dataclass(order=True)
class AnglePatch:
    """Instructions to patch one angle"""

    ai: str
    aj: str
    ak: str
    id: str
    params: dict[str, ParamPatch]

    def __init__(self, elem: Element):
        props = elem.attrib
        ai = props.pop("ai", None)
        aj = props.pop("aj", None)
        ak = props.pop("ak", None)
        if ai is None or aj is None or ak is None:
            raise ValueError("Angle patch must have an ai, aj and ak attribute")

        self.params = props_to_patches(props)
        self.ai = ai
        self.aj = aj
        self.ak = ak
        self.id = "---".join([ai, aj, ak])
        self.id_sym = "---".join(reversed([ai, aj, ak]))


@dataclass(order=True)
class DihedralPatch:
    """Instructions to patch one dihedral"""

    ai: str
    aj: str
    ak: str
    al: str
    func: str
    periodicity: str
    id: str
    params: dict[str, ParamPatch]

    def __init__(self, elem: Element):
        props = elem.attrib
        ai = props.pop("ai", None)
        aj = props.pop("aj", None)
        ak = props.pop("ak", None)
        al = props.pop("al", None)
        func = props.pop("func", None)
        periodicity = props.pop("periodicity", None)
        if (
            ai is None
            or aj is None
            or ak is None
            or al is None
            or func is None
            or periodicity is None
        ):
            raise ValueError(
                "Angle patch must have an ai, aj, ak, al and periodicity attribute"
            )

        self.params = props_to_patches(props)
        self.ai = ai
        self.aj = aj
        self.ak = ak
        self.ak = al
        self.func = func
        self.periodicity = periodicity
        self.id = "---".join([ai, aj, ak, al]) + ":::" + func + ":::" + periodicity
        self.id_sym = (
            "---".join(reversed([ai, aj, ak, al])) + ":::" + func + ":::" + periodicity
        )


class FFPatches:
    """A container for forcefield patches"""

    atompatches: dict[str, AtomPatch]
    bondpatches: dict[str, BondPatch]
    pairpatches: dict[str, PairPatch]
    anglepatches: dict[str, AnglePatch]
    dihedralpatches: dict[str, DihedralPatch]

    def __init__(self, path: Path) -> None:
        xml = read_xml_ff(path)
        self.atompatches = {}
        if elems := xml.findall("Atoms/Atom"):
            for elem in elems:
                atompatch = AtomPatch(elem)
                self.atompatches[atompatch.id] = atompatch

        self.bondpatches = {}
        if elems := xml.findall("Bonds/Bond"):
            for elem in elems:
                bondpatch = BondPatch(elem)
                self.bondpatches[bondpatch.id] = bondpatch

        self.pairpatches = {}
        if elems := xml.findall("Pairs/Pair"):
            for elem in elems:
                pairpatch = PairPatch(elem)
                self.pairpatches[pairpatch.id] = pairpatch

        self.anglepatches = {}
        if elems := xml.findall("Angles/Angle"):
            for elem in elems:
                anglepatch = AnglePatch(elem)
                self.anglepatches[anglepatch.id] = anglepatch

        # note... <https://manual.gromacs.org/current/reference-manual/functions/bonded-interactions.html#proper-dihedrals-periodic-type>
        # periodicity is also an identifier, not a parameter!
        self.dihedralpatches = {}
        if elems := xml.findall("Dihedrals/Dihedral"):
            for elem in elems:
                dihedralpatch = DihedralPatch(elem)
                self.dihedralpatches[dihedralpatch.id] = dihedralpatch

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


Patch = Union[AtomPatch, BondPatch, PairPatch, AnglePatch, DihedralPatch]
Patches = Union[
    dict[str, AtomPatch],
    dict[str, BondPatch],
    dict[str, PairPatch],
    dict[str, AnglePatch],
    dict[str, DihedralPatch],
]
