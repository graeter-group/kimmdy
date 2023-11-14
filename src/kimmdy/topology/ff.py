from __future__ import annotations
import textwrap
import logging
from kimmdy.topology.atomic import (
    Atom,
    Bond,
    Pair,
    Angle,
    Dihedral,
    MultipleDihedrals,
    PositionRestraint,
    DihedralRestraint,
    ResidueImproperSpec,
    ResidueProperSpec,
    Settle,
    Exclusion,
    AtomType,
    BondType,
    AngleType,
    DihedralType,
    ResidueType,
)
from kimmdy.parsing import read_top
from kimmdy.topology.utils import get_top_section

logger = logging.getLogger(__name__)


class FF:
    """Container for parsed forcefield data."""

    def __init__(self, top: dict):
        self.atomtypes: dict[str, AtomType] = {}
        self.bondtypes: dict[tuple[str, str], BondType] = {}
        self.angletypes: dict[tuple[str, str, str], AngleType] = {}
        self.proper_dihedraltypes: dict[
            tuple[str, str, str, str, str], DihedralType
        ] = {}
        self.improper_dihedraltypes: dict[tuple[str, str, str, str], DihedralType] = {}
        self.residuetypes: dict[str, ResidueType] = {}

        ffdir = top["ffdir"]

        atomtypes = get_top_section(top, "atomtypes")
        if atomtypes is not None:
            for l in atomtypes:
                atomtype = AtomType.from_top_line(l)
                self.atomtypes[atomtype.type] = atomtype

        bondtypes = get_top_section(top, "bondtypes")
        if bondtypes is not None:
            for l in bondtypes:
                bondtype = BondType.from_top_line(l)
                self.bondtypes[(bondtype.i, bondtype.j)] = bondtype

        angletypes = get_top_section(top, "angletypes")
        if angletypes is not None:
            for l in angletypes:
                angletype = AngleType.from_top_line(l)
                self.angletypes[(angletype.i, angletype.j, angletype.k)] = angletype

        dihedraltypes = get_top_section(top, "dihedraltypes")
        if dihedraltypes is not None:
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

        if ffdir is None:
            logger.warning("ffdir is None. No residuetypes will be parsed.")
            return
        aminoacids_path = ffdir / "aminoacids.rtp"
        if not aminoacids_path.exists():
            logger.warning(
                "aminoacids.rtp not found in ffdir. No residuetypes will be parsed."
            )
            return
        aminoacids = read_top(aminoacids_path, use_gmx_dir=False)
        for k, v in aminoacids.items():
            if k.startswith("BLOCK") or k in ["bondedtypes", "ffdir", "define"]:
                continue
            if not v.get("subsections"):
                raise AssertionError(f"key {k} has no subsections, only {v}.")
            self.residuetypes[k] = ResidueType.from_section(k, v["subsections"])

    def __str__(self) -> str:
        return textwrap.dedent(
            f"""
        ForceField parameters with
        {len(self.atomtypes)} atomtypes,
        {len(self.bondtypes)} bondtypes,
        {len(self.angletypes)} angletypes,
        {len(self.proper_dihedraltypes)} proper dihedraltypes
        {len(self.improper_dihedraltypes)} improper dihedraltypes
        {len(self.residuetypes)} residuetypes
        """
        )

    def __repr__(self) -> str:
        return f"FF({self.__dict__})"

    def _repr_pretty_(self, p, _):
        """A __repr__ for ipython.

        This whill be used if just the name of the object is entered in the ipython shell
        or a jupyter notebook.

        p is an instance of [IPython.lib.pretty.RepresentationPrinter](https://ipython.org/ipython-doc/3/api/generated/IPython.lib.pretty.html#IPython.lib.pretty.PrettyPrinter)
        """
        p.text(str(self))
