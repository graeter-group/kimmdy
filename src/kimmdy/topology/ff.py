from __future__ import annotations

import logging
import textwrap
from pathlib import Path
from typing import Optional

from kimmdy.constants import FFFUNC
from kimmdy.parsing import read_top
from kimmdy.topology.atomic import (
    AngleId,
    AngleType,
    AtomId,
    AtomType,
    BondId,
    BondType,
    DihedralType,
    ImproperDihedralId,
    NonbondParamType,
    ProperDihedralId,
    ResidueType,
)
from kimmdy.topology.utils import get_top_section

logger = logging.getLogger(__name__)


class FF:
    """Container for parsed forcefield data."""

    def __init__(self, top: dict, residuetypes_path: Optional[Path] = None):
        self.atomtypes: dict[AtomId, AtomType] = {}
        self.bondtypes: dict[BondId, BondType] = {}
        self.angletypes: dict[AngleId, AngleType] = {}
        self.proper_dihedraltypes: dict[ProperDihedralId, DihedralType] = {}
        self.improper_dihedraltypes: dict[ImproperDihedralId, DihedralType] = {}
        self.residuetypes: dict[str, ResidueType] = {}
        self.nonbond_params: dict[BondId, NonbondParamType] = {}

        ffdir: Optional[Path] = top["ffdir"]

        if defaults := get_top_section(top, "defaults"):
            self.defaults = defaults

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

        nonbond_params = get_top_section(top, "nonbond_params")
        if nonbond_params is not None:
            for l in nonbond_params:
                nonbond_param = NonbondParamType.from_top_line(l)
                self.nonbond_params[(nonbond_param.i, nonbond_param.j)] = nonbond_param

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
                if dihedraltype.funct == FFFUNC["mult_improper_dihedral"]:
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

        if residuetypes_path:
            logger.debug(f"Using specified residuetypes file: {residuetypes_path}")
        else:
            logger.debug("Trying to use default amber protein residuetypes file.")
            if ffdir is None:
                logger.warning("ffdir is None. No residuetypes will be parsed.")
                return
            residuetypes_path = ffdir / "aminoacids.rtp"
            if not residuetypes_path.exists():
                logger.warning(
                    "aminoacids.rtp not found in ffdir. No residuetypes will be parsed."
                )
                return
        residuetypes_dict = read_top(residuetypes_path, use_gmx_dir=False)
        for k, v in residuetypes_dict.items():
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
