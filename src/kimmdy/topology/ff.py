from __future__ import annotations
import textwrap
from pathlib import Path
from xml.etree.ElementTree import Element
from kimmdy.topology.atomic import *
from kimmdy.parsing import read_top, read_xml_ff, read_rtp
from typing import Union

from kimmdy.topology.utils import get_top_section


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
        self.residuetypes: dict[str, ResidueType]

        ffdir = top["ffdir"]

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

        p is an instance of IPython.lib.pretty.RepresentationPrinter
        <https://ipython.org/ipython-doc/3/api/generated/IPython.lib.pretty.html#IPython.lib.pretty.PrettyPrinter>
        """
        p.text(str(self))

