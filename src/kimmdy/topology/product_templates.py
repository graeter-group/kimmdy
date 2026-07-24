from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Optional

from kimmdy.parsing import read_top
from kimmdy.topology.utils import (
    get_moleculetype_atomics,
    get_moleculetype_header,
)

if TYPE_CHECKING:
    from kimmdy.topology.atomic import Atom, AtomType
    from kimmdy.topology.topology import Topology

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class ProductTemplateAtom:
    local_nr: str
    name: str
    atom_type: str
    residue: str
    local_resnr: str
    charge_group: str
    charge: str
    mass: Optional[str]


@dataclass(frozen=True)
class ProductTemplate:
    source: Path
    moleculetype: str
    local_resnr: str
    residue: str
    atoms: dict[str, ProductTemplateAtom]
    bonds: frozenset[tuple[str, str]]
    boundary_bonds: frozenset[tuple[str, str, str]]

    def signature(self) -> tuple:
        """Return the chemistry-relevant signature used to remove duplicates."""
        atoms = tuple(
            sorted(
                (
                    atom.name,
                    atom.atom_type,
                    atom.charge,
                    atom.mass,
                )
                for atom in self.atoms.values()
            )
        )
        return (
            self.residue,
            atoms,
            tuple(sorted(self.bonds)),
            tuple(sorted(self.boundary_bonds)),
        )


@dataclass(frozen=True)
class ProductTemplateMatch:
    template: ProductTemplate
    transferred_h_id: str
    acceptor_id: str
    transferred_h_template_name: str
    atom_mapping: dict[str, ProductTemplateAtom]


class ProductTemplateLibrary:
    """Closed-shell residue templates derived from GROMACS molecule ITPs."""

    def __init__(self) -> None:
        self.templates: dict[str, list[ProductTemplate]] = {}
        self._signatures: set[tuple] = set()

    def add(self, template: ProductTemplate) -> None:
        signature = template.signature()
        if signature in self._signatures:
            return
        self._signatures.add(signature)
        self.templates.setdefault(template.residue, []).append(template)

    @classmethod
    def from_itp_directory(
        cls,
        ffdir: Path,
        pattern: str = "*.itp",
    ) -> ProductTemplateLibrary:
        library = cls()
        paths = sorted(ffdir.glob(pattern))
        if not paths:
            raise FileNotFoundError(
                f"No ITP files matching {pattern!r} found in {ffdir}."
            )

        for itp_path in paths:
            parsed = read_top(itp_path, ffdir=ffdir, use_gmx_dir=False)
            for section_name in parsed:
                if not section_name.startswith("moleculetype"):
                    continue
                header = get_moleculetype_header(parsed, section_name)
                atomics = get_moleculetype_atomics(parsed, section_name)
                if header is None or atomics is None:
                    continue
                for template in _templates_from_moleculetype(
                    source=itp_path,
                    moleculetype=header.name,
                    atomics=atomics,
                ):
                    library.add(template)

        logger.info(
            "Loaded %d unique ITP product templates for %d residue names from %s.",
            sum(len(templates) for templates in library.templates.values()),
            len(library.templates),
            ffdir,
        )
        return library

    def validate_atomtypes(self, atomtypes: dict[str, AtomType]) -> None:
        for templates in self.templates.values():
            for template in templates:
                for atom in template.atoms.values():
                    if atom.atom_type not in atomtypes:
                        raise ValueError(
                            f"ITP product template {template.moleculetype}:"
                            f"{template.local_resnr} from {template.source} uses "
                            f"atom type {atom.atom_type!r}, which is missing from "
                            "the parsed force-field [ atomtypes ] section."
                        )

    def match_hydrogenated_residue(
        self,
        topology: Topology,
        acceptor: Atom,
        transferred_h: Atom,
    ) -> Optional[ProductTemplateMatch]:
        matches = []
        for template in self.templates.get(acceptor.residue, []):
            match = _match_template(
                template=template,
                topology=topology,
                acceptor=acceptor,
                transferred_h=transferred_h,
            )
            if match is not None:
                matches.append(match)

        if len(matches) == 1:
            return matches[0]
        if len(matches) > 1:
            descriptions = [
                f"{match.template.source.name}:"
                f"{match.template.moleculetype}:"
                f"{match.template.local_resnr}"
                for match in matches
            ]
            logger.warning(
                "Multiple exact ITP product templates matched hydrogen binding "
                "to %s:%s (atom %s): %s. Using the legacy fallback.",
                acceptor.residue,
                acceptor.atom,
                acceptor.nr,
                descriptions,
            )
        return None


def _templates_from_moleculetype(
    source: Path,
    moleculetype: str,
    atomics: dict,
) -> list[ProductTemplate]:
    atom_lines = atomics.get("atoms") or []
    if not atom_lines:
        return []

    atoms_by_nr = {
        line[0]: ProductTemplateAtom(
            local_nr=line[0],
            atom_type=line[1],
            local_resnr=line[2],
            residue=line[3],
            name=line[4],
            charge_group=line[5],
            charge=line[6],
            mass=line[7] if len(line) > 7 else None,
        )
        for line in atom_lines
    }
    atoms_by_resnr: dict[str, list[ProductTemplateAtom]] = {}
    for atom in atoms_by_nr.values():
        atoms_by_resnr.setdefault(atom.local_resnr, []).append(atom)

    internal_bonds: dict[str, set[tuple[str, str]]] = {
        resnr: set() for resnr in atoms_by_resnr
    }
    boundary_bonds: dict[str, set[tuple[str, str, str]]] = {
        resnr: set() for resnr in atoms_by_resnr
    }
    for line in atomics.get("bonds") or []:
        if len(line) < 2:
            continue
        atom1 = atoms_by_nr.get(line[0])
        atom2 = atoms_by_nr.get(line[1])
        if atom1 is None or atom2 is None:
            raise ValueError(
                f"Bond {line[:2]} in {source} {moleculetype} references an "
                "atom missing from its [ atoms ] section."
            )
        if atom1.local_resnr == atom2.local_resnr:
            internal_bonds[atom1.local_resnr].add(
                _canonical_bond(atom1.name, atom2.name)
            )
        else:
            boundary_bonds[atom1.local_resnr].add(
                (atom1.name, atom2.residue, atom2.name)
            )
            boundary_bonds[atom2.local_resnr].add(
                (atom2.name, atom1.residue, atom1.name)
            )

    templates = []
    for local_resnr, residue_atoms in atoms_by_resnr.items():
        residue_names = {atom.residue for atom in residue_atoms}
        if len(residue_names) != 1:
            raise ValueError(
                f"Molecule {moleculetype} residue number {local_resnr} in "
                f"{source} has multiple residue names: {sorted(residue_names)}."
            )
        names = [atom.name for atom in residue_atoms]
        if len(names) != len(set(names)):
            raise ValueError(
                f"Molecule {moleculetype} residue number {local_resnr} in "
                f"{source} contains duplicate atom names."
            )
        templates.append(
            ProductTemplate(
                source=source,
                moleculetype=moleculetype,
                local_resnr=local_resnr,
                residue=next(iter(residue_names)),
                atoms={atom.name: atom for atom in residue_atoms},
                bonds=frozenset(internal_bonds[local_resnr]),
                boundary_bonds=frozenset(boundary_bonds[local_resnr]),
            )
        )
    return templates


def _match_template(
    template: ProductTemplate,
    topology: Topology,
    acceptor: Atom,
    transferred_h: Atom,
) -> Optional[ProductTemplateMatch]:
    current_atoms = [
        atom
        for atom in topology.atoms.values()
        if atom.resnr == acceptor.resnr and atom.residue == acceptor.residue
    ]
    if len(current_atoms) != len(template.atoms):
        return None

    current_by_name = {
        atom.atom: atom for atom in current_atoms if atom.nr != transferred_h.nr
    }
    if len(current_by_name) != len(current_atoms) - 1:
        return None

    current_names = set(current_by_name)
    template_names = set(template.atoms)
    if not current_names.issubset(template_names):
        return None

    missing_names = template_names - current_names
    if len(missing_names) != 1:
        return None
    missing_name = next(iter(missing_names))
    template_h = template.atoms[missing_name]
    atomtype = topology.ff.atomtypes.get(template_h.atom_type)
    if atomtype is None or atomtype.at_num != "1":
        return None
    if _canonical_bond(acceptor.atom, missing_name) not in template.bonds:
        return None

    atom_mapping = {transferred_h.nr: template_h}
    for name, current_atom in current_by_name.items():
        template_atom = template.atoms[name]
        if current_atom.type != template_atom.atom_type:
            return None
        atom_mapping[current_atom.nr] = template_atom

    mapped_names = {
        atom_nr: template_atom.name
        for atom_nr, template_atom in atom_mapping.items()
    }
    current_bonds: set[tuple[str, str]] = set()
    current_boundary_bonds: set[tuple[str, str, str]] = set()
    current_ids = {atom.nr for atom in current_atoms}
    for atom in current_atoms:
        atom_name = mapped_names[atom.nr]
        for neighbor_nr in atom.bound_to_nrs:
            neighbor = topology.atoms[neighbor_nr]
            if neighbor_nr in current_ids:
                if int(atom.nr) < int(neighbor_nr):
                    current_bonds.add(
                        _canonical_bond(atom_name, mapped_names[neighbor_nr])
                    )
            else:
                current_boundary_bonds.add(
                    (atom_name, neighbor.residue, neighbor.atom)
                )

    # bind_bond calls the matcher before inserting the new H-acceptor bond.
    current_bonds.add(_canonical_bond(acceptor.atom, missing_name))
    if current_bonds != set(template.bonds):
        return None
    if current_boundary_bonds != set(template.boundary_bonds):
        return None

    return ProductTemplateMatch(
        template=template,
        transferred_h_id=transferred_h.nr,
        acceptor_id=acceptor.nr,
        transferred_h_template_name=missing_name,
        atom_mapping=atom_mapping,
    )


def _canonical_bond(atom1: str, atom2: str) -> tuple[str, str]:
    return tuple(sorted((atom1, atom2)))  # type: ignore[return-value]
