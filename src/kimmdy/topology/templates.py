"""Residue templates: the atoms, bonds, charges and types a residue should have.

Two sources describe residues, both read into ResidueTemplate:

- the rtp entries of the force field (aminoacids.rtp etc., FF.residuetypes),
  if changer.topology.product_template_sources contains rtp
- the residues of the moleculetypes in the force field's ITP files, if it
  contains itp

match_residue compares the atoms of a residue (found along its bonds, see
get_residue_by_bonding) with a template by atom names and bonds. It gives
the template name of every atom and the template atoms that are missing
(hydrogens abstracted from a radical). TemplateLibrary holds the templates
of a force field and finds those that may describe a residue: ITP templates
by residue name, rtp entries by the residue name and its protonation and
terminal variants (residue_variants).

The residues recipe steps touched (touched_residues) and the hydrogens
that moved between them (moved_hydrogens) are found here too, for both
users of a recipe segment's templates.

Users:

- Topology.update_names_LJ: names, types and masses of the residues a
  recipe segment touched (assign_recipe_names)
- kimmdy.topology.charges: their charges (Topology.update_charges)
- MoleculeType._get_residue_improper_dihedrals: the rtp entry whose
  impropers a residue gets (TemplateLibrary.rtp_template)
- MoleculeType.find_radicals: atoms whose type allows several numbers of
  bonds (TemplateLibrary.radical_by_template)
"""

from __future__ import annotations

import logging
import math
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Collection, Optional

from kimmdy.parsing import read_top
from kimmdy.recipe import Bind, Break, RecipeStep
from kimmdy.topology.utils import (
    get_moleculetype_atomics,
    get_moleculetype_header,
    get_residue_by_bonding,
)

if TYPE_CHECKING:
    from kimmdy.topology.atomic import Atom, AtomType, ResidueType
    from kimmdy.topology.topology import Topology

logger = logging.getLogger("kimmdy.topology.templates")

DEFAULT_TEMPLATE_TOLERANCE = 1e-3

# how closely a residue has to match a template, tried in this order
MATCH_LEVELS = ("exact", "renamed hydrogen", "radical")

# name of a hydrogen without a free position in its residue's template
UNNAMED_HYDROGEN = "HX"

# Amber types of aromatic carbons, and the hydrogen types by the partner's
# element, for hydrogens that bind where no template has a hydrogen
# (TemplateLibrary.hydrogen_type_on)
AROMATIC_CARBON_TYPES = frozenset({"CA", "CB", "CC", "CN", "CR", "CV", "CW", "C*"})
HYDROGEN_TYPE_BY_ELEMENT = {7: "H", 8: "HO", 16: "HS"}
# by mass, if the force field's atomtypes have no atomic number
ELEMENT_BY_MASS = {1: 1, 12: 6, 14: 7, 16: 8, 32: 16}

# rtp entries of one amino acid in different protonation states (and CYX,
# the disulfide cysteine). Topologies name residues by the amino acid, e.g.
# pdb2gmx and colbuilder write HIS for HID, HIE and HIP, so every entry of
# the group may describe the residue; the atom names and radicals decide.
PROTONATION_STATES = (
    ("HIS", "HID", "HIE", "HIP"),
    ("LYS", "LYN"),
    ("ASP", "ASH"),
    ("GLU", "GLH"),
    ("CYS", "CYX", "CYM"),
    ("DOP", "DO1", "DO2"),
)


def residue_variants(residue: str) -> list[str]:
    """Residue names that may describe a residue and their terminal variants.

    These are the residue itself and the other protonation states of its
    amino acid (PROTONATION_STATES). Assumes the Amber naming N<name> and
    C<name> for the N- and C-terminal variants. This could be improved by
    reading the force field's .r2b files.
    """
    group = next((g for g in PROTONATION_STATES if residue in g), ())
    names = [residue, *(name for name in group if name != residue)]
    return list(dict.fromkeys(v for name in names for v in (name, f"N{name}", f"C{name}")))


def is_hydrogen(atom_type: str) -> bool:
    return atom_type.upper().startswith("H")


def _bond(name_1: str, name_2: str) -> tuple[str, str]:
    return tuple(sorted((name_1, name_2)))  # type: ignore[return-value]


@dataclass(frozen=True)
class TemplateAtom:
    """One atom of a template. mass is None for rtp entries (the atom type's mass)."""

    name: str
    type: str
    charge: str
    mass: Optional[str] = None


@dataclass(frozen=True)
class ResidueTemplate:
    """A residue as an rtp entry or an ITP file describes it.

    source is rtp or itp, name the rtp entry or the ITP moleculetype, residue
    the residue name. bonds are the bonds inside the residue, as sorted
    pairs of atom names. linked are the atoms bonded to a neighboring
    residue. boundary_bonds (ITP only, None for rtp) are those bonds
    exactly: (own atom, other residue, other atom).
    """

    source: str
    name: str
    residue: str
    atoms: dict[str, TemplateAtom]
    bonds: frozenset[tuple[str, str]]
    linked: frozenset[str]
    boundary_bonds: Optional[frozenset[tuple[str, str, str]]] = None
    file: Optional[Path] = None
    local_resnr: Optional[str] = None

    @classmethod
    def from_rtp(cls, residuetype: ResidueType) -> ResidueTemplate:
        """The template of an rtp entry.

        The entry names atoms of the neighboring residues -C, +N etc. Its own
        atom X is linked to a neighbor if X has such a partner, or if the
        entry lists -X or +X (the same atom of the neighbor, e.g. C in -C N).
        """
        bonds = set()
        linked = set()
        for bond in residuetype.bonds:
            own = [name for name in bond if name[0] not in "-+"]
            if len(own) == 2:
                bonds.add(_bond(*own))
            elif len(own) == 1:
                linked.add(own[0])
            linked.update(name[1:] for name in bond if name[0] in "-+")
        return cls(
            source="rtp",
            name=residuetype.residue,
            residue=residuetype.residue,
            atoms={
                name: TemplateAtom(name=name, type=spec.type, charge=spec.charge)
                for name, spec in residuetype.atoms.items()
            },
            bonds=frozenset(bonds),
            linked=frozenset(linked & residuetype.atoms.keys()),
        )

    @property
    def label(self) -> str:
        return f"{self.source}:{self.name}"

    @property
    def charge(self) -> float:
        return sum(float(atom.charge) for atom in self.atoms.values())

    def bonded(self, name: str) -> set[str]:
        """Names of the atoms bonded to this one inside the residue."""
        return {other for bond in self.bonds if name in bond for other in bond if other != name}

    def hydrogens_on(self, name: str) -> list[str]:
        """Names of the hydrogens bonded to this atom, sorted."""
        return sorted(
            other for other in self.bonded(name) if is_hydrogen(self.atoms[other].type)
        )


@dataclass(frozen=True)
class ResidueMatch:
    """How a residue matches a template.

    names gives the template name of every atom of the residue (by atom
    number), missing the template atoms that are not in the residue and
    the number of the atom each was bonded to, a radical.
    """

    template: ResidueTemplate
    level: str
    names: dict[str, str]
    missing: dict[str, str]

    @property
    def radicals(self) -> set[str]:
        return set(self.missing.values())


def match_residue(
    residue: dict[str, Atom],
    template: ResidueTemplate,
    level: str = "exact",
    renamed: Collection[str] = (),
    same_types: bool = False,
    atoms: Optional[dict[str, Atom]] = None,
) -> Optional[ResidueMatch]:
    """How the atoms of a residue match a template at a level; None if they do not.

    exact: the atoms carry the names of the template's atoms, all of them.
    renamed hydrogen: as exact, but the hydrogens in renamed (e.g. those
        that moved into the residue) may carry other names. Each takes a
        free hydrogen position on its partner, its own name if that is one.
    radical: as renamed hydrogen, but hydrogens of the template may be
        missing (abstracted) if each was bonded to one atom of the residue,
        which is then a radical.
    At every level the bonds inside the residue are those of the template
    between the atoms present. With same_types the atoms not in renamed
    have the template's types. A template with boundary_bonds (ITP) also
    needs the same bonds to other residues; atoms are all atoms of the
    topology, to look these up.
    """
    if level not in MATCH_LEVELS:
        raise ValueError(f"Unknown match level {level!r}, expected one of {MATCH_LEVELS}.")
    renamed = set() if level == "exact" else set(renamed) & residue.keys()

    # atoms that keep their names: all of them for an exact match
    names = {nr: atom.atom for nr, atom in residue.items() if nr not in renamed}
    if len(set(names.values())) != len(names):
        return None  # repeated atom names
    if not set(names.values()) <= template.atoms.keys():
        return None
    free = set(template.atoms) - set(names.values())

    # next: renamed hydrogens take free hydrogen positions on their partner
    for nr in sorted(renamed, key=int):
        options = [
            name
            for partner in sorted(
                names[other] for other in residue[nr].bound_to_nrs if other in names
            )
            for name in template.hydrogens_on(partner)
            if name in free
        ]
        if not options:
            return None
        names[nr] = residue[nr].atom if residue[nr].atom in options else options[0]
        free.discard(names[nr])

    # next: atoms of the template still missing; only hydrogens abstracted
    # from a radical, each bonded to one atom of the residue
    if free and level != "radical":
        return None
    numbers = {name: nr for nr, name in names.items()}
    missing = {}
    for name in sorted(free):
        partners = [numbers[other] for other in template.bonded(name) if other in numbers]
        if not is_hydrogen(template.atoms[name].type) or len(partners) != 1:
            return None
        missing[name] = partners[0]

    # next: the bonds inside the residue are the template's
    bonds = {
        _bond(name, names[other])
        for nr, name in names.items()
        for other in residue[nr].bound_to_nrs
        if other in names
    }
    if bonds != {bond for bond in template.bonds if not free & set(bond)}:
        return None

    # next: types and bonds to other residues, where asked for
    if same_types and any(
        residue[nr].type != template.atoms[name].type
        for nr, name in names.items()
        if nr not in renamed
    ):
        return None
    if template.boundary_bonds is not None:
        if atoms is None:
            raise ValueError(f"Matching template {template.label} needs all atoms.")
        boundary = {
            (name, atoms[other].residue, atoms[other].atom)
            for nr, name in names.items()
            for other in residue[nr].bound_to_nrs
            if other not in residue
        }
        if boundary != template.boundary_bonds:
            return None
    return ResidueMatch(template=template, level=level, names=names, missing=missing)


def _is_copy(residue: dict[str, Atom], template: ResidueTemplate, atoms: dict[str, Atom]) -> bool:
    """Whether a residue is exactly the template: names, types, charges, masses and bonds."""
    match = match_residue(residue, template, "exact", same_types=True, atoms=atoms)
    if match is None:
        return False
    for nr, name in match.names.items():
        atom, spec = residue[nr], template.atoms[name]
        if not math.isclose(float(atom.charge), float(spec.charge), rel_tol=0.0, abs_tol=1e-8):
            return False
        if spec.mass is not None and (
            atom.mass is None
            or not math.isclose(float(atom.mass), float(spec.mass), rel_tol=0.0, abs_tol=1e-8)
        ):
            return False
    return True


def radical_mismatch(match: ResidueMatch, radicals: Collection[str]) -> Optional[str]:
    """Why a match does not fit the radicals of the topology (radicals); None if it fits.

    Each hydrogen missing from the template must have been bonded to a
    radical, and each radical of the residue must be one in the template
    too: bonded there to a missing hydrogen or to a neighboring residue (a
    broken bond between residues). Otherwise the template describes another
    protonation state, e.g. HID for a HIP that lost HE2.
    """
    for name, nr in match.missing.items():
        if nr not in radicals:
            return f"{match.names[nr]} lost {name} but is not a radical"
    for nr, name in match.names.items():
        if nr in radicals and nr not in match.radicals and name not in match.template.linked:
            return f"radical {name} is not a radical in the entry"
    return None


class TemplateLibrary:
    """The residue templates of a force field.

    itp holds the ITP templates by residue name, rtp the rtp entries by
    entry name. transitions maps a precursor ITP moleculetype to its
    product (changer.topology.product_template_transitions). atomtypes
    are the force field's (FF.atomtypes), for the types of hydrogens no
    template places (hydrogen_type_on).
    """

    def __init__(self, atomtypes: Optional[dict[str, AtomType]] = None) -> None:
        self.itp: dict[str, list[ResidueTemplate]] = {}
        self.rtp: dict[str, ResidueTemplate] = {}
        self.transitions: dict[str, str] = {}
        self.atomtypes: dict[str, AtomType] = atomtypes if atomtypes is not None else {}
        # moleculetype of the kept ITP template with each signature, and
        # moleculetypes dropped as copies: the moleculetype kept instead
        self._itp_signatures: dict[tuple, str] = {}
        self._itp_copies: dict[str, str] = {}
        # types of the hydrogens the rtp entries bond to an atom of each type
        self._rtp_hydrogen_types: dict[str, set[str]] = {}

    def add_rtp(self, residuetypes: dict[str, ResidueType]) -> None:
        for name, residuetype in residuetypes.items():
            template = ResidueTemplate.from_rtp(residuetype)
            self.rtp[name] = template
            for atom in template.atoms.values():
                for h in template.hydrogens_on(atom.name):
                    self._rtp_hydrogen_types.setdefault(atom.type, set()).add(
                        template.atoms[h].type
                    )

    def add_itp(self, template: ResidueTemplate) -> None:
        """Add an ITP template, unless one with the same atoms and bonds is there.

        Charges and masses are compared as numbers, so copies that write them
        differently (0.09 and 0.090000) are copies too. The moleculetype of a
        dropped copy is remembered, for the transitions.
        """
        signature = (
            template.residue,
            tuple(
                sorted(
                    (
                        a.name,
                        a.type,
                        round(float(a.charge), 6),
                        None if a.mass is None else round(float(a.mass), 6),
                    )
                    for a in template.atoms.values()
                )
            ),
            tuple(sorted(template.bonds)),
            tuple(sorted(template.boundary_bonds or ())),
        )
        kept = self._itp_signatures.get(signature)
        if kept is not None:
            if kept != template.name:
                self._itp_copies.setdefault(template.name, kept)
            return
        self._itp_signatures[signature] = template.name
        self.itp.setdefault(template.residue, []).append(template)

    def load_itp(self, ffdir: Path, pattern: str = "*.itp") -> None:
        """Add the residues of every moleculetype in the ITP files of ffdir."""
        paths = sorted(ffdir.glob(pattern))
        if not paths:
            raise FileNotFoundError(f"No ITP files matching {pattern!r} found in {ffdir}.")
        for itp_path in paths:
            parsed = read_top(itp_path, ffdir=ffdir, use_gmx_dir=False)
            for section_name in parsed:
                if not section_name.startswith("moleculetype"):
                    continue
                header = get_moleculetype_header(parsed, section_name)
                atomics = get_moleculetype_atomics(parsed, section_name)
                if header is None or atomics is None:
                    continue
                for template in _itp_templates(itp_path, header.name, atomics):
                    self.add_itp(template)
        logger.info(
            "Loaded %d unique ITP product templates for %d residue names from %s.",
            sum(len(templates) for templates in self.itp.values()),
            len(self.itp),
            ffdir,
        )

    def validate_atomtypes(self, atomtypes: dict[str, AtomType]) -> None:
        for templates in self.itp.values():
            for template in templates:
                for atom in template.atoms.values():
                    if atom.type not in atomtypes:
                        raise ValueError(
                            f"ITP product template {template.name}:"
                            f"{template.local_resnr} from {template.file} uses "
                            f"atom type {atom.type!r}, which is missing from "
                            "the parsed force-field [ atomtypes ] section."
                        )

    def set_transitions(self, transitions: dict[str, str]) -> None:
        """Configure explicit precursor-to-product moleculetype transitions.

        A moleculetype whose residues were dropped as copies of another's
        (add_itp) stands for that one, with a warning.
        """
        known = {template.name for templates in self.itp.values() for template in templates}

        def resolve(name: str) -> str:
            if name in known:
                return name
            kept = self._itp_copies[name]
            logger.warning(
                f"ITP moleculetype {name} has the same residues as {kept}; "
                f"the product template transitions use {kept} for it."
            )
            return kept

        resolved = {}
        for precursor, product in transitions.items():
            if not isinstance(precursor, str) or not isinstance(product, str):
                raise TypeError(
                    "Product template transition names must be strings, got "
                    f"{precursor!r} -> {product!r}."
                )
            if precursor == product:
                raise ValueError(
                    "Product template transition precursor and product must "
                    f"differ, got {precursor!r}."
                )
            unknown = {precursor, product} - known - self._itp_copies.keys()
            if unknown:
                raise ValueError(
                    "Unknown moleculetype(s) in product template transition "
                    f"{precursor!r} -> {product!r}: {sorted(unknown)}."
                )
            resolved[resolve(precursor)] = resolve(product)
        self.transitions = resolved

    def rtp_variants(self, residue: str) -> list[ResidueTemplate]:
        """The rtp entries that may describe a residue of this name (residue_variants)."""
        return [self.rtp[name] for name in residue_variants(residue) if name in self.rtp]

    def rtp_template(self, residue: dict[str, Atom]) -> Optional[ResidueTemplate]:
        """The rtp entry of a residue, by its atoms and bonds.

        The first entry of the residue's name and its variants (rtp_variants)
        that the residue matches (match_residue), trying the levels of
        MATCH_LEVELS in order: e.g. NSER for an N-terminal SER (H1, H2, H3),
        HIE for a HIS with HE2, the plain entry for a radical without its
        abstracted hydrogen. Hydrogens named HX may take any free hydrogen
        position. If they cannot all take one, the residue is matched again
        without them: an HX has no place in the entry and is in none of its
        impropers, so it does not decide; a position it would have taken
        counts as abstracted (radical level). If several entries match at one
        level, the first is used. None if no entry matches.
        """
        if not residue:
            return None
        name = next(iter(residue.values())).residue
        unnamed = {nr for nr, atom in residue.items() if atom.atom == UNNAMED_HYDROGEN}
        templates = self.rtp_variants(name)
        attempts = [(residue, unnamed)]
        if unnamed:
            without = {nr: atom for nr, atom in residue.items() if nr not in unnamed}
            attempts.append((without, set()))
        for atoms, renamed in attempts:
            for level in MATCH_LEVELS:
                for template in templates:
                    if match_residue(atoms, template, level, renamed=renamed) is not None:
                        return template
        return None

    def _element(self, atom_type: str) -> Optional[int]:
        """Atomic number of an atom type: its at.num, else by its mass; None if unknown."""
        atomtype = self.atomtypes.get(atom_type)
        if atomtype is None:
            return None
        if atomtype.at_num.strip().isdigit():
            return int(atomtype.at_num)
        try:
            return ELEMENT_BY_MASS.get(round(float(atomtype.mass)))
        except ValueError:
            return None

    def hydrogen_type_on(
        self, partner: Atom, hydrogen: Atom, atoms: dict[str, Atom]
    ) -> Optional[str]:
        """Type of a hydrogen bonded to partner where no template has a hydrogen.

        1. the one hydrogen type the rtp entries bond to atoms of partner's
           type (CA: HA, CW: H4, CR: H5, N: H, OH: HO), if there is one
        2. else by partner's element and its other neighbors (atoms):
           aromatic carbon (AROMATIC_CARBON_TYPES) HA, carbon bonded to N,
           O or S H1, other carbon HC, N H, O HO, S HS

        None if the type is not in the force field or partner's element is
        unknown; the caller then keeps the hydrogen's type.
        """
        rtp_types = self._rtp_hydrogen_types.get(partner.type, set())
        if len(rtp_types) == 1:
            h_type = next(iter(rtp_types))
        else:
            element = self._element(partner.type)
            if element == 6:
                neighbors = (
                    self._element(atoms[nr].type)
                    for nr in partner.bound_to_nrs
                    if nr != hydrogen.nr and nr in atoms
                )
                if partner.type in AROMATIC_CARBON_TYPES:
                    h_type = "HA"
                elif any(e in (7, 8, 16) for e in neighbors):
                    h_type = "H1"
                else:
                    h_type = "HC"
            else:
                h_type = HYDROGEN_TYPE_BY_ELEMENT.get(element)  # type: ignore[arg-type]
        if h_type not in self.atomtypes:
            return None
        return h_type

    def hydrogen_identity(
        self,
        residue: dict[str, Atom],
        partner: Atom,
        hydrogen: Atom,
        atoms: Optional[dict[str, Atom]] = None,
    ) -> Optional[tuple[str, str]]:
        """Name and type of a hydrogen that binds to partner, from the rtp entry of partner's residue.

        residue are the atoms of partner's residue without the hydrogen
        (rtp_template). The hydrogen takes a free hydrogen position on
        partner, its own name if that is one, and the position's type. If
        all are taken it is named HX with the type of the first; if the
        entry has none on partner it is named HX with the type from
        hydrogen_type_on (partner's neighbors from atoms, default residue),
        or keeps its type. None if the residue has no rtp entry.
        """
        template = self.rtp_template(residue)
        if template is None:
            return None
        positions = template.hydrogens_on(partner.atom)
        if not positions:
            h_type = self.hydrogen_type_on(
                partner, hydrogen, atoms if atoms is not None else residue
            )
            return UNNAMED_HYDROGEN, h_type or hydrogen.type
        taken = {atom.atom for nr, atom in residue.items() if nr != hydrogen.nr}
        free = [name for name in positions if name not in taken]
        if not free:
            return UNNAMED_HYDROGEN, template.atoms[positions[0]].type
        name = hydrogen.atom if hydrogen.atom in free else free[0]
        return name, template.atoms[name].type

    def match_itp(
        self,
        residue: dict[str, Atom],
        joined: Collection[str],
        atoms: dict[str, Atom],
        atoms_a: Optional[dict[str, Atom]] = None,
    ) -> Optional[ResidueMatch]:
        """The ITP template of a residue, with the residue's types and bonds.

        Hydrogens in joined (they moved into the residue) may take other
        names and types (match_residue, renamed hydrogen). If no template
        fits, a residue that was exactly a precursor template before (in
        atoms_a, without the joined hydrogens) may become its product
        (transitions), whose types may differ from the residue's. None if no
        template fits. If several fit, e.g. two charge sets of one molecule,
        the first (by ITP file name, then order in the file) is used, with a
        warning.
        """
        name = next((atom.residue for nr, atom in residue.items() if nr not in joined), None)
        templates = self.itp.get(name, []) if name is not None else []
        matches = [
            match
            for template in templates
            if (
                match := match_residue(
                    residue, template, "renamed hydrogen", joined, same_types=True, atoms=atoms
                )
            )
            is not None
        ]
        if not matches and self.transitions and atoms_a is not None:
            residue_a = {nr: atoms_a[nr] for nr in residue if nr not in joined and nr in atoms_a}
            products = {
                self.transitions[template.name]
                for template in templates
                if template.name in self.transitions and _is_copy(residue_a, template, atoms_a)
            }
            matches = [
                match
                for template in templates
                if template.name in products
                and (
                    match := match_residue(
                        residue, template, "renamed hydrogen", joined, atoms=atoms
                    )
                )
                is not None
            ]
        if len(matches) > 1:
            resnr = next(iter(residue.values())).resnr
            fitting = ", ".join(
                f"{m.template.label} ({m.template.file.name if m.template.file else '?'})"
                for m in matches
            )
            logger.warning(
                f"Several ITP templates fit residue {resnr} {name}: {fitting}; "
                f"the first, {matches[0].template.label}, is used."
            )
        return matches[0] if matches else None

    def radical_by_template(
        self,
        atoms: dict[str, Atom],
        atom: Atom,
        tolerance: float = DEFAULT_TEMPLATE_TOLERANCE,
    ) -> Optional[bool]:
        """Whether the rtp entry of an atom's residue makes the atom a radical; None if no entry fits.

        The residue fits an entry (rtp_variants) if it has the entry's atoms,
        types and bonds, maybe without hydrogens (abstracted from a radical),
        and the entry's charge. Types and charge decide between states with
        the same atoms: CYS without HG has the atoms of CYX (SG of type S)
        and CYM, LYS without HZ3 those of LYN, but a radical keeps the charge
        of its entry (radical convention).
        """
        residue = get_residue_by_bonding(atom, atoms)
        charge = sum(float(a.charge) for a in residue.values())
        templates = self.rtp_variants(atom.residue)
        for level in ("exact", "radical"):
            for template in templates:
                match = match_residue(residue, template, level, same_types=True)
                if match is not None and abs(template.charge - charge) <= tolerance:
                    return atom.nr in match.radicals
        return None


def _itp_templates(file: Path, moleculetype: str, atomics: dict) -> list[ResidueTemplate]:
    """The templates of the residues of one moleculetype in an ITP file."""
    # atom number: (residue number, residue name, atom)
    atoms_by_nr = {
        line[0]: (
            line[2],
            line[3],
            TemplateAtom(
                name=line[4],
                type=line[1],
                charge=line[6],
                mass=line[7] if len(line) > 7 else None,
            ),
        )
        for line in atomics.get("atoms") or []
    }
    residues: dict[str, list[tuple[str, TemplateAtom]]] = {}
    for resnr, residue, atom in atoms_by_nr.values():
        residues.setdefault(resnr, []).append((residue, atom))

    # next: bonds inside each residue and to other residues
    bonds: dict[str, set[tuple[str, str]]] = {resnr: set() for resnr in residues}
    boundary: dict[str, set[tuple[str, str, str]]] = {resnr: set() for resnr in residues}
    for line in atomics.get("bonds") or []:
        if len(line) < 2:
            continue
        if line[0] not in atoms_by_nr or line[1] not in atoms_by_nr:
            raise ValueError(
                f"Bond {line[:2]} in {file} {moleculetype} references an "
                "atom missing from its [ atoms ] section."
            )
        resnr_1, residue_1, atom_1 = atoms_by_nr[line[0]]
        resnr_2, residue_2, atom_2 = atoms_by_nr[line[1]]
        if resnr_1 == resnr_2:
            bonds[resnr_1].add(_bond(atom_1.name, atom_2.name))
        else:
            boundary[resnr_1].add((atom_1.name, residue_2, atom_2.name))
            boundary[resnr_2].add((atom_2.name, residue_1, atom_1.name))

    templates = []
    for resnr, members in residues.items():
        residue_names = {residue for residue, _ in members}
        if len(residue_names) != 1:
            raise ValueError(
                f"Molecule {moleculetype} residue number {resnr} in "
                f"{file} has multiple residue names: {sorted(residue_names)}."
            )
        names = [atom.name for _, atom in members]
        if len(names) != len(set(names)):
            raise ValueError(
                f"Molecule {moleculetype} residue number {resnr} in "
                f"{file} contains duplicate atom names."
            )
        templates.append(
            ResidueTemplate(
                source="itp",
                name=moleculetype,
                residue=residue_names.pop(),
                atoms={atom.name: atom for _, atom in members},
                bonds=frozenset(bonds[resnr]),
                linked=frozenset(name for name, _, _ in boundary[resnr]),
                boundary_bonds=frozenset(boundary[resnr]),
                file=file,
                local_resnr=resnr,
            )
        )
    return templates


# --- the residues recipe steps touched, and their names -------------------------------------


@dataclass
class TouchedResidue:
    """A residue that recipe steps touched.

    atoms are its atoms in B, found along the bonds (touched_residues); the
    fragments of a split residue count together. charge_a is the charge of
    the same residue in A.
    """

    resnr: str
    atoms: dict[str, Atom]
    charge_a: float


def recipe_bond_atoms(recipe_steps: list[RecipeStep]) -> tuple[set[str], set[str]]:
    """The atoms of the Break and Bind steps, and those of the Break steps."""
    bond_steps = [step for step in recipe_steps if isinstance(step, (Break, Bind))]
    atom_nrs = {nr for step in bond_steps for nr in (step.atom_id_1, step.atom_id_2)}
    broken = {
        nr
        for step in bond_steps
        if isinstance(step, Break)
        for nr in (step.atom_id_1, step.atom_id_2)
    }
    return atom_nrs, broken


def _residue_walks(top: Topology, atom_nrs: set[str]) -> list[dict[str, Atom]]:
    """The residues containing these atoms, found along the bonds, each once."""
    walks: list[dict[str, Atom]] = []
    for nr in sorted(atom_nrs, key=int):
        atom = top.atoms.get(nr)
        if atom is not None and not any(nr in walk for walk in walks):
            walks.append(get_residue_by_bonding(atom, top.atoms))
    return walks


def touched_residues(top_a: Topology, top_b: Topology, atom_nrs: set[str]) -> list[TouchedResidue]:
    """The residues containing these atoms, with their atoms in B and their charge in A.

    Residues are found along the bonds, in A and in B. Residue numbers are
    not unique in a merged topology, so they do not group atoms. A residue
    of A and one of B are the same residue if the atoms of one are a subset
    of the other's and at least one of their atoms keeps its residue number:
    a residue split by a broken bond is one residue in A and several in B, a
    residue that lost or gained a hydrogen keeps all its other atoms, and
    pieces joined by a new bond are several residues in A and one in B.
    Atoms that only moved into a residue, such as those of a water without
    bonds (settles) in a hydrolysis, are residues of their own in A; their
    charge in A does not count for the residue they join. Matches are
    joined transitively.
    """
    walks_a = _residue_walks(top_a, atom_nrs)
    walks_b = _residue_walks(top_b, atom_nrs)

    # union-find over the walks, those of A first, then those of B
    parent = list(range(len(walks_a) + len(walks_b)))

    def find(i: int) -> int:
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    for i, walk_a in enumerate(walks_a):
        for j, walk_b in enumerate(walks_b):
            nested = walk_a.keys() <= walk_b.keys() or walk_b.keys() <= walk_a.keys()
            # an atom that keeps its residue number links them; atoms that only
            # moved into the residue (the O of a water without bonds) do not
            stays = any(
                nr in walk_b and walk_b[nr].resnr == atom.resnr for nr, atom in walk_a.items()
            )
            if nested and stays:
                parent[find(len(walks_a) + j)] = find(i)

    groups: dict[int, tuple[list[dict[str, Atom]], list[dict[str, Atom]]]] = {}
    for i, walk in enumerate(walks_a + walks_b):
        groups.setdefault(find(i), ([], []))[i >= len(walks_a)].append(walk)

    residues = []
    for group_a, group_b in groups.values():
        if not group_b:
            continue  # none of its atoms are in B
        atoms = {nr: atom for walk in group_b for nr, atom in walk.items()}
        resnr = atoms[min(atoms, key=int)].resnr
        if group_a:
            charge_a = sum(float(atom.charge) for walk in group_a for atom in walk.values())
        else:
            logger.warning(
                f"No residue of the topology before the recipe matches residue {resnr} "
                "after it; its charge before is taken from the same atoms."
            )
            charge_a = sum(float(top_a.atoms[nr].charge) for nr in atoms if nr in top_a.atoms)
        residues.append(TouchedResidue(resnr=resnr, atoms=atoms, charge_a=charge_a))
    return residues


def moved_hydrogens(
    top_a: Topology, top_b: Topology, atom_nrs: set[str]
) -> dict[str, tuple[Optional[str], Optional[str]]]:
    """Hydrogens among these atoms whose bonding partner changed: {hydrogen: (partner in A, partner in B)}."""
    moved = {}
    for nr in atom_nrs:
        if nr not in top_a.atoms or nr not in top_b.atoms:
            continue
        if not is_hydrogen(top_b.atoms[nr].type):
            continue
        partner_a = next(iter(top_a.atoms[nr].bound_to_nrs), None)
        partner_b = next(iter(top_b.atoms[nr].bound_to_nrs), None)
        if partner_a != partner_b:
            moved[nr] = (partner_a, partner_b)
    return moved


def joined_hydrogens(
    residue: TouchedResidue, moved: dict[str, tuple[Optional[str], Optional[str]]]
) -> set[str]:
    """The moved hydrogens whose new partner is in the residue."""
    return {h for h, (_, partner_b) in moved.items() if partner_b in residue.atoms}


def _residue_identities(
    top_a: Topology,
    top_b: Topology,
    residue: TouchedResidue,
    joined: set[str],
    moved: dict[str, tuple[Optional[str], Optional[str]]],
) -> tuple[str, dict[str, tuple[str, str, Optional[str]]]]:
    """Where the names of a residue come from, and {atom: (name, type, mass)}; mass None is the type's."""
    templates = top_b.ff.templates

    # 1. the ITP template: all atoms
    match = templates.match_itp(residue.atoms, joined, top_b.atoms, top_a.atoms)
    if match is not None:
        atoms = match.template.atoms
        return match.template.label, {
            nr: (name, atoms[name].type, atoms[name].mass) for nr, name in match.names.items()
        }

    # next: 2. the first rtp entry that fits: the hydrogens that moved in
    if not joined:
        return "", {}
    name = next(atom.residue for nr, atom in residue.atoms.items() if nr not in joined)
    for level in MATCH_LEVELS:
        for template in templates.rtp_variants(name):
            match = match_residue(residue.atoms, template, level, renamed=joined)
            if match is not None and radical_mismatch(match, top_b.radicals) is None:
                return f"{template.label} ({level})", {
                    nr: (match.names[nr], template.atoms[match.names[nr]].type, None)
                    for nr in joined
                }

    # next: 3. each hydrogen that moved in by its position on its partner;
    # if no rtp entry matches (e.g. a residue split by a broken bond) it is
    # HX with the type for its partner, or its own type if there is none
    others = {nr: atom for nr, atom in residue.atoms.items() if nr not in joined}
    identities = {}
    for h in sorted(joined, key=int):
        partner = top_b.atoms[moved[h][1]]  # type: ignore[index]
        hydrogen = top_b.atoms[h]
        identity = templates.hydrogen_identity(others, partner, hydrogen, top_b.atoms)
        if identity is None:
            h_type = templates.hydrogen_type_on(partner, hydrogen, top_b.atoms)
            identity = (UNNAMED_HYDROGEN, h_type or hydrogen.type)
        identities[h] = (*identity, None)
    return "rtp positions of the hydrogens", identities


def _set_identities(top: Topology, identities: dict[str, tuple[str, str, Optional[str]]]) -> None:
    """Set names, types and masses (None: the type's mass)."""
    for nr, (name, atom_type, mass) in identities.items():
        atom = top.atoms[nr]
        atom.atom = name
        atom.type = atom_type
        if mass is not None:
            atom.mass = mass
        elif atom_type in top.ff.atomtypes:
            atom.mass = top.ff.atomtypes[atom_type].mass


def assign_recipe_names(
    top_a: Topology, top_b: Topology, recipe_steps: list[RecipeStep]
) -> list[TouchedResidue]:
    """Set names, types and masses of the residues recipe steps touched, from their templates.

    top_a is the topology before the steps and stays unchanged, top_b the
    one after them. Each residue containing an atom of a Break or Bind gets:

    1. from the ITP template that fits (TemplateLibrary.match_itp), the
       names, types and masses of all its atoms
    2. else from the first rtp entry that fits (rtp_variants; levels exact,
       renamed hydrogen, radical; radical_mismatch), the names and types of
       the hydrogens that moved into the residue
    3. else each hydrogen that moved in takes its position on its partner in
       the residue's rtp entry, or HX (TemplateLibrary.hydrogen_identity);
       if no entry matches the residue it is HX too. Where neither names a
       type, the type comes from its partner (TemplateLibrary.hydrogen_type_on),
       or it stays if that type is not in the force field

    The templates decide by the atoms and bonds alone, not by the charges.
    Returns the touched residues (touched_residues), with their atoms in B.
    """
    atom_nrs, _ = recipe_bond_atoms(recipe_steps)
    moved = moved_hydrogens(top_a, top_b, atom_nrs)
    residues = touched_residues(top_a, top_b, atom_nrs)
    for residue in residues:
        joined = joined_hydrogens(residue, moved)
        source, identities = _residue_identities(top_a, top_b, residue, joined, moved)
        if identities:
            _set_identities(top_b, identities)
            logger.info(f"Names and types of residue {residue.resnr} from {source}.")
    return residues
