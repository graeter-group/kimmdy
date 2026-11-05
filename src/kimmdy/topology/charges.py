"""Partial charges of the product of a reaction recipe.

KIMMDY sets the charges once per recipe segment, the steps up to a Relax or
to the end of the recipe (see Topology.update_charges_bonded): after its
topology changes (Break, Bind, CustomTopMod) and before the parameterization,
from the topology before the segment (A) and after it (B). Every residue the
segment touched goes
through a list of charge schemes, most specific first, and takes the first
one that changes its charge by an integer compared to A:

1. an exact ITP template from bind_bond
2. the rtp entries of `residue_variants`: first exact matches, then matches
   with renamed joined hydrogens, then radicals after H abstraction (the
   charge of a missing hydrogen goes to the radical atom it was bonded to)
3. the local rules (`local_charges`)

Whether the recipe kept the net charge is checked separately
(`check_charge_conservation`).

A CustomTopMod that sets charges runs before this, and the schemes compute
the touched residues anew from A. Charges a CustomTopMod set on a touched
residue are therefore replaced, unless no scheme fits that residue. Example:
the hydrolysis plugin sets the charges of the new termini from the rtp,
which the rtp scheme then sets again.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Callable, Optional

from kimmdy.recipe import Bind, Break, RecipeStep
from kimmdy.topology.utils import get_residue_by_bonding

if TYPE_CHECKING:
    from kimmdy.topology.atomic import Atom, ResidueType
    from kimmdy.topology.product_templates import ProductTemplateMatch
    from kimmdy.topology.topology import Topology

logger = logging.getLogger("kimmdy.topology.charges")

CONSERVATION_CHECK_ACTIONS = ("warn", "error", "off")
DEFAULT_TEMPLATE_TOLERANCE = 1e-3


@dataclass
class _Residue:
    """A residue the recipe touched.

    atoms are its atoms in B: all atoms bonded to a recipe atom within its
    residue number, so the fragments of a split residue count together.
    charge_a is the charge of this residue number in A.
    """

    resnr: str
    atoms: dict[str, Atom]
    fragments: list[set[str]]
    charge_a: float


@dataclass(frozen=True)
class _Scheme:
    """Charges a scheme gives a residue.

    charges are written to the topology as they are: template charges and
    unchanged charges verbatim, only computed ones formatted. identities
    sets name, type and mass of atoms: of moved hydrogens for templates, of
    the whole residue for an ITP template.
    """

    source: str
    charges: dict[str, str]
    identities: dict[str, tuple[str, str, Optional[str]]] = field(default_factory=dict)


def _format_charge(charge: float) -> str:
    return f"{charge:7.4f}"


def _is_hydrogen(atom_type: str) -> bool:
    return atom_type.upper().startswith("H")


def _fragments(atoms: dict[str, Atom]) -> list[set[str]]:
    """Bonded parts of a set of atoms."""
    fragments = []
    unvisited = set(atoms)
    while unvisited:
        frontier = {unvisited.pop()}
        fragment = set(frontier)
        while frontier:
            frontier = {
                other
                for nr in frontier
                for other in atoms[nr].bound_to_nrs
                if other in unvisited
            }
            unvisited -= frontier
            fragment |= frontier
        fragments.append(fragment)
    return fragments


def _residues_by_number(top: Topology, atom_nrs: set[str]) -> dict[str, dict[str, Atom]]:
    """Atoms bonded to these atoms within their residue number, grouped by residue number."""
    residues: dict[str, dict[str, Atom]] = {}
    for nr in sorted(atom_nrs, key=int):
        atom = top.atoms.get(nr)
        if atom is not None:
            residues.setdefault(atom.resnr, {}).update(
                get_residue_by_bonding(atom, top.atoms)
            )
    return residues


def touched_residues(top_a: Topology, top_b: Topology, atom_nrs: set[str]) -> list[_Residue]:
    """The residues containing these atoms, with their atoms in B and their charge in A."""
    residues_a = _residues_by_number(top_a, atom_nrs)
    return [
        _Residue(
            resnr=resnr,
            atoms=atoms,
            fragments=_fragments(atoms),
            charge_a=sum(float(atom.charge) for atom in residues_a.get(resnr, {}).values()),
        )
        for resnr, atoms in _residues_by_number(top_b, atom_nrs).items()
    ]


def moved_hydrogens(
    top_a: Topology, top_b: Topology, atom_nrs: set[str]
) -> dict[str, tuple[Optional[str], Optional[str]]]:
    """Hydrogens among these atoms whose bonding partner changed: {hydrogen: (partner in A, partner in B)}."""
    moved = {}
    for nr in atom_nrs:
        if nr not in top_a.atoms or nr not in top_b.atoms:
            continue
        if not _is_hydrogen(top_b.atoms[nr].type):
            continue
        partner_a = next(iter(top_a.atoms[nr].bound_to_nrs), None)
        partner_b = next(iter(top_b.atoms[nr].bound_to_nrs), None)
        if partner_a != partner_b:
            moved[nr] = (partner_a, partner_b)
    return moved


def local_charges(
    top_a: Topology,
    residue: _Residue,
    moved: dict[str, tuple[Optional[str], Optional[str]]],
    broken: set[str],
) -> dict[str, str]:
    """Charges of a residue from its charges in A by local rules.

    - A moved hydrogen leaves its charge on its old partner (radical
      convention). At its new partner it keeps that charge and the new
      partner gives up the same amount. A hydrogen without a new partner
      gets 0.
    - If the residue is split, each fragment is made integer by
      correcting its atom of a broken bond (broken).

    Atoms the rules do not change keep their charge from A as written there.
    """
    charges = {
        nr: top_a.atoms[nr].charge if nr in top_a.atoms else atom.charge
        for nr, atom in residue.atoms.items()
    }
    # computed charges, formatted once at the end
    changed: dict[str, float] = {}

    def add(nr: str, charge: float) -> None:
        changed[nr] = changed.get(nr, float(charges[nr])) + charge

    for hydrogen, (partner_a, partner_b) in moved.items():
        h_charge = float(top_a.atoms[hydrogen].charge)
        if partner_a in charges:
            add(partner_a, h_charge)
        if hydrogen in charges and partner_b is None:
            changed[hydrogen] = 0.0
        if partner_b in charges:
            add(partner_b, -h_charge)
    if len(residue.fragments) > 1:
        for fragment in residue.fragments:
            broken_atoms = sorted(fragment & broken, key=int)
            if broken_atoms:
                charge = sum(changed.get(nr, float(charges[nr])) for nr in fragment)
                if abs(charge - round(charge)) > 1e-9:
                    add(broken_atoms[0], round(charge) - charge)
    charges.update({nr: _format_charge(charge) for nr, charge in changed.items()})
    return charges


def residue_variants(residue: str) -> list[str]:
    """Residue names that may describe a residue: itself and its terminal variants.

    Assumes the Amber naming N<name> and C<name> for the N- and C-terminal
    variants. This could be improved by reading the force field's .r2b files.
    """
    return list(dict.fromkeys([residue, f"N{residue}", f"C{residue}"]))


def _bonded_in_template(residuetype: ResidueType, name: str) -> set[str]:
    return {other for bond in residuetype.bonds if name in bond for other in bond} - {name}


# how closely a residue has to match an rtp entry, tried in this order
RTP_MATCH_LEVELS = ("exact", "renamed hydrogen", "radical")


def _rtp_scheme(
    top_b: Topology, residue: _Residue, joined: set[str], variant: str, level: str
) -> Optional[_Scheme]:
    """Charges from one rtp entry, if the residue matches it at this level.

    exact: the atom names of the residue are those of the entry.
    renamed hydrogen: as exact, but joined hydrogens may carry other names
        (bind_bond names them after the plain residue name, or HX); each
        takes a free hydrogen of the entry bonded to its partner.
    radical: as renamed hydrogen, but hydrogens of the entry may be missing
        (H abstraction) if each was bonded to one radical atom, which gets
        its charge (radical convention).
    Joined hydrogens get the name and type of their position in the entry.
    """
    residuetype = top_b.ff.residuetypes.get(variant)
    if residuetype is None:
        return None

    # atoms that must carry names of the entry: all of them for an exact
    # match, all but the joined hydrogens otherwise
    renamed = set() if level == "exact" else joined
    named = {atom.atom: atom for nr, atom in residue.atoms.items() if nr not in renamed}
    if len(named) != len(residue.atoms) - len(renamed):
        return None  # repeated atom names
    if not named.keys() <= residuetype.atoms.keys():
        return None
    missing = set(residuetype.atoms) - named.keys()
    charges = {atom.nr: residuetype.atoms[name].charge for name, atom in named.items()}
    identities = {
        nr: (residuetype.atoms[atom.atom].name, residuetype.atoms[atom.atom].type, None)
        for nr, atom in residue.atoms.items()
        if nr in joined and nr not in renamed
    }

    # next: joined hydrogens with other names take free hydrogen positions
    # bonded to their partner
    for hydrogen in sorted(renamed, key=int):
        partners = {top_b.atoms[nr].atom for nr in top_b.atoms[hydrogen].bound_to_nrs}
        options = sorted(
            name
            for name in missing
            if _is_hydrogen(residuetype.atoms[name].type)
            and _bonded_in_template(residuetype, name) & partners
        )
        if not options:
            return None
        spec = residuetype.atoms[options[0]]
        missing.discard(spec.name)
        charges[hydrogen] = spec.charge
        identities[hydrogen] = (spec.name, spec.type, None)

    # next: atoms of the entry still missing; only abstracted hydrogens of a
    # radical are allowed, their charge goes to the radical atom
    if missing and level != "radical":
        return None
    for name in missing:
        radicals = [
            named[other]
            for other in _bonded_in_template(residuetype, name)
            if other in named and named[other].nr in top_b.radicals
        ]
        if not _is_hydrogen(residuetype.atoms[name].type) or len(radicals) != 1:
            return None
        radical = radicals[0].nr
        charges[radical] = _format_charge(
            float(charges[radical]) + float(residuetype.atoms[name].charge)
        )
    return _Scheme(f"rtp:{variant} ({level})", charges, identities)


def _itp_scheme(residue: _Residue, match: ProductTemplateMatch) -> Optional[_Scheme]:
    """Charges, names, types and masses from an exact ITP match covering the residue."""
    if set(match.atom_mapping) != set(residue.atoms):
        return None
    return _Scheme(
        f"itp:{match.template.moleculetype}",
        {nr: atom.charge for nr, atom in match.atom_mapping.items()},
        {nr: (atom.name, atom.atom_type, atom.mass) for nr, atom in match.atom_mapping.items()},
    )


def _schemes(
    top_a: Topology,
    top_b: Topology,
    residue: _Residue,
    moved: dict[str, tuple[Optional[str], Optional[str]]],
    broken: set[str],
    matches: list[ProductTemplateMatch],
) -> list[Callable[[], list[_Scheme]]]:
    """The charge schemes for a residue, most specific first, each computed on demand."""
    joined = {h for h, (_, partner_b) in moved.items() if partner_b in residue.atoms}

    # step 1: the exact ITP match from bind_bond
    def itp() -> list[_Scheme]:
        if len(residue.fragments) > 1:
            return []  # templates describe whole residues, not fragments
        return [s for s in (_itp_scheme(residue, match) for match in matches) if s]

    # step 2: rtp entries; first exact matches, then matches with renamed
    # joined hydrogens, then radicals after H abstraction, each for all
    # residue variants
    def rtp() -> list[_Scheme]:
        if len(residue.fragments) > 1:
            return []
        name = next(atom.residue for nr, atom in residue.atoms.items() if nr not in joined)
        found = []
        for level in RTP_MATCH_LEVELS:
            for variant in residue_variants(name):
                scheme = _rtp_scheme(top_b, residue, joined, variant, level)
                if scheme is not None:
                    found.append(scheme)
        return found

    # step 3: local rules
    def local() -> list[_Scheme]:
        return [_Scheme("local", local_charges(top_a, residue, moved, broken))]

    return [itp, rtp, local]


def _apply(top: Topology, scheme: _Scheme) -> None:
    for nr, charge in scheme.charges.items():
        top.atoms[nr].charge = charge
    for nr, (name, atom_type, mass) in scheme.identities.items():
        atom = top.atoms[nr]
        atom.atom = name
        atom.type = atom_type
        if mass is not None:
            atom.mass = mass
        elif atom_type in top.ff.atomtypes:
            atom.mass = top.ff.atomtypes[atom_type].mass


def assign_recipe_charges(
    top_a: Topology,
    top_b: Topology,
    recipe_steps: list[RecipeStep],
    tolerance: float = DEFAULT_TEMPLATE_TOLERANCE,
) -> None:
    """Set the charges of the residues a recipe touched, once, before parameterization.

    top_a is the topology before the recipe and stays unchanged, top_b the
    topology after all its steps. Each residue containing an atom of a Break
    or Bind takes the first of its schemes (see the module docstring) that
    changes its charge by an integer, within tolerance. A residue no scheme
    fits keeps its charges, with a warning.
    """
    # what the recipe did: atoms of its bonds, atoms of broken bonds, ITP
    # matches found by bind_bond and hydrogens that changed their partner
    bond_steps = [step for step in recipe_steps if isinstance(step, (Break, Bind))]
    atom_nrs = {nr for step in bond_steps for nr in (step.atom_id_1, step.atom_id_2)}
    broken = {
        nr
        for step in bond_steps
        if isinstance(step, Break)
        for nr in (step.atom_id_1, step.atom_id_2)
    }
    matches = []
    for step in bond_steps:
        if isinstance(step, Bind):
            key = tuple(sorted((step.atom_id_1, step.atom_id_2), key=int))
            match = top_b.pending_product_template_matches.pop(key, None)
            if match is not None:
                matches.append(match)
    moved = moved_hydrogens(top_a, top_b, atom_nrs)

    # next: each touched residue takes the first scheme that changes its
    # charge by an integer (steps 1 to 3, see _schemes)
    for residue in touched_residues(top_a, top_b, atom_nrs):
        scheme = _first_fitting_scheme(
            residue, _schemes(top_a, top_b, residue, moved, broken, matches), tolerance
        )
        if scheme is None:
            logger.warning(
                f"No charge scheme changes the charge of residue {residue.resnr} by "
                "an integer; its charges are kept."
            )
            continue
        _apply(top_b, scheme)
        logger.info(f"Charges of residue {residue.resnr} from {scheme.source}.")


def _first_fitting_scheme(
    residue: _Residue,
    scheme_lists: list[Callable[[], list[_Scheme]]],
    tolerance: float,
) -> Optional[_Scheme]:
    """The first scheme that changes the residue's charge by an integer; a failing scheme list is skipped."""
    for scheme_list in scheme_lists:
        try:
            schemes = scheme_list()
        except Exception:
            logger.exception(
                f"Charge scheme {scheme_list.__name__} failed for residue "
                f"{residue.resnr}; trying the next one."
            )
            continue
        for scheme in schemes:
            change = sum(float(charge) for charge in scheme.charges.values()) - residue.charge_a
            if abs(change - round(change)) <= tolerance:
                return scheme
            logger.debug(
                f"{scheme.source} changes the charge of residue {residue.resnr} "
                f"by {change:.4f} e."
            )
    return None


def check_charge_conservation(
    top_a: Topology,
    top_b: Topology,
    tolerance: float = DEFAULT_TEMPLATE_TOLERANCE,
    action: str = "warn",
) -> float:
    """Compare the net charge of the reactive molecule before and after a recipe.

    Returns the change. With action warn, a change larger than tolerance is
    logged as a warning, with error it raises, with off nothing is checked.
    """
    if action not in CONSERVATION_CHECK_ACTIONS:
        raise ValueError(
            f"Unknown charge conservation action {action!r}, expected one of "
            f"{', '.join(CONSERVATION_CHECK_ACTIONS)}."
        )
    if action == "off":
        return 0.0
    charge_a = sum(float(atom.charge) for atom in top_a.atoms.values())
    charge_b = sum(float(atom.charge) for atom in top_b.atoms.values())
    delta = charge_b - charge_a
    if abs(delta) <= tolerance:
        return delta
    m = (
        f"The recipe changed the net charge of the reactive molecule by "
        f"{delta:.4f} e (from {charge_a:.4f} to {charge_b:.4f} e)."
    )
    if action == "error":
        raise ValueError(m)
    logger.warning(m)
    return delta
