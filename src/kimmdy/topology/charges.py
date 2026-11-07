"""Partial charges of the residues a reaction recipe touched.

KIMMDY sets the charges once per recipe segment, the steps up to a Relax or
to the end of the recipe (see Topology.update_all): after its topology
changes (Break, Bind, CustomTopMod) and the names and types
(Topology.update_names_LJ), and before the parameterization, from the
topology before the segment (A) and after it (B). Every residue the
segment touched goes
through a list of charge schemes, most specific first, and takes the first
one that changes its charge by an integer compared to A:

1. an ITP template (TemplateLibrary.match_itp)
2. the rtp entries of `residue_variants`: first exact matches, then matches
   with renamed joined hydrogens, then radicals after H abstraction (the
   charge of a missing hydrogen goes to the radical atom it was bonded to)
3. the local rules (`local_charges`)

The templates, how a residue matches one and which residues the recipe
touched are in kimmdy.topology.templates.

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
from dataclasses import dataclass
from typing import TYPE_CHECKING, Callable, Optional

from kimmdy.recipe import RecipeStep
from kimmdy.topology.templates import (
    DEFAULT_TEMPLATE_TOLERANCE,
    MATCH_LEVELS,
    ResidueMatch,
    ResidueTemplate,
    TouchedResidue,
    joined_hydrogens,
    match_residue,
    moved_hydrogens,
    radical_mismatch,
    recipe_bond_atoms,
    residue_variants,  # noqa: F401, its old place
    touched_residues,
)

if TYPE_CHECKING:
    from kimmdy.topology.atomic import Atom
    from kimmdy.topology.topology import Topology

logger = logging.getLogger("kimmdy.topology.charges")

CONSERVATION_CHECK_ACTIONS = ("warn", "error", "off")


@dataclass(frozen=True)
class _Scheme:
    """Charges a scheme gives a residue.

    charges are written to the topology as they are: template charges and
    unchanged charges verbatim, only computed ones formatted.
    """

    source: str
    charges: dict[str, str]


def _format_charge(charge: float) -> str:
    return f"{charge:7.4f}"


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


def local_charges(
    top_a: Topology,
    residue: TouchedResidue,
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
    fragments = _fragments(residue.atoms)
    if len(fragments) > 1:
        for fragment in fragments:
            broken_atoms = sorted(fragment & broken, key=int)
            if broken_atoms:
                charge = sum(changed.get(nr, float(charges[nr])) for nr in fragment)
                if abs(charge - round(charge)) > 1e-9:
                    add(broken_atoms[0], round(charge) - charge)
    charges.update({nr: _format_charge(charge) for nr, charge in changed.items()})
    return charges


def _template_charges(match: ResidueMatch) -> dict[str, str]:
    """Charges of the template's atoms; a missing hydrogen's charge goes to the radical it was bonded to."""
    atoms = match.template.atoms
    charges = {nr: atoms[name].charge for nr, name in match.names.items()}
    added: dict[str, float] = {}
    for name, radical in match.missing.items():
        added[radical] = added.get(radical, 0.0) + float(atoms[name].charge)
    charges.update({nr: _format_charge(float(charges[nr]) + q) for nr, q in added.items()})
    return charges


def _rtp_scheme(
    top_b: Topology,
    residue: TouchedResidue,
    joined: set[str],
    template: ResidueTemplate,
    level: str,
) -> Optional[_Scheme]:
    """Charges from one rtp entry, if the residue matches it at this level (see match_residue).

    Joined hydrogens may carry other names at the levels renamed hydrogen
    and radical (HX, if no template named them). The radicals of the
    topology must fit the entry (radical_mismatch).
    """
    match = match_residue(residue.atoms, template, level, renamed=joined)
    if match is None:
        return None
    source = f"{template.label} ({level})"
    mismatch = radical_mismatch(match, top_b.radicals)
    if mismatch is not None:
        logger.debug(f"{source} does not fit residue {residue.resnr}: {mismatch}.")
        return None
    return _Scheme(source, _template_charges(match))


def _itp_scheme(
    top_a: Topology, top_b: Topology, residue: TouchedResidue, joined: set[str]
) -> Optional[_Scheme]:
    """Charges from the ITP template that fits the residue (TemplateLibrary.match_itp)."""
    match = top_b.ff.templates.match_itp(residue.atoms, joined, top_b.atoms, top_a.atoms)
    if match is None:
        return None
    return _Scheme(match.template.label, _template_charges(match))


def _schemes(
    top_a: Topology,
    top_b: Topology,
    residue: TouchedResidue,
    moved: dict[str, tuple[Optional[str], Optional[str]]],
    broken: set[str],
    tolerance: float = DEFAULT_TEMPLATE_TOLERANCE,
) -> list[Callable[[], list[_Scheme]]]:
    """The charge schemes for a residue, most specific first, each computed on demand.

    A split residue fits no template (its bonds are not the template's), so
    it gets the local rules.
    """
    joined = joined_hydrogens(residue, moved)

    # step 1: the ITP template
    def itp() -> list[_Scheme]:
        scheme = _itp_scheme(top_a, top_b, residue, joined)
        return [] if scheme is None else [scheme]

    # step 2: rtp entries; first exact matches, then matches with renamed
    # joined hydrogens, then radicals after H abstraction, each for all
    # residue variants
    def rtp() -> list[_Scheme]:
        name = next(atom.residue for nr, atom in residue.atoms.items() if nr not in joined)
        templates = top_b.ff.templates.rtp_variants(name)
        found = []
        for level in MATCH_LEVELS:
            schemes = [
                scheme
                for template in templates
                if (scheme := _rtp_scheme(top_b, residue, joined, template, level))
                is not None
            ]
            # protonation states can fit the same atoms with different
            # charges (CYX and CYM); those that keep the residue charge first
            schemes.sort(key=lambda s: abs(_charge_change(s, residue)) > tolerance)
            found.extend(schemes)
        return found

    # step 3: local rules
    def local() -> list[_Scheme]:
        return [_Scheme("local", local_charges(top_a, residue, moved, broken))]

    return [itp, rtp, local]


def _apply(top: Topology, scheme: _Scheme) -> None:
    for nr, charge in scheme.charges.items():
        top.atoms[nr].charge = charge


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
    # what the recipe did: atoms of its bonds, atoms of broken bonds and
    # hydrogens that changed their partner
    atom_nrs, broken = recipe_bond_atoms(recipe_steps)
    moved = moved_hydrogens(top_a, top_b, atom_nrs)

    # next: each touched residue takes the first scheme that changes its
    # charge by an integer (steps 1 to 3, see _schemes)
    for residue in touched_residues(top_a, top_b, atom_nrs):
        scheme = _first_fitting_scheme(
            residue,
            _schemes(top_a, top_b, residue, moved, broken, tolerance),
            tolerance,
        )
        if scheme is None:
            logger.warning(
                f"No charge scheme changes the charge of residue {residue.resnr} by "
                "an integer; its charges are kept."
            )
            continue
        _apply(top_b, scheme)
        logger.info(f"Charges of residue {residue.resnr} from {scheme.source}.")


def _charge_change(scheme: _Scheme, residue: TouchedResidue) -> float:
    """How much a scheme changes the charge of a residue compared to A."""
    return sum(float(charge) for charge in scheme.charges.values()) - residue.charge_a


def _first_fitting_scheme(
    residue: TouchedResidue,
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
            change = _charge_change(scheme, residue)
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
