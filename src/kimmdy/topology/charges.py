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


def _residue_walks(top: Topology, atom_nrs: set[str]) -> list[dict[str, Atom]]:
    """The residues containing these atoms, found along the bonds, each once."""
    walks: list[dict[str, Atom]] = []
    for nr in sorted(atom_nrs, key=int):
        atom = top.atoms.get(nr)
        if atom is not None and not any(nr in walk for walk in walks):
            walks.append(get_residue_by_bonding(atom, top.atoms))
    return walks


def touched_residues(top_a: Topology, top_b: Topology, atom_nrs: set[str]) -> list[_Residue]:
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
        residues.append(
            _Residue(resnr=resnr, atoms=atoms, fragments=_fragments(atoms), charge_a=charge_a)
        )
    return residues


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


def _bonded_in_template(residuetype: ResidueType, name: str) -> set[str]:
    return {other for bond in residuetype.bonds if name in bond for other in bond} - {name}


def _bonded_to_neighbor_in_template(residuetype: ResidueType, name: str) -> bool:
    """Whether the entry bonds this atom to the previous or next residue.

    The rtp names atoms of the neighbors -C, +N etc. A residue's own atom X
    bonds to a neighbor if X has such a partner, or if the entry lists -X or
    +X (the same atom of the neighbor, e.g. C in -C N).
    """
    names = {n for bond in residuetype.bonds for n in bond}
    return (
        any(other[0] in "-+" for other in _bonded_in_template(residuetype, name))
        or f"-{name}" in names
        or f"+{name}" in names
    )


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
    At every level, each radical of the residue must be one in the entry
    too: bonded there to a missing hydrogen or to a neighboring residue
    (a broken bond between residues). Otherwise the entry describes another
    protonation state, e.g. HID for a HIP that lost HE2.
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

    # next: every radical must be explained by the entry (missing hydrogen
    # or bond to a neighboring residue)
    for name, atom in named.items():
        if atom.nr not in top_b.radicals:
            continue
        if not (
            _bonded_in_template(residuetype, name) & missing
            or _bonded_to_neighbor_in_template(residuetype, name)
        ):
            logger.debug(
                f"rtp:{variant} ({level}) does not fit residue {residue.resnr}: "
                f"radical {name} is not a radical in the entry."
            )
            return None
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
    tolerance: float = DEFAULT_TEMPLATE_TOLERANCE,
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
            schemes = [
                scheme
                for variant in residue_variants(name)
                if (scheme := _rtp_scheme(top_b, residue, joined, variant, level))
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
            residue,
            _schemes(top_a, top_b, residue, moved, broken, matches, tolerance),
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


def _charge_change(scheme: _Scheme, residue: _Residue) -> float:
    """How much a scheme changes the charge of a residue compared to A."""
    return sum(float(charge) for charge in scheme.charges.values()) - residue.charge_a


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
