import logging
import textwrap
from copy import copy, deepcopy
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Optional, Union

from kimmdy.constants import (
    ATOM_ID_FIELDS,
    ATOMTYPE_BONDORDER_FLAT,
    FFFUNC,
    REACTIVE_MOLECULEYPE,
    RESNR_ID_FIELDS,
)
from kimmdy.parsing import TopologyDict, empty_section, read_top, write_top
from kimmdy.plugins import BasicParameterizer, Parameterizer
from kimmdy.recipe import Bind, Break, RecipeStep
from kimmdy.topology.atomic import (
    Angle,
    AngleType,
    Atom,
    Bond,
    BondType,
    Dihedral,
    DihedralType,
    DihedralRestraint,
    Exclusion,
    ImproperDihedralId,
    MoleculeTypeHeader,
    MultipleDihedrals,
    Pair,
    PositionRestraint,
    ResidueImproperSpec,
    combine_periodic_terms,
    ResidueType,
    Settle,
)
from kimmdy.topology.ff import FF
from kimmdy.topology.product_templates import ProductTemplateMatch
from kimmdy.topology.utils import (
    attributes_to_list,
    get_moleculetype_atomics,
    get_moleculetype_header,
    get_residue_by_bonding,
    get_residue_fragments,
    get_top_section,
    increment_field,
    is_not_solvent_or_ion,
    match_atomic_item_to_atomic_type,
    set_top_section,
)
from kimmdy.utils import TopologyAtomAddress

logger = logging.getLogger("kimmdy.topology")


@dataclass(frozen=True)
class ParameterizationReport:
    event_focus_ids: set[str]
    active_radical_ids: set[str]
    restored_templates: int


class MoleculeType:
    """One moleculetype in the topology

    Attributes
    ----------
    atoms : dict[str, Atom]
    bonds : dict[tuple[str, str], Bond]
    pairs : dict[tuple[str, str], Pair]
    angles : dict[tuple[str, str, str], Angle]
    proper_dihedrals : dict[tuple[str, str, str, str], MultipleDihedrals]
    improper_dihedrals : dict[tuple[str, str, str, str], Dihedral]
    position_restraints : dict[str, PositionRestraint]
    dihedral_restraints : dict[tuple[str, str, str, str], DihedralRestraint]
    radicals : dict[str, Atom]
        dict mapping atom indices to atom objects for storing all radical atoms

    """

    def __init__(
        self,
        header: MoleculeTypeHeader,
        atomics: dict,
        radicals: Optional[str] = None,
        combine_same_periodicity: bool = True,
        handle_duplicate_lines: bool = True,
    ) -> None:
        # merge funct-4 terms that share atoms and periodicity while parsing
        self.combine_same_periodicity = combine_same_periodicity
        # merge repeated explicit proper dihedral lines and stop on any other
        # repeated line instead of silently keeping only the last one
        self.handle_duplicate_lines = handle_duplicate_lines
        self.name = header.name
        self.nrexcl = header.nrexcl
        self.atomics = atomics

        logger.debug(f"parsing molecule {self.name}")

        self.atoms: dict[str, Atom] = {}
        self.bonds: dict[tuple[str, str], Bond] = {}
        self.pairs: dict[tuple[str, str], Pair] = {}
        self.angles: dict[tuple[str, str, str], Angle] = {}
        self.proper_dihedrals: dict[tuple[str, str, str, str], MultipleDihedrals] = {}
        self.improper_dihedrals: dict[tuple[str, str, str, str], MultipleDihedrals] = {}
        self.position_restraints: dict[str, PositionRestraint] = {}
        self.dihedral_restraints: dict[tuple[str, str, str, str], DihedralRestraint] = (
            {}
        )
        self.settles: dict[str, Settle] = {}
        self.exclusions: dict[tuple, Exclusion] = {}

        self._parse_atoms()
        self._parse_bonds()
        self._parse_pairs()
        self._parse_angles()
        self._parse_dihedrals()
        self._parse_restraints()
        self._parse_settles()
        self._parse_exclusions()
        self._initialize_graph()

        # must be after self._parse_atoms
        self.radicals: dict[str, Atom] = {}
        if radicals is not None:
            if self.name == "Reactive" and len(radicals) > 0:
                logger.debug(
                    f"Using 'radicals' section from config file with entries: '{radicals}'."
                )
                for radical in radicals.strip().split(sep=" "):
                    atom = self.atoms.get(radical)
                    if atom:
                        self.radicals[radical] = atom
                        atom.is_radical = True
                    else:
                        logger.debug(
                            f"Atom nr {radical} from 'radicals' section in "
                            "config file not in topology. Ignoring this entry."
                        )
        else:
            self.find_radicals()

    def __str__(self) -> str:
        return textwrap.dedent(
            f"""\
        Moleculetype {self.name} with:
        {len(self.atoms)} atoms,
        {len(self.bonds)} bonds,
        {len(self.angles)} angles,
        {len(self.pairs)} pairs,
        {len(self.settles)} settles,
        {len(self.exclusions)} exclusions,
        {len(self.proper_dihedrals)} proper dihedrals
        {len(self.improper_dihedrals)} improper dihedrals
        {len(self.position_restraints)} position restraints
        {len(self.dihedral_restraints)} dihedral restraints
        """
        )

    def __repr__(self) -> str:
        return f"Molecule({(self.name, self.nrexcl)}, {self.atomics})"

    def _repr_pretty_(self, p, _):
        """A __repr__ for ipython.

        This whill be used if just the name of the object is entered in the ipython shell
        or a jupyter notebook.

        p is an instance of
        [IPython.lib.pretty.PrettyPrinter](https://ipython.org/ipython-doc/3/api/generated/IPython.lib.pretty.html)
        """
        p.text(str(self))

    def _parse_atoms(self):
        """Parse atoms from topology dictionary."""
        ls = self.atomics.get("atoms")
        if ls is None:
            raise ValueError(f"No atoms found in topology for molecule {self.name}.")
        for l in ls:
            atom = Atom.from_top_line(l)
            self.atoms[atom.nr] = atom

    def _parse_bonds(self):
        """Parse bond from topology dictionary."""
        ls = self.atomics.get("bonds")
        if ls is None:
            return
        for l in ls:
            bond = Bond.from_top_line(l)
            key = (bond.ai, bond.aj)
            # a repeated line stops with an error (handle_duplicate_lines)
            if key in self.bonds:
                self._reject_duplicate_line("bonds", key)
            self.bonds[key] = bond

    def _parse_pairs(self):
        """Parse pairs from topology dictionary."""
        ls = self.atomics.get("pairs")
        if ls is None:
            return
        for l in ls:
            pair = Pair.from_top_line(l)
            key = (pair.ai, pair.aj)
            # a repeated line stops with an error (handle_duplicate_lines)
            if key in self.pairs:
                self._reject_duplicate_line("pairs", key)
            self.pairs[key] = pair

    def _parse_angles(self):
        """Parse angles from topology dictionary."""
        ls = self.atomics.get("angles")
        if ls is None:
            return
        for l in ls:
            angle = Angle.from_top_line(l)
            key = (angle.ai, angle.aj, angle.ak)
            # a repeated line stops with an error (handle_duplicate_lines)
            if key in self.angles:
                self._reject_duplicate_line("angles", key)
            self.angles[key] = angle

    def _parse_dihedrals(self):
        """Read the [ dihedrals ] lines of this molecule into proper and improper dihedrals.

        Each dihedral is stored under its four atoms, with one term per
        periodicity. If a line repeats atoms and periodicity of an earlier
        line, the two terms are merged into one equivalent term where that is
        exact (see combine_periodic_terms); any other repeat stops with an
        error (see _reject_duplicate_line), because GROMACS would use both
        lines but only one can be stored.
        """
        ls = self.atomics.get("dihedrals")
        if ls is None:
            return
        for l in ls:
            dihedral = Dihedral.from_top_line(l)
            key = (dihedral.ai, dihedral.aj, dihedral.ak, dihedral.al)
            # periodic improper dihedrals (function 4)
            if dihedral.funct == FFFUNC["mult_improper_dihedral"]:
                if self.improper_dihedrals.get(key) is None:
                    self.improper_dihedrals[key] = MultipleDihedrals(
                        *key, dihedral.funct, dihedrals={}
                    )
                existing = self.improper_dihedrals[key].dihedrals.get(
                    dihedral.periodicity
                )
                if (
                    self.combine_same_periodicity
                    and existing is not None
                    and None not in (existing.c0, existing.c1, dihedral.c0, dihedral.c1)
                ):
                    # same atoms and periodicity listed twice (a99SB-disp psi):
                    # keep both contributions instead of the later one only
                    dihedral.c0, dihedral.c1 = combine_periodic_terms(
                        [(existing.c0, existing.c1), (dihedral.c0, dihedral.c1)]
                    )
                elif existing is not None:
                    self._reject_duplicate_line(
                        "dihedrals (improper)", key, dihedral.periodicity
                    )
                self.improper_dihedrals[key].dihedrals[dihedral.periodicity] = dihedral
            # all other dihedrals (proper dihedrals, function 1 or 9)
            else:
                if self.proper_dihedrals.get(key) is None:
                    self.proper_dihedrals[key] = MultipleDihedrals(
                        *key, dihedral.funct, dihedrals={}
                    )
                existing = self.proper_dihedrals[key].dihedrals.get(
                    dihedral.periodicity
                )
                if existing is not None:
                    if (
                        self.handle_duplicate_lines
                        and dihedral.funct in ("1", "9")
                        and existing.funct == dihedral.funct
                        and None not in (existing.c0, existing.c1, dihedral.c0, dihedral.c1)
                    ):
                        # same atoms and periodicity listed twice with numbers:
                        # GROMACS adds both terms, so merge them into one
                        # equivalent term instead of keeping the later one only
                        dihedral.c0, dihedral.c1 = combine_periodic_terms(
                            [(existing.c0, existing.c1), (dihedral.c0, dihedral.c1)]
                        )
                    else:
                        self._reject_duplicate_line(
                            "dihedrals (proper)", key, dihedral.periodicity
                        )
                self.proper_dihedrals[key].dihedrals[dihedral.periodicity] = dihedral

    def _reject_duplicate_line(self, section: str, key: tuple, periodicity: str = ""):
        """Stop with an error when a topology line repeats an earlier one that cannot be merged.

        KIMMDY keeps one term per atom set (and per periodicity for
        dihedrals), while GROMACS uses every line, so silently keeping only
        the last line would change the simulation. The error names the
        molecule, section and atoms. With handle_duplicate_lines off it does
        nothing, so the later line replaces the earlier as before.
        """
        if not self.handle_duplicate_lines:
            return
        where = f" with periodicity {periodicity}" if periodicity else ""
        m = (
            f"Molecule {self.name}: [ {section} ] lists atoms {' '.join(key)}{where} "
            "more than once. KIMMDY can keep only one such term, but GROMACS would "
            "apply every line, so the topology cannot be represented without "
            "changing it. Merge or remove the repeated line, or set "
            "changer.topology.handle_duplicate_topology_lines to false to keep only the last one."
        )
        logger.error(m)
        raise ValueError(m)

    def _parse_restraints(self):
        """Parse restraints from topology dictionary."""
        ls = self.atomics.get("position_restraints")
        if not ls:
            ls = []
        condition = None
        for l in ls:
            restraint = PositionRestraint.from_top_line(l, condition=condition)
            self.position_restraints[restraint.ai] = restraint
        ls = self.atomics.get("dihedral_restraints")
        if not ls:
            ls = []
        for l in ls:
            restraint = DihedralRestraint.from_top_line(l)
            self.dihedral_restraints[
                (restraint.ai, restraint.aj, restraint.ak, restraint.al)
            ] = restraint

    def _parse_settles(self):
        """Parse settles from topology dictionary."""
        ls = self.atomics.get("settles")
        if not ls:
            return
        for l in ls:
            settle = Settle.from_top_line(l)
            self.settles[settle.nr] = settle

    def _parse_exclusions(self):
        """Parse exclusions from topology dictionary."""
        ls = self.atomics.get("exclusions")
        if not ls:
            return
        for _, l in enumerate(ls):
            exclusion = Exclusion.from_top_line(l)
            self.exclusions[exclusion.key()] = exclusion

    def _initialize_graph(self):
        """Add a list of atom nrs bound to an atom to each atom."""
        for bond in self.bonds.values():
            i = bond.ai
            j = bond.aj
            self.atoms[i].bound_to_nrs.append(j)
            self.atoms[j].bound_to_nrs.append(i)

    def find_radicals(self):
        """Find atoms that are radicals and update self.radicals.

        Iterate over all atoms and designate them as radicals if they have
        fewer bounds than their natural bond order.
        """
        logger.debug(
            "Trying to infer radical status based on AMBER(!!) atomtype bond "
            f"order in molecule {self.name}."
        )
        self.radicals = {}
        error_atoms = set()
        for atom in self.atoms.values():
            bo = ATOMTYPE_BONDORDER_FLAT.get(atom.type)
            if bo is None:
                error_atoms.add(atom.type)
            if bo and bo > len(atom.bound_to_nrs):
                self.radicals[atom.nr] = atom
                atom.is_radical = True
            else:
                atom.is_radical = False
        if len(error_atoms) > 0:
            logger.warning(
                "Some atomtypes not found in AMBER atomtypes. Cannot infer radical status.\n"
                "Missing atom types:\n"
                f"{error_atoms}"
            )
        return None

    def _update_atomics_dict(self):
        self.atomics["atoms"] = [attributes_to_list(x) for x in self.atoms.values()]
        self.atomics["bonds"] = [attributes_to_list(x) for x in self.bonds.values()]
        self.atomics["angles"] = [attributes_to_list(x) for x in self.angles.values()]
        self.atomics["pairs"] = [attributes_to_list(x) for x in self.pairs.values()]
        self.atomics["dihedrals"] = [
            attributes_to_list(x)
            for dihedrals in self.proper_dihedrals.values()
            for x in dihedrals.dihedrals.values()
        ] + [
            attributes_to_list(x)
            for dihedrals in self.improper_dihedrals.values()
            for x in dihedrals.dihedrals.values()
        ]
        self.atomics["position_restraints"] = [
            attributes_to_list(x) for x in self.position_restraints.values()
        ]
        self.atomics["dihedral_restraints"] = [
            attributes_to_list(x) for x in self.dihedral_restraints.values()
        ]
        self.atomics["settles"] = [attributes_to_list(x) for x in self.settles.values()]
        self.atomics["exclusions"] = [
            attributes_to_list(x) for x in self.exclusions.values()
        ]

    def _get_atom_bonds(self, atom_nr: str) -> list[tuple[str, str]]:
        """Get all bonds a particular atom is involved in."""
        ai = atom_nr
        bonds = []
        for aj in self.atoms[ai].bound_to_nrs:
            if int(ai) < int(aj):
                bonds.append((ai, aj))
            else:
                bonds.append((aj, ai))
        return bonds

    def _get_atom_pairs(self, _: str) -> list[tuple[str, str]]:
        raise NotImplementedError(
            "get_atom_pairs is not implemented. Get the pairs as the endpoints of dihedrals instead."
        )

    def _get_atom_angles(self, atom_nr: str) -> list[tuple[str, str, str]]:
        """
        Gets a list of angles that one atom is involved in based on the bond
        information stored in the atom.

        Angles between atoms ai, aj, ak satisfy ai < ak
        """
        return self._get_center_atom_angles(atom_nr) + self._get_margin_atom_angles(
            atom_nr
        )

    def _get_center_atom_angles(self, atom_nr: str) -> list[tuple[str, str, str]]:
        # atom_nr in the middle of an angle
        angles = []
        aj = atom_nr
        partners = self.atoms[aj].bound_to_nrs
        for ai in partners:
            for ak in partners:
                if int(ai) < int(ak):
                    angles.append((ai, aj, ak))
        return angles

    def _get_margin_atom_angles(self, atom_nr: str) -> list[tuple[str, str, str]]:
        # atom_nr at the outer corner of angle
        angles = []
        ai = atom_nr
        for aj in self.atoms[ai].bound_to_nrs:
            for ak in self.atoms[aj].bound_to_nrs:
                if ai == ak:
                    continue
                if int(ai) < int(ak):
                    angles.append((ai, aj, ak))
                else:
                    angles.append((ak, aj, ai))
        return angles

    def _get_atom_proper_dihedrals(
        self, atom_nr: str
    ) -> list[tuple[str, str, str, str]]:
        """
        Gets a list of dihedrals that one atom is involved in based on the bond
        information stored in the atom.
        """
        return self._get_center_atom_dihedrals(
            atom_nr
        ) + self._get_margin_atom_dihedrals(atom_nr)

    def _get_center_atom_dihedrals(
        self, atom_nr: str, improper: bool = False
    ) -> list[tuple[str, str, str, str]]:
        """
        propers are ordered by the middle atoms j,k
        impropers are not ordered
        """
        dihedrals = []
        aj = atom_nr
        partners = self.atoms[aj].bound_to_nrs
        for ai in partners:
            for ak in partners:
                if ai == ak:
                    continue
                for al in self.atoms[ak].bound_to_nrs:
                    if int(aj) < int(ak) or improper:
                        key = (ai, aj, ak, al)
                    else:
                        key = (al, ak, aj, ai)
                    if len(key) == len(set(key)):
                        # one atom can not be in two positions
                        dihedrals.append(key)
        return dihedrals

    def _get_margin_atom_dihedrals(
        self, atom_nr: str, improper: bool = False
    ) -> list[tuple[str, str, str, str]]:
        """
        propers are ordered by the middle atoms j,k
        impropers are not ordered
        """
        dihedrals = []
        ai = atom_nr
        for aj in self.atoms[ai].bound_to_nrs:
            for ak in self.atoms[aj].bound_to_nrs:
                if ai == ak:
                    continue
                for al in self.atoms[ak].bound_to_nrs:
                    if int(aj) < int(ak) or improper:
                        key = (ai, aj, ak, al)
                    else:
                        key = (al, ak, aj, ai)
                    if len(key) == len(set(key)):
                        # one atom can not be in two positions
                        dihedrals.append(key)
        return dihedrals

    def _get_atom_improper_dihedrals(
        self, atom_nr: str, ff: FF
    ) -> dict[ImproperDihedralId, list[ResidueImproperSpec]]:
        """Find the residue-template impropers that contain one atom, with real atom numbers.

        Looks up the atom's residue template, translates the template's atom
        names (including -C, the previous residue's C, and +N, the next
        residue's N) into atom numbers, and returns every improper term that
        contains the atom: {four atom numbers: list of terms, one per
        periodicity}.
        <https://manual.gromacs.org/current/reference-manual/functions/bonded-interactions.html#improper-dihedrals>
        """

        atom = self.atoms[atom_nr]
        residuetype = ff.residuetypes.get(atom.residue)
        if not isinstance(residuetype, ResidueType):
            return {}

        # find the previous and next residue (needed for -C and +N)
        residue = get_residue_by_bonding(atom=atom, atoms=self.atoms)
        previous_atom, next_atom = self._get_neighbor_residue_atoms(residue)
        if previous_atom is not None:
            previous_residue = get_residue_by_bonding(
                atom=previous_atom, atoms=self.atoms
            )
            previous_residuetype = ff.residuetypes.get(previous_atom.residue)
        else:
            previous_residue = {}
            previous_residuetype = None
        if next_atom is not None:
            next_residue = get_residue_by_bonding(atom=next_atom, atoms=self.atoms)
            next_residuetype = ff.residuetypes.get(next_atom.residue)
        else:
            next_residue = {}
            next_residuetype = None

        atomname_to_nr = {atom.atom: atom_nr for atom_nr, atom in residue.items()}
        # residues on aminoacids.rtp specify a dihedral to the next or previous
        # AA with -C and +N as the atomname
        prev_atomname_to_nr = {
            atom.atom: atom_nr for atom_nr, atom in previous_residue.items()
        }
        next_atomname_to_nr = {
            atom.atom: atom_nr for atom_nr, atom in next_residue.items()
        }

        dihedrals: dict[ImproperDihedralId, list[ResidueImproperSpec]] = {}
        if previous_residuetype is not None and atom.atom == "N":
            # the atom is N, so it may have an improper defined
            # in the previous residue, where it would be referred to as +N
            for key, impropers in previous_residuetype.improper_dihedrals.items():
                # find the dihedral with +N
                if any(["+" in k for k in key]):
                    nr_key: ImproperDihedralId = tuple(
                        atom.nr if "+" in x else prev_atomname_to_nr.get(x) for x in key
                    )  # pyright: ignore
                    if None in nr_key:
                        m = f"Improper dihedral {key} with {nr_key} not found in residue {atom.residue}."
                        logger.debug(m)
                        continue
                    dihedrals[nr_key] = impropers

        if next_residuetype is not None and atom.atom == "C":
            # the atom is C, so it may have an improper defined
            # in the previous residue, where it would be referred to as -C
            for key, impropers in next_residuetype.improper_dihedrals.items():
                # find the dihedral with -C
                if any(["-" in k for k in key]):
                    nr_key = tuple(
                        atom.nr if "-" in x else next_atomname_to_nr.get(x) for x in key
                    )  # pyright: ignore
                    if None in nr_key:
                        m = f"Improper dihedral {key} with {nr_key} not found in residue {atom.residue}."
                        logger.debug(m)
                        continue
                    dihedrals[nr_key] = impropers

        def get_number(name):
            """Get the atom number for an atom name from a dihedral specification"""
            if "+" in name:
                # + refers to the next residue
                return next_atomname_to_nr.get(name.replace("+", ""), None)
            elif "-" in name:
                # - refers to the previous residue
                return prev_atomname_to_nr.get(name.replace("-", ""), None)
            else:
                # regular atom name
                return atomname_to_nr.get(name, None)

        # regularly defined improper dihedrals of the same residue
        for key, impropers in residuetype.improper_dihedrals.items():
            nr_key = tuple(get_number(x) for x in key)  # pyright: ignore
            if atom_nr not in nr_key:
                continue
            if None in nr_key:
                m = f"Improper dihedral {key} {nr_key} not found in residue {atom.residue}."
                logger.debug(m)
                continue
            dihedrals[nr_key] = impropers

        # gromacs sorts the keys such that the first atom is the lowest numbered
        new_dihedrals: dict[ImproperDihedralId, list[ResidueImproperSpec]] = {}
        for key, impropers in dihedrals.items():
            if int(key[0]) > int(key[3]):
                key = (key[3], key[2], key[1], key[0])
                new_dihedrals[key] = [improper.reversed() for improper in impropers]
            else:
                new_dihedrals[key] = list(impropers)

        return new_dihedrals

    def _get_neighbor_residue_atoms(
        self, residue: dict[str, Atom]
    ) -> tuple[Optional[Atom], Optional[Atom]]:
        """Find one atom of the previous residue and one of the next residue, through bonds.

        The previous residue is the one whose C is bonded to this residue's N;
        the next residue is the one whose N is bonded to this residue's C.
        Atom numbers cannot be used for this: an H that moved into this residue
        keeps its old, possibly distant number, so "highest number + 1" can
        point into another chain. Returns (atom of previous residue, atom of
        next residue), None where there is none. Residues without atoms named
        N or C (not protein) keep the old atom-number rule.

        ToDo: Generalise the protein recognition, maybe using residuetypes.dat
        """

        def across(own: str, other: str) -> Optional[Atom]:
            """Follow a bond from this residue's atom named own to an atom named other outside the residue; None if there is no such bond."""
            for atom in residue.values():
                if atom.atom != own:
                    continue
                for nr in atom.bound_to_nrs:
                    if nr not in residue and self.atoms[nr].atom == other:
                        return self.atoms[nr]
            return None

        # protein residue: previous through N-C, next through C-N
        names = {atom.atom for atom in residue.values()}
        if "N" in names or "C" in names:
            return across("N", "C"), across("C", "N")

        # anything else: the atoms numbered just below and above the residue
        lowest_atomnr = min(residue.keys(), key=int)
        highest_atomnr = max(residue.keys(), key=int)
        return (
            self.atoms.get(str(int(lowest_atomnr) - 1)),
            self.atoms.get(str(int(highest_atomnr) + 1)),
        )

    def _regenerate_topology_from_bound_to(self, ff):
        # clear all bonds, angles, dihedrals
        self.bonds.clear()
        self.angles.clear()
        self.proper_dihedrals.clear()
        self.proper_dihedrals.clear()
        self.improper_dihedrals.clear()

        # bonds
        keys = []
        for atom in self.atoms.values():
            keys = self._get_atom_bonds(atom.nr)
            for key in keys:
                self.bonds[key] = Bond(key[0], key[1], FFFUNC["harmonic_bond"])

        # angles
        for atom in self.atoms.values():
            keys = self._get_atom_angles(atom.nr)
            for key in keys:
                self.angles[key] = Angle(
                    key[0], key[1], key[2], FFFUNC["harmonic_angle"]
                )

        # dihedrals and pass
        for atom in self.atoms.values():
            keys = self._get_atom_proper_dihedrals(atom.nr)
            for key in keys:
                self.proper_dihedrals[key] = MultipleDihedrals(
                    *key,
                    FFFUNC["mult_proper_dihedral"],
                    dihedrals={"": Dihedral(*key, FFFUNC["mult_proper_dihedral"])},
                )
                pairkey: tuple[str, str] = tuple(
                    str(x) for x in sorted([key[0], key[3]], key=int)
                )  # type: ignore
                if self.pairs.get(pairkey) is None:
                    self.pairs[pairkey] = Pair(pairkey[0], pairkey[1], FFFUNC["pair"])

        # impropers from the residue templates: all their terms, or only the
        # last one when keep_all_improper_periodicities is off
        for atom in self.atoms.values():
            impropers = self._get_atom_improper_dihedrals(atom.nr, ff)
            keep_all = getattr(ff, "keep_all_improper_periodicities", True)
            for key, specs in impropers.items():
                dihedrals = {}
                for improper in specs if keep_all else specs[-1:]:
                    periodicity = "" if improper.c2 is None else improper.c2
                    # legacy behaviour stored the single term under key ""
                    dihedrals[periodicity if keep_all else ""] = Dihedral(
                        *key,
                        FFFUNC["mult_improper_dihedral"],
                        c0=improper.c0,
                        c1=improper.c1,
                        periodicity=periodicity,
                    )
                self.improper_dihedrals[key] = MultipleDihedrals(
                    *key,
                    FFFUNC["mult_improper_dihedral"],
                    dihedrals=dihedrals,
                )

    def reindex_atomnrs(self) -> dict[str, str]:
        """Reindex atom numbers in topology.

        Starts at index 1.
        This also updates the numbers for bonds, angles, dihedrals and pairs.
        Returns a dict, mapping of old atom number strings to new ones
        """

        update_map = {
            atom_nr: str(i + 1)
            for i, atom_nr in enumerate(sorted(list(self.atoms.keys()), key=int))
        }

        new_atoms = {}
        for atom in self.atoms.values():
            atom.nr = update_map[atom.nr]
            bound_to = []
            for nr in atom.bound_to_nrs:
                new_nr = update_map.get(nr)
                if new_nr is not None:
                    bound_to.append(new_nr)
            atom.bound_to_nrs = bound_to
            new_atoms[atom.nr] = atom
        self.atoms = new_atoms

        new_bonds = {}
        for bond in self.bonds.values():
            ai = update_map.get(bond.ai)
            aj = update_map.get(bond.aj)
            # drop bonds to a deleted atom
            if ai is None or aj is None:
                continue
            bond.ai = ai
            bond.aj = aj
            new_bonds[(ai, aj)] = bond
        self.bonds = new_bonds

        new_angles = {}
        for angle in self.angles.values():
            ai = update_map.get(angle.ai)
            aj = update_map.get(angle.aj)
            ak = update_map.get(angle.ak)
            # drop angles to a deleted atom
            if None in (ai, aj, ak):
                continue
            # pyright does not grok `None in ...`
            angle.ai = ai  # type: ignore
            angle.aj = aj  # type: ignore
            angle.ak = ak  # type: ignore
            new_angles[(angle.ai, angle.aj, angle.ak)] = angle
        self.angles = new_angles

        new_pairs = {}
        new_multiple_dihedrals = {}
        for dihedrals in self.proper_dihedrals.values():
            new_dihedrals = {}
            ai = update_map.get(dihedrals.ai)
            aj = update_map.get(dihedrals.aj)
            ak = update_map.get(dihedrals.ak)
            al = update_map.get(dihedrals.al)
            # drop dihedrals to a deleted atom
            if None in (ai, aj, ak, al):
                continue

            # do pairs before the dihedrals are updated
            pair = self.pairs.pop((dihedrals.ai, dihedrals.al), None)
            if pair is not None:
                pair_ai = update_map.get(pair.ai)
                pair_aj = update_map.get(pair.aj)

                if None not in (pair_ai, pair_aj):
                    pair.ai = pair_ai  # type: ignore (pyright bug)
                    pair.aj = pair_aj  # type: ignore
                    new_pairs[(pair_ai, pair_aj)] = pair

            dihedrals.ai = ai  # type: ignore
            dihedrals.aj = aj  # type: ignore
            dihedrals.ak = ak  # type: ignore
            dihedrals.al = al  # type: ignore

            for dihedral in dihedrals.dihedrals.values():
                dihedral.ai = ai  # type: ignore
                dihedral.aj = aj  # type: ignore
                dihedral.ak = ak  # type: ignore
                dihedral.al = al  # type: ignore
                new_dihedrals[dihedral.periodicity] = dihedral
            dihedrals.dihedrals = new_dihedrals

            new_multiple_dihedrals[
                (
                    dihedrals.ai,
                    dihedrals.aj,
                    dihedrals.ak,
                    dihedrals.al,
                )
            ] = dihedrals

        self.proper_dihedrals = new_multiple_dihedrals
        self.pairs = new_pairs

        new_multiple_dihedrals = {}
        for dihedrals in self.improper_dihedrals.values():
            new_dihedrals = {}
            ai = update_map.get(dihedrals.ai)
            aj = update_map.get(dihedrals.aj)
            ak = update_map.get(dihedrals.ak)
            al = update_map.get(dihedrals.al)
            # drop dihedrals to a deleted atom
            if None in (ai, aj, ak, al):
                continue
            dihedrals.ai = ai  # type: ignore
            dihedrals.aj = aj  # type: ignore
            dihedrals.ak = ak  # type: ignore
            dihedrals.al = al  # type: ignore

            for dihedral in dihedrals.dihedrals.values():
                dihedral.ai = ai  # type: ignore
                dihedral.aj = aj  # type: ignore
                dihedral.ak = ak  # type: ignore
                dihedral.al = al  # type: ignore
                new_dihedrals[dihedral.periodicity] = dihedral
            dihedrals.dihedrals = new_dihedrals

            new_multiple_dihedrals[(ai, aj, ak, al)] = dihedrals
        self.improper_dihedrals = new_multiple_dihedrals

        return update_map


class Topology:
    """Smart container for parsed topology data.

    A topology keeps track of connections when bonds are broken or formed.
    Reparametrization is triggerd automatically if `to_dict` is called
    after bonds have changed.

    Also see <https://manual.gromacs.org/current/reference-manual/topologies/topology-file-formats.html#topology-file>

    Assumptions:

    - the topology of interest (the Reactive moleculetype) consists of the first moleculetypes (non-solvent).

    Parameters
    ----------
    top
        A dictionary containing the parsed topology data, produced by
        [](`kimmdy.parsing.read_top`)
    parametrizer
        The parametrizer to use when reparametrizing the topology.
    is_reactive_predicate_f
        A function that takes a moleculetype name and returns True if the moleculetype
        should be merged into the reactive moleculetype.
    radicals
        A string of atom numbers that are radicals.
    residuetypes_path
        Path to the residue types file.
    reactive_nrexcl
        Overwrite nrexcl value for the reactive moleculetype. Otherwise takes the nrexcl of the first reactive moleculetype.
    needs_parameterization
        Does the topology currently need to be re-parameterized due to changes?
    parameterization_focus_ids
        list of atoms ids around which the parameterization happens (to avoid re-parameterizing the whole Reactive moleculetype)
    """

    @classmethod
    def from_path(cls, path: str | Path, ffdir: str | Path | None = None, **kwargs):
        """Load a Topology object from a file path.

        Simplifies the `Topology(read_top(...), ...)` pattern.

        Parameters
        ----------
        path
            Path to the topology file.
        kwargs
            Additional keyword arguments to pass to the Topology constructor.

        Returns
        -------
        Topology
        """
        if isinstance(path, str):
            path = Path(path)
        if isinstance(ffdir, str):
            ffdir = Path(ffdir)
        return cls(read_top(path=path, ffdir=ffdir), **kwargs)

    def to_path(self, path: str | Path) -> None:
        """Write the topology to a file.

        Parameters
        ----------
        path
            Path to the topology file.
        """
        if isinstance(path, str):
            path = Path(path)
        write_top(self.to_dict(), path)

    def __init__(
        self,
        top: TopologyDict,
        parametrizer: Parameterizer = BasicParameterizer(),
        is_reactive_predicate_f: Callable[[str], bool] = is_not_solvent_or_ion,
        radicals: Optional[str] = None,
        residuetypes_path: Optional[Path] = None,
        reactive_nrexcl: Optional[str] = None,
        product_template_sources: tuple[str, ...] = ("rtp",),
        product_template_itp_glob: str = "*.itp",
        product_template_transitions: Optional[dict[str, str]] = None,
        keep_all_improper_periodicities: bool = True,
        combine_same_periodicity_dihedrals: bool = True,
        conserve_transferred_h_charge: bool = True,
        handle_duplicate_topology_lines: bool = True,
    ):
        """Topology fixes, all on by default; set to False for the 16ce888 behaviour:

        keep_all_improper_periodicities
            Keep every periodic term a residue template lists for one atom
            quadruple when impropers are rebuilt after a reaction.
        combine_same_periodicity_dihedrals
            Merge funct-4 terms with the same atoms and periodicity into one
            equivalent term (forces unchanged, energy shifted by a constant)
            instead of keeping only the last one.
        conserve_transferred_h_charge
            A transferred H keeps its charge and the acceptor gives up the
            same amount, instead of both taking residue-template charges.
        handle_duplicate_topology_lines
            Merge repeated explicit proper dihedral lines (function 1 or 9)
            with the same atoms and periodicity, and stop with an error on any
            other repeated line, instead of keeping only the last one.
        """
        self.combine_same_periodicity_dihedrals = combine_same_periodicity_dihedrals
        self.handle_duplicate_topology_lines = handle_duplicate_topology_lines
        if top == {}:
            raise NotImplementedError(
                "Generating an empty Topology from an empty TopologyDict is not implemented."
            )

        self.top = top
        self.moleculetypes: dict[str, MoleculeType] = {}
        molecules = get_top_section(top, "molecules")
        if molecules is None:
            raise ValueError("molecules not found in top file")
        self.molecules = [(l[0], l[1]) for l in molecules]

        self.ff = FF(
            top,
            residuetypes_path,
            product_template_sources=product_template_sources,
            product_template_itp_glob=product_template_itp_glob,
            product_template_transitions=product_template_transitions,
            keep_all_improper_periodicities=keep_all_improper_periodicities,
            combine_same_periodicity_dihedrals=combine_same_periodicity_dihedrals,
        )
        self.pending_product_template_matches: dict[
            tuple[str, str], ProductTemplateMatch
        ] = {}
        self.pending_transferred_h_charges: dict[str, str] = {}
        self.conserve_transferred_h_charge = conserve_transferred_h_charge
        self.pending_product_template_restorations: list[
            ProductTemplateMatch
        ] = []
        if self.handle_duplicate_topology_lines:
            self._expand_parameter_macros()
        self._parse_molecules()
        self.parametrizer = parametrizer
        self._check_is_reactive_molecule = is_reactive_predicate_f

        self.needs_parameterization: bool = False
        self.parameterization_focus_ids: set[str] = set()

        self._merge_moleculetypes(radicals=radicals, nrexcl=reactive_nrexcl)
        self._link_atomics()

    def _link_atomics(self):
        """link atoms, bonds etc. properties of self (=top) to the reactive moleculeype.

        Call this any time a large change is made to the reactive moleculetype
        that re-assigns a property and not just modifies parts of it in place.
        """
        self.reactive_molecule = self.moleculetypes[REACTIVE_MOLECULEYPE]
        self.atoms = self.reactive_molecule.atoms
        self.bonds = self.reactive_molecule.bonds
        self.angles = self.reactive_molecule.angles
        self.exclusions = self.reactive_molecule.exclusions
        self.settles = self.reactive_molecule.settles
        self.proper_dihedrals = self.reactive_molecule.proper_dihedrals
        self.improper_dihedrals = self.reactive_molecule.improper_dihedrals
        self.pairs = self.reactive_molecule.pairs
        self.position_restraints = self.reactive_molecule.position_restraints
        self.dihedral_restraints = self.reactive_molecule.dihedral_restraints
        self.radicals = self.reactive_molecule.radicals
        self.nrexcl = self.reactive_molecule.nrexcl

    def _extract_mergable_molecules(self) -> dict[str, int]:
        """Extract all molecules that are to be merged into one moleculetype.
        And replaces them with a single moleculetype with the name REACTIVE_MOLECULEYPE
        """
        reactive_molecules = {}
        new_molecules = []
        started_merging = False
        stopped_merging = False
        added_reactive_molecule = False
        non_reactive = None
        for m, n in self.molecules:
            if self._check_is_reactive_molecule(m):
                if stopped_merging:
                    m = f"""Attempting to merge a moleculetype {m} interspersed with non-reactive moleculetypes.
            Please make sure that all moleculetypes to be merged (all non-solvent molecules or ions by default)
            are listed consecutively in the [molecules] section of the topology and in the coordinates file (.gro).
            The latest non-reactive molecule was {non_reactive}.
            """
                    logger.error(m)
                    raise ValueError(m)
                started_merging = True
                reactive_molecules[m] = int(n)
                if not added_reactive_molecule:
                    new_molecules += [(REACTIVE_MOLECULEYPE, "1")]
                    added_reactive_molecule = True
            else:
                if started_merging:
                    stopped_merging = True
                    non_reactive = m
                new_molecules += [(m, n)]

        self.molecules = new_molecules
        logger.debug(
            "Merging the following molecules into the Reactive moleculetype and making their multiples explicit:"
        )
        for m, n in reactive_molecules.items():
            logger.debug(f"\t{m} {n}")
        return reactive_molecules

    def _merge_moleculetypes(
        self, radicals: Optional[str] = None, nrexcl: Optional[str] = None
    ):
        """
        Merge all moleculetypes within which reactions can happen into one moleculetype.
        This also makes multiples explicit.
        Reactive molecules have to be in a continuous block in the [molecules] section of the topology
        and the gro file, preferably the start.
        Atom ids (gromacs, 1-based) are updated accordingly and thus correspond directly
        to the atom numbers in the coordinates file (.gro) (ignoring gromacs overflow
        problems in the atomnr column, the correct internal atom nr is always "gro file line number" - 2).
        """
        molecules = self._extract_mergable_molecules()
        first_reactive_molecule = list(molecules.keys())[0]
        if nrexcl is None:
            nrexcl = self.moleculetypes[first_reactive_molecule].nrexcl
        reactive_atomics = {}
        atomnr_offset = 0
        resnr_offset = 0
        for name, n in molecules.items():
            add_atomics = self.moleculetypes[name].atomics
            n_atoms = len(add_atomics["atoms"])
            highest_resnr = int(add_atomics["atoms"][-1][RESNR_ID_FIELDS["atoms"][0]])
            for _ in range(n):
                for section_name, section in add_atomics.items():
                    section = deepcopy(section)
                    atomnr_fields = ATOM_ID_FIELDS.get(section_name, [])
                    resnr_fields = RESNR_ID_FIELDS.get(section_name, [])
                    for line in section:
                        if atomnr_fields is True:
                            # if atomnr_fields is True, then all fields are atomnr fields
                            # NOTE: This is why this is testing for actually being True, not just truthiness!
                            for field, _ in enumerate(line):
                                increment_field(line, field, atomnr_offset)
                            continue
                        for field in atomnr_fields:
                            increment_field(line, field, atomnr_offset)
                        for field in resnr_fields:
                            increment_field(line, field, resnr_offset)
                    if reactive_atomics.get(section_name) is None:
                        reactive_atomics[section_name] = []
                    reactive_atomics[section_name] += section

                atomnr_offset += n_atoms
                resnr_offset += highest_resnr
            # remove now merged moleculetype from topology
            # and the topology dict
            del self.moleculetypes[name]
            del self.top[f"moleculetype_{name}"]

        # add merged moleculetype to topology
        reactive_moleculetype = MoleculeType(
            MoleculeTypeHeader(name=REACTIVE_MOLECULEYPE, nrexcl=nrexcl),
            reactive_atomics,
            radicals,
            combine_same_periodicity=self.combine_same_periodicity_dihedrals,
            handle_duplicate_lines=self.handle_duplicate_topology_lines,
        )
        self.moleculetypes[REACTIVE_MOLECULEYPE] = reactive_moleculetype
        # update topology dict
        self._update_dict()

    def validate_bond(self, atm1: Atom, atm2: Atom) -> bool:
        """Validates bond consistency between both atoms and top
        Returns True if bond exists, False if not.
        Raises RuntimeError if bond is not consistent.
        """

        counter = 0
        if (atm1.nr, atm2.nr) in self.bonds.keys():
            counter += 1
        if atm1.nr in atm2.bound_to_nrs:
            counter += 1
        if atm2.nr in atm1.bound_to_nrs:
            counter += 1
        if counter == 3:
            return True
        elif counter == 0:
            return False
        raise RuntimeError(
            f"Bond between atom {atm1.nr} {atm1.type} and {atm2.nr} {atm2.type} "
            "(1-based) is ill-defined!\n"
            f"\tatm1.bound_to_nrs {atm1.bound_to_nrs}\n"
            f"\tatm2.bound_to_nrs {atm2.bound_to_nrs}\n"
        )

    def _expand_parameter_macros(self):
        """Replace #define shortcut names in bonded lines by the numbers they stand for.

        Some force fields (for example Amber99SB-ILDN) write torsions as named
        shortcuts: `180 165 167 170 9 torsion_LEU_C_CA_CB_CG_mult1`, with
        `_mult2` and `_mult3` on the same four atoms. The periodicity is hidden
        inside the name, so KIMMDY stored all three in one slot and kept only
        the last. This swaps each name for its numbers (phase, strength,
        periodicity), exactly as GROMACS does when it reads the file, so every
        line is kept. Changes the topology dict in place.
        """
        # shortcuts that stand for numbers (defines without a value are flags)
        defines = {k: v for k, v in (self.top.get("define") or {}).items() if v}
        if not defines:
            return
        # in each section, the numbers start after the function-type column
        funct_index = {"bonds": 2, "pairs": 2, "angles": 3, "dihedrals": 4}

        # replace the shortcut names in every bonded line of every molecule
        expanded = 0
        for name, section in self.top.items():
            if not name.startswith("moleculetype"):
                continue
            subsections = section.get("subsections") or {}
            for sec, index in funct_index.items():
                block = subsections.get(sec)
                if not block:
                    continue
                for key in ("content", "else_content"):
                    lines = block.get(key) or []
                    for i, line in enumerate(lines):
                        if not any(tok in defines for tok in line[index + 1 :]):
                            continue
                        new = list(line[: index + 1])
                        for tok in line[index + 1 :]:
                            new.extend(defines.get(tok, [tok]))
                        lines[i] = new
                        expanded += 1
        if expanded:
            logger.info(f"Expanded {expanded} parameter macro line(s) into explicit numbers")

    def _parse_molecules(self):
        moleculetypes = [k for k in self.top.keys() if k.startswith("moleculetype")]
        for moleculetype in moleculetypes:
            header = get_moleculetype_header(self.top, moleculetype)
            if header is None:
                logger.warning(f"moleculetype {moleculetype} has no header. Skipping.")
                continue
            atomics = get_moleculetype_atomics(self.top, moleculetype)
            if atomics is None:
                logger.warning(
                    f"moleculetype {moleculetype} has no atoms, bonds, angles etc. Skipping."
                )
                continue
            name = header.name
            # Only search for radicals in reactive MoleculeType during merging
            self.moleculetypes[name] = MoleculeType(
                header,
                atomics,
                radicals="",
                combine_same_periodicity=self.combine_same_periodicity_dihedrals,
                handle_duplicate_lines=self.handle_duplicate_topology_lines,
            )

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Topology):
            return False
        return self.to_dict() == other.to_dict()

    def _update_dict_from_moleculetype_atomics(
        self, moleculetype: MoleculeType, create: bool = False
    ):
        """Set content of the atomics (atoms/bonds/angles etc.) in the a topology dict from of a moleculetype.

        Resolves any `#ifdef` statements by check in the top['define'] dict
        and chooses the 'content' or 'else_content' depending on the result.

        If create is True, a new section is created if it does not exist.
        """
        name = moleculetype.name
        atomics = moleculetype.atomics
        top = self.top
        moleculetype_name = f"moleculetype_{name}"
        section = top.get(moleculetype_name)
        if section is None:
            if create:
                logger.debug(
                    f"topology does not contain {moleculetype_name}. Creating new section."
                )
                section = empty_section()
                section["content"] = [[name, moleculetype.nrexcl]]
                section["subsections"] = {k: empty_section() for k in atomics.keys()}
                top[moleculetype_name] = section
            else:
                logger.warning(
                    f"topology does not contain {moleculetype_name} and create=False. Not creating new section."
                )
                return None

        subsections = section["subsections"]
        section["content"] = [[name, moleculetype.nrexcl]]
        for k in list(subsections.keys()):
            v = subsections[k]
            condition = v.get("condition")
            if condition is not None:
                condition_type = condition.get("type")
                condition_value = condition.get("value")
                if condition_type == "ifdef":
                    if condition_value in top["define"].keys():
                        v["content"] = atomics[k]
                    else:
                        v["else_content"] = atomics[k]
                elif condition_type == "ifndef":
                    if condition_value not in top["define"].keys():
                        v["content"] = atomics[k]
                    else:
                        v["else_content"] = atomics[k]
                else:
                    raise NotImplementedError(
                        f"condition type {condition_type} is not supported"
                    )
            else:
                v["content"] = atomics[k]

    def _update_dict(self):
        """Update the topology dictionary with the current state of the topology."""
        for name, moleculetype in self.moleculetypes.items():
            moleculetype._update_atomics_dict()
            if name == REACTIVE_MOLECULEYPE:
                self._update_dict_from_moleculetype_atomics(moleculetype, create=True)
            else:
                self._update_dict_from_moleculetype_atomics(moleculetype)

        set_top_section(
            self.top,
            "molecules",
            [[name, n] for name, n in self.molecules],
        )

        set_top_section(
            self.top,
            "atomtypes",
            [attributes_to_list(x) for x in self.ff.atomtypes.values()],
        )

        set_top_section(
            self.top,
            "bondtypes",
            [attributes_to_list(x) for x in self.ff.bondtypes.values()],
        )
        if self.ff.nonbond_params != {}:
            set_top_section(
                self.top,
                "nonbond_params",
                [attributes_to_list(x) for x in self.ff.nonbond_params.values()],
            )
        set_top_section(
            self.top,
            "angletypes",
            [attributes_to_list(x) for x in self.ff.angletypes.values()],
        )

    def _parameterize_focus(self, focus_ids: set[str]) -> None:
        if self.parametrizer is None:
            raise RuntimeError("No Parametrizer was initialized in this topology!")
        logger.info(
            "Starting parametrization using %s with focus on %s",
            self.parametrizer.__class__.__name__,
            focus_ids,
        )
        self.parametrizer.parameterize_topology(self, focus_nrs=focus_ids)

    def update_parameters(self):
        if self.needs_parameterization:
            self._parameterize_focus(set(self.parameterization_focus_ids))
            self.needs_parameterization = False

    def prepare_product_parameters(self) -> ParameterizationReport:
        """Prepare topology B and retain only still-active radical focus IDs."""
        event_focus_ids = set(self.parameterization_focus_ids)
        if self.needs_parameterization:
            self._parameterize_focus(event_focus_ids)

        restored_templates = self.restore_pending_product_templates()
        active_radical_ids = set(self.radicals)
        if restored_templates and active_radical_ids:
            # Active radical regions own overlapping interactions.
            self._parameterize_focus(active_radical_ids)

        # TODO: Restore exact classical RTP bonded templates before compacting
        # focus IDs; RTP products currently retain their final Grappa terms.
        self.parameterization_focus_ids = active_radical_ids
        self.needs_parameterization = False
        return ParameterizationReport(
            event_focus_ids=event_focus_ids,
            active_radical_ids=active_radical_ids,
            restored_templates=restored_templates,
        )

    def restore_pending_product_templates(self) -> int:
        """Restore exact single-residue ITP products after Grappa parameterization."""
        restored = 0
        pending = self.pending_product_template_restorations
        self.pending_product_template_restorations = []
        for match in pending:
            if not match.template.supports_bonded_restoration:
                logger.info(
                    "ITP template %s:%s spans multiple residues; retaining "
                    "Grappa bonded parameters.",
                    match.template.moleculetype,
                    match.template.local_resnr,
                )
                continue
            self._restore_product_template_interactions(match)
            restored += 1
        return restored

    def _restore_product_template_interactions(
        self, match: ProductTemplateMatch
    ) -> None:
        template = match.template
        mapped_ids = set(match.local_to_global.values())

        bonds = [
            self._resolve_template_bond(match, line)
            for line in template.bond_lines
        ]
        pairs = [self._map_template_pair(match, line) for line in template.pair_lines]
        angles = [
            self._resolve_template_angle(match, line) for line in template.angle_lines
        ]
        proper_dihedrals, improper_dihedrals = self._resolve_template_dihedrals(
            match
        )

        for key in list(self.bonds):
            if set(key) <= mapped_ids:
                self.bonds.pop(key)
        self.bonds.update({(bond.ai, bond.aj): bond for bond in bonds})

        for key in list(self.pairs):
            if set(key) <= mapped_ids:
                self.pairs.pop(key)
        self.pairs.update({(pair.ai, pair.aj): pair for pair in pairs})

        for key in list(self.angles):
            if set(key) <= mapped_ids:
                self.angles.pop(key)
        self.angles.update(
            {(angle.ai, angle.aj, angle.ak): angle for angle in angles}
        )

        for key, interaction in list(self.proper_dihedrals.items()):
            if (
                set(key) <= mapped_ids
                and interaction.funct == FFFUNC["mult_proper_dihedral"]
            ):
                self.proper_dihedrals.pop(key)
        self.proper_dihedrals.update(proper_dihedrals)

        # Only replace periodic impropers KIMMDY can interpolate safely. CHARMM
        # harmonic function-2 impropers remain under Grappa ownership for now.
        for key in improper_dihedrals:
            self.improper_dihedrals.pop(key, None)
        self.improper_dihedrals.update(improper_dihedrals)

        logger.info(
            "Restored bonded parameters from exact ITP product template %s:%s "
            "(%d bonds, %d pairs, %d angles, %d proper and %d improper dihedrals).",
            template.moleculetype,
            template.local_resnr,
            len(bonds),
            len(pairs),
            len(angles),
            len(proper_dihedrals),
            len(improper_dihedrals),
        )

    def _mapped_template_line(
        self,
        match: ProductTemplateMatch,
        line: tuple[str, ...],
        atom_count: int,
    ) -> list[str]:
        mapped = list(line)
        try:
            mapped[:atom_count] = [
                match.local_to_global[atom_nr] for atom_nr in line[:atom_count]
            ]
        except KeyError as error:
            raise RuntimeError(
                "Exact ITP product-template interaction references an unmapped "
                f"atom {error.args[0]!r} in {match.template.source}."
            ) from error
        return mapped

    def _resolve_template_bond(
        self, match: ProductTemplateMatch, line: tuple[str, ...]
    ) -> Bond:
        mapped = self._mapped_template_line(match, line, 2)
        key = tuple(sorted(mapped[:2], key=int))
        bond = Bond.from_top_line([*key, *mapped[2:]])
        if bond.c0 is not None and bond.c1 is not None:
            return bond
        atom_types = [self.atoms[atom_nr].type for atom_nr in key]
        bond_type = match_atomic_item_to_atomic_type(atom_types, self.ff.bondtypes)
        if not isinstance(bond_type, BondType):
            raise ValueError(
                f"Could not resolve ITP bond parameters for {atom_types} in "
                f"{match.template.source}."
            )
        return Bond(*key, funct=bond_type.funct, c0=bond_type.c0, c1=bond_type.c1)

    def _map_template_pair(
        self, match: ProductTemplateMatch, line: tuple[str, ...]
    ) -> Pair:
        mapped = self._mapped_template_line(match, line, 2)
        key = tuple(sorted(mapped[:2], key=int))
        return Pair.from_top_line([*key, *mapped[2:]])

    def _resolve_template_angle(
        self, match: ProductTemplateMatch, line: tuple[str, ...]
    ) -> Angle:
        mapped = self._mapped_template_line(match, line, 3)
        if int(mapped[0]) > int(mapped[2]):
            mapped[0], mapped[2] = mapped[2], mapped[0]
        angle = Angle.from_top_line(mapped)
        if angle.c0 is not None and angle.c1 is not None:
            return angle
        key = (angle.ai, angle.aj, angle.ak)
        atom_types = [self.atoms[atom_nr].type for atom_nr in key]
        angle_type = match_atomic_item_to_atomic_type(atom_types, self.ff.angletypes)
        if not isinstance(angle_type, AngleType):
            raise ValueError(
                f"Could not resolve ITP angle parameters for {atom_types} in "
                f"{match.template.source}."
            )
        return Angle(
            *key,
            funct=angle_type.funct,
            c0=angle_type.c0,
            c1=angle_type.c1,
            c2=angle_type.c2,
            c3=angle_type.c3,
        )

    def _resolve_template_dihedrals(
        self, match: ProductTemplateMatch
    ) -> tuple[
        dict[tuple[str, str, str, str], MultipleDihedrals],
        dict[tuple[str, str, str, str], MultipleDihedrals],
    ]:
        proper: dict[tuple[str, str, str, str], MultipleDihedrals] = {}
        improper: dict[tuple[str, str, str, str], MultipleDihedrals] = {}
        unsupported_functions: set[str] = set()

        for line in match.template.dihedral_lines:
            mapped = self._mapped_template_line(match, line, 4)
            funct = mapped[4]
            if funct not in {
                FFFUNC["mult_proper_dihedral"],
                FFFUNC["mult_improper_dihedral"],
            }:
                unsupported_functions.add(funct)
                continue
            if (
                funct == FFFUNC["mult_proper_dihedral"]
                and int(mapped[1]) > int(mapped[2])
            ):
                mapped[:4] = reversed(mapped[:4])

            dihedral = Dihedral.from_top_line(mapped)
            key = (dihedral.ai, dihedral.aj, dihedral.ak, dihedral.al)
            target = proper if funct == FFFUNC["mult_proper_dihedral"] else improper
            types = (
                self.ff.proper_dihedraltypes
                if target is proper
                else self.ff.improper_dihedraltypes
            )

            resolved: list[Dihedral] = []
            if dihedral.c0 is not None and dihedral.c1 is not None:
                resolved.append(dihedral)
            else:
                atom_types = [self.atoms[atom_nr].type for atom_nr in key]
                for periodicity in map(str, range(1, 7)):
                    dihedral_type = match_atomic_item_to_atomic_type(
                        atom_types, types, periodicity
                    )
                    if isinstance(dihedral_type, DihedralType):
                        resolved.append(
                            Dihedral(
                                *key,
                                funct=funct,
                                c0=dihedral_type.c0,
                                c1=dihedral_type.c1,
                                periodicity=dihedral_type.periodicity,
                            )
                        )

            if not resolved:
                raise ValueError(
                    "Could not resolve ITP dihedral parameters for "
                    f"{[self.atoms[atom_nr].type for atom_nr in key]} in "
                    f"{match.template.source}."
                )
            multiple = target.setdefault(
                key, MultipleDihedrals(*key, funct=funct, dihedrals={})
            )
            for item in resolved:
                multiple.dihedrals[item.periodicity] = item

        if unsupported_functions:
            logger.warning(
                "ITP template %s contains unsupported dihedral functions %s; "
                "retaining Grappa parameters for those interactions.",
                match.template.source,
                sorted(unsupported_functions),
            )
        return proper, improper

    def update_partial_charges(self, recipe_steps: list[RecipeStep]) -> None:
        """Update the topology atom partial charges.

        Calls must reflect the topology state represented by their steps. A
        standalone Break can be processed after the recipe, preserving legacy
        slow-growth behavior. A Break associated with a later Bind must be
        processed before ``bind_bond()`` changes the shared atom's residue;
        that Bind is then processed afterward, when product-template
        information is available.

        Update rules follow a simple assignment scheme that works well with
        grappa. If fragments are created, their partial charges are kept
        integers. If previously broken bonds are formed again, the original
        partial charges are restored.
        If no product template matches a transferred H, its original charge
        is restored and balanced locally by the acceptor atom. With
        ``conserve_transferred_h_charge`` this local balancing is used for
        every transferred H, so templates only provide names and atom types.
        """

        # build updated list of atoms by residue for topology so that this does
        # not need to be repeated
        # Warning: resnr not unique for each residue, 'residues' can map to
        # more than one real residue
        residues = {}
        for atom in self.atoms.values():
            if residues.get(atom.resnr) is None:
                residues[atom.resnr] = []
            residues[atom.resnr].append(atom)

        for step in recipe_steps:
            if isinstance(step, Break):
                # make partial charges on either side of the break integer
                if self.atoms[step.atom_id_1].resnr == self.atoms[step.atom_id_2].resnr:
                    atom1 = self.atoms[step.atom_id_1]
                    atom2 = self.atoms[step.atom_id_2]
                    fragment1, fragment2 = get_residue_fragments(
                        self, residues[atom1.resnr], atom1, atom2
                    )
                    if len(fragment2) in (0, 1):
                        # intra-residue HAT case (or similar)
                        if atom1.type.upper().startswith("H"):
                            self.pending_transferred_h_charges[atom1.nr] = (
                                atom1.charge
                            )
                            atom2.charge = (
                                f"{float(atom2.charge) + float(atom1.charge):7.4f}"
                            )
                            atom1.charge = "0.0"
                        elif atom2.type.upper().startswith("H"):
                            self.pending_transferred_h_charges[atom2.nr] = (
                                atom2.charge
                            )
                            atom1.charge = (
                                f"{float(atom1.charge) + float(atom2.charge):7.4f}"
                            )
                            atom2.charge = "0.0"
                        else:
                            logger.warning(
                                "Residue is not fragmented and no break atom is hydrogen. Doing nothing!"
                            )
                            continue
                    else:
                        # homolysis case (or similar)

                        charge_fragment1 = [
                            float(self.atoms[nr].charge) for nr in fragment1
                        ]
                        charge_fragment2 = [
                            float(self.atoms[nr].charge) for nr in fragment2
                        ]
                        diff1 = sum(charge_fragment1) - round(sum(charge_fragment1))
                        diff2 = sum(charge_fragment2) - round(sum(charge_fragment2))
                        atom1.charge = f"{float(atom1.charge) - diff1:7.4f}"
                        atom2.charge = f"{float(atom2.charge) - diff2:7.4f}"
                else:
                    # no change in partial charges necessary for break between residues
                    continue
            elif isinstance(step, Bind):
                if self.atoms[step.atom_id_1].resnr == self.atoms[step.atom_id_2].resnr:
                    atom1 = self.atoms[step.atom_id_1]
                    atom2 = self.atoms[step.atom_id_2]
                    if self.conserve_transferred_h_charge and (
                        atom1.nr in self.pending_transferred_h_charges
                        or atom2.nr in self.pending_transferred_h_charges
                    ):
                        # charge-conserving transfer (conserve_transferred_h_charge):
                        # the H keeps its old charge and the acceptor gives up the
                        # same amount, instead of taking the template charges
                        self.pending_product_template_matches.pop(
                            tuple(sorted((step.atom_id_1, step.atom_id_2), key=int)),
                            None,
                        )
                        self._apply_transferred_h_charge_fallback(atom1, atom2)
                        continue
                    match_key = tuple(
                        sorted((step.atom_id_1, step.atom_id_2), key=int)
                    )
                    template_match = self.pending_product_template_matches.pop(
                        match_key,
                        None,
                    )
                    if template_match is not None:
                        for atom_nr, template_atom in (
                            template_match.atom_mapping.items()
                        ):
                            current_atom = self.atoms.get(atom_nr)
                            if current_atom is None:
                                raise RuntimeError(
                                    "Exact ITP product template match became "
                                    f"invalid: atom {atom_nr} no longer exists."
                                )
                            if current_atom.atom != template_atom.name:
                                raise RuntimeError(
                                    "Exact ITP product template match became "
                                    f"invalid for atom {atom_nr}: expected name "
                                    f"{template_atom.name!r}, found "
                                    f"{current_atom.atom!r}."
                                )
                            current_atom.type = template_atom.atom_type
                            current_atom.charge = template_atom.charge
                            current_atom.mass = (
                                template_atom.mass
                                if template_atom.mass is not None
                                else self.ff.atomtypes[
                                    template_atom.atom_type
                                ].mass
                            )
                        logger.info(
                            "Restored atom types, masses, and charges from exact "
                            "ITP product template %s:%s (%s).",
                            template_match.template.moleculetype,
                            template_match.template.local_resnr,
                            template_match.template.residue,
                        )
                        self.pending_product_template_restorations.append(
                            template_match
                        )
                        self._discard_pending_transferred_h_charge(atom1, atom2)
                        continue

                    # check whether bond exists in topology
                    residue = atom1.residue
                    bondtype = (atom1.atom, atom2.atom)
                    residuetype = self.ff.residuetypes.get(residue)
                    if residuetype is None:
                        m = f"Residue {residue} not found in residuetypes. Can't update partial charges!"
                        logger.warning(m)
                        self._apply_transferred_h_charge_fallback(atom1, atom2)
                        continue
                    residue_bond_spec = residuetype.bonds.get(
                        bondtype,
                        residuetype.bonds.get((bondtype[-1], bondtype[-2])),
                    )
                    if residue_bond_spec:
                        atom1.charge = residuetype.atoms[atom1.atom].charge
                        atom2.charge = residuetype.atoms[atom2.atom].charge
                        self._discard_pending_transferred_h_charge(atom1, atom2)
                    else:
                        logger.warning(
                            f"New bond defined in {step} but can't be found in residuetypes definition. Not changing the partial charges!"
                        )
                        self._apply_transferred_h_charge_fallback(atom1, atom2)
                else:
                    # no change in partial charges necessary
                    continue
            else:
                continue

    def _discard_pending_transferred_h_charge(
        self, atom1: Atom, atom2: Atom
    ) -> None:
        """Discard saved H charge after a product template restored charges."""
        self.pending_transferred_h_charges.pop(atom1.nr, None)
        self.pending_transferred_h_charges.pop(atom2.nr, None)

    def _apply_transferred_h_charge_fallback(
        self, atom1: Atom, atom2: Atom
    ) -> bool:
        """Conserve residue charges when no product charge template is available."""
        if atom1.nr in self.pending_transferred_h_charges:
            hydrogen, acceptor = atom1, atom2
        elif atom2.nr in self.pending_transferred_h_charges:
            hydrogen, acceptor = atom2, atom1
        else:
            return False

        hydrogen_charge = self.pending_transferred_h_charges.pop(hydrogen.nr)
        acceptor.charge = (
            f"{float(acceptor.charge) - float(hydrogen_charge):7.4f}"
        )
        hydrogen.charge = hydrogen_charge
        logger.warning(
            "No exact product charge template was available for transferred "
            "hydrogen %s. Restored its pre-transfer charge and subtracted that "
            "charge from acceptor atom %s.",
            hydrogen.nr,
            acceptor.nr,
        )
        return True

    def to_dict(self) -> TopologyDict:
        self.update_parameters()
        self._update_dict()
        return self.top

    def del_atom(
        self, atom_nr: Union[list[str], str], parameterize: bool = True
    ) -> dict[str, str]:
        """Deletes atom

        Deletes atom and all attached bonds. Reindexes the top and updates the
        parameters if requested. Also moves charges to first bound_nrs atom.

        Parameters
        ----------
        atom_nr
            1-based atom number as string to delete
        parameterize
            If true and bonds are removed triggers reparameterization,
            by default True

        Returns
        -------
        update_map
            Dict, mapping of old atom number strings to new ones.
        """
        if not isinstance(atom_nr, list):
            atom_nr = [atom_nr]
        for _atom_nr in atom_nr:
            atom = self.atoms[_atom_nr]
            logger.debug(
                f"Deleting Atom nr {atom.nr}, type {atom.type}, res {atom.residue}"
            )
            if len(atom.bound_to_nrs) > 0:
                # move charge to first neighbor
                self.atoms[atom.bound_to_nrs[0]].charge = (
                    f"{float(self.atoms[atom.bound_to_nrs[0]].charge) + float(atom.charge):7.4f}"
                )

                # break all bonds and delete all pairs, diheadrals with these bonds
                for bound_nr in copy(atom.bound_to_nrs):
                    self.break_bond((bound_nr, _atom_nr))

                for an in tuple(self.angles.keys()):
                    if _atom_nr in an:
                        self.angles.pop(an)

                for pd in tuple(self.proper_dihedrals.keys()):
                    if _atom_nr in pd:
                        self.proper_dihedrals.pop(pd)

                for id in tuple(self.improper_dihedrals.keys()):
                    if _atom_nr in id:
                        self.improper_dihedrals.pop(id)

            self.radicals.pop(_atom_nr)
            self.atoms.pop(_atom_nr)

        update_map_all = self.reindex_atomnrs()
        update_map = update_map_all[REACTIVE_MOLECULEYPE]

        if parameterize:
            self.update_parameters()
            # Overwriting in case of no parameterization wanted
            self.needs_parameterization = False
            self.parameterization_focus_ids = set()

        return update_map

    def reindex_atomnrs(self) -> dict[str, dict[str, str]]:
        """Reindex atom numbers in topology.

        Starts at index 1.
        This also updates the numbers for bonds, angles, dihedrals and pairs.
        Returns a dict of all moleculetypes to their update maps (old -> new).
        """
        update_map_all = {}
        for id, moleculetype in self.moleculetypes.items():
            update_map_all[id] = moleculetype.reindex_atomnrs()
        self._link_atomics()  # necessary because reindex_atomnrs creates a new dict
        return update_map_all

    def _regenerate_topology_from_bound_to(self):
        """Regenerate the topology from the bound_to lists of the atoms."""
        for moleculetype in self.moleculetypes.values():
            moleculetype._regenerate_topology_from_bound_to(self.ff)

    def __str__(self) -> str:
        molecules = "\n".join([name + ": " + n for name, n in self.molecules])
        reactive_moleculetype = self.moleculetypes.get(REACTIVE_MOLECULEYPE)
        text = (
            textwrap.dedent(f"Topology with the following molecules: \n{molecules}\n\n")
            + textwrap.dedent(f"With {reactive_moleculetype}")
            + textwrap.dedent(f"{self.ff}")
        )
        return text

    def __repr__(self) -> str:
        self._update_dict()
        return f"Topology({self.top})"

    def _repr_pretty_(self, p, _):
        """A __repr__ for ipython.

        This whill be used if just the name of the object is entered in the ipython shell
        or a jupyter notebook.

        p is an instance of
        [IPython.lib.pretty.PrettyPrinter](https://ipython.org/ipython-doc/3/api/generated/IPython.lib.pretty.html)
        """
        p.text(str(self))

    def get_neighbors(self, initial_atoms: set[str], order: int):
        explored_atoms = initial_atoms
        frontier = initial_atoms
        for _ in range(order):
            new_frontier = set()
            for atom in frontier:
                new_frontier.update(self.atoms[atom].bound_to_nrs)
            frontier = new_frontier - explored_atoms
            explored_atoms |= new_frontier
        return explored_atoms

    def break_bond(
        self, atompair_addresses: tuple[TopologyAtomAddress, TopologyAtomAddress]
    ):
        """Break bonds in topology homolytically.

        Removes bond, angles and dihedrals where atompair was involved.
        Modifies the topology dictionary in place.
        Atom pairs become radicals.

        Parameters
        ----------
        atompair_addresses
            Between which atoms to break the bond.
        """
        logger.debug(f"Breaking bond between {atompair_addresses}")
        reactive_moleculetype = self.reactive_molecule

        # tuple -> list -> sorted -> tuple still makes it a tuple of two strings
        # so pyright can chill.
        atompair_nrs: tuple[str, str] = tuple(sorted(atompair_addresses, key=int))  # type: ignore

        atompair = [
            reactive_moleculetype.atoms[atompair_nrs[0]],
            reactive_moleculetype.atoms[atompair_nrs[1]],
        ]

        if not self.validate_bond(*atompair):
            raise ValueError(
                "Trying to break non-existing bond!"
                f"\n\tatom 1, nr {atompair[0].nr}, type {atompair[0].type}, res {atompair[0].residue}"
                f"\n\tatom 2, nr {atompair[1].nr}, type {atompair[1].type}, res {atompair[1].residue}"
            )

        self.needs_parameterization = True
        self.parameterization_focus_ids.update(atompair_nrs)

        # mark atoms as radicals
        for atom in atompair:
            atom.is_radical = True
            reactive_moleculetype.radicals[atom.nr] = atom

        # bonds
        # remove bond
        removed_bond = reactive_moleculetype.bonds.pop(atompair_nrs)
        logging.debug(f"removed bond: {removed_bond}")

        # remove angles
        angle_keys = reactive_moleculetype._get_center_atom_angles(
            atompair_nrs[0]
        ) + reactive_moleculetype._get_center_atom_angles(atompair_nrs[1])
        for key in angle_keys:
            if all([x in key for x in atompair_nrs]):
                # angle contained a now deleted bond because
                # it had both atoms of the broken bond
                reactive_moleculetype.angles.pop(key, None)

        # remove proper dihedrals
        # and pairs
        dihedral_keys = reactive_moleculetype._get_atom_proper_dihedrals(
            atompair_nrs[0]
        ) + reactive_moleculetype._get_atom_proper_dihedrals(atompair_nrs[1])
        for key in dihedral_keys:
            # don't use periodicity in key for checking atompair_nrs
            if all([x in key for x in atompair_nrs]):
                # dihedral contained a now deleted bond because
                # it had both atoms of the broken bond
                reactive_moleculetype.proper_dihedrals.pop(key, None)
                pairkey: tuple[str, str] = tuple(sorted((key[0], key[3]), key=int))  # type: ignore
                reactive_moleculetype.pairs.pop(pairkey, None)

        # and improper dihedrals
        keys_remove = []
        for key in reactive_moleculetype.improper_dihedrals:
            if all([x in key for x in atompair_nrs]):
                keys_remove.append(key)
        for key in keys_remove:
            reactive_moleculetype.improper_dihedrals.pop(key, None)

        # warn if atom is part of restraint
        for atom in atompair:
            posres = atom.nr in self.position_restraints.keys()
            dihres = any([atom.nr in atms for atms in self.position_restraints.keys()])
            if posres or dihres:
                logging.warning(
                    "Deleting bond with atom involved in a restraint!"
                    "This may causes unwanted behaviour.\n"
                    f"\n\tatom 1, nr {atompair[0].nr}, type {atompair[0].type}, res {atompair[0].residue}"
                    f"\n\tatom 2, nr {atompair[1].nr}, type {atompair[1].type}, res {atompair[1].residue}"
                )

        # update bound_to
        atompair[0].bound_to_nrs.remove(atompair[1].nr)
        atompair[1].bound_to_nrs.remove(atompair[0].nr)

    def bind_bond(
        self, atompair_addresses: tuple[TopologyAtomAddress, TopologyAtomAddress]
    ):
        """Add a bond in topology.

        Modifies the topology dictionary in place.
        It keeps track of affected terms in the topology via a graph representation of the topology
        and applies the necessary changes to bonds, angles and dihedrals (proper and improper).
        Furthermore, it modifies to function types in the topology to account for radicals.

        Parameters
        ----------
        atompair_addresses
            Atoms to bind together.
        """
        logger.debug(f"Creating bond between {atompair_addresses}")

        reactive_moleculetype = self.reactive_molecule

        atompair_nrs: tuple[str, str] = tuple(sorted(atompair_addresses, key=int))  # type: ignore
        atompair = [
            reactive_moleculetype.atoms[atompair_nrs[0]],
            reactive_moleculetype.atoms[atompair_nrs[1]],
        ]

        # check whether they are bound already
        if self.validate_bond(*atompair):
            raise ValueError(
                "Trying to bind to atoms already bound!"
                f"\n\tatom 1, nr {atompair[0].nr}, type {atompair[0].type}, res {atompair[0].residue}"
                f"\n\tatom 2, nr {atompair[1].nr}, type {atompair[1].type}, res {atompair[1].residue}"
            )

        self.needs_parameterization = True
        self.parameterization_focus_ids.update(atompair_nrs)

        # de-radialize if re-combining two radicals
        if all(map(lambda x: x.is_radical, atompair)):
            atompair[0].is_radical = False
            atompair[1].is_radical = False
            for a in atompair:
                if reactive_moleculetype.radicals.pop(a.nr) is None:
                    logging.warning(
                        f"Atom {a} should have been a radical but was "
                        "not registed in moleculetype.radicals"
                    )

        # TODO: Do we still need this with grappa?
        # quickfix for jumping hydrogens
        # to make them adopt the correct residuetype and atomtype
        # when bound to a new heavy atom
        for i, atom in enumerate(atompair):
            if atom.type.upper().startswith("H"):
                other_i = abs(i - 1)
                other_atom = atompair[other_i]
                other_res = other_atom.residue
                atom.residue = other_res
                atom.resnr = other_atom.resnr

                template_match = None
                if "itp" in self.ff.product_template_sources:
                    template_match = (
                        self.ff.product_templates.match_hydrogenated_residue(
                            topology=self,
                            acceptor=other_atom,
                            transferred_h=atom,
                        )
                    )
                if template_match is not None:
                    template_h = template_match.template.atoms[
                        template_match.transferred_h_template_name
                    ]
                    atom.atom = template_h.name
                    atom.type = template_h.atom_type
                    atom.mass = (
                        template_h.mass
                        if template_h.mass is not None
                        else self.ff.atomtypes[template_h.atom_type].mass
                    )
                    self.pending_product_template_matches[
                        atompair_nrs
                    ] = template_match
                    logger.info(
                        "Matched hydrogenated product to ITP template "
                        "%s:%s (%s); assigned transferred H as %s with type %s.",
                        template_match.template.moleculetype,
                        template_match.template.local_resnr,
                        template_match.template.residue,
                        template_h.name,
                        template_h.atom_type,
                    )
                    continue

                if "rtp" not in self.ff.product_template_sources:
                    logger.warning(
                        "No exact ITP product template matched hydrogen binding "
                        "to %s:%s (atom %s). Keeping the legacy atom identity "
                        "and using Grappa bonded-parameter fallback.",
                        other_atom.residue,
                        other_atom.atom,
                        other_atom.nr,
                    )
                    continue

                # TODO: Extend exact RTP product matching to restore resolved
                # bonded interactions, including cross-residue boundary terms.
                aa = self.ff.residuetypes.get(other_res)
                if not aa:
                    logging.warning(
                        f"No residuetype found in ff for {other_res}, to which a H would bind"
                    )
                    continue

                other_residue = get_residue_by_bonding(
                    other_atom, reactive_moleculetype.atoms
                )
                residue_atomnames_current = [a.atom for a in other_residue.values()]
                name_set = False
                for key, bond in aa.bonds.items():
                    if other_atom.atom in key and any(
                        k.upper().startswith("H") for k in key
                    ):
                        if key[0] == other_atom.atom:
                            h_name = bond.atom2
                        else:
                            h_name = bond.atom1

                        atom.type = aa.atoms[h_name].type
                        if (
                            not h_name in residue_atomnames_current
                            or h_name == atom.atom
                        ):
                            atom.atom = h_name
                            name_set = True
                        else:
                            atom.atom = "HX"
                            name_set = True
                            continue
                        logger.debug(f"Hydrogen will be bound to {other_atom}.")
                        break
                else:
                    if name_set:
                        logger.warning(
                            "Found new atomtype but not atomname for HAT hydrogen!"
                        )
                    else:
                        logger.warning(
                            "Found neither new atomtype nor atomname for HAT hydrogen!"
                        )
                if not name_set:
                    atom.atom = "HX"
                    logger.debug(f"Named newly bonded hydrogen 'HX'")

        # update bound_to
        atompair[0].bound_to_nrs.append(atompair[1].nr)
        atompair[1].bound_to_nrs.append(atompair[0].nr)

        # bonds
        # add bond
        bond = Bond(atompair_nrs[0], atompair_nrs[1], FFFUNC["harmonic_bond"])
        reactive_moleculetype.bonds[atompair_nrs] = bond
        logging.info(f"added bond: {bond}")

        # make sure that there are no additional pairs
        # where bonds or angles are already between the atoms
        bond_pairkey: tuple[str, str] = tuple(
            str(x) for x in sorted([atompair[0].nr, atompair[1].nr], key=int)
        )  # type: ignore
        reactive_moleculetype.pairs.pop(bond_pairkey, None)

        # and remove improper dihedrals connected to the bonds,
        # as they will be re-added based on potentially changed residuetypes
        to_pop = []
        for key in reactive_moleculetype.improper_dihedrals.keys():
            if atompair_nrs[0] in key or atompair_nrs[1] in key:
                to_pop.append(key)

        for key in to_pop:
            reactive_moleculetype.improper_dihedrals.pop(key, None)

        # add angles
        angle_keys = reactive_moleculetype._get_atom_angles(
            atompair_nrs[0]
        ) + reactive_moleculetype._get_atom_angles(atompair_nrs[1])
        for key in angle_keys:
            if reactive_moleculetype.angles.get(key) is None:
                reactive_moleculetype.angles[key] = Angle(
                    key[0], key[1], key[2], FFFUNC["harmonic_angle"]
                )
            # make sure that there are no additional pairs
            # where bonds or angles are already between the atoms
            angle_pairkey: tuple[str, str] = tuple(
                str(x) for x in sorted([key[0], key[2]], key=int)
            )  # type: ignore
            reactive_moleculetype.pairs.pop(angle_pairkey, None)

        # add proper and improper dihedrals
        # add proper dihedrals and pairs
        # this assumes function type 9 for all proper dihedrals
        # and adds all possible dihedrals (and pairs)
        dihedral_keys = reactive_moleculetype._get_atom_proper_dihedrals(
            atompair_nrs[0]
        ) + reactive_moleculetype._get_atom_proper_dihedrals(atompair_nrs[1])
        for key in dihedral_keys:
            if reactive_moleculetype.proper_dihedrals.get(key) is None:
                reactive_moleculetype.proper_dihedrals[key] = MultipleDihedrals(
                    *key,
                    FFFUNC["mult_proper_dihedral"],
                    dihedrals={"": Dihedral(*key, FFFUNC["mult_proper_dihedral"])},
                )
            pairkey: tuple[str, str] = tuple(
                str(x) for x in sorted([key[0], key[3]], key=int)
            )  # type: ignore
            if reactive_moleculetype.pairs.get(pairkey) is None:
                # only add pair if the atoms are not otherwise bound (in a circle)
                if (
                    pairkey[0]
                    not in reactive_moleculetype.atoms[pairkey[1]].bound_to_nrs
                ):
                    reactive_moleculetype.pairs[pairkey] = Pair(
                        pairkey[0], pairkey[1], FFFUNC["pair"]
                    )

        # improper dihedral: re-add the template impropers of both atoms,
        # every listed term (one per periodicity)
        atom_impropers = reactive_moleculetype._get_atom_improper_dihedrals(
            atompair_nrs[0], self.ff
        )
        atom_impropers.update(
            reactive_moleculetype._get_atom_improper_dihedrals(atompair_nrs[1], self.ff)
        )
        for key, specs in atom_impropers.items():
            if reactive_moleculetype.improper_dihedrals.get(key) is None:
                reactive_moleculetype.improper_dihedrals[key] = MultipleDihedrals(
                    *key, FFFUNC["mult_improper_dihedral"], dihedrals={}
                )
            for value in specs:
                periodicity = "" if value.c2 is None else value.c2
                reactive_moleculetype.improper_dihedrals[key].dihedrals[periodicity] = (
                    Dihedral(
                        *key,
                        FFFUNC["mult_improper_dihedral"],
                        c0=value.c0,
                        c1=value.c1,
                        periodicity=periodicity,
                    )
                )

        # remove settles and explicit exclusions
        # those are used by solvent molecules
        # but if a solvent molecule get's bound to other parts
        # of the reactive molecule it is no longer a solvent molecule
        for ai in atompair_nrs:
            to_delete = []
            for exclusion_key in reactive_moleculetype.exclusions.keys():
                if ai in exclusion_key:
                    to_delete.append(exclusion_key)

            if len(to_delete) > 0:
                logger.debug(f"Removing exclusions {to_delete}")
            for key in to_delete:
                reactive_moleculetype.exclusions.pop(key)
            settles = reactive_moleculetype.settles.get(ai)
            if settles is not None:
                logger.debug(f"Removing settles {ai}")
                reactive_moleculetype.settles.pop(ai)
