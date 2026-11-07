import logging
import textwrap
from copy import copy, deepcopy
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Collection, Optional, Union

from kimmdy.constants import (
    ATOM_ID_FIELDS,
    ATOMTYPE_BONDORDER_FLAT,
    FFFUNC,
    REACTIVE_MOLECULEYPE,
    RESNR_ID_FIELDS,
)
from kimmdy.parsing import TopologyDict, empty_section, read_top, write_top
from kimmdy.plugins import BasicParameterizer, Parameterizer
from kimmdy.recipe import Bind, RecipeStep
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
from kimmdy.topology.charges import (
    DEFAULT_TEMPLATE_TOLERANCE,
    assign_recipe_charges,
    check_charge_conservation,
)
from kimmdy.topology.ff import FF
from kimmdy.topology.templates import assign_recipe_names, is_hydrogen
from kimmdy.topology.utils import (
    attributes_to_list,
    get_moleculetype_atomics,
    get_moleculetype_header,
    get_residue_by_bonding,
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
        repeated_lines: str = "merge",
        ff: Optional[FF] = None,
    ) -> None:
        # what to do with repeated topology lines: merge, error or last
        # (changer.topology.repeated_topology_lines)
        self.repeated_lines = repeated_lines
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
            self.find_radicals(ff)

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
            # a repeated line stops with an error (see _reject_duplicate_line)
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
            # a repeated line stops with an error (see _reject_duplicate_line)
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
            # a repeated line stops with an error (see _reject_duplicate_line)
            if key in self.angles:
                self._reject_duplicate_line("angles", key)
            self.angles[key] = angle

    def _parse_dihedrals(self):
        """Read the [ dihedrals ] lines of this molecule into proper and improper dihedrals.

        Each dihedral is stored under its four atoms, with one term per
        periodicity. If a line repeats atoms and periodicity of an earlier
        line, the two terms are merged into one equivalent term where that is
        exact (see _may_combine and combine_periodic_terms); any other repeat
        stops with an error (see _reject_duplicate_line), because GROMACS
        would use both lines but only one can be stored.
        """
        ls = self.atomics.get("dihedrals")
        if ls is None:
            return
        merged = 0
        for l in ls:
            dihedral = Dihedral.from_top_line(l)
            key = (dihedral.ai, dihedral.aj, dihedral.ak, dihedral.al)
            # periodic impropers (function 4) and everything else (propers) are kept apart
            is_improper = dihedral.funct == FFFUNC["mult_improper_dihedral"]
            kind = "improper" if is_improper else "proper"
            store = self.improper_dihedrals if is_improper else self.proper_dihedrals
            if key not in store:
                store[key] = MultipleDihedrals(*key, dihedral.funct, dihedrals={})
            terms = store[key].dihedrals

            existing = terms.get(dihedral.periodicity)
            if existing is not None:
                if self._may_combine(existing, dihedral, is_improper):
                    # same atoms and periodicity listed twice: GROMACS adds both,
                    # so store one equivalent term instead of only the later one
                    old_terms = [(existing.c0, existing.c1), (dihedral.c0, dihedral.c1)]
                    dihedral.c0, dihedral.c1 = combine_periodic_terms(old_terms)
                    # written at the end of the line, so the file shows why one line is missing
                    dihedral.comment = "; merged from 2 lines with the same atoms and periodicity"
                    inputs = " + ".join(f"(phase {p}, k {k})" for p, k in old_terms)
                    logger.debug(
                        f"Molecule {self.name}: merged repeated {kind} dihedral {' '.join(key)} "
                        f"(funct {dihedral.funct}, periodicity {dihedral.periodicity}): "
                        f"{inputs} -> phase {dihedral.c0}, k {dihedral.c1} "
                        "(same forces, energy shifted by a constant)"
                    )
                    merged += 1
                else:
                    self._reject_duplicate_line(f"dihedrals ({kind})", key, dihedral.periodicity)
            terms[dihedral.periodicity] = dihedral

        # one summary warning per molecule; the single merges are logged at debug level
        if merged:
            logger.warning(
                f"Molecule {self.name}: merged {merged} repeated dihedral term(s) "
                "with the same atoms and periodicity (same forces, energy shifted by a constant)"
            )

    def _may_combine(self, existing: Dihedral, new: Dihedral, is_improper: bool) -> bool:
        """Tell whether two dihedral lines with the same atoms and periodicity can be merged exactly.

        Only in merge mode, only when both lines have numbers and the same
        function, and for propers only for the periodic functions 1 and 9.
        """
        if self.repeated_lines != "merge":
            return False
        if None in (existing.c0, existing.c1, new.c0, new.c1):
            return False  # without numbers there is nothing to add up
        if existing.funct != new.funct:
            return False
        return is_improper or new.funct in ("1", "9")

    def _reject_duplicate_line(self, section: str, key: tuple, periodicity: str = ""):
        """Stop with an error when a topology line repeats an earlier one that cannot be merged.

        KIMMDY keeps one term per atom set (and per periodicity for
        dihedrals), while GROMACS uses every line, so silently keeping only
        the last line would change the simulation. The error names the
        molecule, section and atoms. In last mode it does nothing, so the
        later line replaces the earlier as older KIMMDY did.
        """
        if self.repeated_lines == "last":
            return
        where = f" with periodicity {periodicity}" if periodicity else ""
        m = (
            f"Molecule {self.name}: [ {section} ] lists atoms {' '.join(key)}{where} "
            "more than once. KIMMDY can keep only one such term, but GROMACS would "
            "apply every line, so the topology cannot be represented without "
            "changing it. Merge or remove the repeated line, or set "
            "changer.topology.repeated_topology_lines to last to keep only the last one "
            "(this can drop interactions)."
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

    def find_radicals(self, ff: Optional[FF] = None):
        """Find atoms that are radicals and update self.radicals.

        Iterate over all atoms and designate them as radicals if they have
        fewer bonds than the smallest natural bond order of their atom type.
        Atom types can have several natural bond orders
        (see ATOMTYPE_BONDORDER_FLAT), e.g. N3 with 3 or 4 bonds. An atom
        with fewer bonds than the largest, e.g. an N3 with 3 bonds (a neutral
        amine or an NH3+ that lost a hydrogen), is a radical if the rtp entry
        of its residue makes it one (TemplateLibrary.radical_by_template),
        given ff. Without ff, or without an entry that fits, it is not.
        """
        logger.debug(
            "Trying to infer radical status based on AMBER(!!) atomtype bond "
            f"order in molecule {self.name}."
        )
        self.radicals = {}
        error_atoms = set()
        unexpected_bonds = set()
        # atoms whose number of bonds fits several states, decided by templates
        ambiguous = []
        for atom in self.atoms.values():
            bond_orders = ATOMTYPE_BONDORDER_FLAT.get(atom.type)
            if bond_orders is None:
                error_atoms.add(atom.type)
                atom.is_radical = False
                continue
            n_bonds = len(atom.bound_to_nrs)
            if n_bonds < min(bond_orders):
                self.radicals[atom.nr] = atom
                atom.is_radical = True
            else:
                atom.is_radical = False
                if n_bonds not in bond_orders:
                    unexpected_bonds.add((atom.type, n_bonds))
                elif 0 not in bond_orders and n_bonds < max(bond_orders):
                    ambiguous.append(atom)
        if ff is not None:
            for atom in ambiguous:
                if ff.templates.radical_by_template(self.atoms, atom):
                    self.radicals[atom.nr] = atom
                    atom.is_radical = True
        if len(error_atoms) > 0:
            logger.warning(
                "Some atomtypes not found in AMBER atomtypes. Cannot infer radical status.\n"
                "Missing atom types:\n"
                f"{error_atoms}"
            )
        if len(unexpected_bonds) > 0:
            logger.warning(
                f"Some atoms in molecule {self.name} have a number of bonds "
                "that does not match the bond orders of their atom type. "
                "They are not treated as radicals.\n"
                "(atom type, number of bonds):\n"
                f"{sorted(unexpected_bonds)}"
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

        Looks up the atom's residue template (by residue name, else by atom
        names, e.g. HIE for a HIS), translates the template's atom
        names (including -C, the previous residue's C, and +N, the next
        residue's N) into atom numbers, and returns every improper term that
        contains the atom: {four atom numbers: list of terms, one per
        periodicity}.
        <https://manual.gromacs.org/current/reference-manual/functions/bonded-interactions.html#improper-dihedrals>
        """

        def residuetype_of(residue: dict[str, Atom]) -> Optional[ResidueType]:
            """The rtp entry of a residue, by its name and atom names (TemplateLibrary.rtp_template)."""
            template = ff.templates.rtp_template(residue)
            return None if template is None else ff.residuetypes.get(template.name)

        atom = self.atoms[atom_nr]
        residue = get_residue_by_bonding(atom=atom, atoms=self.atoms)
        residuetype = residuetype_of(residue)
        if residuetype is None:
            return {}

        # find the previous and next residue (needed for -C and +N)
        previous_atom, next_atom = self._get_neighbor_residue_atoms(residue)
        if previous_atom is not None:
            previous_residue = get_residue_by_bonding(
                atom=previous_atom, atoms=self.atoms
            )
            previous_residuetype = residuetype_of(previous_residue)
        else:
            previous_residue = {}
            previous_residuetype = None
        if next_atom is not None:
            next_residue = get_residue_by_bonding(atom=next_atom, atoms=self.atoms)
            next_residuetype = residuetype_of(next_residue)
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

    def update_impropers(self, atom_nrs: Collection[str], ff: FF) -> None:
        """Set the improper dihedrals of these atoms again from their residue templates.

        Removes every improper dihedral that contains one of the atoms and
        adds the template impropers that contain one of them, every listed
        term (one per periodicity), by the current atom names.
        """
        atom_nrs = set(atom_nrs)
        if not atom_nrs:
            return
        for key in [key for key in self.improper_dihedrals if atom_nrs & set(key)]:
            self.improper_dihedrals.pop(key)

        impropers: dict[ImproperDihedralId, list[ResidueImproperSpec]] = {}
        for atom_nr in sorted(atom_nrs, key=int):
            impropers.update(self._get_atom_improper_dihedrals(atom_nr, ff))
        for key, specs in impropers.items():
            if self.improper_dihedrals.get(key) is None:
                self.improper_dihedrals[key] = MultipleDihedrals(
                    *key, FFFUNC["mult_improper_dihedral"], dihedrals={}
                )
            for value in specs:
                periodicity = "" if value.c2 is None else value.c2
                self.improper_dihedrals[key].dihedrals[periodicity] = Dihedral(
                    *key,
                    FFFUNC["mult_improper_dihedral"],
                    c0=value.c0,
                    c1=value.c1,
                    periodicity=periodicity,
                )

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

        # impropers from the residue templates, with all their periodic terms
        # TODO: terminal residues are rebuilt from the plain residue template.
        # pdb2gmx applies its terminus rules (.n.tdb/.c.tdb) and, for example,
        # writes no a99SB-disp psi correction (N CA C +N) for the N-terminal
        # residue, so a rebuild there adds a psi term pdb2gmx does not have.
        for atom in self.atoms.values():
            impropers = self._get_atom_improper_dihedrals(atom.nr, ff)
            for key, specs in impropers.items():
                dihedrals = {}
                for improper in specs:
                    periodicity = "" if improper.c2 is None else improper.c2
                    dihedrals[periodicity] = Dihedral(
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
        repeated_topology_lines: str = "merge",
    ):
        """Topology options (changer.topology in kimmdy.yml):

        repeated_topology_lines
            What to do when a topology or residue template repeats the same
            atoms (and periodicity): "merge" adds up repeated dihedral terms
            with numbers into one equivalent term (forces unchanged, energy
            shifted by a constant) and warns, and stops on any other repeat;
            "error" stops on any repeat; "last" keeps only the last line, as
            older KIMMDY did. #define names are replaced first in every mode.
        """
        if repeated_topology_lines not in ("merge", "error", "last"):
            raise ValueError(
                "repeated_topology_lines must be merge, error or last; "
                f"got {repeated_topology_lines!r}"
            )
        self.repeated_topology_lines = repeated_topology_lines
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
            repeated_topology_lines=repeated_topology_lines,
        )
        # #define names first, so periodicities are known whatever the mode
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
            repeated_lines=self.repeated_topology_lines,
            ff=self.ff,
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

        Names are replaced in one pass: a #define whose value is another
        #define name (a nested alias) is not supported and stays unresolved.
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
                repeated_lines=self.repeated_topology_lines,
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

    def update_bonded(self) -> ParameterizationReport:
        """Parameterize the bonded interactions around the focus atoms, if the topology changed.

        Returns the atoms in focus and the radicals. Charges and atom types
        must be final before, because the parameterizer may use them
        (update_nonbonded).
        """
        focus_ids = set(self.parameterization_focus_ids)
        if self.needs_parameterization:
            self._parameterize_focus(focus_ids)
            self.needs_parameterization = False
        return ParameterizationReport(
            event_focus_ids=focus_ids,
            active_radical_ids=set(self.radicals),
        )

    def update_parameters(self) -> ParameterizationReport:
        """Deprecated name of update_bonded."""
        return self.update_bonded()

    def update_names_LJ(self, top_a: "Topology", recipe_steps: list[RecipeStep]) -> None:
        """Set names, atom types (Lennard-Jones) and masses of the residues recipe steps touched.

        This topology is the one after the steps (B), top_a the one before
        them (A). The residue templates decide (assign_recipe_names); atoms
        that get another name get their improper dihedrals again, and so do
        the atoms of the Bind steps: bind_bond built theirs while a moved
        hydrogen still had its old name.
        """
        logger.info("Updating atom names and types")
        renamed = assign_recipe_names(top_a, self, recipe_steps)
        if renamed:
            bound = {
                nr
                for step in recipe_steps
                if isinstance(step, Bind)
                for nr in (step.atom_id_1, step.atom_id_2)
            }
            self.reactive_molecule.update_impropers(renamed | bound, self.ff)

    def update_charges(
        self,
        top_a: "Topology",
        recipe_steps: list[RecipeStep],
        template_tolerance: float = DEFAULT_TEMPLATE_TOLERANCE,
        conservation_tolerance: float = DEFAULT_TEMPLATE_TOLERANCE,
        conservation_check: str = "warn",
    ) -> None:
        """Set the charges of the residues recipe steps touched and check the net charge.

        This topology is the one after the steps (B), top_a the one before
        them (A). See assign_recipe_charges and check_charge_conservation.
        """
        logger.info("Updating partial charges")
        assign_recipe_charges(top_a, self, recipe_steps, tolerance=template_tolerance)
        check_charge_conservation(
            top_a,
            self,
            tolerance=conservation_tolerance,
            action=conservation_check,
        )

    def update_nonbonded(
        self,
        top_a: "Topology",
        recipe_steps: list[RecipeStep],
        template_tolerance: float = DEFAULT_TEMPLATE_TOLERANCE,
        conservation_tolerance: float = DEFAULT_TEMPLATE_TOLERANCE,
        conservation_check: str = "warn",
    ) -> None:
        """Set names, atom types, masses and charges after recipe steps.

        Names and types come first, so the charges match the residue
        templates by the final names.
        """
        self.update_names_LJ(top_a, recipe_steps)
        self.update_charges(
            top_a,
            recipe_steps,
            template_tolerance=template_tolerance,
            conservation_tolerance=conservation_tolerance,
            conservation_check=conservation_check,
        )

    def update_all(
        self,
        top_a: Optional["Topology"],
        recipe_steps: list[RecipeStep],
        template_tolerance: float = DEFAULT_TEMPLATE_TOLERANCE,
        conservation_tolerance: float = DEFAULT_TEMPLATE_TOLERANCE,
        conservation_check: str = "warn",
    ) -> ParameterizationReport:
        """Update the topology after recipe steps: nonbonded, then bonded.

        This topology is the one after the steps (B), top_a the one before
        them (A), which stays unchanged. If top_a is None the steps did not
        change the topology, and only the bonded parameters are updated. The
        parameterizer may use charges and types, so they are final before it
        runs. Afterwards only the radicals stay in the parameterization focus,
        for the next event. Returns what update_bonded returns.
        """
        if top_a is not None:
            self.update_nonbonded(
                top_a,
                recipe_steps,
                template_tolerance=template_tolerance,
                conservation_tolerance=conservation_tolerance,
                conservation_check=conservation_check,
            )
        logger.info("Updating bonded parameters")
        report = self.update_bonded()
        self.parameterization_focus_ids = set(self.radicals)
        return report

    def to_dict(self) -> TopologyDict:
        self.update_bonded()
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
            self.update_bonded()
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

        # a hydrogen joins the residue of its new partner; update_names_LJ
        # gives it its name and type from the residue templates
        for i, atom in enumerate(atompair):
            if is_hydrogen(atom.type):
                partner = atompair[1 - i]
                atom.residue = partner.residue
                atom.resnr = partner.resnr

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

        # improper dihedrals of both atoms from their residue templates
        reactive_moleculetype.update_impropers(atompair_nrs, self.ff)

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
