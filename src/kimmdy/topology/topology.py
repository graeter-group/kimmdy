from pathlib import Path
from typing import Optional
from kimmdy.constants import ATOMTYPE_BONDORDER_FLAT
from kimmdy.parsing import TopologyDict
from kimmdy.topology.atomic import *
from kimmdy.topology.utils import (
    match_id_to_patch,
    attributes_to_list,
    match_atomic_item_to_atomic_type,
    set_protein_section,
    set_top_section,
)
from kimmdy.topology.ff import FF, FFPatches, Patch
from itertools import permutations, combinations
import textwrap
import logging
import re


PROTEIN_SECTION = "moleculetype_0"


class Topology:
    """Smart container for parsed topology data.

    A topology keeps track of connections and applies patches to parameters when bonds are broken or formed.
    Assumptions:
    - the topology of interest (the protein) is in section 'moleculetype_0'.

    Parameters
    ----------
    top :
        A dictionary containing the parsed topology data.
    ffpatch : Optional[Path]
        Path to a force field patch file. If None, no patching is applied.
    """

    def __init__(
        self,
        top: TopologyDict,
        ffpatch: Optional[Path] = None,
    ) -> None:
        if top == {}:
            raise NotImplementedError(
                "Generating an empty Topology from an empty TopologyDict is not implemented."
            )

        if not top.get(PROTEIN_SECTION) or not top[PROTEIN_SECTION].get("subsections"):
            raise ValueError(
                "The topology does not contain a protein section."
                "Please make sure the topology contains a section"
                "called [ moleculetype ]."
                "The first of which is assumed to be the protein of interest."
            )
        self.protein = top[PROTEIN_SECTION]["subsections"]
        self.top = top
        self.atoms: dict[str, Atom] = {}
        self.bonds: dict[tuple[str, str], Bond] = {}
        self.pairs: dict[tuple[str, str], Pair] = {}
        self.angles: dict[tuple[str, str, str], Angle] = {}
        self.proper_dihedrals: dict[tuple[str, str, str, str], Dihedral] = {}
        self.improper_dihedrals: dict[tuple[str, str, str, str], Dihedral] = {}
        self.position_restraints: dict[str, PositionRestraint] = {}
        self.dihedral_restraints: dict[
            tuple[str, str, str, str], DihedralRestraint
        ] = {}
        self.radicals: dict[str, Atom] = {}

        self.ff = FF(top)

        self.ffpatches = None
        if ffpatch:
            self.ffpatches = FFPatches(ffpatch)

        self._parse_atoms()
        self._parse_bonds()
        self._parse_pairs()
        self._parse_angles()
        self._parse_dihedrals()
        self._parse_restraints()
        self._initialize_graph()
        self._test_for_radicals()

    def _update_dict(self):
        set_protein_section(
            self.top, "atoms", [attributes_to_list(x) for x in self.atoms.values()]
        )

        set_protein_section(
            self.top, "bonds", [attributes_to_list(x) for x in self.bonds.values()]
        )

        set_protein_section(
            self.top, "pairs", [attributes_to_list(x) for x in self.pairs.values()]
        )

        set_protein_section(
            self.top, "angles", [attributes_to_list(x) for x in self.angles.values()]
        )

        set_protein_section(
            self.top,
            "dihedrals",
            [attributes_to_list(x) for x in self.proper_dihedrals.values()]
            + [attributes_to_list(x) for x in self.improper_dihedrals.values()],
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

        set_top_section(
            self.top,
            "angletypes",
            [attributes_to_list(x) for x in self.ff.angletypes.values()],
        )

        # activate again once this works
        # set_top_section(
        #     self.top,
        #     "dihedraltypes",
        #     [attributes_to_list(x) for x in self.ff.improper_dihedraltypes.values()]
        #     + [attributes_to_list(x) for x in self.ff.proper_dihedraltypes.values()],
        # )

    def to_dict(self) -> TopologyDict:
        self._update_dict()
        return self.top

    def __repr__(self) -> str:
        return textwrap.dedent(
            f"""\
        Topology with
        {len(self.atoms)} atoms,
        {len(self.bonds)} bonds,
        {len(self.angles)} angles,
        {len(self.pairs)} pairs,
        {len(self.proper_dihedrals)} proper dihedrals
        {len(self.improper_dihedrals)} improper dihedrals
        """
        )

    def reindex_atomnrs(self):
        """Reindex atom numbers in topology.

        Starts at index 1.
        This also updates the numbers for bonds, angles, dihedrals and pairs.
        """
        update_map = {
            atom_nr: str(i + 1) for i, atom_nr in enumerate(self.atoms.keys())
        }
        print(update_map)

        new_atoms = {}
        for atom in self.atoms.values():
            atom.nr = update_map[atom.nr]
            new_atoms[atom.nr] = atom
        self.atoms = new_atoms

        new_bonds = {}
        for bond in self.bonds.values():
            bond.ai = update_map[bond.ai]
            bond.aj = update_map[bond.aj]
            new_bonds[(update_map[bond.ai], update_map[bond.aj])] = bond
        self.bonds = new_bonds

        new_angles = {}
        for angle in self.angles.values():
            angle.ai = update_map[angle.ai]
            angle.aj = update_map[angle.aj]
            angle.ak = update_map[angle.ak]
            new_angles[
                (update_map[angle.ai], update_map[angle.aj], update_map[angle.ak])
            ] = angle
        self.angles = new_angles

        new_dihedrals = {}
        for dihedral in self.proper_dihedrals.values():
            dihedral.ai = update_map[dihedral.ai]
            dihedral.aj = update_map[dihedral.aj]
            dihedral.ak = update_map[dihedral.ak]
            dihedral.al = update_map[dihedral.al]
            new_dihedrals[
                (
                    update_map[dihedral.ai],
                    update_map[dihedral.aj],
                    update_map[dihedral.ak],
                    update_map[dihedral.al],
                )
            ] = dihedral
        self.dihedrals = new_dihedrals

        new_pairs = {}
        for pair in self.pairs.values():
            pair.ai = update_map[pair.ai]
            pair.aj = update_map[pair.aj]
            new_pairs[(update_map[pair.ai], update_map[pair.aj])] = pair
        self.pairs = new_pairs

    def __str__(self) -> str:
        return str(self.atoms)

    def _parse_atoms(self):
        """Parse atoms from topology dictionary."""
        ls = self.protein["atoms"]["content"]
        for l in ls:
            atom = Atom.from_top_line(l)
            self.atoms[atom.nr] = atom

    def _parse_bonds(self):
        """Parse bond from topology dictionary."""
        ls = self.protein["bonds"]["content"]
        for l in ls:
            bond = Bond.from_top_line(l)
            self.bonds[(bond.ai, bond.aj)] = bond

    def _parse_pairs(self):
        """Parse pairs from topology dictionary."""
        ls = self.protein["pairs"]["content"]
        for l in ls:
            pair = Pair.from_top_line(l)
            self.pairs[(pair.ai, pair.aj)] = pair

    def _parse_angles(self):
        """Parse angles from topology dictionary."""
        ls = self.protein["angles"]["content"]
        for l in ls:
            angle = Angle.from_top_line(l)
            self.angles[(angle.ai, angle.aj, angle.ak)] = angle

    def _parse_dihedrals(self):
        """Parse improper and proper dihedrals from topology dictionary."""
        ls = self.protein["dihedrals"]["content"]
        for l in ls:
            dihedral = Dihedral.from_top_line(l)
            if dihedral.funct == "4":
                self.improper_dihedrals[
                    (dihedral.ai, dihedral.aj, dihedral.ak, dihedral.al)
                ] = dihedral
            else:
                self.proper_dihedrals[
                    (dihedral.ai, dihedral.aj, dihedral.ak, dihedral.al)
                ] = dihedral

    def _parse_restraints(self):
        """Parse restraints from topology dictionary."""
        ls = self.protein.get("position_restraints")
        if ls is None:
            return
        condition = None
        for l in ls.get("content"):
            if l[0] == "#ifdef":
                condition = l[1]
                continue
            elif l[0] == "#endif":
                condition = None
                continue
            restraint = PositionRestraint.from_top_line(l, condition=condition)
            self.position_restraints[restraint.ai] = restraint
        ls = self.protein.get("dihedral_restraints")
        if ls is None:
            return
        for l in ls.get("content"):
            restraint = DihedralRestraint.from_top_line(l)
            self.dihedral_restraints[
                (restraint.ai, restraint.aj, restraint.ak, restraint.al)
            ] = restraint

    def _initialize_graph(self):
        """Add a list of atom nrs bound to an atom to each atom."""
        for bond in self.bonds.values():
            i = bond.ai
            j = bond.aj
            self.atoms[i].bound_to_nrs.append(j)
            self.atoms[j].bound_to_nrs.append(i)

    def _test_for_radicals(self):
        """Iterate over all atoms and designate them as radicals if they have fewer bounds than their natural bond order"""
        for atom in self.atoms.values():
            bo = ATOMTYPE_BONDORDER_FLAT.get(atom.type)
            if bo and bo > len(atom.bound_to_nrs):
                atom.is_radical = True
                self.radicals[atom.nr] = atom
            else:
                atom.is_radical = False

        return None

    def _apply_param_patch(
        self,
        atomic_item: Atomic,
        atomic_id: list[str],
        patch: Patch,
        types: AtomicTypes,
    ):
        """Apply a patch to an atomic item (Atom, Bond, Angle etc.).

        Initial values are taken from the topology or the force field (supplied via `types`).
        """
        item_type = match_atomic_item_to_atomic_type(atomic_id, types)
        for param, param_patch in patch.params.items():
            initial = atomic_item.__dict__.get(param)
            if initial is None:
                # get initial value from the FF
                initial = item_type.__dict__.get(param)
            if initial is None:
                logging.warning(
                    f"Can't patch parameter because no initial parameter was found in the topology or the FF for: {atomic_id}, {atomic_item}."
                )
                continue

            try:
                initial = float(initial)
            except ValueError as _:
                logging.warning(
                    f"Malformed patchfile. Some parameter patches couldn't be converted to a number: {param} with value {initial} of {atomic_item}."
                )
                continue
            except TypeError as _:
                continue

            result = param_patch.apply(initial)
            atomic_item.__dict__[param] = result

    def _revert_param_patch(
        self,
        atomic_item: Atomic,
        atomic_id: list[str],
        types: AtomicTypes,
    ):
        """Revert a patch to an atomic item (Atom, Bond, Angle etc.).

        Values reset to `None` if found in the forcefield (supplied via `types`).
        """
        item_type = match_atomic_item_to_atomic_type(atomic_id, types)
        if item_type is None:
            logging.warning(
                f"Won't revert patch because no initial parameter was found in the FF for: {atomic_id}, {atomic_item}."
            )
            return

        print(f"item_type: {item_type}")
        for param in atomic_item.__dict__.keys():
            if re.match(r"^c\d", param):
                atomic_item.__dict__[param] = None

    def patch_parameters(self, focus_nr: list[str]):
        if self.ffpatches is None:
            return
        focus = [self.atoms[nr] for nr in focus_nr]
        logging.info(f"Applying parameter patches around these atoms: {focus}.")

        # atoms
        for atom in focus:
            patch = match_id_to_patch([atom.radical_type()], self.ffpatches.atompatches)
            if patch is None:
                continue
            self._apply_param_patch(atom, [atom.type], patch, self.ff.atomtypes)

        # bonds
        for atom in focus:
            for bond_key in self._get_atom_bonds(atom.nr):
                bond = self.bonds.get(bond_key)
                if (
                    bond is None
                    or self.ffpatches is None
                    or self.ffpatches.bondpatches is None
                ):
                    continue
                id = [self.atoms[i].radical_type() for i in [bond.ai, bond.aj]]
                patch = match_id_to_patch(id, self.ffpatches.bondpatches)
                if patch is None:
                    continue
                id_base = [self.atoms[i].type for i in [bond.ai, bond.aj]]
                self._apply_param_patch(bond, id_base, patch, self.ff.bondtypes)

        angle_keys = self._get_atom_angles(focus_nr[0]) + self._get_atom_angles(
            focus_nr[1]
        )
        for key in angle_keys:
            angle = self.angles.get(key)
            if (
                angle is None
                or self.ffpatches is None
                or self.ffpatches.anglepatches is None
            ):
                continue
            id = [self.atoms[i].radical_type() for i in key]
            patch = match_id_to_patch(id, self.ffpatches.anglepatches)
            if patch is None:
                continue
            id_base = [self.atoms[i].radical_type() for i in key]
            self._apply_param_patch(angle, id_base, patch, self.ff.angletypes)

        # proper dihedrals and pairs
        dihedral_keys = self._get_atom_proper_dihedrals(
            focus_nr[0]
        ) + self._get_atom_proper_dihedrals(focus_nr[1])
        for key in dihedral_keys:
            dihedral = self.proper_dihedrals.get(key)
            if (
                dihedral is None
                or self.ffpatches is None
                or self.ffpatches.anglepatches is None
            ):
                continue
            id = [self.atoms[i].radical_type() for i in key]
            patch = match_id_to_patch(id, self.ffpatches.dihedralpatches)
            if patch is None:
                continue
            id_base = [self.atoms[i].radical_type() for i in key]
            self._apply_param_patch(
                dihedral, id_base, patch, self.ff.proper_dihedraltypes
            )

        # TODO: improper dihedrals

    def break_bond(self, atompair_nrs: tuple[str, str]):
        """Break bonds in topology.

        removes bond, angles and dihedrals where atompair was involved.
        Modifies the topology dictionary in place.
        """
        atompair_nrs = tuple(sorted(atompair_nrs, key=int))
        atompair = [self.atoms[atompair_nrs[0]], self.atoms[atompair_nrs[1]]]

        # mark atoms as radicals
        for atom in atompair:
            atom.is_radical = True
            self.radicals[atom.nr] = atom

        # bonds
        # remove bond
        removed_bond = self.bonds.pop(atompair_nrs, None)
        logging.info(f"removed bond: {removed_bond}")

        # remove angles
        angle_keys = self._get_atom_angles(atompair_nrs[0]) + self._get_atom_angles(
            atompair_nrs[1]
        )
        for key in angle_keys:
            if all([x in key for x in atompair_nrs]):
                # angle contained a now deleted bond because
                # it had both atoms of the broken bond
                self.angles.pop(key, None)

        # remove proper dihedrals
        # and pairs
        dihedral_keys = self._get_atom_proper_dihedrals(
            atompair_nrs[0]
        ) + self._get_atom_proper_dihedrals(atompair_nrs[1])
        for key in dihedral_keys:
            if all([x in key for x in atompair_nrs]):
                # dihedral contained a now deleted bond because
                # it had both atoms of the broken bond
                self.proper_dihedrals.pop(key, None)
                pairkey = tuple(sorted((key[0], key[3]), key=int))
                self.pairs.pop(pairkey, None)

        # and improper dihedrals
        dihedral_k_v = self._get_atom_improper_dihedrals(
            atompair_nrs[0]
        ) + self._get_atom_improper_dihedrals(atompair_nrs[1])
        for key, _ in dihedral_k_v:
            if all([x in key for x in atompair_nrs]):
                self.improper_dihedrals.pop(key, None)

        # update bound_to
        try:
            atompair[0].bound_to_nrs.remove(atompair[1].nr)
            atompair[1].bound_to_nrs.remove(atompair[0].nr)
        except ValueError as _:
            m = f"tried to remove bond between already disconnected atoms: {atompair}."
            logging.warning(m)

    def bind_bond(self, atompair_nrs: tuple[str, str]):
        """Add a bond in topology.

        Modifies the topology dictionary in place.
        It keeps track of affected terms in the topology via a graph representation of the topology
        and applies the necessary changes to bonds, angles and dihedrals (proper and improper).
        Furthermore, it modifies to function types in the topology to account for radicals.

        Parameters
        ----------
        atompair:
            A tuple of integers with the atoms indices (id, starting at 1)
            with `from`, the atom being moved and
            `to`, the atom to which the `from` atom will be bound
        """

        atompair_nrs = tuple(sorted(atompair_nrs, key=int))
        atompair = [self.atoms[atompair_nrs[0]], self.atoms[atompair_nrs[1]]]

        # de-radialize if re-combining two radicals
        if all(map(lambda x: x.is_radical, atompair)):
            atompair[0].is_radical = False
            atompair[1].is_radical = False
            for a in atompair:
                self.radicals.pop(a.nr)

        # quickfix for jumping hydrogens
        # to make them adopt the correct residuetype and atomtype
        # when bound to a new heavy atom
        for i, atom in enumerate(atompair):
            if atom.type.startswith("H"):
                other_i = abs(i - 1)
                other_atom = atompair[other_i]
                other_res = other_atom.residue
                atom.residue = other_res
                aa = self.ff.residuetypes.get(other_res)
                if not aa:
                    logging.warning(
                        f"No AA found for {other_res}, to which a H would jump"
                    )
                    continue

                for key, bond in aa.bonds.items():
                    if other_atom.atom in key and any(k.startswith("H") for k in key):
                        if key[0] == other_atom.atom:
                            h_name = bond.atom2
                        else:
                            h_name = bond.atom1

                        atom.type = aa.atoms[h_name].type
                        logging.info(f"HAT hydrogen will be bound to {other_atom}.")
                        break
                else:
                    logging.warn(f"Found no new atomtype for HAT hydrogen!")

        # update bound_to
        atompair[0].bound_to_nrs.append(atompair[1].nr)
        atompair[1].bound_to_nrs.append(atompair[0].nr)

        # bonds
        # add bond
        bond = Bond(atompair_nrs[0], atompair_nrs[1], "1")
        self.bonds[atompair_nrs] = bond
        logging.info(f"added bond: {bond}")

        # add angles
        angle_keys = self._get_atom_angles(atompair_nrs[0]) + self._get_atom_angles(
            atompair_nrs[1]
        )
        for key in angle_keys:
            if self.angles.get(key) is None:
                self.angles[key] = Angle(key[0], key[1], key[2], "1")

        # add proper and improper dihedrals
        # add proper dihedrals and pairs
        dihedral_keys = self._get_atom_proper_dihedrals(
            atompair_nrs[0]
        ) + self._get_atom_proper_dihedrals(atompair_nrs[1])
        for key in dihedral_keys:
            if self.proper_dihedrals.get(key) is None:
                self.proper_dihedrals[key] = Dihedral(
                    key[0], key[1], key[2], key[3], "9"
                )
            pairkey = tuple(str(x) for x in sorted([key[0], key[3]], key=int))
            if self.pairs.get(pairkey) is None:
                self.pairs[pairkey] = Pair(pairkey[0], pairkey[1], "1")

        # improper dihedral
        dihedral_k_v = self._get_atom_improper_dihedrals(
            atompair_nrs[0]
        ) + self._get_atom_improper_dihedrals(atompair_nrs[1])
        for key, value in dihedral_k_v:
            if self.improper_dihedrals.get(key) is None:
                # TODO: fix this
                c2 = None
                if value.q0 is not None:
                    c2 = "1"
                self.improper_dihedrals[key] = Dihedral(
                    key[0], key[1], key[2], key[3], "4", value.q0, value.cq, c2
                )

    def move_hydrogen(self, from_to: tuple[str, str]):
        """Move a singly bound atom to a new location.

        This is typically H for Hydrogen Atom Transfer (HAT).
        """
        f, t = list(map(str, from_to))
        assert (
            self.atoms[f].type[0] == "H"
        ), f"move_hydrogen called for non-hydrogen! type: {self.atoms[f].type}"
        heavy = self.atoms[f].bound_to_nrs.pop()
        if heavy is None:
            logging.error(f"Atom {f} is not bound to anything.")
            return
        self.break_bond((f, heavy))
        self.bind_bond((f, t))

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
            "get_atom_pairs is not implementes. Get the pairs as the endpoints of dihedrals instead."
        )

    def _get_atom_angles(self, atom_nr: str) -> list[tuple[str, str, str]]:
        """
        each atom has a list of atoms it is bound to
        get a list of angles that one atom is involved in
        based in these lists.
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
        each atom has a list of atoms it is bound to.
        get a list of dihedrals that one atom is involved in
        based in these lists.
        """
        return self._get_center_atom_dihedrals(
            atom_nr
        ) + self._get_margin_atom_dihedrals(atom_nr)

    def _get_center_atom_dihedrals(
        self, atom_nr: str
    ) -> list[tuple[str, str, str, str]]:
        dihedrals = []
        aj = atom_nr
        partners = self.atoms[aj].bound_to_nrs
        for ai in partners:
            for ak in partners:
                if ai == ak:
                    continue
                for al in self.atoms[ak].bound_to_nrs:
                    if al == ak or aj == al:
                        continue
                    if int(aj) < int(ak):
                        dihedrals.append((ai, aj, ak, al))
                    else:
                        dihedrals.append((al, ak, aj, ai))
        return dihedrals

    def _get_margin_atom_dihedrals(
        self, atom_nr: str
    ) -> list[tuple[str, str, str, str]]:
        dihedrals = []
        ai = atom_nr
        for aj in self.atoms[ai].bound_to_nrs:
            for ak in self.atoms[aj].bound_to_nrs:
                if ai == ak:
                    continue
                for al in self.atoms[ak].bound_to_nrs:
                    if al == ak or aj == al:
                        continue
                    if int(aj) < int(ak):
                        dihedrals.append((ai, aj, ak, al))
                    else:
                        dihedrals.append((al, ak, aj, ai))
        return dihedrals

    def _get_atom_improper_dihedrals(
        self, atom_nr: str
    ) -> list[tuple[tuple[str, str, str, str], ResidueImproperSpec]]:
        # TODO: cleanup and make more efficient
        # which improper dihedrals are used is defined for each residue
        # in aminoacids.rtp
        # get improper diheldrals from FF based on residue
        # TODO: handle impropers defined for the residue that
        # belongs to an adjacent atom, not just the the specied one
        atom = self.atoms[atom_nr]
        residue = self.ff.residuetypes.get(atom.residue)
        if residue is None:
            return []

        # <https://manual.gromacs.org/current/reference-manual/functions/bonded-interactions.html#improper-dihedrals>
        # atom in a line, like a regular dihedral:
        dihedrals = []
        dihedral_candidate_keys = self._get_margin_atom_dihedrals(
            atom_nr
        ) + self._get_center_atom_dihedrals(atom_nr)

        # atom in the center of a star/tetrahedron:
        ai = atom_nr
        partners = self.atoms[ai].bound_to_nrs
        if len(partners) >= 3:
            combs = combinations(partners, 3)
            for comb in combs:
                aj, ak, al = comb
                dihedral_candidate_keys.append((ai, aj, ak, al))

        # atom in corner of a star/tetrahedron:
        for a in self.atoms[atom_nr].bound_to_nrs:
            partners = self.atoms[a].bound_to_nrs
            if len(partners) >= 3:
                combs = permutations(partners + [a], 4)
                for comb in combs:
                    ai, aj, ak, al = comb
                    dihedral_candidate_keys.append((ai, aj, ak, al))

        # residues on aminoacids.rtp specify a dihedral to the next or previous
        # AA with -C and +N as the atomname
        for candidate in dihedral_candidate_keys:
            candidate_key = [self.atoms[atom_nr].atom for atom_nr in candidate]
            for i, nr in enumerate(candidate):
                if self.atoms[nr].resnr != atom.resnr:
                    if candidate_key[i] == "C":
                        candidate_key[i] = "-C"
                    elif candidate_key[i] == "N":
                        candidate_key[i] = "+N"

            candidate_key = tuple(candidate_key)
            dihedral = residue.improper_dihedrals.get(candidate_key)
            if dihedral:
                dihedrals.append((candidate, dihedral))

        return dihedrals

    def _regenerate_topology_from_bound_to(self):
        # clear all bonds, angles, dihedrals
        self.bonds = {}
        self.angles = {}
        self.dihedrals = {}
        self.proper_dihedrals = {}
        self.improper_dihedrals = {}

        # bonds
        keys = []
        for atom in self.atoms.values():
            keys = self._get_atom_bonds(atom.nr)
            for key in keys:
                self.bonds[key] = Bond(key[0], key[1], "1")

        # angles
        for atom in self.atoms.values():
            keys = self._get_atom_angles(atom.nr)
            for key in keys:
                self.angles[key] = Angle(key[0], key[1], key[2], "1")

        # dihedrals and pass
        for atom in self.atoms.values():
            keys = self._get_atom_proper_dihedrals(atom.nr)
            for key in keys:
                self.proper_dihedrals[key] = Dihedral(
                    key[0], key[1], key[2], key[3], "9"
                )
                pairkey = tuple(str(x) for x in sorted([key[0], key[3]], key=int))
                if self.pairs.get(pairkey) is None:
                    self.pairs[pairkey] = Pair(pairkey[0], pairkey[1], "1")

        for atom in self.atoms.values():
            impropers = self._get_atom_improper_dihedrals(atom.nr)
            for key, improper in impropers:
                self.improper_dihedrals[key] = Dihedral(
                    improper.atom1,
                    improper.atom2,
                    improper.atom3,
                    improper.atom4,
                    "4",
                    improper.cq,
                )
