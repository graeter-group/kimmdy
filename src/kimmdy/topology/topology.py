from pathlib import Path
from typing import Optional
from kimmdy.parsing import TopologyDict
from kimmdy.topology.atomic import *
from kimmdy.topology.utils import (
    match_id_to_patch,
    attributes_to_list,
    match_atomic_item_to_atomic_type,
)
from kimmdy.topology.ff import FF, FFPatches, Patch
from itertools import permutations, combinations
import textwrap
import logging


class Topology:
    """Smart container for parsed topology data.

    A topology keeps track of connections and applies patches to parameters when bonds are broken or formed.
    """

    def __init__(
        self,
        top: TopologyDict,
        ffdir: Optional[Path] = None,
        ffpatch: Optional[Path] = None,
    ) -> None:
        self.top = top
        self.forcefield_directory = ffdir
        self.atoms: dict[str, Atom] = {}
        self.bonds: dict[tuple[str, str], Bond] = {}
        self.pairs: dict[tuple[str, str], Pair] = {}
        self.angles: dict[tuple[str, str, str], Angle] = {}
        self.proper_dihedrals: dict[tuple[str, str, str, str], Dihedral] = {}
        self.improper_dihedrals: dict[tuple[str, str, str, str], Dihedral] = {}
        self.radicals: dict[str, Atom] = {}

        if ffdir:
            self.ff = FF(ffdir)
        self.ffpatches = None
        if ffpatch:
            self.ffpatches = FFPatches(ffpatch)

        # generate empty Topology if empty TopologyDict
        if self.top == {}:
            return

        self._parse_atoms()
        self._parse_bonds()
        self._parse_pairs()
        self._parse_angles()
        self._parse_dihedrals()
        self._initialize_graph()

    def _update_dict(self):

        self.top["atoms"] = [attributes_to_list(x) for x in self.atoms.values()]
        self.top["bonds"] = [attributes_to_list(x) for x in self.bonds.values()]
        self.top["pairs"] = [attributes_to_list(x) for x in self.pairs.values()]
        self.top["angles"] = [attributes_to_list(x) for x in self.angles.values()]
        self.top["dihedrals"] = [
            attributes_to_list(x) for x in self.proper_dihedrals.values()
        ]

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

    def __str__(self) -> str:
        return str(self.atoms)

    def _parse_atoms(self):
        ls = self.top["atoms"]
        for l in ls:
            atom = Atom.from_top_line(l)
            self.atoms[atom.nr] = atom

    def _parse_bonds(self):
        ls = self.top["bonds"]
        for l in ls:
            bond = Bond.from_top_line(l)
            self.bonds[(bond.ai, bond.aj)] = bond

    def _parse_pairs(self):
        ls = self.top["pairs"]
        for l in ls:
            pair = Pair.from_top_line(l)
            self.pairs[(pair.ai, pair.aj)] = pair

    def _parse_angles(self):
        ls = self.top["angles"]
        for l in ls:
            angle = Angle.from_top_line(l)
            self.angles[(angle.ai, angle.aj, angle.ak)] = angle

    def _parse_dihedrals(self):
        ls = self.top["dihedrals"]
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

    def _initialize_graph(self):
        for bond in self.bonds.values():
            i = bond.ai
            j = bond.aj
            self.atoms[i].bound_to_nrs.append(j)
            self.atoms[j].bound_to_nrs.append(i)

    def _apply_param_patch(
        self,
        atomic_item: Atomic,
        atomic_id: list[str],
        patch: Patch,
        types: AtomicTypes,
    ):
        # atomic_id_base = [item.removesuffix("_R") for item in atomic_id]
        item_type = match_atomic_item_to_atomic_type(atomic_id, types)
        for param, param_patch in patch.params.items():
            initial = atomic_item.__dict__.get(param)
            if initial is None:
                # get initial value from the FF
                # TODO:
                # what about matching up differently named parameters?
                # in atomtypes and the topology?
                initial = item_type.__dict__.get(param)

            try:
                initial = float(initial)
            except ValueError as _:
                logging.warning(
                    "Malformed patchfile. Some parameter pachtes couldn't be converted to a number:"
                )
                continue
            except TypeError as _:
                logging.warning(
                    "Can't patch parameter because no initial parameter was found in the topology or the FF: "
                )
                continue

            result = param_patch.apply(initial)
            atomic_item.__dict__[param] = result

    def break_bond(self, atompair_nrs: tuple[str, str]):
        """Break bonds in topology.

        removes bond, angles and dihedrals where atompair was involved.
        Modifies the topology dictionary in place.
        It modifies the function types and parameters in the topology to account for radicals.
        """
        atompair_nrs = tuple(sorted(atompair_nrs, key=int))
        atompair = [self.atoms[atompair_nrs[0]], self.atoms[atompair_nrs[1]]]

        # mark atoms as radicals
        for atom in atompair:
            atom.is_radical = True
            self.radicals[atom.nr] = atom

            # patch parameters
            if self.ffpatches is None or self.ffpatches.atompatches is None:
                continue
            atom_id = atom.type + "_R"

            patch = match_id_to_patch([atom_id], self.ffpatches.atompatches)
            if patch is None:
                continue
            self._apply_param_patch(atom, [atom.type], patch, self.ff.atomtypes)

        # bonds
        # remove bond
        removed_bond = self.bonds.pop(atompair_nrs, None)
        logging.info(f"removed bond: {removed_bond}")

        # get other possibly affected bonds
        for atom in atompair:
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

        # remove angles
        angle_keys = self._get_atom_angles(atompair_nrs[0]) + self._get_atom_angles(
            atompair_nrs[1]
        )
        for key in angle_keys:
            if all([x in key for x in atompair_nrs]):
                # angle contained a now deleted bond because
                # it had both atoms of the broken bond
                self.angles.pop(key, None)
            else:
                # angle only contains one of the affecte atoms
                # angle is not removed but might need to be patched
                # patch parameters
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

        # remove proper dihedrals
        # and pairs
        dihedral_keys = self._get_atom_proper_dihedrals(
            atompair_nrs[0]
        ) + self._get_atom_proper_dihedrals(atompair_nrs[1])
        for key in dihedral_keys:
            if all([x in key for x in atompair_nrs]):
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

        Move an atom (typically H for Hydrogen Atom Transfer) to a new location.
        Modifies the topology dictionary in place.
        It keeps track of affected terms in the topology via a graph representation of the topology
        and applies the necessary changes to bonds, angles and dihedrals (proper and improper).
        Furthermore, it modifies to function types in the topology to account for radicals.

        Parameters
        ----------
        atompair: a tuple of integers with the atoms indices
            `from`, the atom being moved and
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
                # TODO: fix this after the demonstration
                c2 = None
                if value.q0 is not None:
                    c2 = "1"
                self.improper_dihedrals[key] = Dihedral(
                    key[0], key[1], key[2], key[3], "4", value.q0, value.cq, c2
                )

        # if there are no changed parameters for radicals, exit here
        if self.ffpatches is None:
            self._update_dict()
            return

    def _get_atom_bonds(self, atom_nr: str) -> list[tuple[str, str]]:
        ai = atom_nr
        bonds = []
        for aj in self.atoms[ai].bound_to_nrs:
            # TODO: maybe sort here and filter later instead
            if int(ai) < int(aj):
                bonds.append((ai, aj))
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
    ) -> list[tuple[tuple[str, str, str, str], ResidueImroperSpec]]:
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


def generate_topology_from_bound_to(
    atoms: list[Atom], ffdir: Path, ffpatch: Path
) -> Topology:
    top = Topology({}, ffdir, ffpatch)
    for atom in atoms:
        top.atoms[atom.nr] = atom

    # bonds
    keys = []
    for atom in top.atoms.values():
        keys = top._get_atom_bonds(atom.nr)
        for key in keys:
            top.bonds[key] = Bond(key[0], key[1], "1")

    # angles
    for atom in top.atoms.values():
        keys = top._get_atom_angles(atom.nr)
        for key in keys:
            top.angles[key] = Angle(key[0], key[1], key[2], "1")

    # dihedrals and pass
    for atom in top.atoms.values():
        keys = top._get_atom_proper_dihedrals(atom.nr)
        for key in keys:
            top.proper_dihedrals[key] = Dihedral(key[0], key[1], key[2], key[3], "9")
            pairkey = tuple(str(x) for x in sorted([key[0], key[3]], key=int))
            if top.pairs.get(pairkey) is None:
                top.pairs[pairkey] = Pair(pairkey[0], pairkey[1], "1")

    for atom in top.atoms.values():
        impropers = top._get_atom_improper_dihedrals(atom.nr)
        for key, improper in impropers:
            top.improper_dihedrals[key] = Dihedral(
                improper.atom1,
                improper.atom2,
                improper.atom3,
                improper.atom4,
                "4",
                improper.cq,
            )

    return top
