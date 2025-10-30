import logging
from copy import copy
from pathlib import Path
from typing import Callable, Optional, Union

from gmx_top4py.parsing import TopologyDict
from gmx_top4py.parameterizing import Parameterizer, BasicParameterizer
from gmx_top4py.constants import FFFUNC
from gmx_top4py.topology.atomic import (
    Angle,
    Atom,
    Bond,
    Dihedral,
    MultipleDihedrals,
    Pair,
)
from gmx_top4py.topology.utils import is_not_solvent_or_ion, get_residue_by_bonding
from gmx_top4py.topology.topology import (
    MoleculeType,
    Topology as BasicTopology,
)

from kimmdy.constants import REACTIVE_MOLECULEYPE
from kimmdy.recipe import Bind, Break, RecipeStep
from kimmdy.topology.ff import FF
from kimmdy.topology.utils import get_residue_fragments

from kimmdy.utils import TopologyAtomAddress

logger = logging.getLogger("kimmdy.topology")


class Topology(BasicTopology):
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
        [](`gmx_top4py.parsing.read_top`) or deprecated by [](`kimmdy.parsing.read_top`)
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

    def __init__(
        self,
        top: TopologyDict,
        parametrizer: Parameterizer = BasicParameterizer(),
        is_reactive_predicate_f: Callable[[str], bool] = is_not_solvent_or_ion,
        radicals: Optional[str] = None,
        residuetypes_path: Optional[Path] = None,
        reactive_nrexcl: Optional[str] = None,
    ):
        self.selected_moleculetype = REACTIVE_MOLECULEYPE  # In kimmdy: selected_moleculetype = REACTIVE_MOLECULEYPE = "reactive".
        super().__init__(
            top,
            parametrizer=parametrizer,
            is_selected_moleculetype_f=is_reactive_predicate_f,
            radicals=radicals,
            residuetypes_path=residuetypes_path,
            nrexcl=reactive_nrexcl,
        )

    ################## Wrapper for backward compatibility ##################
    @property
    def _check_is_reactive_molecule(self) -> Callable[[str], bool]:
        return self._check_is_selected_molecule

    @_check_is_reactive_molecule.setter  # Is the setter necessary? Not sure but doesn't hurt to have.
    def _check_is_reactive_molecule(self, f: Callable[[str], bool]) -> None:
        self._check_is_selected_molecule = f

    @property
    def reactive_molecule(self) -> MoleculeType:
        return self.selected_molecule

    @reactive_molecule.setter  # Is the setter necessary? Not sure but doesn't hurt to have.
    def reactive_molecule(self, molecule: MoleculeType) -> None:
        self.selected_molecule = molecule

    ##########################################################################

    def _parse_ff(self, top: TopologyDict, residuetypes_path: Optional[Path] = None):
        """Parse force field data from topology dict into a container."""
        self.ff = FF(top, residuetypes_path)

    def _link_atomics(self):
        """Link atoms, bonds etc. properties of self (=top) to the reactive moleculeype.

        Call this any time a large change is made to the reactive moleculetype
        that re-assigns a property and not just modifies parts of it in place.
        """
        super()._link_atomics()

    def _extract_mergable_molecules(self) -> dict[str, int]:
        """Extract all molecules that are to be merged into one moleculetype.
        And replaces them with a single moleculetype with the name REACTIVE_MOLECULEYPE
        """
        reactive_molecules = {}
        new_molecules = []
        started_merging = False
        stopped_merging = False
        added_reactive_molecule = False
        for m, n in self.molecules:
            if self._check_is_reactive_molecule(m):
                if stopped_merging:
                    m = f"""Attempting to merge a moleculetype {m} interspersed with non-merging moleculetypes.
            Please make sure that all moleculetypes to be merged (all non-solvent molecules or ions by default)
            are listed consecutively in the [molecules] section of the topology and in the coordinates file (.gro)
            """
                    logger.error(m)
                    raise ValueError(m)
                started_merging = True
                reactive_molecules[m] = int(n)
                if not added_reactive_molecule:
                    new_molecules += [(self.selected_moleculetype, "1")]
                    added_reactive_molecule = True
            else:
                if started_merging:
                    stopped_merging = True
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
        super()._merge_moleculetypes(radicals=radicals, nrexcl=nrexcl)

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

    def update_partial_charges(self, recipe_steps: list[RecipeStep]) -> None:
        """Update the topology atom partial charges.

        This function must be called after the recipe_steps are applied.
        Changes are based on the recipe_steps. Update rules follow a simple
        assignment scheme that works well with grappa. If fragments are created,
        their partial charges are kept integers. If previously broken bonds are
        formed again, the original partial charges are restored.
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
                            atom2.charge = (
                                f"{float(atom2.charge) + float(atom1.charge):7.4f}"
                            )
                            atom1.charge = "0.0"
                        elif atom2.type.upper().startswith("H"):
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
                    # check whether bond exists in topology
                    residue = atom1.residue
                    bondtype = (atom1.atom, atom2.atom)
                    residuetype = self.ff.residuetypes.get(residue)
                    if residuetype is None:
                        m = f"Residue {residue} not found in residuetypes. Can't update partial charges!"
                        logger.warning(m)
                        continue
                    residue_bond_spec = residuetype.bonds.get(
                        bondtype,
                        residuetype.bonds.get((bondtype[-1], bondtype[-2])),
                    )
                    if residue_bond_spec:
                        atom1.charge = residuetype.atoms[atom1.atom].charge
                        atom2.charge = residuetype.atoms[atom2.atom].charge
                    else:
                        logger.warning(
                            f"New bond defined in {step} but can't be found in residuetypes definition. Not changing the partial charges!"
                        )
                else:
                    # no change in partial charges necessary
                    continue
            else:
                continue

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
        update_map = update_map_all[self.selected_moleculetype]

        if parameterize:
            self.update_parameters()
            # Overwriting in case of no parameterization wanted
            self.needs_parameterization = False
            self.parameterization_focus_ids = set()

        return update_map

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

        # improper dihedral
        atom_impropers = reactive_moleculetype._get_atom_improper_dihedrals(
            atompair_nrs[0], self.ff
        )
        atom_impropers.update(
            reactive_moleculetype._get_atom_improper_dihedrals(atompair_nrs[1], self.ff)
        )
        for key, value in atom_impropers.items():
            if value.c2 is None:
                value.c2 = ""
            if reactive_moleculetype.improper_dihedrals.get(key) is None:
                reactive_moleculetype.improper_dihedrals[key] = MultipleDihedrals(
                    *key, FFFUNC["mult_improper_dihedral"], dihedrals={}
                )
            reactive_moleculetype.improper_dihedrals[key].dihedrals[value.c2] = (
                Dihedral(
                    *key,
                    FFFUNC["mult_improper_dihedral"],
                    c0=value.c0,
                    c1=value.c1,
                    periodicity=value.c2,
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
