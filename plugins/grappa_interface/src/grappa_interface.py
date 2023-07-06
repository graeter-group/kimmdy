import logging

from kimmdy.topology.topology import Topology
from kimmdy.topology.atomic import Bond, Angle, Dihedral, MultipleDihedrals
from kimmdy.parameterize import Parameterizer


def generate_input(top: Topology) -> tuple[list, list, list]:
    at_map = top.ff.atomtypes
    atoms = [
        [
            int(atom.nr),
            atom.residue,
            atom.type,
            [float(at_map[atom.type].sigma), float(at_map[atom.type].epsilon)],
            int(at_map[atom.type].at_num),
        ]
        for atom in top.atoms.values()
    ]
    bonds = [(int(bond[0]), int(bond[1])) for bond in top.bonds.keys()]
    radicals = [int(radical) for radical in top.radicals.keys()]

    return {"atoms": atoms, "bonds": bonds, "radicals": radicals}


def apply_parameters(top: Topology, parameters: dict) -> Topology:
    ## atoms
    for atom_id, atom_prm in parameters["atoms"].items():
        if not top.atoms.get(atom_id):
            logging.warning(f"Ignored parameters with invalid ids: {atom_id} for atoms")
            continue
        # can anything but charge change??
        top.atoms[atom_id].charge = atom_prm[0]
        top.atoms[atom_id].chargeB = None

    ## bonds
    for bond_ids, bond_prm in parameters["bonds"].items():
        if not top.bonds.get(bond_ids):
            logging.warning(
                f"Ignored parameters with invalid ids: {bond_ids} for bonds"
            )
            continue
        # bond_prm order should be [b0, kb]
        top.bonds[bond_ids] = Bond(*bond_ids, "1", *bond_prm)

    ## angles
    for angle_ids, angle_prm in parameters["angles"].items():
        if not top.angles.get(angle_ids):
            logging.warning(
                f"Ignored parameters with invalid ids: {angle_ids} for angles"
            )
            continue
        # angle_prm order should be [th0, cth]
        top.angles[angle_ids] = Angle(*angle_ids, "1", *angle_prm)

    ## proper dihedrals
    for dihedral_ids, dihedral in parameters["proper_dihedrals"].items():
        if not top.proper_dihedrals.get(dihedral_ids):
            logging.warning(
                f"Ignored parameters with invalid ids: {dihedral_ids} for proper dihedrals"
            )
            continue
        dihedral_dict = {}
        for periodicity, dihedral_prm in dihedral.items():
            # dihedral_prm order should be [phase, kd]
            dihedral_dict[periodicity] = Dihedral(
                *dihedral_ids, "9", *dihedral_prm, periodicity
            )

        top.proper_dihedrals[dihedral_ids] = MultipleDihedrals(
            *dihedral_ids, "9", dihedral_dict
        )

    ## improper dihedrals
    for dihedral_ids, dihedral_prm in parameters["improper_dihedrals"].items():
        if not top.improper_dihedrals.get(dihedral_ids):
            logging.warning(
                f"Ignored parameters with invalid ids: {dihedral_ids} for improper dihedrals"
            )
            continue
        # dihedral_prm order should be [phase, kd]
        top.improper_dihedrals[dihedral_ids] = Dihedral(
            *dihedral_ids, "4", dihedral_prm, "2"
        )

    return top


class GrappaInterface(Parameterizer):
    def parameterize_topology(
        self, current_topology: Topology, focus_nr: list[str] = []
    ) -> Topology:
        ## get atoms, bonds, radicals in required format
        input_dict = generate_input(current_topology)

        parameters: dict = predict_parameters(input_dict)
        # assume parameters is {'atoms':{},'bonds':{},'angles':{},'proper_dihedrals:{},'improper_dihedrals:{}}
        # assume units are according to https://manual.gromacs.org/current/reference-manual/definitions.html
        # namely: length [nm], mass [kg], time [ps], energy [kJ/mol], force [kJ mol-1 nm-1], angle [deg]

        assert list(parameters.keys()) == [
            "atoms",
            "bonds",
            "angles",
            "proper_dihedrals",
            "improper_dihedrals",
        ], f"GrAPPa returned unexpected parameters with sections {list(parameters.keys())}."
        apply_parameters(current_topology, parameters)

        return current_topology
