from kimmdy.topology.topology import Topology
from kimmdy.topology.atomic import Bond, Angle, Dihedral, MultipleDihedrals
from kimmdy.parameterize import Parameterizer


class GrappaInterface(Parameterizer):
    def parameterize_topology(
        self, current_topology: Topology, focus_nr: list[str] = []
    ) -> Topology:
        at_num_to_ele = {
            "1": "H",
            "3": "Li",
            "6": "C",
            "7": "N",
            "8": "O",
            "9": "F",
            "11": "Na",
            "12": "Mg",
            "15": "P",
            "16": "S",
            "17": "Cl",
            "19": "K",
            "20": "Ca",
            "26": "Fe",
            "29": "Cu",
            "30": "Zn",
            "35": "Br",
            "37": "Rb",
            "53": "I",
            "55": "Cs",
            "0": "Du",
        }
        at_map = current_topology.ff.atomtypes

        ## get atoms, bonds, radicals in required format
        atoms = [
            [
                atom.residue,
                [at_map[atom.type].sigma, at_map[atom.type].epsilon],
                at_num_to_ele[at_map[atom.type].at_num],
            ]
            for atom in current_topology.atoms.values()
        ]
        bonds = [(int(bond[0]), int(bond[1])) for bond in current_topology.bonds.keys()]
        radicals = [int(radical) for radical in current_topology.radicals.keys()]

        parameters: dict = predict_parameters(atoms, bonds, radicals)
        # assume parameters is {'atoms':{},'bonds':{},'angles':{},'proper_dihedrals:{},'improper_dihedrals:{}}
        # assume units are according to https://manual.gromacs.org/current/reference-manual/definitions.html
        # namely: length [nm], mass [kg], time [ps], energy [kJ/mol], force [kJ mol-1 nm-1], angle [deg]

        ## atoms
        for atom_id, atom_prm in parameters["atoms"].items():
            # can anything but charge change??
            current_topology.atoms[atom_id].charge = atom_prm[0]
            current_topology.atoms[atom_id].chargeB = None

        ## bonds
        for bond_ids, bond_prm in parameters["bonds"].items():
            # bond_prm order should be [b0, kb]
            current_topology.bonds[bond_ids] = Bond(*bond_ids, "1", *bond_prm)

        ## angles
        for angle_ids, angle_prm in parameters["angles"].items():
            # angle_prm order should be [th0, cth]
            current_topology.angles[angle_ids] = Angle(*angle_ids, "1", *angle_prm)

        ## proper dihedrals
        for dihedral_ids, dihedral in parameters["proper_dihedrals"].items():
            dihedral_dict = {}
            for periodicity, dihedral_prm in dihedral.items():
                # dihedral_prm order should be [phase, kd]
                dihedral_dict[periodicity] = Dihedral(
                    *dihedral_ids, "9", *dihedral_prm, periodicity
                )

            current_topology.proper_dihedrals[dihedral_ids] = MultipleDihedrals(
                *dihedral_ids, "9", dihedral_dict
            )

        ## improper dihedrals
        for dihedral_ids, dihedral_prm in parameters["improper_dihedrals"].items():
            # dihedral_prm order should be [phase, kd]
            current_topology.proper_dihedrals[dihedral_ids] = Dihedral(
                *dihedral_ids, "4", dihedral_prm, "2"
            )

        return current_topology
