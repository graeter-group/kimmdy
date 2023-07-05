from kimmdy.topology.topology import Topology
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

        # predict_parameters(atoms,bonds,radicals)

        return current_topology
