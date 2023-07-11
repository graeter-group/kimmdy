import logging
import numpy as np
import math

from kimmdy.topology.topology import Topology
from kimmdy.topology.atomic import Bond, Angle, Dihedral, MultipleDihedrals
from kimmdy.parameterize import Parameterizer

import grappa.ff
import openmm.unit


def in_degree(vals: np.ndarray):
    if max(vals.flat) <= 2 * math.pi:
        return vals * 180.0 / math.pi
    else:
        return vals


def order_proper(idxs: np.ndarray):
    # center atoms of dihedral must have ascending value
    # idx_list is array(array(list(i)),array(list(j)),array(list(k)),array(list(l)))
    if idxs[1] < idxs[2]:
        return idxs
    else:
        return np.flip(idxs)


def generate_input(top: Topology) -> dict:
    at_map = top.ff.atomtypes
    atoms = [
        [
            int(atom.nr),
            atom.atom,
            atom.residue,
            int(atom.resnr),
            [float(at_map[atom.type].sigma), float(at_map[atom.type].epsilon)],
            int(at_map[atom.type].at_num),
        ]
        for atom in top.atoms.values()
    ]
    bonds = [(int(bond[0]), int(bond[1])) for bond in top.bonds.keys()]
    radicals = [int(radical) for radical in top.radicals.keys()]

    return {"atoms": atoms, "bonds": bonds, "radicals": radicals}


def apply_parameters(top: Topology, parameters: dict) -> Topology:
    # assume parameters keys are ("atom_idxs","atom_q","atom_sigma","atom_epsilon","atom_mass",
    # "bond_idxs","bond_k","bond_eq","angle_idxs","angle_k","angle_eq","proper_idxs","proper_ks",
    # "proper_ns","proper_phases","improper_idxs","improper_ks","improper_ns","improper_phases",)
    # assume units are according to https://manual.gromacs.org/current/reference-manual/definitions.html
    # namely: length [nm], mass [kg], time [ps], energy [kJ/mol], force [kJ mol-1 nm-1], angle [deg]

    assert list(parameters.keys()) == [
        "atom_idxs",
        "atom_q",
        "atom_sigma",
        "atom_epsilon",
        "atom_mass",
        "bond_idxs",
        "bond_k",
        "bond_eq",
        "angle_idxs",
        "angle_k",
        "angle_eq",
        "proper_idxs",
        "proper_ks",
        "proper_ns",
        "proper_phases",
        "improper_idxs",
        "improper_ks",
        "improper_ns",
        "improper_phases",
    ], f"GrAPPa returned unexpected parameters with sections {list(parameters.keys())}."
    assert all(
        [
            len(parameters[y]) == len(parameters["atom_idxs"])
            if y.startswith("atom")
            else True
            for y in parameters.keys()
        ]
    ), f"Different length of atom parameters: { {k:len(v) for k, v in parameters.items()} }"
    assert all(
        [
            len(parameters[y]) == len(parameters["bond_idxs"])
            if y.startswith("bond")
            else True
            for y in parameters.keys()
        ]
    ), f"Different length of bond parameters: { {k:len(v) for k, v in parameters.items()} }"
    assert all(
        [
            len(parameters[y]) == len(parameters["angle_idxs"])
            if y.startswith("angle")
            else True
            for y in parameters.keys()
        ]
    ), f"Different length of angle parameters: { {k:len(v) for k, v in parameters.items()} }"
    assert all(
        [
            len(parameters[y]) == len(parameters["proper_idxs"])
            if y.startswith("proper")
            else True
            for y in parameters.keys()
        ]
    ), f"Different length of proper parameters: { {k:len(v) for k, v in parameters.items()} }"
    assert all(
        [
            len(parameters[y]) == len(parameters["improper_idxs"])
            if y.startswith("imporper")
            else True
            for y in parameters.keys()
        ]
    ), f"Different length of improper parameters: { {k:len(v) for k, v in parameters.items()} }"

    for tup in parameters["proper_idxs"]:
        p1, p2, p3, p4 = tup
        assert p1 == [], f"{p1},{p2},{p3},{p4}"

    # parameters are way too nested at the moment. some functions need to change if this gets fixed
    # maybe write clean_parameters function to convert to tuples and string

    parameters["proper_phases"] = in_degree(parameters["proper_phases"])
    parameters["improper_phases"] = in_degree(parameters["improper_phases"])
    parameters["proper_idxs"] = np.asarray(
        [order_proper(x) for x in parameters["proper_idxs"]]
    )

    ## atoms
    for idx in range(len(parameters["atom_idxs"])):
        if not (atom := top.atoms.get(str(*parameters["atom_idxs"][idx]))):
            raise KeyError(
                f"bad index {str(*parameters['atom_idxs'][idx])} in {list(top.atoms.keys())}"
            )
            logging.warning(
                f"Ignored parameters with invalid ids: {parameters['atom_idxs'][idx]} for atoms"
            )
            continue
        # can anything but charge change??
        atom.charge = str(parameters["atom_q"][idx])
        atom.chargeB = None

    ## bonds
    for idx in range(len(parameters["bond_idxs"])):
        if not top.bonds.get(tuple(str(*x) for x in parameters["bond_idxs"][idx])):
            raise KeyError(
                f"bad index {tuple(str(*x) for x in parameters['bond_idxs'][idx])} in {list(top.bonds.keys())}"
            )
            logging.warning(
                f"Ignored parameters with invalid ids: {parameters['bond_idxs'][idx]} for bonds"
            )
            continue
        top.bonds[
            tuple(str(*x) for x in parameters["bond_idxs"][idx])
        ] = Bond.from_top_line(
            [
                str(x)
                for x in [
                    *parameters["bond_idxs"][idx],
                    "1",
                    parameters["bond_eq"][idx],
                    parameters["bond_k"][idx],
                ]
            ]
        )

    ## angles
    for idx in range(len(parameters["angle_idxs"])):
        if not top.angles.get(tuple(str(*x) for x in parameters["angle_idxs"][idx])):
            raise KeyError(
                f"bad index {tuple(str(*x) for x in parameters['angle_idxs'][idx])} in {list(top.angles.keys())}"
            )
            logging.warning(
                f"Ignored parameters with invalid ids: {parameters['angle_idxs'][idx]} for angles"
            )
            continue
        top.angles[
            tuple(str(*x) for x in parameters["angle_idxs"][idx])
        ] = Angle.from_top_line(
            [
                str(x)
                for x in [
                    *parameters["angle_idxs"][idx],
                    "1",
                    parameters["angle_eq"][idx],
                    parameters["angle_k"][idx],
                ]
            ]
        )

    ## proper dihedrals
    for idx in range(len(parameters["proper_idxs"])):
        if not top.proper_dihedrals.get(
            tuple(str(*x) for x in parameters["proper_idxs"][idx])
        ):
            raise KeyError(
                f"bad index {tuple(str(*x) for x in parameters['proper_idxs'][idx])} in {list(top.proper_dihedrals.keys())}"
            )
            logging.warning(
                f"Ignored parameters with invalid ids: {parameters['proper_idxs'][idx]} for proper dihedrals"
            )
            continue
        dihedral_dict = {}
        for i, n in enumerate(parameters["proper_ns"][0]):
            dihedral_dict[str(n)] = Dihedral.from_top_line(
                [
                    str(x)
                    for x in [
                        *parameters["proper_idxs"][idx],
                        "9",
                        parameters["proper_phases"][i],
                        parameters["proper_ks"][i],
                        n,
                    ]
                ]
            )

        top.proper_dihedrals[
            tuple(str(*x) for x in parameters["proper_idxs"][idx])
        ] = MultipleDihedrals(
            *[str(*x) for x in parameters["proper_idxs"][idx]], "9", dihedral_dict
        )
    assert top.proper_dihedrals == [], f"{top.proper_dihedrals}"

    ## improper dihedrals
    for idx in range(len(parameters["improper_idxs"])):
        if not top.improper_dihedrals.get(
            tuple(str(*x) for x in parameters["improper_idxs"][idx])
        ):
            raise KeyError(
                f"bad index {tuple(str(*x) for x in parameters['improper_idxs'][idx])} in {list(top.improper_dihedrals.keys())}"
            )
            logging.warning(
                f"Ignored parameters with invalid ids: {parameters['improper_idxs'][idx]} for improper dihedrals"
            )
            continue
        # dihedral_prm order should be [phase, kd]
        for i, n in enumerate(parameters["improper_ns"]):
            if str(n) != "2":
                logging.warning(
                    f"Ignored improper with periodicity of {n} for idxs {parameters['improper_idxs'][idx]}. This term should not be part of an amber style force field."
                )
            else:
                top.improper_dihedrals[
                    tuple(str(*x) for x in parameters["improper_idxs"][idx])
                ] = Dihedral.from_top_line(
                    [
                        *[str(*x) for x in parameters["improper_idxs"][idx]],
                        "4",
                        parameters["improper_phases"][idx][i],
                        parameters["improper_ks"][idx][i],
                        "2",
                    ]
                )

    return top


class GrappaInterface(Parameterizer):
    def parameterize_topology(
        self, current_topology: Topology, focus_nr: list[str] = []
    ) -> Topology:
        ## get atoms, bonds, radicals in required format
        input_dict = generate_input(current_topology)

        mpath = "/hits/fast/mbm/seutelf/grappa/mains/runs/param_search/versions/5__dotted_512_6_6_full_ds/best_model.pt"
        ff = grappa.ff.ForceField(model_path=mpath)
        ff.units["angle"] = openmm.unit.degree
        # gromacs angle force constant are already in kJ/mol/rad-2]
        parameters = ff.params_from_topology_dict(input_dict)

        apply_parameters(current_topology, parameters)

        return current_topology
