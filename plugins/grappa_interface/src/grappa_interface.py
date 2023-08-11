import logging
import numpy as np
import math
from typing import Union

from kimmdy.topology.topology import Topology
from kimmdy.topology.atomic import Bond, Angle, Dihedral, MultipleDihedrals
from kimmdy.topology.utils import get_by_permutations
from kimmdy.parameterize import Parameterizer

import grappa.ff
import openmm.unit


def order_proper(idxs: np.ndarray):
    # center atoms of dihedral must have ascending value
    # idx_list is array(array(list(i)),array(list(j)),array(list(k)),array(list(l)))
    if idxs[1] < idxs[2]:
        return idxs
    else:
        return np.flip(idxs)


def check_equal_length(d: dict, name: str):
    lengths = [len(y) for y in d.values()]
    assert (
        len(set(lengths)) == 1
    ), f"Different length of {name} parameters: { {k:len(v) for k, v in d.items()} }"


def convert_to_python_types(array: Union[list, np.ndarray]) -> list:
    return getattr(array, "tolist", lambda: array)()


def order_proper(idxs: np.ndarray):
    # center atoms of dihedral must have ascending value
    # idx_list is array(array(list(i)),array(list(j)),array(list(k)),array(list(l)))
    if idxs[1] < idxs[2]:
        return idxs
    else:
        return np.flip(idxs)


def clean_parameters(parameters: dict) -> None:
    harmonic_keys = {"idxs", "eq", "k"}
    dihedral_keys = {"idxs", "phases", "ks", "ns"}
    parameters_clean = {
        "atom": {"idxs": [], "q": []},
        "bond": {k: [] for k in harmonic_keys},
        "angle": {k: [] for k in harmonic_keys},
        "proper": {k: [] for k in dihedral_keys},
        "improper": {k: [] for k in dihedral_keys},
    }
    parameters["proper_idxs"] = [order_proper(x) for x in parameters["proper_idxs"]]
    try:
        for atomic in parameters_clean.keys():
            for parameter in parameters_clean[atomic].keys():
                key = atomic + "_" + parameter
                parameters_clean[atomic][parameter] = convert_to_python_types(
                    parameters[key]
                )
    except KeyError:
        raise KeyError(
            f"GrAPPa returned parameters {list(parameters.keys())}, which do not contain the required sections to fill {parameters_clean}"
        )

    for name, atomic in parameters_clean.items():
        check_equal_length(atomic, name)

    ## sample type check
    assert isinstance(
        parameters_clean["atom"]["idxs"][0], int
    ), f"atom idxs element has wrong type {type(parameters_clean['atom']['idxs'][0])}, should be int."
    assert isinstance(
        parameters_clean["bond"]["k"][0], float
    ), f"b k element has wrong type {type(parameters_clean['bond']['k'][0])}, should be float."
    assert isinstance(
        parameters_clean["proper"]["ns"][0], list
    ), f"proper ns element has wrong type {type(parameters_clean['proper']['ns'][0])}, should be list."
    assert isinstance(
        parameters_clean["improper"]["phases"][0][0], float
    ), f"improper phases element has wrong type {type(parameters_clean['improper']['phases'][0][0])}, should be float."
    return parameters_clean


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
    # parameter structure is defined in clean_parameters()
    # assume units are according to https://manual.gromacs.org/current/reference-manual/definitions.html
    # namely: length [nm], mass [kg], time [ps], energy [kJ/mol], force [kJ mol-1 nm-1], angle [deg]

    ## atoms
    for idx in range(len(parameters["atom"]["idxs"])):
        if not (atom := top.atoms.get(str(parameters["atom"]["idxs"][idx]))):
            raise KeyError(
                f"bad index {str(parameters['atom']['idxs'][idx])} in {list(top.atoms.keys())}"
            )
            logging.warning(
                f"Ignored parameters with invalid ids: {parameters['atom']['idxs'][idx]} for atoms"
            )
            continue
        # can anything but charge change??
        atom.charge = str(parameters["atom"]["q"][idx])
        atom.chargeB = None

    ## bonds
    for idx in range(len(parameters["bond"]["idxs"])):
        if not top.bonds.get(tuple(str(x) for x in parameters["bond"]["idxs"][idx])):
            raise KeyError(
                f"bad index {tuple(str(x) for x in parameters['bond']['idxs'][idx])} in {list(top.bonds.keys())}"
            )
            logging.warning(
                f"Ignored parameters with invalid ids: {parameters['bond']['idxs'][idx]} for bonds"
            )
            continue
        top.bonds[
            tuple(str(x) for x in parameters["bond"]["idxs"][idx])
        ] = Bond.from_top_line(
            [
                str(x)
                for x in [
                    *parameters["bond"]["idxs"][idx],
                    "1",
                    parameters["bond"]["eq"][idx],
                    parameters["bond"]["k"][idx],
                ]
            ]
        )

    ## angles
    for idx in range(len(parameters["angle"]["idxs"])):
        if not top.angles.get(tuple(str(x) for x in parameters["angle"]["idxs"][idx])):
            raise KeyError(
                f"bad index {tuple(str(x) for x in parameters['angle']['idxs'][idx])} in {list(top.angles.keys())}"
            )
            logging.warning(
                f"Ignored parameters with invalid ids: {parameters['angle']['idxs'][idx]} for angles"
            )
            continue
        top.angles[
            tuple(str(x) for x in parameters["angle"]["idxs"][idx])
        ] = Angle.from_top_line(
            [
                str(x)
                for x in [
                    *parameters["angle"]["idxs"][idx],
                    "1",
                    parameters["angle"]["eq"][idx],
                    parameters["angle"]["k"][idx],
                ]
            ]
        )

    ## proper dihedrals
    for idx in range(len(parameters["proper"]["idxs"])):
        if not top.proper_dihedrals.get(
            tuple(str(x) for x in parameters["proper"]["idxs"][idx])
        ):
            raise KeyError(
                f"bad index {tuple(str(x) for x in parameters['proper']['idxs'][idx])} in {list(top.proper_dihedrals.keys())}"
            )
            logging.warning(
                f"Ignored parameters with invalid ids: {parameters['proper']['idxs'][idx]} for proper dihedrals"
            )
            continue
        dihedral_dict = {}
        for i, n in enumerate(parameters["proper"]["ns"][idx]):
            dihedral_dict[str(n)] = Dihedral.from_top_line(
                [
                    str(x)
                    for x in [
                        *parameters["proper"]["idxs"][idx],
                        "9",
                        parameters["proper"]["phases"][i],
                        parameters["proper"]["ks"][i],
                        n,
                    ]
                ]
            )

        top.proper_dihedrals[
            tuple(str(x) for x in parameters["proper"]["idxs"][idx])
        ] = MultipleDihedrals(
            *[str(x) for x in parameters["proper"]["idxs"][idx]], "9", dihedral_dict
        )

    ## improper dihedrals
    for idx in range(len(parameters["improper"]["idxs"])):
        if not (
            term := get_by_permutations(
                top.improper_dihedrals,
                tuple(str(x) for x in parameters["improper"]["idxs"][idx]),
            )
        ):
            raise KeyError(
                f"bad index {tuple(str(x) for x in parameters['improper']['idxs'][idx])} in {list(top.improper_dihedrals.keys())}"
            )
            logging.warning(
                f"Ignored parameters with invalid ids: {parameters['improper']['idxs'][idx]} for improper dihedrals"
            )
            continue
        # dihedral_prm order should be [phase, kd]
        for i, n in enumerate(parameters["improper"]["ns"][idx]):
            term.periodicity = "2" if term.periodicity == "" else term.periodicity
            if str(n) != term.periodicity and not math.isclose(
                parameters["improper"]["ks"][idx][i], 0.0
            ):
                logging.warning(
                    f"Ignored improper with periodicity of {n} for idxs {parameters['improper']['idxs'][idx]}. This term should not be part of an amber style force field."
                )
            else:
                term.c0 = str(parameters["improper"]["phases"][idx][i])
                term.c1 = str(parameters["improper"]["ks"][idx][i])

    return top


class GrappaInterface(Parameterizer):
    def parameterize_topology(
        self, current_topology: Topology, focus_nr: list[str] = []
    ) -> Topology:
        ## get atoms, bonds, radicals in required format
        input_dict = generate_input(current_topology)

        ff = grappa.ff.ForceField.from_tag("radical_example")
        ff.units["angle"] = openmm.unit.degree
        # gromacs angle force constant are already in kJ/mol/rad-2]
        parameters = ff.params_from_topology_dict(input_dict)

        parameters = clean_parameters(parameters)

        apply_parameters(current_topology, parameters)
        return current_topology
