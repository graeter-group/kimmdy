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


def elements_to_string(l: list):
    for i, e in enumerate(l):
        if isinstance(e, list):
            elements_to_string(e)
        else:
            l[i] = str(e)


def clean_parameters(parameters: dict) -> dict:
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
                elements_to_string(parameters_clean[atomic][parameter])
    except KeyError:
        raise KeyError(
            f"GrAPPa returned parameters {list(parameters.keys())}, which do not contain the required sections to fill {parameters_clean}"
        )

    for name, atomic in parameters_clean.items():
        check_equal_length(atomic, name)

    ## sample type check
    assert parameters_clean["atom"]["idxs"][
        0
    ].isdigit(), f"atom idxs element does not look like int {parameters_clean['atom']['idxs'][0]}."
    assert (
        parameters_clean["bond"]["k"][0]
        .strip()
        .lstrip("-")
        .replace(".", "", 1)
        .isdigit()
    ), f"b k element does not look like float {parameters_clean['bond']['k'][0]}."
    assert isinstance(
        parameters_clean["proper"]["ns"][0], list
    ), f"proper ns element has wrong type {type(parameters_clean['proper']['ns'][0])}, should be list."
    assert (
        parameters_clean["improper"]["phases"][0][0]
        .strip()
        .lstrip("-")
        .replace(".", "", 1)
        .isdigit()
    ), f"improper phases element does not look like float {type(parameters_clean['improper']['phases'][0][0])}"
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


def apply_parameters(top: Topology, parameters: dict):
    # parameter structure is defined in clean_parameters()
    # assume units are according to https://manual.gromacs.org/current/reference-manual/definitions.html
    # namely: length [nm], mass [kg], time [ps], energy [kJ/mol], force [kJ mol-1 nm-1], angle [deg]

    ## atoms
    for i, idx in enumerate(parameters["atom"]["idxs"]):
        if not (atom := top.atoms.get(idx)):
            raise KeyError(f"bad index {idx} in {list(top.atoms.keys())}")
            logging.warning(f"Ignored parameters with invalid ids: {idx} for atoms")
            continue
        # can anything but charge change??
        atom.charge = parameters["atom"]["q"][i]
        atom.chargeB = None

    ## bonds
    for i, idx in enumerate(parameters["bond"]["idxs"]):
        tup = tuple(idx)
        if not top.bonds.get(tup):
            raise KeyError(f"bad index {tup} in {list(top.bonds.keys())}")
            logging.warning(f"Ignored parameters with invalid ids: {tup} for bonds")
            continue
        top.bonds[tup] = Bond(
            *parameters["bond"]["idxs"][i],
            funct="1",
            c0=parameters["bond"]["eq"][i],
            c1=parameters["bond"]["k"][i],
        )

    ## angles
    for i, idx in enumerate(parameters["angle"]["idxs"]):
        tup = tuple(idx)
        if not top.angles.get(tup):
            raise KeyError(f"bad index {tup} in {list(top.angles.keys())}")
            logging.warning(f"Ignored parameters with invalid ids: {tup} for angles")
            continue
        top.angles[tup] = Angle(
            *parameters["angle"]["idxs"][i],
            funct="1",
            c0=parameters["angle"]["eq"][i],
            c1=parameters["angle"]["k"][i],
        )

    ## proper dihedrals
    for i, idx in enumerate(parameters["proper"]["idxs"]):
        tup = tuple(idx)
        if not top.proper_dihedrals.get(tup):
            raise KeyError(f"bad index {tup} in {list(top.proper_dihedrals.keys())}")
            logging.warning(
                f"Ignored parameters with invalid ids: {tup} for proper dihedrals"
            )
            continue
        dihedral_dict = {}
        for ii, n in enumerate(parameters["proper"]["ns"][i]):
            dihedral_dict[n] = Dihedral(
                *tup,
                funct="9",
                c0=parameters["proper"]["phases"][i][ii],
                c1=parameters["proper"]["ks"][i][ii],
                periodicity=n,
            )
        top.proper_dihedrals[tup] = MultipleDihedrals(
            *tup, funct="9", dihedrals=dihedral_dict
        )

    ## improper dihedrals
    for i, idx in enumerate(parameters["improper"]["idxs"]):
        tup = tuple(idx)
        if not (term := get_by_permutations(top.improper_dihedrals, tup)):
            raise KeyError(f"bad index {tup} in {list(top.improper_dihedrals.keys())}")
            logging.warning(
                f"Ignored parameters with invalid ids: {tup} for improper dihedrals"
            )
            continue
        # dihedral_prm order should be [phase, kd]
        for ii, n in enumerate(parameters["improper"]["ns"][i]):
            term.periodicity = "2" if term.periodicity == "" else term.periodicity
            if n != term.periodicity and not math.isclose(
                parameters["improper"]["ks"][i][ii], 0.0
            ):
                logging.warning(
                    f"Ignored improper with periodicity of {n} for idxs {tup}. This term should not be part of an amber style force field."
                )
            else:
                term.c0 = parameters["improper"]["phases"][i][ii]
                term.c1 = parameters["improper"]["ks"][i][ii]

    return


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
