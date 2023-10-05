import pytest

pytest.importorskip("grappa_interface")
pytest.importorskip("grappa")

import os
from pathlib import Path
from copy import deepcopy
import numpy as np
import json

from grappa_interface import (
    clean_parameters,
    generate_input,
    apply_parameters,
    GrappaInterface,
)
from kimmdy.topology.topology import Topology
from kimmdy.constants import AA3
from kimmdy.parsing import write_json, read_json, read_top


## fixtures ##
@pytest.fixture
def grappa_input():
    return read_json(Path(__file__).parent / "GrAPPa_input_alanine.json")


@pytest.fixture
def grappa_output():
    return read_json(Path(__file__).parent / "GrAPPa_output_alanine.json")


## test scripts ##
def test_generate_input():
    top = Topology(read_top(Path(__file__).parent / "Ala_out.top"))
    input_dict = generate_input(top)

    # with open("GrAPPa_input_alanine.json", "w") as f:
    #     json.dump(input_dict, f)
    assert len(input_dict["atoms"]) == 21
    assert len(input_dict["bonds"]) == 20
    assert input_dict["radicals"] == [9]
    assert input_dict["atoms"][0] == [1, "CH3", "ACE", 1, [0.339967, 0.45773], 6]

    assert all([len(x) == 6 for x in input_dict["atoms"]])
    assert all([isinstance(x[s], str) for x in input_dict["atoms"] for s in [1, 2]])
    assert all([isinstance(x[i], int) for x in input_dict["atoms"] for i in [0, 3, 5]])
    assert all([isinstance(x[4], list) for x in input_dict["atoms"]])
    assert all([len(x[4]) == 2 for x in input_dict["atoms"]])
    assert all([x[2] in AA3 for x in input_dict["atoms"]])


def test_clean_parameters(grappa_output):
    parameter_prepared = deepcopy(grappa_output)
    parameter_prepared["angle_k"] = np.array(parameter_prepared["angle_k"])
    clean_parameters(parameter_prepared)
    pass


def test_generate_parameters(grappa_input):
    import grappa.ff

    ff = grappa.ff.ForceField.from_tag("radical_example")
    # in rad to avoid import openmm
    parameters = ff.params_from_topology_dict(grappa_input)
    clean_parameters(parameters)
    # write_json(parameters, "GrAPPa_output_alanine.json")
    # with open("GrAPPa_output_alanine.json", "r") as f:
    #     parameters = json.load(f,parse_float=lambda x: round(float(x), 3))
    # write_json(parameters, "GrAPPa_output_alanine.json")


def test_apply_parameters():
    top = Topology(read_top(Path(__file__).parent / "Ala_out.top"))

    parameters_clean = {
        "atom": {"idxs": ["1"], "q": ["0.1"]},
        "bond": {"idxs": [["1", "3"]], "eq": ["400.0"], "k": ["0.3"]},
        "angle": {"idxs": [["2", "1", "5"]], "eq": ["1000.0"], "k": ["40.0"]},
        "proper": {
            "idxs": [["17", "16", "18", "21"]],
            "phases": [["20.0"]],
            "ks": [["0.0"]],
            "ns": [["1"]],
        },
        "improper": {
            "idxs": [["18", "14", "16", "17"]],
            "phases": [["0.0"]],
            "ks": [["2.0"]],
            "ns": [["2"]],
        },
    }

    apply_parameters(top, parameters_clean)

    assert isinstance(top, Topology)
    assert (
        top.atoms[parameters_clean["atom"]["idxs"][0]].charge
        == parameters_clean["atom"]["q"][0]
    )
    assert (
        top.bonds[tuple(parameters_clean["bond"]["idxs"][0])].c1
        == parameters_clean["bond"]["k"][0]
    )
    assert (
        top.angles[tuple(parameters_clean["angle"]["idxs"][0])].c0
        == parameters_clean["angle"]["eq"][0]
    )
    assert (
        top.proper_dihedrals[tuple(parameters_clean["proper"]["idxs"][0])]
        .dihedrals["1"]
        .c0
        == parameters_clean["proper"]["phases"][0][0]
    )
    assert (
        top.improper_dihedrals[tuple(parameters_clean["improper"]["idxs"][0])].c1
        == parameters_clean["improper"]["ks"][0][0]
    )


def test_parameterize_topology(tmp_path):
    os.chdir(tmp_path.resolve())
    parameterizer = GrappaInterface()
    top = Topology(read_top(Path(__file__).parent / "Ala_out.top"), parameterizer)

    curr_top = deepcopy(top)
    curr_top.needs_parameterization = True
    curr_top.update_parameters()
    assert top != curr_top
