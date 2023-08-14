from grappa_interface import (
    clean_parameters,
    generate_input,
    apply_parameters,
    GrappaInterface,
)
from kimmdy.topology.topology import Topology
from copy import deepcopy

# from kimmdy.parsing import write_json, read_json
# import json


def test_generate_input(generic_topology):
    AA3 = [
        "ALA",
        "CYS",
        "ASP",
        "GLU",
        "PHE",
        "GLY",
        "HIE",
        "ILE",
        "LYS",
        "LEU",
        "MET",
        "ASN",
        "PRO",
        "HYP",
        "GLN",
        "ARG",
        "SER",
        "THR",
        "VAL",
        "TRP",
        "TYR",
        "DOP",
        "ACE",
        "NME",
    ]

    input_dict = generate_input(generic_topology)

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


def test_clean_parameters(generic_parameter_output):
    clean_parameters(generic_parameter_output)
    pass


def test_generate_parameters(generic_parameter_input):
    import grappa.ff

    ff = grappa.ff.ForceField.from_tag("radical_example")
    # in rad to avoid import openmm
    parameters = ff.params_from_topology_dict(generic_parameter_input)
    clean_parameters(parameters)
    # write_json(parameters, "GrAPPa_output_alanine.json")
    # with open("GrAPPa_output_alanine.json", "r") as f:
    #     parameters = json.load(f,parse_float=lambda x: round(float(x), 3))
    # write_json(parameters, "GrAPPa_output_alanine.json")


def test_apply_parameters(generic_topology):
    parameters_clean = {
        "atom": {"idxs": ["1"], "q": ["0.1"]},
        "bond": {"idxs": [["1", "3"]], "eq": ["400.0"], "k": ["0.3"]},
        "angle": {"idxs": [["2", "1", "5"]], "eq": ["1000.0"], "k": ["40.0"]},
        "proper": {
            "idxs": [["17", "16", "18", "21"]],
            "phases": [["20.0"]],
            "ks": [["0.0"]],
            "ns": [["2"]],
        },
        "improper": {
            "idxs": [["18", "14", "16", "17"]],
            "phases": [["0.0"]],
            "ks": [["2.0"]],
            "ns": [["2"]],
        },
    }

    apply_parameters(generic_topology, parameters_clean)

    assert isinstance(generic_topology, Topology)
    assert generic_topology.atoms[
        str(parameters_clean["atom"]["idxs"][0])
    ].charge == str(parameters_clean["atom"]["q"][0])
    assert generic_topology.bonds[
        tuple(str(x) for x in parameters_clean["bond"]["idxs"][0])
    ].c1 == str(parameters_clean["bond"]["k"][0])
    assert generic_topology.angles[
        tuple(str(x) for x in parameters_clean["angle"]["idxs"][0])
    ].c0 == str(parameters_clean["angle"]["eq"][0])
    assert generic_topology.proper_dihedrals[
        tuple(str(x) for x in parameters_clean["proper"]["idxs"][0])
    ].dihedrals["2"].c0 == str(parameters_clean["proper"]["phases"][0][0])
    assert generic_topology.improper_dihedrals[
        tuple(["14", "18", "16", "17"])
    ].c0 == str(parameters_clean["improper"]["phases"][0][0])


def test_parameterize_topology(generic_topology):
    parameterizer = GrappaInterface()
    curr_top = deepcopy(generic_topology)
    parameterizer.parameterize_topology(curr_top)
    assert generic_topology != curr_top
