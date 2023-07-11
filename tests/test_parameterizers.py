from grappa_interface import generate_input, apply_parameters, GrappaInterface
from kimmdy.topology.topology import Topology
import numpy as np
import json


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

    with open("GrAPPa_input_tripelhelix.json", "w") as f:
        json.dump(input_dict, f)

    assert len(input_dict["atoms"]) == 72
    assert len(input_dict["bonds"]) == 70
    assert input_dict["radicals"] == [29, 35]
    assert input_dict["atoms"][0] == [1, "CH3", "ACE", 1, [0.339967, 0.45773], 6]

    assert all([len(x) == 6 for x in input_dict["atoms"]])
    assert all([isinstance(x[s], str) for x in input_dict["atoms"] for s in [1, 2]])
    assert all([isinstance(x[i], int) for x in input_dict["atoms"] for i in [0, 3, 5]])
    assert all([isinstance(x[4], list) for x in input_dict["atoms"]])
    assert all([len(x[4]) == 2 for x in input_dict["atoms"]])
    assert all([x[2] in AA3 for x in input_dict["atoms"]])


def test_apply_parameters_assertE(generic_topology, caplog):
    parameters = {
        k: {("100"): []}
        for k in ["atoms", "bonds", "angles", "proper_dihedrals", "improper_dihedrals"]
    }

    apply_parameters(generic_topology, parameters)

    warnings = 0
    for record in caplog.records:
        if record.levelname == "WARNING":
            warnings += 1
    assert warnings == 5


def test_apply_parameters(generic_topology):
    idxs = {
        "atoms": ("1"),
        "bonds": ("1", "3"),
        "angles": ("2", "1", "5"),
        "proper_dihedrals": ("68", "67", "69", "72"),
    }
    vals = {
        "atoms": ["0.1"],
        "bonds": ["0.3", "400.0"],
        "angles": ["40.0", "1000.0"],
        "proper_dihedrals": {"2": ["0.0", "20.0"]},
    }
    parameters = {
        "atoms": {idxs["atoms"]: vals["atoms"]},
        "bonds": {idxs["bonds"]: vals["bonds"]},
        "angles": {idxs["angles"]: vals["angles"]},
        "proper_dihedrals": {idxs["proper_dihedrals"]: vals["proper_dihedrals"]},
        "improper_dihedrals": {},
    }

    apply_parameters(generic_topology, parameters)

    assert isinstance(generic_topology, Topology)
    assert generic_topology.atoms[idxs["atoms"]].charge == vals["atoms"][0]
    assert generic_topology.bonds[idxs["bonds"]].c0 == vals["bonds"][0]
    assert generic_topology.angles[idxs["angles"]].c1 == vals["angles"][1]
    assert (
        generic_topology.proper_dihedrals[idxs["proper_dihedrals"]].dihedrals["2"].c0
        == vals["proper_dihedrals"]["2"][0]
    )


def test_parameterize_topology(generic_topology):
    parameterizer = GrappaInterface()
    new_topology = parameterizer.parameterize_topology(generic_topology)
    assert new_topology == []
