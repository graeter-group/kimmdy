import os
import subprocess as sp
from pathlib import Path

import pytest

from kimmdy.cmd import kimmdy_run
from kimmdy.parsing import read_top, write_top
from kimmdy.topology.topology import Topology
from kimmdy.plugins import parameterization_plugins


def read_last_line(file):
    with open(file, "rb") as f:
        try:  # catch OSError in case of a one line file
            f.seek(-2, os.SEEK_END)
            while f.read(1) != b"\n":
                f.seek(-2, os.SEEK_CUR)
        except OSError:
            f.seek(0)
        return f.readline().decode()


@pytest.mark.parametrize(
    "arranged_tmp_path", (["test_integration/emptyrun"]), indirect=True
)
def test_integration_emptyrun(arranged_tmp_path):
    # not expecting this to run
    # because the topology is empty
    Path("emptyrun.txt").touch()
    with pytest.raises(ValueError):
        kimmdy_run()
    assert len(list(Path.cwd().glob("emptyrun_001/*"))) == 2


@pytest.mark.parametrize(
    "arranged_tmp_path", (["test_integration/minimal_input_files"]), indirect=True
)
def test_integration_valid_input_files(arranged_tmp_path):
    kimmdy_run()
    assert "Finished running tasks" in read_last_line(Path("kimmdy.log"))
    assert len(list(Path.cwd().glob("minimal/*"))) == 3


@pytest.mark.parametrize(
    "arranged_tmp_path", (["test_integration/minimal_input_files"]), indirect=True
)
def test_grompp_with_kimmdy_topology(arranged_tmp_path):
    raw_top = read_top(Path("minimal.top"))
    top = Topology(raw_top)
    top_dict = top.to_dict()
    write_top(top.to_dict(), Path("output.top"))
    assert sp.run(
        [
            "gmx",
            "grompp",
            "-f",
            "minimal.mdp",
            "-c",
            "minimal.gro",
            "-p",
            "output.top",
            "-o",
            "minial.tpr",
        ],
        check=True,
    )


@pytest.mark.require_grappa
@pytest.mark.parametrize(
    "arranged_tmp_path", (["test_integration/triplehelix_pull"]), indirect=True
)
def test_grappa_partial_parameterization(arranged_tmp_path):
    parametrizer = parameterization_plugins["grappa"]()

    # parameterize everything
    top = Topology(read_top(Path("topol.top")), parametrizer=parametrizer)
    top.needs_parameterization = True
    top.update_parameters()

    # parameterize around delete
    top_del = Topology(read_top(Path("topol.top")), parametrizer=parametrizer)
    top_del.needs_parameterization = True
    top_del.update_parameters()
    top_del.break_bond(("1", "5"))
    top_del.needs_parameterization = True
    top_del.update_parameters(focus_nrs=set(["1", "5"]))

    # check values are changed around homolysis atoms
    assert top.bonds[("5", "6")] != top_del.bonds[("5", "6")]
    assert top.angles[("6", "5", "7")] != top_del.angles[("6", "5", "7")]
    assert (
        top.proper_dihedrals[("13", "11", "14", "15")]
        != top_del.proper_dihedrals[("13", "11", "14", "15")]
    )
    assert (
        top.improper_dihedrals[("8", "5", "7", "9")]
        != top_del.improper_dihedrals[("8", "5", "7", "9")]
    )

    # check values are the same further away
    assert top.bonds[("200", "203")] == top_del.bonds[("200", "203")]
    assert top.angles[("228", "230", "232")] == top_del.angles[("228", "230", "232")]
    assert (
        top.proper_dihedrals[("1714", "1713", "1715", "1717")]
        == top_del.proper_dihedrals[("1714", "1713", "1715", "1717")]
    )
    assert (
        top.improper_dihedrals[("349", "357", "355", "356")]
        == top_del.improper_dihedrals[("349", "357", "355", "356")]
    )


@pytest.mark.parametrize(
    "arranged_tmp_path",
    (["test_integration/hexalanine_single_reaction"]),
    indirect=True,
)
def test_integration_single_reaction(arranged_tmp_path):
    kimmdy_run()
    assert "Finished running tasks" in read_last_line(Path("kimmdy.log"))
    assert len(list(Path.cwd().glob("single_reaction_000/*"))) == 8


@pytest.mark.slow
@pytest.mark.parametrize(
    "arranged_tmp_path", (["test_integration/alanine_hat_naive"]), indirect=True
)
def test_integration_hat_naive_reaction(arranged_tmp_path):
    kimmdy_run()
    assert "Finished running tasks" in read_last_line(Path("kimmdy.log"))
    assert len(list(Path.cwd().glob("alanine_hat_000/*"))) == 16


@pytest.mark.slow
@pytest.mark.parametrize(
    "arranged_tmp_path", (["test_integration/hexalanine_homolysis"]), indirect=True
)
def test_integration_homolysis_reaction(arranged_tmp_path):
    kimmdy_run()
    assert "Finished running tasks" in read_last_line(Path("kimmdy.log"))
    assert len(list(Path.cwd().glob("hexalanine_homolysis_000/*"))) == 13


@pytest.mark.slow
@pytest.mark.parametrize(
    "arranged_tmp_path", (["test_integration/triplehelix_pull"]), indirect=True
)
def test_integration_pull(arranged_tmp_path):
    kimmdy_run()
    assert "Finished running tasks" in read_last_line(Path("kimmdy.log"))
    assert len(list(Path.cwd().glob("kimmdy_001/*"))) == 12


@pytest.mark.require_grappa
@pytest.mark.slow
@pytest.mark.parametrize(
    "arranged_tmp_path",
    (["test_integration/charged_peptide_homolysis_hat_naive"]),
    indirect=True,
)
def test_integration_whole_run(arranged_tmp_path):
    kimmdy_run()
    assert "Finished running tasks" in read_last_line(Path("kimmdy.log"))
    assert len(list(Path.cwd().glob("kimmdy_001/*"))) == 25
