import os
import shutil
import subprocess as sp
from pathlib import Path

import pytest

from kimmdy.cmd import kimmdy_run
from kimmdy.parsing import read_top, write_top
from kimmdy.topology.topology import Topology


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
    with pytest.raises(ValueError) as e:
        kimmdy_run()


@pytest.mark.parametrize(
    "arranged_tmp_path", (["test_integration/minimal_input_files"]), indirect=True
)
def test_integration_valid_input_files(arranged_tmp_path):
    kimmdy_run()
    assert "Finished running tasks" in read_last_line(Path("kimmdy.log"))


@pytest.mark.parametrize(
    "arranged_tmp_path", (["test_integration/minimal_input_files"]), indirect=True
)
def test_grompp_with_kimmdy_topology(arranged_tmp_path):
    raw_top = read_top(Path("minimal.top"))
    top = Topology(raw_top)
    top._update_dict()
    write_top(top.top, Path("output.top"))
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


@pytest.mark.parametrize(
    "arranged_tmp_path",
    (["test_integration/hexalanine_single_reaction"]),
    indirect=True,
)
def test_integration_single_reaction(arranged_tmp_path):
    kimmdy_run()
    assert "Finished running tasks" in read_last_line(Path("kimmdy.log"))


@pytest.mark.slow
@pytest.mark.parametrize(
    "arranged_tmp_path", (["test_integration/alanine_hat_naive"]), indirect=True
)
def test_integration_hat_naive_reaction(arranged_tmp_path):
    kimmdy_run()
    assert "Finished running tasks" in read_last_line(Path("kimmdy.log"))


@pytest.mark.slow
@pytest.mark.parametrize(
    "arranged_tmp_path", (["test_integration/hexalanine_homolysis"]), indirect=True
)
def test_integration_homolysis_reaction(arranged_tmp_path):
    kimmdy_run()
    assert "Finished running tasks" in read_last_line(Path("kimmdy.log"))


@pytest.mark.slow
@pytest.mark.parametrize(
    "arranged_tmp_path", (["test_integration/triplehelix_pull"]), indirect=True
)
def test_integration_pull(arranged_tmp_path):
    kimmdy_run()
    assert "Finished running tasks" in read_last_line(Path("kimmdy.log"))


@pytest.mark.slow
@pytest.mark.parametrize(
    "arranged_tmp_path",
    (["test_integration/charged_peptide_homolysis_hat_naive"]),
    indirect=True,
)
def test_integration_whole_run(arranged_tmp_path):
    kimmdy_run()
    assert "Finished running tasks" in read_last_line(Path("kimmdy.log"))
