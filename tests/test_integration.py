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


def setup_testdir(tmp_path, dirname) -> Path:
    try:
        filedir = Path(__file__).parent / "test_files" / "test_integration" / dirname
        assetsdir = Path(__file__).parent / "test_files" / "assets"
    except NameError:
        filedir = Path("./tests/test_files") / "test_integration" / dirname
        assetsdir = Path("./tests/test_files") / "assets"
    testdir = tmp_path / "test_integration" / dirname
    shutil.copytree(filedir, testdir)
    if not Path(testdir / "amber99sb-star-ildnp.ff").exists():
        Path(testdir / "amber99sb-star-ildnp.ff").symlink_to(
            assetsdir / "amber99sb-star-ildnp.ff",
            target_is_directory=True,
        )
    os.chdir(testdir.resolve())
    return testdir


def test_integration_emptyrun(tmp_path):
    testdir = setup_testdir(tmp_path, "emptyrun")
    (testdir / "emptyrun.txt").touch()

    # not expecting this to run
    # because the topology is empty
    with pytest.raises(ValueError) as e:
        kimmdy_run()


def test_integration_valid_input_files(tmp_path):
    testdir = setup_testdir(tmp_path, "minimal_input_files")
    (testdir / "emptyrun.txt").touch()

    kimmdy_run()
    assert "Finished running tasks" in read_last_line(testdir / "kimmdy.log")


def test_grompp_with_kimmdy_topology(tmp_path):
    testdir = setup_testdir(tmp_path, "minimal_input_files")
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


def test_integration_single_reaction(tmp_path):
    testdir = setup_testdir(tmp_path, "single_reaction")

    kimmdy_run()
    assert "Finished running tasks" in read_last_line(testdir / "kimmdy.log")


@pytest.mark.slow
def test_integration_hat_naive_reaction(tmp_path):
    testdir = setup_testdir(tmp_path, "hat_naive")

    kimmdy_run()
    assert "Finished running tasks" in read_last_line(testdir / "kimmdy.log")


@pytest.mark.slow
def test_integration_homolysis_reaction(tmp_path):
    testdir = setup_testdir(tmp_path, "homolysis")

    kimmdy_run()
    assert "Finished running tasks" in read_last_line(testdir / "kimmdy.log")


@pytest.mark.slow
def test_integration_pull(tmp_path):
    testdir = setup_testdir(tmp_path, "pull")

    kimmdy_run()
    assert "Finished running tasks" in read_last_line(testdir / "kimmdy.log")


@pytest.mark.slow
def test_integration_whole_run(tmp_path):
    testdir = setup_testdir(tmp_path, "whole_run")

    kimmdy_run()
    assert "Finished running tasks" in read_last_line(testdir / "kimmdy.log")
