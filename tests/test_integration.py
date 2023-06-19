from kimmdy.parsing import (
    read_topol,
    read_plumed,
)
from kimmdy.cmd import kimmdy_run
from kimmdy.utils import run_shell_cmd

import os
import shutil
from pathlib import Path
import pytest

from kimmdy.topology.topology import Topology
from kimmdy.changemanager import break_bond_plumed


def check_succesful_run(logfile):
    with open(logfile, "r") as f:
        log = f.readlines()
    if len(log) > 1:
        if set(["Finished", "running", "tasks,"]).issubset(set(log[-1].split(sep=" "))):
            return
    raise ValueError("Log file indicates an unplanned exit of KIMMDY.")


def setup_testdir(tmp_path, dirname) -> Path:
    try:
        filedir = Path(__file__).parent / "test_files" / "test_integration" / dirname
        assetsdir = Path(__file__).parent / "test_files" / "assets"
    except NameError:
        filedir = Path("./tests/test_files" / "test_integration" / dirname)
        assetsdir = Path("./tests/test_files" / "assets")
    testdir = tmp_path / "test_integration" / dirname
    shutil.copytree(filedir, testdir)
    shutil.copy2(assetsdir / "amber99sb_patches.xml", testdir)
    Path(testdir / "amber99sb-star-ildnp.ff").symlink_to(
        assetsdir / "amber99sb-star-ildnp.ff",
        target_is_directory=True,
    )
    os.chdir(testdir.resolve())
    return testdir


def test_integration_emptyrun(tmp_path):
    testdir = setup_testdir(tmp_path, "emptyrun")

    (testdir / "emptyrun.txt").touch()
    run_shell_cmd("kimmdy")
    check_succesful_run(testdir / "kimmdy.log")


def test_integration_hat_reaction(tmp_path):
    testdir = setup_testdir(tmp_path, "HAT_reaction")

    (testdir / "emptyrun.txt").touch()
    run_shell_cmd("kimmdy")
    check_succesful_run(testdir / "kimmdy.log")


def test_integration_homolysis_reaction(tmp_path):
    tmpdir = tmp_path / "homolysis"
    shutil.copytree(
        Path(__file__).parent / "test_files/test_integration/homolysis", tmpdir
    )
    os.chdir(tmpdir)
    Path(tmpdir / "amber99sb-star-ildnp.ff").symlink_to(
        Path(__file__).parent / "test_files/assets/amber99sb-star-ildnp.ff",
        target_is_directory=True,
    )
    kimmdy_run(tmpdir / "kimmdy.yml")


def test_integration_whole_run(tmp_path):
    tmpdir = tmp_path / "whole_run"
    shutil.copytree(
        Path(__file__).parent / "test_files/test_integration/whole_run", tmpdir
    )
    os.chdir(tmpdir)
    Path(tmpdir / "amber99sb-star-ildnp.ff").symlink_to(
        Path(__file__).parent / "test_files/assets/amber99sb-star-ildnp.ff",
        target_is_directory=True,
    )
    kimmdy_run(tmpdir / "kimmdy.yml")
