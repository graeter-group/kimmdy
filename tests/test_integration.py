from kimmdy.cmd import kimmdy_run

import subprocess as sp
import os
import shutil
import logging
from pathlib import Path
from kimmdy.topology.topology import Topology

import pytest

from kimmdy.parsing import read_top, write_top


def setup_testdir(tmp_path, dirname) -> Path:
    try:
        filedir = Path(__file__).parent / "test_files" / "test_integration" / dirname
        assetsdir = Path(__file__).parent / "test_files" / "assets"
    except NameError:
        filedir = Path("./tests/test_files") / "test_integration" / dirname
        assetsdir = Path("./tests/test_files") / "assets"
    testdir = tmp_path / "test_integration" / dirname
    shutil.copytree(filedir, testdir)
    shutil.copy2(assetsdir / "amber99sb_patches.xml", testdir)
    Path(testdir / "amber99sb-star-ildnp.ff").symlink_to(
        assetsdir / "amber99sb-star-ildnp.ff",
        target_is_directory=True,
    )
    os.chdir(testdir.resolve())
    return testdir


def test_integration_emptyrun(tmp_path, caplog):
    testdir = setup_testdir(tmp_path, "emptyrun")
    caplog.set_level(logging.INFO)
    (testdir / "emptyrun.txt").touch()

    # not expecting this to run
    # because the topology is empty
    with pytest.raises(ValueError) as e:
        kimmdy_run()


def test_integration_valid_input_files(tmp_path, caplog):
    testdir = setup_testdir(tmp_path, "minimal_input_files")
    caplog.set_level(logging.INFO)
    (testdir / "emptyrun.txt").touch()

    kimmdy_run()
    for record in caplog.records:
        # assert record.levelname != "WARNING"
        assert record.levelname != "CRITICAL"
    assert set(["Finished", "running", "tasks,"]).issubset(
        set(caplog.records[-1].message.split(sep=" "))
    )


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


def test_integration_hat_reaction(tmp_path, caplog):
    testdir = setup_testdir(tmp_path, "hat_naive")
    caplog.set_level(logging.INFO)

    kimmdy_run()
    for record in caplog.records:
        # assert record.levelname != "WARNING"
        assert record.levelname != "CRITICAL"
    assert set(["Finished", "running", "tasks,"]).issubset(
        set(caplog.records[-1].message.split(sep=" "))
    )


def test_integration_homolysis_reaction(tmp_path, caplog):
    testdir = setup_testdir(tmp_path, "homolysis")
    caplog.set_level(logging.INFO)

    kimmdy_run()

    for record in caplog.records:
        # assert record.levelname != "WARNING"
        assert record.levelname != "CRITICAL"
    assert set(["Finished", "running", "tasks,"]).issubset(
        set(caplog.records[-1].message.split(sep=" "))
    )


def test_integration_whole_run(tmp_path, caplog):
    testdir = setup_testdir(tmp_path, "whole_run")
    caplog.set_level(logging.INFO)

    kimmdy_run()

    for record in caplog.records:
        # assert record.levelname != "WARNING"
        assert record.levelname != "CRITICAL"
    assert set(["Finished", "running", "tasks,"]).issubset(
        set(caplog.records[-1].message.split(sep=" "))
    )
