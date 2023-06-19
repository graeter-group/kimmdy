from kimmdy.parsing import (
    read_topol,
    read_plumed,
)
from kimmdy.cmd import kimmdy_run

import os
import shutil
from pathlib import Path
import pytest

from kimmdy.topology.topology import Topology
from kimmdy.changemanager import break_bond_plumed


@pytest.fixture(scope="function")
def testdir(tmp_path) -> Path:
    dirname = "test_integration"
    try:
        filedir = Path(__file__).parent / "test_files" / dirname
        assetsdir = Path(__file__).parent / "test_files" / "assets"
    except NameError:
        filedir = Path("./tests/test_files" / dirname)
        assetsdir = Path("./tests/test_files" / "assets")
    testdir = tmp_path / dirname
    shutil.copytree(filedir, testdir)
    shutil.copy2(assetsdir / "amber99sb_patches.xml", testdir)
    Path(testdir / "amber99sb-star-ildnp.ff").symlink_to(
        assetsdir / "amber99sb-star-ildnp.ff",
        target_is_directory=True,
    )
    return testdir


def test_integration_emptyrun(tmp_path):
    tmpdir = tmp_path / "emptyrun"
    shutil.copytree(
        Path(__file__).parent / "test_files/test_integration/emptyrun", tmpdir
    )
    os.chdir(tmpdir)
    Path(tmpdir / "emptyrun.txt").touch()
    Path(tmpdir / "amber99sb-star-ildnp.ff").symlink_to(
        Path(__file__).parent / "test_files/assets/amber99sb-star-ildnp.ff",
        target_is_directory=True,
    )
    kimmdy_run(tmpdir / "kimmdy_emptyrun.yml")


def test_integration_hat_reaction(tmp_path):
    tmpdir = tmp_path / "HAT_reaction"
    shutil.copytree(
        Path(__file__).parent / "test_files/test_integration/HAT_reaction", tmpdir
    )
    os.chdir(tmpdir)
    Path(tmpdir / "amber99sb-star-ildnp.ff").symlink_to(
        Path(__file__).parent / "test_files/assets/amber99sb-star-ildnp.ff",
        target_is_directory=True,
    )
    kimmdy_run(tmpdir / "kimmdy.yml")


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
