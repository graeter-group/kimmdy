from kimmdy.parsing import (
    read_topol,
    read_plumed,
)
from kimmdy.cmd import kimmdy_run

import os
import shutil
from pathlib import Path

from kimmdy.topology.topology import Topology


def test_integration_break_bond_top():
    ffdir = Path("../assets/amber99sb-star-ildnp.ff")
    ffpatch = Path("amber99sb_patches.xml")
    toppath = Path(__file__).parent / "test_files/test_integration/hexala_nat.top"
    toppath_compare = (
        Path(__file__).parent / "test_files/test_integration/hexala_break29-35.top"
    )

    top = Topology(read_topol(toppath), ffdir, ffpatch)
    top_broken = Topology(read_topol(toppath_compare), ffdir, ffpatch)

    pair = ("29", "35")
    top.break_bond(pair)

    assert top.bonds == top_broken.bonds
    assert top.pairs == top_broken.pairs
    assert top.angles == top_broken.angles
    assert top.proper_dihedrals == top_broken.proper_dihedrals
    assert top.improper_dihedrals == top_broken.improper_dihedrals


def test_integration_break_plumed():
    # FIXME
    pass


def test_integration_move_top():
    ffdir = Path("../assets/amber99sb-star-ildnp.ff")
    ffpatch = Path("amber99sb_patches.xml")
    toppath = (
        Path(__file__).parent / "test_files/test_integration/hexala_break29-35.top"
    )
    toppath_compare = (
        Path(__file__).parent / "test_files/test_integration/hexala_move34-29.top"
    )
    ffdir = Path(__file__).parent / "test_files/assets/amber99sb-star-ildnp.ff"

    top = Topology(read_topol(toppath), ffdir, ffpatch)
    top_moved = Topology(read_topol(toppath_compare), ffdir, ffpatch)

    break_bond = ("31", "34")
    bind_bond = ("34", "29")
    top.break_bond(break_bond)
    top.bind_bond(bind_bond)

    assert top.bonds == top_moved.bonds
    assert top.pairs == top_moved.pairs
    assert top.angles == top_moved.angles
    assert top.proper_dihedrals == top_moved.proper_dihedrals
    assert top.improper_dihedrals == top_moved.improper_dihedrals


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
