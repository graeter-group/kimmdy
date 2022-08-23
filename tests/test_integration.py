from kimmdy.parsing import read_topol,read_plumed, topol_split_dihedrals, topol_merge_propers_impropers, write_topol
from kimmdy.changemanager import break_bond_top, break_bond_plumed, move_bond_top
from kimmdy.reaction import ConversionRecipe, ConversionType
from kimmdy.cmd import kimmdy_run

import pytest
import os
import shutil
from pathlib import Path

def test_integration_break_bond_top():
    toppath = Path(__file__).parent / "test_files/test_integration/hexala_nat.top"
    toppath_compare = Path(__file__).parent / "test_files/test_integration/hexala_break29-35.top"

    topology = read_topol(toppath)
    topology_compare = read_topol(toppath_compare)

    pair = (29,35)

    topology = break_bond_top(topology, pair)

    for key in ['bonds','pairs','angles','dihedrals']:
        assert topology[key] == topology_compare[key]

def test_integration_break_plumed():
    plumedpath = Path(__file__).parent / "test_files/test_integration/plumed_nat.dat"
    plumedpath_compare = Path(__file__).parent / "test_files/test_integration/plumed_break29-35.dat"

    plumed = read_plumed(plumedpath)
    plumed_compare = read_plumed(plumedpath_compare)

    pair = (29,35)

    plumed = break_bond_plumed(plumed,pair,Path("distances.dat"))

    for key in ["distances","prints"]:
        assert plumed[key] == plumed_compare[key]
    pass

def test_integration_move_top():
    toppath = Path(__file__).parent / "test_files/test_integration/hexala_break29-35.top"
    toppath_compare = Path(__file__).parent / "test_files/test_integration/hexala_move34-29.top"
    ffdir = Path(__file__).parent / "test_files/test_integration/amber99sb-star-ildnp.ff"

    topology = read_topol(toppath)
    topology_compare = read_topol(toppath_compare)

    pair = (34,29)

    topology = move_bond_top(topology, pair,ffdir)

    for key in ['bonds','pairs','angles','dihedrals']:
        assert topology[key] == topology_compare[key]

def test_integration_emptyrun(tmp_path):   
    tmpdir = tmp_path / "emptyrun"
    shutil.copytree(Path(__file__).parent / "test_files/test_integration/emptyrun", tmpdir)
    os.chdir(tmpdir)
    Path(tmpdir / "emptyrun.txt").touch()
    kimmdy_run(tmpdir / "kimmdy_emptyrun.yml")

def test_integration_hat_reaction(tmp_path):
    tmpdir = tmp_path / "HAT_reaction"
    shutil.copytree(Path(__file__).parent / "test_files/test_integration/HAT_reaction", tmpdir)
    os.chdir(tmpdir)
    kimmdy_run(tmpdir / "kimmdy.yml")

def test_integration_whole_run(tmp_path):
    tmpdir = tmp_path / "whole_run"
    shutil.copytree(Path(__file__).parent / "test_files/test_integration/whole_run", tmpdir)
    os.chdir(tmpdir)
    kimmdy_run(tmpdir / "kimmdy.yml")


