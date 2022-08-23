from kimmdy.parsing import read_topol,read_plumed, topol_split_dihedrals, topol_merge_propers_impropers, write_topol
from kimmdy.changemanager import break_bond_top, break_bond_plumed, move_bond_top
from kimmdy.reaction import ConversionRecipe, ConversionType

import pytest
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

