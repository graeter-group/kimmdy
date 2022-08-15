from kimmdy.parsing import read_topol,read_plumed
from kimmdy.changemanager import break_bond_top,break_bond_plumed

from pathlib import Path
from copy import deepcopy
import pytest

def test_break_bond_top():
    input_f = Path(__file__).parent / "test_files/test_changemanager/hexala_out.top"
    topology = read_topol(input_f)
    breakpair = (9,15)
    breakpair = (str(breakpair[0]),str(breakpair[1]))
    
    topology_new = break_bond_top(deepcopy(topology),breakpair)
    
    diffdict = {}
    for key in topology.keys():
        oldset = set(tuple(x) for x in topology[key])
        newset = set(tuple(x) for x in topology_new[key])
        diffdict[key] = list(oldset - newset)

    assert len(diffdict['bonds']) == 1 and all([x in diffdict['bonds'][0] for x in breakpair])
    assert len(diffdict['pairs']) == 13
    assert len(diffdict['angles']) == 5 and all([x in angle for angle in diffdict['angles'] for x in breakpair])
    assert len(diffdict['dihedrals']) == 15 and all([x in dih for dih in diffdict['dihedrals'] for x in breakpair])  #includes impropers which might change

def test_break_bond_plumed():
    input_f = Path(__file__).parent / "test_files/test_changemanager/plumed.dat"
    plumeddat = read_plumed(input_f)
    breakpair = (9,15)
    breakpair = (str(breakpair[0]),str(breakpair[1]))

    newplumeddat = break_bond_plumed(deepcopy(plumeddat),breakpair,'distances.dat')

    oldset = set(tuple(x['atoms']) for x in plumeddat['distances'])
    newset = set(tuple(x['atoms']) for x in newplumeddat['distances'])
    diffs = list(oldset - newset)
    assert len(diffs) == 1 and diffs[0] == breakpair




