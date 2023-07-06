from grappa_interface import generate_input, apply_parameters, GrappaInterface
from kimmdy.topology.topology import Topology
import pytest

def test_generate_input(generic_topology):
    AA3  = ['ALA','CYS','ASP','GLU','PHE','GLY','HIE','ILE','LYS','LEU','MET','ASN','PRO','HYP','GLN','ARG','SER','THR','VAL','TRP','TYR','DOP','ACE','NME']

    atoms, bonds, radicals = generate_input(generic_topology)
    
    assert len(atoms) == 72
    assert len(bonds) == 70
    assert radicals == [29, 35]
    assert atoms[0] == ['ACE', 'CT', ['3.39967e-01', '4.57730e-01'], '6']

    assert all([len(x) == 4 for x in atoms])
    assert all([isinstance(x[y],str) for x in atoms for y in [0,1,3]])
    assert all([isinstance(x[2],list) for x in atoms])
    assert all([len(x[2]) == 2 for x in atoms])
    assert all([x[0] in AA3 for x in atoms])

def test_apply_parameters_assertE(generic_topology, caplog):
    parameters = {k:{('100'):[]} for k in ['atoms','bonds','angles','proper_dihedrals','improper_dihedrals']}

    apply_parameters(generic_topology,parameters)

    warnings = 0
    for record in caplog.records:
        if record.levelname == 'WARNING':
            warnings += 1
    assert warnings == 5


def test_apply_parameters(generic_topology):
    idxs = {'atoms':('1'),'bonds':('1','3'),'angles':('2','1','5'),'proper_dihedrals':('68','67','69','72')}
    vals = {'atoms':['0.1'],'bonds':['0.3','400.0'],'angles':['40.0','1000.0'],'proper_dihedrals':{'2':['0.0','20.0']}}
    parameters = {'atoms':{idxs['atoms']:vals['atoms']},'bonds':{idxs['bonds']:vals['bonds']},'angles':{idxs['angles']:vals['angles']},'proper_dihedrals':{idxs['proper_dihedrals']:vals['proper_dihedrals']},'improper_dihedrals':{}}

    apply_parameters(generic_topology,parameters)

    assert isinstance(generic_topology, Topology)
    assert generic_topology.atoms[idxs['atoms']].charge == vals['atoms'][0]
    assert generic_topology.bonds[idxs['bonds']].c0 == vals['bonds'][0]
    assert generic_topology.angles[idxs['angles']].c1 == vals['angles'][1]
    assert generic_topology.proper_dihedrals[idxs['proper_dihedrals']].dihedrals['2'].c0 == vals['proper_dihedrals']['2'][0]

def test_parameterize_topology(generic_topology):
    raise NotImplementedError