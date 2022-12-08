#%%
from copy import deepcopy
from itertools import takewhile, permutations
from pathlib import Path
import os
from xml.etree.ElementTree import Element
from kimmdy.parsing import is_not_comment, read_topol, read_xml_ff
from kimmdy.changemanager import LocalGraph
from hypothesis import given, strategies as st
from kimmdy.topology import FF, Angle, Topology, Atom, Bond, Dihedral, Pair, generate_topology_from_bound_to, get_by_permutations, get_element_id
import string
import logging

#%%
def set_dir():
    try:
        test_dir = Path(__file__).parent / "test_files/test_topology"
    except NameError:
        test_dir = Path("./tests/test_files/test_topology")
    os.chdir(test_dir)


set_dir()

#%%
allowed_text = st.text(
    "COHT1" + "*+", min_size=1, max_size=5
)

@st.composite
def random_atomlist(draw):
    n = draw(st.integers(min_value=5, max_value=10))
    atom_nrs = [str(x + 1) for x in range(n)]
    atoms = []
    for i in atom_nrs:
        type = draw(allowed_text)
        resnr = draw(allowed_text)
        residue = draw(allowed_text)
        a = draw(allowed_text)
        cgnr = draw(allowed_text)
        charge = draw(allowed_text)
        mass = draw(allowed_text)
        atom = Atom(str(i), type, resnr, residue, a, cgnr, charge, mass)
        atom.bound_to_nrs = [str(x) for x in draw(bound_to) if str(x) != i]
        atoms.append(atom)
    
    return atoms

@st.composite
def random_topology_and_break(draw):
    atoms = draw(random_atomlist())
    top = generate_topology_from_bound_to(atoms)
    break_this = draw(st.sampled_from(list(top.bonds.values())))
    return (top, (break_this.ai, break_this.aj))

class TestTopology:
    hexala_top = read_topol(Path('hexala.top'))
    ffdir = Path("../assets/amber99sb-star-ildnp.ff")
    ffpatch = Path('amber99sb_patches.xml')
    top = Topology(hexala_top, ffdir, ffpatch)
    oldtop = deepcopy(top)


    def test_break_bind_bond_hexala(self):
        top = Topology(self.hexala_top, self.ffdir, self.ffpatch)
        og_top = deepcopy(top)
        # bond = ('9', '10')
        bond = ('1', '2')
        top.break_bond(bond)
        top.bind_bond(bond)
        assert top.bonds == og_top.bonds
        assert top.pairs == og_top.pairs
        assert top.angles == og_top.angles
        assert top.proper_dihedrals == og_top.proper_dihedrals

    @given(bondindex = st.integers(min_value=0, max_value=70))
    def test_break_bind_random_bond_hexala(self, bondindex):
        top = Topology(self.hexala_top, self.ffdir, self.ffpatch)
        og_top = deepcopy(top)
        bond = list(top.bonds.keys())[bondindex]
        top.break_bond(bond)
        top.bind_bond(bond)
        assert top.bonds == og_top.bonds
        assert top.pairs == og_top.pairs
        assert top.angles == og_top.angles
        assert top.proper_dihedrals == og_top.proper_dihedrals


    def test_generate_topology_from_bound_to(self):
        hexala_top = read_topol(Path('hexala.top'))
        ffdir = Path("../assets/amber99sb-star-ildnp.ff")
        ffpatch = Path('amber99sb_patches.xml')
        og_top = Topology(hexala_top, ffdir, ffpatch)

        atoms = list(og_top.atoms.values())
        newtop = generate_topology_from_bound_to(atoms)
        assert newtop.bonds == og_top.bonds
        assert newtop.pairs == og_top.pairs
        assert newtop.angles == og_top.angles
        assert newtop.proper_dihedrals == og_top.proper_dihedrals


    # @given(
    #     top_break = random_topology_and_break()
    # )
    # def test_break_bind_bond_invertible(self, top_break):
    #     top, to_break = top_break
    #     og_top = deepcopy(top)
    #     top.break_bond(to_break)
    #     top.bind_bond(to_break)
    #     assert top.bonds == og_top.bonds
    #     assert top.pairs == og_top.pairs
    #     assert top.angles == og_top.angles
    #     assert top.proper_dihedrals == og_top.proper_dihedrals


