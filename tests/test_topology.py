#%%
from copy import deepcopy
from pathlib import Path
import os

from kimmdy.parsing import read_topol
from hypothesis import Phase, given, settings, strategies as st
from kimmdy.topology.topology import Topology, generate_topology_from_bound_to
from kimmdy.topology.atomic import *

#%%
def set_dir():
    try:
        test_dir = Path(__file__).parent / "test_files/test_topology"
    except NameError:
        test_dir = Path("./tests/test_files/test_topology")
    os.chdir(test_dir)


set_dir()

#%%
ffdir = Path("../assets/amber99sb-star-ildnp.ff")
ffpatch = Path("amber99sb_patches.xml")

allowed_text = st.text("COHT1" + "*+", min_size=1, max_size=5)


@st.composite
def random_atomlist(draw):
    n = draw(st.integers(min_value=3, max_value=5))
    elements = st.integers(min_value=1, max_value=n)
    bound_to = st.lists(elements, min_size=2, max_size=4, unique=True)
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
    top = generate_topology_from_bound_to(atoms, ffdir, ffpatch)
    break_this = draw(st.sampled_from(list(top.bonds.keys())))
    return (top, break_this)


class TestTopology:
    hexala_top = read_topol(Path("hexala.top"))
    ffdir = Path("../assets/amber99sb-star-ildnp.ff")
    ffpatch = Path("amber99sb_patches.xml")
    top = Topology(hexala_top, ffdir, ffpatch)

    def test_break_bind_bond_hexala(self):
        top = deepcopy(self.top)
        og_top = deepcopy(top)

        bondindex = 24
        bond_key = list(top.bonds.keys())[bondindex]
        # bond_key = ('9', '10')
        # bond = ('1', '2')
        # bond = ('5', '6')
        top.break_bond(bond_key)
        top.bind_bond(bond_key)
        assert top.bonds == og_top.bonds
        assert top.pairs == og_top.pairs
        assert top.angles == og_top.angles
        assert top.proper_dihedrals == og_top.proper_dihedrals
        assert top.improper_dihedrals == og_top.improper_dihedrals

    @given(bondindex=st.integers(min_value=0, max_value=70))
    def test_break_bind_random_bond_hexala(self, bondindex):
        top = deepcopy(self.top)
        og_top = deepcopy(top)
        bond_key = list(top.bonds.keys())[bondindex]
        top.break_bond(bond_key)
        top.bind_bond(bond_key)
        assert top.bonds == og_top.bonds
        assert top.pairs == og_top.pairs
        assert top.angles == og_top.angles
        assert top.proper_dihedrals == og_top.proper_dihedrals
        assert top.improper_dihedrals == og_top.improper_dihedrals

    def test_generate_topology_from_bound_to(self):
        og_top = deepcopy(self.top)
        atoms = list(og_top.atoms.values())
        newtop = generate_topology_from_bound_to(atoms, ffdir, ffpatch)
        assert newtop.bonds == og_top.bonds
        assert newtop.pairs == og_top.pairs
        assert newtop.angles == og_top.angles
        assert newtop.proper_dihedrals == og_top.proper_dihedrals

    @settings(max_examples=1, phases=[Phase.generate])
    @given(top_break=random_topology_and_break())
    def test_break_bind_bond_invertible(self, top_break):
        top, to_break = top_break
        og_top = deepcopy(top)
        top.break_bond(to_break)
        top.bind_bond(to_break)
        assert top.bonds == og_top.bonds
        assert top.pairs == og_top.pairs
        assert top.angles == og_top.angles
        assert top.proper_dihedrals == og_top.proper_dihedrals
