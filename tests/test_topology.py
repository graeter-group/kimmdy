# %%
from copy import deepcopy
from pathlib import Path
import os

from kimmdy.parsing import read_topol
from hypothesis import Phase, given, settings, strategies as st
from kimmdy.topology.topology import Topology, generate_topology_from_bound_to
from kimmdy.topology.atomic import *
from kimmdy.topology.utils import match_atomic_item_to_atomic_type
import logging


# %%
def set_dir():
    try:
        test_dir = Path(__file__).parent / "test_files/test_topology"
    except NameError:
        test_dir = Path("./tests/test_files/test_topology")
    os.chdir(test_dir)


set_dir()

# %%
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

class TestFFPatches:
    hexala_top = read_topol(Path("hexala.top"))
    top = Topology(hexala_top, ffdir, ffpatch)

    def test_match_atomic_item_to_atomic_type(self):
        types = self.top.ff.angletypes

        atomic_id = ['CT', 'C_R', 'N']
        want = ('CT', 'C', 'N')
        types_wanted = {want: types[want]}
        item_type = match_atomic_item_to_atomic_type(atomic_id, types_wanted)
        expected = AngleType(i='CT', j='C', k='N', id='CT---C---N', id_sym='N---C---CT', funct='1', c0='116.600', c1='585.760', c2=None, c3=None)
        assert item_type == expected

        atomic_id = ['C_R', 'CA', 'HA']
        want = ('C', 'CA', 'HA')
        types_wanted = {want: types[want]}
        item_type = match_atomic_item_to_atomic_type(atomic_id, types_wanted)
        expected = AngleType(i='C', j='CA', k='HA', id='C---CA---HA', id_sym='HA---CA---C', funct='1', c0='120.000', c1='418.400', c2=None, c3=None)
        assert item_type == expected

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
        logging.info(f"bond_key: {bond_key}")
        # 25, 27
        # C, N
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



class TestHexalaTopology:
    hexala_top = read_topol(Path("hexala.top"))
    hexala_break_29_35 = read_topol(Path("hexala_break29-35.top"))
    hexala_move_34_39 = read_topol(Path("hexala_move34-29.top"))
    ffdir = Path("../assets/amber99sb-star-ildnp.ff")
    ffpatch = Path("amber99sb_patches.xml")
    top = Topology(hexala_top, ffdir, ffpatch)
    top_break_29_35 = Topology(hexala_break_29_35, ffdir, ffpatch)
    top_move_34_29 = Topology(hexala_move_34_39, ffdir, ffpatch)

    def all_terms_accounted_for(self):
        top = self.top
        hexala_top = self.hexala_top
        assert len(top.atoms) == len(hexala_top["atoms"])
        assert len(top.bonds) == len(hexala_top["bonds"])
        assert len(top.pairs) == len(hexala_top["pairs"])
        assert len(top.angles) == len(hexala_top["angles"])
        assert len(top.proper_dihedrals) == len(hexala_top["propers"])
        assert len(top.improper_dihedrals) == len(hexala_top["impropers"]) 

    def break_bond_29_35(self):
        top = deepcopy(self.top)
        top_broken = deepcopy(self.top_break_29_35)
        top.break_bond(('29', '35'))
        assert top.bonds == top_broken.bonds
        assert top.pairs == top_broken.pairs
        assert top.angles == top_broken.angles
        assert top.proper_dihedrals == top_broken.proper_dihedrals
        assert top.improper_dihedrals == top_broken.improper_dihedrals

    def test_break_bond_9_15(self):
        top = deepcopy(self.top)
        og_top = deepcopy(self.top)
        breakpair = ('9', '15')

        top.break_bond(breakpair)
        top._update_dict()

        topology = og_top.top
        topology_new = top.top

        diffdict = {}
        for key in topology.keys():
            oldset = set(tuple(x) for x in topology[key])
            newset = set(tuple(x) for x in topology_new[key])
            diffdict[key] = list(oldset - newset)

        assert len(diffdict["bonds"]) == 1 and all(
            [x in diffdict["bonds"][0] for x in breakpair]
        )
        assert len(diffdict["pairs"]) == 13
        assert len(diffdict["angles"]) == 5 and all(
            [x in angle for angle in diffdict["angles"] for x in breakpair]
        )
        assert len(diffdict["dihedrals"]) == 15 and all(
            [x in dih for dih in diffdict["dihedrals"] for x in breakpair]
        )
        # includes impropers which might change

    def test_top_properties(self):
        top = deepcopy(self.top)
        og_top = deepcopy(self.top)

        # initial correct number of atoms and bonds
        assert len(top.atoms) == 16
        assert len(top.bonds) == 15

        # order is correct
        val = 0
        for atom_idx in top.atoms.keys():
            assert int(atom_idx) > val
            val = int(atom_idx)

        val = 0
        for bond in top.bonds.keys():
            assert int(bond[0]) >= val
            val = int(bond[0])

        # specific atoms
        for atom in top.atoms.values():
            assert atom.type == "CH3"
            assert atom.residue == "ACE"
            break

        # everything is bound to something
        for atom in top.atoms.values():
            assert len(atom.bound_to_nrs) > 0



    def test_move_34_29_after_break(self):
        """Move H at 34 to C at 39
        """
        top = deepcopy(self.top)
        top_moved = deepcopy(self.top_move_34_29)
        top.break_bond(('29', '25'))
        top.break_bond(('31', '34'))
        top.bind_bond(('34', '29'))

        # compare topologies
        assert top.bonds == top_moved.bonds
        assert top.pairs == top_moved.pairs
        assert top.angles == top_moved.angles
        assert top.proper_dihedrals == top_moved.proper_dihedrals
        assert top.improper_dihedrals == top_moved.improper_dihedrals

        # inspect HAT hydrogen
        h = top.atoms['34']
        assert h.is_radical == False
        assert h.atom == "H9"
        assert h.type == "HX"
        assert h.residue == "ALA"
        assert h.bound_to_nrs == ['29']

    def test_ff(self):
        top = deepcopy(self.top)

        residues = list(top.ff.residuetypes.keys())
        assert len(residues) == 125

        res = ResidueType('HOH',
                          atoms={
                              'OW': ResidueAtomSpec('OW', 'OW', '-0.834', '0'),
                              'HW1': ResidueAtomSpec('HW1', 'HW', '0.417', '0'),
                              'HW2': ResidueAtomSpec('HW2', 'HW', '0.417', '0')
                          },
                          bonds={
                              ('OW', 'HW1'): ResidueBondSpec('OW', 'HW1'),
                              ('OW', 'HW2'): ResidueBondSpec('OW', 'HW2')
                          },
                          improper_dihedrals={}
                          )

        assert residues[0] == res




