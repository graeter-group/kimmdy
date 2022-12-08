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

#%%
def set_dir():
    try:
        test_dir = Path(__file__).parent / "test_files/test_topology"
    except NameError:
        test_dir = Path("./tests/test_files/test_topology")
    os.chdir(test_dir)


set_dir()

#%%
hexala_top = read_topol(Path('hexala.top'))
ffdir = Path("../assets/amber99sb-star-ildnp.ff")
ffpatch = Path('amber99sb_patches.xml')
top = Topology(hexala_top, ffdir, ffpatch)
oldtop = deepcopy(top)

allowed_text = st.text(
    string.ascii_letters + "*+", min_size=1, max_size=5
)

@st.composite
def random_atomlist(draw):
    n = draw(st.integers(min_value=10, max_value=50))
    atom_nrs = [str(x + 1) for x in range(n)]
    bound_to = st.lists(st.integers(min_value=1, max_value=n), min_size=1, max_size=n)
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


def test_generate_topology_from_bound_to():
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
#     atoms = random_atomlist()
# )
# def test_break_bind_bond_invertible(atoms):
#     top = generate_topology_from_bound_to(atoms)
#     print(top)
#     assert False


