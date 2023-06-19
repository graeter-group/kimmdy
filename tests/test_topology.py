# %%
from copy import deepcopy
from pathlib import Path
import os
import pytest

from kimmdy.parsing import read_topol, write_topol,TopologyDict
from hypothesis import Phase, given, settings, strategies as st
from kimmdy.topology.topology import Topology, generate_topology_from_bound_to
from kimmdy.topology.atomic import *
from kimmdy.topology.utils import match_atomic_item_to_atomic_type
import logging

# %%

@pytest.fixture
def filedir() -> Path:
    dirname = "test_topology"
    try:
        file_dir = Path(__file__).parent / "test_files" / dirname
    except NameError:
        file_dir = Path("./tests/test_files" / dirname)
    return file_dir


@pytest.fixture
def assetsdir() -> Path:
    return Path(__file__).parent / "test_files" / "assets"


# %%


@st.composite
def random_atomlist(draw):
    n = draw(st.integers(min_value=3, max_value=5))
    elements = st.integers(min_value=1, max_value=n)
    bound_to = st.lists(elements, min_size=2, max_size=4, unique=True)
    atom_nrs = [str(x + 1) for x in range(n)]
    atoms = []
    allowed_text = st.text("COHT1" + "*+", min_size=1, max_size=5)
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
    try:
        assets_dir = Path(__file__).parent / "test_files" / "assets"
    except NameError:
        assets_dir = Path("./tests/test_files" / "assets")
    ffdir = assets_dir / "amber99sb-star-ildnp.ff"
    ffpatch = assets_dir / "amber99sb_patches.xml"
    atoms = draw(random_atomlist())
    top = generate_topology_from_bound_to(atoms, ffdir, ffpatch)
    break_this = draw(st.sampled_from(list(top.bonds.keys())))
    return (top, break_this)


class TestFFPatches:
    @pytest.fixture
    def top_fix(self, assetsdir, filedir) -> Topology:
        ffdir = assetsdir / "amber99sb-star-ildnp.ff"
        ffpatch = assetsdir / "amber99sb_patches.xml"
        hexala_top = read_topol(filedir / "hexala.top")
        return Topology(hexala_top, ffdir, ffpatch)

    def test_match_atomic_item_to_atomic_type(self, top_fix):
        types = top_fix.ff.angletypes

        atomic_id = ["CT", "C_R", "N"]
        want = ("CT", "C", "N")
        types_wanted = {want: types[want]}
        item_type = match_atomic_item_to_atomic_type(atomic_id, types_wanted)
        expected = AngleType(
            i="CT",
            j="C",
            k="N",
            id="CT---C---N",
            id_sym="N---C---CT",
            funct="1",
            c0="116.600",
            c1="585.760",
            c2=None,
            c3=None,
        )
        assert item_type == expected

        atomic_id = ["C_R", "CA", "HA"]
        want = ("C", "CA", "HA")
        types_wanted = {want: types[want]}
        item_type = match_atomic_item_to_atomic_type(atomic_id, types_wanted)
        expected = AngleType(
            i="C",
            j="CA",
            k="HA",
            id="C---CA---HA",
            id_sym="HA---CA---C",
            funct="1",
            c0="120.000",
            c1="418.400",
            c2=None,
            c3=None,
        )
        assert item_type == expected


class TestTopology:
    @pytest.fixture
    def top_fix(self, assetsdir, filedir) -> Topology:
        ffdir = assetsdir / "amber99sb-star-ildnp.ff"
        ffpatch = assetsdir / "amber99sb_patches.xml"
        hexala_top = read_topol(filedir / "hexala.top")
        return Topology(hexala_top, ffdir, ffpatch)

    def test_break_bind_bond_hexala(self, top_fix):
        top = deepcopy(top_fix)
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
    def test_break_bind_random_bond_hexala(self, top_fix, bondindex):
        top = deepcopy(top_fix)
        og_top = deepcopy(top)
        bond_key = list(top.bonds.keys())[bondindex]
        top.break_bond(bond_key)
        top.bind_bond(bond_key)
        assert top.bonds == og_top.bonds
        assert top.pairs == og_top.pairs
        assert top.angles == og_top.angles
        assert top.proper_dihedrals == og_top.proper_dihedrals
        assert top.improper_dihedrals == og_top.improper_dihedrals

    def test_generate_topology_from_bound_to(self, top_fix, assetsdir):
        ffdir = assetsdir / "amber99sb-star-ildnp.ff"
        ffpatch = assetsdir / "amber99sb_patches.xml"
        og_top = deepcopy(top_fix)
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
    @pytest.fixture
    def hexala_top(self, filedir) -> TopologyDict:
        return read_topol(filedir / "hexala.top")

    @pytest.fixture
    def top_fix(self, assetsdir, filedir) -> Topology:
        ffdir = assetsdir / "amber99sb-star-ildnp.ff"
        ffpatch = assetsdir / "amber99sb_patches.xml"
        hexala_top = read_topol(filedir / "hexala.top")
        return Topology(hexala_top, ffdir, ffpatch)

    @pytest.fixture
    def top_break_29_35_fix(self, assetsdir, filedir) -> Topology:
        ffdir = assetsdir / "amber99sb-star-ildnp.ff"
        ffpatch = assetsdir / "amber99sb_patches.xml"
        hexala_top = read_topol(filedir / "hexala_break29-35.top")
        return Topology(hexala_top, ffdir, ffpatch)

    @pytest.fixture
    def top_move_34_29_fix(self, assetsdir, filedir) -> Topology:
        ffdir = assetsdir / "amber99sb-star-ildnp.ff"
        ffpatch = assetsdir / "amber99sb_patches.xml"
        hexala_top = read_topol(filedir / "hexala_move34-29.top")
        return Topology(hexala_top, ffdir, ffpatch)

    def test_all_terms_accounted_for(self, top_fix, hexala_top):
        assert len(top_fix.atoms) == len(hexala_top["atoms"])
        assert len(top_fix.bonds) == len(hexala_top["bonds"])
        assert len(top_fix.pairs) == len(hexala_top["pairs"])
        assert len(top_fix.angles) == len(hexala_top["angles"])
        assert len(top_fix.proper_dihedrals) + len(top_fix.improper_dihedrals) == len(
            hexala_top["dihedrals"]
        )

    def test_find_bondtypes(self, top_fix):
        top_cp = deepcopy(top_fix)

        id = ["C", "CT"]
        result = match_atomic_item_to_atomic_type(id, top_cp.ff.bondtypes)
        id = ["CT", "C"]
        result = match_atomic_item_to_atomic_type(id, top_cp.ff.bondtypes)
        assert result is not None

    def test_break_bond_29_35(self, top_fix, top_break_29_35_fix):
        top = deepcopy(top_fix)
        top_broken = deepcopy(top_break_29_35_fix)
        top.break_bond(("29", "35"))
        assert len(top.bonds) == len(top_broken.bonds)
        assert len(top.pairs) == len(top_broken.pairs)
        assert len(top.angles) == len(top_broken.angles)
        assert len(top.proper_dihedrals) == len(top_broken.proper_dihedrals)
        assert len(top.improper_dihedrals) == len(top_broken.improper_dihedrals)

    def test_break_bond_9_15(self, top_fix):
        top = deepcopy(top_fix)
        og_top = deepcopy(top_fix)
        breakpair = ("9", "15")

        top.break_bond(breakpair)
        top._update_dict()

        topology = og_top.top
        topology_new = top.top

        bonddiff = set([(x[0], x[1]) for x in topology["bonds"]]) - set(
            [(x[0], x[1]) for x in topology_new["bonds"]]
        )
        pairdiff = set([tuple(x[:2]) for x in topology["pairs"]]) - set(
            [tuple(x[:2]) for x in topology_new["pairs"]]
        )
        anglediff = set([tuple(x[:3]) for x in topology["angles"]]) - set(
            [tuple(x[:3]) for x in topology_new["angles"]]
        )
        dihedraldiff = set([tuple(x[:4]) for x in topology["dihedrals"]]) - set(
            [tuple(x[:4]) for x in topology_new["dihedrals"]]
        )

        assert bonddiff == set([breakpair])
        assert pairdiff == set(
            [
                ("5", "15"),
                ("7", "16"),
                ("7", "17"),
                ("8", "15"),
                ("9", "18"),
                ("9", "19"),
                ("10", "16"),
                ("10", "17"),
                ("11", "16"),
                ("11", "17"),
                ("12", "15"),
                ("13", "15"),
                ("14", "15"),
            ]
        )
        assert anglediff == set(
            [
                ("7", "9", "15"),
                ("10", "9", "15"),
                ("11", "9", "15"),
                ("9", "15", "16"),
                ("9", "15", "17"),
            ]
        )
        assert dihedraldiff == set(
            [
                ("5", "7", "9", "15"),
                ("8", "7", "9", "15"),
                ("15", "9", "11", "12"),
                ("15", "9", "11", "13"),
                ("15", "9", "11", "14"),
                ("7", "9", "15", "16"),
                ("7", "9", "15", "17"),
                ("10", "9", "15", "16"),
                ("10", "9", "15", "17"),
                ("11", "9", "15", "16"),
                ("11", "9", "15", "17"),
                ("9", "15", "17", "18"),
                ("9", "15", "17", "19"),
                ("9", "17", "15", "16"),  # improper
            ]
        )

        assert len(bonddiff) == 1
        assert len(pairdiff) == 13
        assert len(anglediff) == 5
        assert len(dihedraldiff) == 14

    def test_top_properties(self, top_fix):
        top = deepcopy(top_fix)

        # initial correct number of atoms and bonds
        assert len(top.atoms) == 72
        assert len(top.bonds) == 71

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
            assert atom.type == "CT"
            assert atom.atom == "CH3"
            assert atom.residue == "ACE"
            break

        # everything is bound to something
        for atom in top.atoms.values():
            assert len(atom.bound_to_nrs) > 0

    def test_find_terms_around_atom(self, top_fix):
        top = deepcopy(top_fix)
        atomnr = "29"
        # 29       CT       4        ALA      CA       29       0.0337   12.01

        bonds = top._get_atom_bonds(atomnr)
        angles = top._get_atom_angles(atomnr)
        proper_dihedrals = top._get_atom_proper_dihedrals(atomnr)
        improper_dihedrals = top._get_atom_improper_dihedrals(atomnr)

        assert len(bonds) == 4
        assert len(angles) == 13
        assert len(proper_dihedrals) == 25
        assert len(improper_dihedrals) == 3

        atomnr = "9"
        bonds = top._get_atom_bonds(atomnr)
        angles = top._get_atom_angles(atomnr)
        proper_dihedrals_center = top._get_center_atom_dihedrals(atomnr)

        assert bonds == [("7", "9"), ("9", "10"), ("9", "11"), ("9", "15")]
        assert set(proper_dihedrals_center) == set(
            [
                ("5", "7", "9", "10"),
                ("5", "7", "9", "11"),
                ("5", "7", "9", "15"),
                ("8", "7", "9", "10"),
                ("8", "7", "9", "11"),
                ("8", "7", "9", "15"),
                ("7", "9", "11", "12"),
                ("7", "9", "11", "13"),
                ("7", "9", "11", "14"),
                ("10", "9", "11", "12"),
                ("10", "9", "11", "13"),
                ("10", "9", "11", "14"),
                ("15", "9", "11", "12"),
                ("15", "9", "11", "13"),
                ("15", "9", "11", "14"),
                ("7", "9", "15", "16"),
                ("7", "9", "15", "17"),
                ("10", "9", "15", "16"),
                ("10", "9", "15", "17"),
                ("11", "9", "15", "16"),
                ("11", "9", "15", "17"),
            ]
        )

    def test_move_34_29_after_break(self, top_fix, top_move_34_29_fix):
        """Move H at 34 to C at 29"""
        top = deepcopy(top_fix)
        top_moved = deepcopy(top_move_34_29_fix)
        top.break_bond(("29", "35"))
        top.break_bond(("31", "34"))
        top.bind_bond(("34", "29"))
        # patches are applied here, but parameters don't match, only number of parameters
        # the reference is shit anyway
        focus = set(["29", "31", "34", "35"])
        top.patch_parameters(list(focus))

        # compare topologies
        assert len(top.bonds) == len(top_moved.bonds)
        assert len(top.pairs) == len(top_moved.pairs)
        assert len(top.angles) == len(top_moved.angles)
        assert len(top.proper_dihedrals) == len(top_moved.proper_dihedrals)
        # improper dihedral patches not implemented yet!
        # assert len(top.improper_dihedrals) == len(top_moved.improper_dihedrals)

        # inspect HAT hydrogen
        h = top.atoms["34"]
        assert h.is_radical == False
        # keeps name
        assert h.atom == "HB3"
        # changes type
        assert h.type == "H1"
        assert h.residue == "ALA"
        assert h.bound_to_nrs == ["29"]

    def test_ff(self, top_fix):
        top = deepcopy(top_fix)

        residues = list(top.ff.residuetypes.keys())
        assert len(residues) == 121

        res = ResidueType(
            "HOH",
            atoms={
                "OW": ResidueAtomSpec("OW", "OW", "-0.834", "0"),
                "HW1": ResidueAtomSpec("HW1", "HW", "0.417", "0"),
                "HW2": ResidueAtomSpec("HW2", "HW", "0.417", "0"),
            },
            bonds={
                ("OW", "HW1"): ResidueBondSpec("OW", "HW1"),
                ("OW", "HW2"): ResidueBondSpec("OW", "HW2"),
            },
            proper_dihedrals={},
            improper_dihedrals={},
        )

        assert top.ff.residuetypes["HOH"] == res


class TestRadicalAla:
    @pytest.fixture
    def top_noprm_fix(self, assetsdir, filedir) -> Topology:
        ffdir = assetsdir / "amber99sb-star-ildnp.ff"
        ffpatch = assetsdir / "amber99sb_patches.xml"
        hexala_top = read_topol(filedir / "Ala_R_noprm.top")
        return Topology(hexala_top, ffdir, ffpatch)

    @pytest.fixture
    def top_prm_fix(self, assetsdir, filedir) -> Topology:
        ffdir = assetsdir / "amber99sb-star-ildnp.ff"
        ffpatch = assetsdir / "amber99sb_patches.xml"
        hexala_top = read_topol(filedir / "Ala_R_prm.top")
        return Topology(hexala_top, ffdir, ffpatch)

    def test_is_radical(self, top_noprm_fix):
        assert top_noprm_fix.atoms["9"].is_radical == True
        assert top_noprm_fix.atoms["10"].is_radical == False

    # use this test when parameter assignments from graph are working
    # def test_parameters_applied(self, top_noprm_fix, top_prm_fix):
    #     top = deepcopy(top_noprm_fix)
    #     focus = set(["9","10"])
    #     top.patch_parameters(list(focus))
    #     assert top_noprm_fix.atoms == top_prm_fix.atoms
    #     top_dict = top.to_dict()
    #     write_topol(top_dict,Path("/hits/fast/mbm/hartmaec/kimmdys/kimmdy/tests/test_files/test_topology/Ala_R_prm_curr.top"))
    #     assert top_noprm_fix.bonds == top_prm_fix.bonds
    #     assert top_noprm_fix.pairs == top_prm_fix.pairs
    #     assert top_noprm_fix.angles == top_prm_fix.angles
    #     assert top_noprm_fix.proper_dihedrals == top_prm_fix.proper_dihedrals
    #     assert top_noprm_fix.improper_dihedrals == top_prm_fix.improper_dihedrals


# {('9', '10'): Bond(ai='9', aj='10', funct='1', c0='0.14955', c1='259408.000000', c2=None, c3=None)} != {('9', '10'): Bond(ai='9', aj='10', funct='1', c0=None, c1=None, c2=None, c3=None)}
# {('10', '13'): Bond(ai='10', aj='13', funct='1', c0='0.10900', c1='284512.0', c2=None, c3=None)} != {('10', '13'): Bond(ai='10', aj='13', funct='1', c0=None, c1=None, c2=None, c3=None)}
# {('7', '9'): Bond(ai='7', aj='9', funct='1', c0='0.13600', c1='282001.600000', c2=None, c3=None)} != {('7', '9'): Bond(ai='7', aj='9', funct='1', c0=None, c1=None, c2=None, c3=None)}
# {('9', '14'): Bond(ai='9', aj='14', funct='1', c0='0.14916', c1='265265.600000', c2=None, c3=None)} != {('9', '14'): Bond(ai='9', aj='14', funct='1', c0=None, c1=None, c2=None, c3=None)}
# {('10', '11'): Bond(ai='10', aj='11', funct='1', c0='0.10900', c1='284512.0', c2=None, c3=None)} != {('10', '11'): Bond(ai='10', aj='11', funct='1', c0=None, c1=None, c2=None, c3=None)}
# {('10', '12'): Bond(ai='10', aj='12', funct='1', c0='0.10900', c1='284512.0', c2=None, c3=None)} != {('10', '12'): Bond(ai='10', aj='12', funct='1', c0=None, c1=None, c2=None, c3=None)}

# %%
