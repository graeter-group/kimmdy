# %%
from copy import deepcopy
from pathlib import Path
import pytest

from kimmdy.parsing import read_top, TopologyDict
from hypothesis import Phase, given, settings, HealthCheck, strategies as st
from kimmdy.topology.topology import Topology
from kimmdy.topology.atomic import *
from kimmdy.topology.utils import (
    get_top_section,
    match_atomic_item_to_atomic_type,
    get_protein_section,
)
import logging


@pytest.fixture(scope="module")
def filedir() -> Path:
    dirname = "test_topology"
    try:
        file_dir = Path(__file__).parent / "test_files" / dirname
    except NameError:
        file_dir = Path("./tests/test_files") / dirname
    return file_dir


@pytest.fixture(scope="module")
def assetsdir() -> Path:
    return Path(__file__).parent / "test_files" / "assets"


@pytest.fixture()
def raw_hexala_top_fix(filedir) -> TopologyDict:
    return read_top(filedir / "hexala.top")


@pytest.fixture()
def raw_top_a_fix(filedir) -> TopologyDict:
    return read_top(filedir / "topol_stateA.top")


@pytest.fixture()
def raw_top_b_fix(filedir) -> TopologyDict:
    return read_top(filedir / "topol_stateB.top")


@pytest.fixture()
def raw_urea_top_fix(filedir) -> TopologyDict:
    return read_top(filedir / "urea.top")


@pytest.fixture()
def hexala_top_fix(assetsdir, filedir) -> Topology:
    ffpatch = assetsdir / "amber99sb_patches.xml"
    hexala_top = read_top(filedir / "hexala.top")
    return Topology(hexala_top, ffpatch)


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
        dir = Path(__file__).parent / "test_files" / "test_topology"
    except NameError:
        dir = Path("./tests/test_files") / "test_topology"
    ffdir = dir / "amber99sb-star-ildnp.ff"
    ffpatch = dir / "amber99sb_patches.xml"
    hexala_top = read_top(dir / "hexala.top")
    top = Topology(hexala_top, ffpatch)
    atomlist = draw(random_atomlist())
    top.atoms = {atom.nr: atom for atom in atomlist}
    top._regenerate_topology_from_bound_to()
    break_this = draw(st.sampled_from(list(top.bonds.keys())))
    return (top, break_this)


class TestFFPatches:
    @pytest.fixture
    def top_fix(self, assetsdir, filedir) -> Topology:
        ffpatch = assetsdir / "amber99sb_patches.xml"
        hexala_top = read_top(filedir / "hexala.top")
        return Topology(hexala_top, ffpatch)

    def test_match_atomic_item_to_atomic_type(self, top_fix):
        types = top_fix.ff.angletypes

        atomic_id = ["CT", "C_R", "N"]
        want = ("CT", "C", "N")
        types_wanted: dict[AngleId, AngleType] = {want: types[want]}
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


class TestUrea:
    def test_urea(self, raw_urea_top_fix):
        raw = deepcopy(raw_urea_top_fix)
        top = Topology(raw)
        assert len(top.atoms) == 8
        assert len(top.bonds) == 7
        assert len(top.pairs) == 0
        assert len(top.angles) == 0
        assert len(top.proper_dihedrals) == 8
        assert len(top.improper_dihedrals) == 3

    def test_reindex_atomnumbers_for_already_ordered(self, raw_urea_top_fix):
        raw = deepcopy(raw_urea_top_fix)
        top = Topology(raw)
        og_top = deepcopy(top)
        top.reindex_atomnrs()
        assert top == og_top

    def test_reindex_atomnumbers_after_deletion(self, raw_urea_top_fix):
        raw = deepcopy(raw_urea_top_fix)
        top = Topology(raw)
        og_top = deepcopy(top)
        top.atoms.pop('3')
        top.reindex_atomnrs()

        og_atoms = list(og_top.atoms.keys())
        atoms = list(top.atoms.keys())

        og_bonds = list(og_top.bonds.keys())
        bonds = list(top.bonds.keys())

        og_dihedrals = list(og_top.proper_dihedrals.keys())
        dihedrals = list(top.proper_dihedrals.keys())

        og_impropers = list(og_top.improper_dihedrals.keys())
        impropers = list(top.improper_dihedrals.keys())

        assert og_atoms == [str(x) for x in range(1,9)]
        assert atoms == [str(x) for x in range(1,8)]

        assert top.atoms['3'] == Atom('3', 'H', '1', 'URE', 'H11', '4', '0.395055', '1.00800')

        assert og_bonds == [('1', '2'), ('1', '3'), ('1', '6'), ('3', '4'), ('3', '5'), ('6', '7'), ('6', '8')]
        assert bonds == [('1', '2'), ('1', '5'), ('5', '6'), ('5', '7')]

        assert og_dihedrals == [('2', '1', '3', '4'), ('2', '1', '3', '5'), ('2', '1', '6', '7'), ('2', '1', '6', '8'), ('3', '1', '6', '7'), ('3', '1', '6', '8'), ('6', '1', '3', '4'), ('6', '1', '3', '5')]
        assert dihedrals == [('2', '1', '5', '6'), ('2', '1', '5', '7')]

        assert og_impropers == []
        assert impropers == []


class TestTopAB:
    def test_top_ab(self, raw_top_a_fix, raw_top_b_fix):
        topA = Topology(raw_top_a_fix)
        topB = Topology(raw_top_b_fix)
        assert topA
        assert topB
        assert len(topA.atoms) == 41
        assert len(topA.proper_dihedrals) == 84
        proper_dihedrals_counts = 0
        for dihedral in topA.proper_dihedrals.values():
            proper_dihedrals_counts += len(dihedral.dihedrals)
        assert proper_dihedrals_counts == 88


class TestTopology:
    def test_break_bind_bond_hexala(self, hexala_top_fix):
        top = deepcopy(hexala_top_fix)
        og_top = deepcopy(top)

        bondindex = 24
        bond_key = list(top.bonds.keys())[bondindex]
        logging.info(f"bond_key: {bond_key}")
        top.break_bond(bond_key)
        top.bind_bond(bond_key)
        assert top.bonds == og_top.bonds
        assert top.pairs == og_top.pairs
        assert top.angles == og_top.angles
        assert top.proper_dihedrals == og_top.proper_dihedrals
        assert top.improper_dihedrals == og_top.improper_dihedrals

    @given(bondindex=st.integers(min_value=0, max_value=70))
    @settings(suppress_health_check=[HealthCheck.function_scoped_fixture])
    def test_break_bind_random_bond_hexala(self, hexala_top_fix, bondindex):
        top = deepcopy(hexala_top_fix)
        og_top = deepcopy(top)
        bond_key = list(top.bonds.keys())[bondindex]
        top.break_bond(bond_key)
        top.bind_bond(bond_key)
        assert top.bonds == og_top.bonds
        assert top.pairs == og_top.pairs
        assert top.angles == og_top.angles
        assert top.proper_dihedrals == og_top.proper_dihedrals
        assert top.improper_dihedrals == og_top.improper_dihedrals

    def test_generate_topology_from_bound_to(self, hexala_top_fix):
        og_top = deepcopy(hexala_top_fix)
        newtop = deepcopy(hexala_top_fix)
        newtop._regenerate_topology_from_bound_to()
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
    def top_break_29_35_fix(self, assetsdir, filedir) -> Topology:
        ffpatch = assetsdir / "amber99sb_patches.xml"
        hexala_top = read_top(filedir / "hexala_break29-35.top")
        return Topology(hexala_top, ffpatch)

    @pytest.fixture
    def top_move_34_29_fix(self, assetsdir, filedir) -> Topology:
        ffpatch = assetsdir / "amber99sb_patches.xml"
        hexala_top = read_top(filedir / "hexala_move34-29.top")
        return Topology(hexala_top, ffpatch)

    def test_all_terms_accounted_for(self, raw_hexala_top_fix, hexala_top_fix):
        atoms = get_protein_section(raw_hexala_top_fix, "atoms")
        bonds = get_protein_section(raw_hexala_top_fix, "bonds")
        pairs = get_protein_section(raw_hexala_top_fix, "pairs")
        angles = get_protein_section(raw_hexala_top_fix, "angles")
        dihedrals = get_protein_section(raw_hexala_top_fix, "dihedrals")

        assert atoms
        assert bonds
        assert pairs
        assert angles
        assert dihedrals
        assert len(hexala_top_fix.atoms) == len(atoms)
        assert len(hexala_top_fix.bonds) == len(bonds)
        assert len(hexala_top_fix.pairs) == len(pairs)
        assert len(hexala_top_fix.angles) == len(angles)
        assert len(hexala_top_fix.proper_dihedrals) + len(
            hexala_top_fix.improper_dihedrals
        ) == len(dihedrals)

    def test_find_bondtypes(self, hexala_top_fix):
        top_cp = deepcopy(hexala_top_fix)

        id = ["C", "CT"]
        result = match_atomic_item_to_atomic_type(id, top_cp.ff.bondtypes)
        id = ["CT", "C"]
        result = match_atomic_item_to_atomic_type(id, top_cp.ff.bondtypes)
        assert result is not None

    def test_break_bond_29_35(self, hexala_top_fix, top_break_29_35_fix):
        top = deepcopy(hexala_top_fix)
        top_broken = deepcopy(top_break_29_35_fix)
        top.break_bond(("29", "35"))

        assert top.bonds == top_broken.bonds
        assert top.pairs == top_broken.pairs
        assert top.angles == top_broken.angles
        assert top.proper_dihedrals == top_broken.proper_dihedrals
        assert top.improper_dihedrals == top_broken.improper_dihedrals

    def test_break_bond_9_15(self, hexala_top_fix):
        top = deepcopy(hexala_top_fix)
        og_top = deepcopy(hexala_top_fix)
        breakpair = ("9", "15")

        top.break_bond(breakpair)
        top._update_dict()

        topology = og_top.top
        topology_new = top.top

        atoms = get_protein_section(topology, "atoms")
        bonds = get_protein_section(topology, "bonds")
        pairs = get_protein_section(topology, "pairs")
        angles = get_protein_section(topology, "angles")
        dihedrals = get_protein_section(topology, "dihedrals")

        atoms_new = get_protein_section(topology_new, "atoms")
        bonds_new = get_protein_section(topology_new, "bonds")
        pairs_new = get_protein_section(topology_new, "pairs")
        angles_new = get_protein_section(topology_new, "angles")
        dihedrals_new = get_protein_section(topology_new, "dihedrals")
        assert atoms is not None
        assert bonds is not None
        assert pairs is not None
        assert angles is not None
        assert dihedrals is not None
        assert atoms_new is not None
        assert bonds_new is not None
        assert pairs_new is not None
        assert angles_new is not None
        assert dihedrals_new is not None

        bonddiff = set([(x[0], x[1]) for x in bonds]) - set(
            [(x[0], x[1]) for x in bonds_new]
        )
        pairdiff = set([tuple(x[:2]) for x in pairs]) - set(
            [tuple(x[:2]) for x in pairs_new]
        )
        anglediff = set([tuple(x[:3]) for x in angles]) - set(
            [tuple(x[:3]) for x in angles_new]
        )
        dihedraldiff = set([tuple(x[:4]) for x in dihedrals]) - set(
            [tuple(x[:4]) for x in dihedrals_new]
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

    def test_top_update_dict(self, raw_hexala_top_fix):
        raw = raw_hexala_top_fix
        raw_copy = deepcopy(raw)
        top = Topology(raw_copy)
        top._update_dict()
        assert top.top["dihedraltypes"]["content"] == raw["dihedraltypes"]["content"]
        assert (
            top.top["moleculetype_0"]["subsections"]["dihedrals"]["content"]
            == raw["moleculetype_0"]["subsections"]["dihedrals"]["content"]
        )
        assert top.top == raw

    def test_top_properties(self, hexala_top_fix):
        top = deepcopy(hexala_top_fix)

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

    def test_find_terms_around_atom(self, hexala_top_fix):
        top = deepcopy(hexala_top_fix)
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

    def test_move_34_29_after_break(self, hexala_top_fix, top_move_34_29_fix):
        """Move H at 34 to C at 29"""
        top = deepcopy(hexala_top_fix)
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

    def test_ff(self, hexala_top_fix):
        top = deepcopy(hexala_top_fix)

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
        ffpatch = assetsdir / "amber99sb_patches.xml"
        hexala_top = read_top(filedir / "Ala_R_noprm.top")
        return Topology(hexala_top, ffpatch)

    @pytest.fixture
    def top_prm_fix(self, assetsdir, filedir) -> Topology:
        ffpatch = assetsdir / "amber99sb_patches.xml"
        hexala_top = read_top(filedir / "Ala_R_prm.top")
        return Topology(hexala_top, ffpatch)

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
