import logging
from copy import deepcopy
from pathlib import Path

import numpy as np
import pytest
from hypothesis import HealthCheck, Phase, given, settings
from hypothesis import strategies as st

from kimmdy.parsing import TopologyDict, read_top, write_top
from kimmdy.recipe import Bind, Break
from kimmdy.topology.atomic import *
from kimmdy.topology.topology import REACTIVE_MOLECULEYPE, Topology
from kimmdy.topology.utils import (
    get_protein_section,
    get_reactive_section,
    get_residue_by_bonding,
    match_atomic_item_to_atomic_type,
)
from kimmdy.utils import get_gmx_dir


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
def raw_urea_times_two_top_fix(filedir) -> TopologyDict:
    return read_top(filedir / "urea-times-2.top")


@pytest.fixture()
def raw_two_different_ureas_top_fix(filedir) -> TopologyDict:
    return read_top(filedir / "two-different-ureas.top")


@pytest.fixture()
def hexala_rawtop_fix(assetsdir, filedir) -> TopologyDict:
    return read_top(filedir / "hexala.top")


@pytest.fixture()
def hexala_top_fix(assetsdir, filedir) -> Topology:
    hexala_top = read_top(filedir / "hexala.top")
    return Topology(hexala_top)


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
    hexala_top = read_top(dir / "hexala.top")
    top = Topology(hexala_top)
    atomlist = draw(random_atomlist())
    top.atoms = {atom.nr: atom for atom in atomlist}
    top._regenerate_topology_from_bound_to()
    break_this = draw(st.sampled_from(list(top.bonds.keys())))
    return (top, break_this)


class TestGMX:
    def test_gmx_dir_is_found(self):
        gmx = get_gmx_dir()
        assert gmx
        assert Path(gmx).is_dir()


class TestMatch:
    @pytest.fixture
    def top_fix(self, filedir) -> Topology:
        hexala_top = read_top(filedir / "hexala.top")
        return Topology(hexala_top)

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
        assert top.moleculetypes["Reactive"].nrexcl == "3"

    def test_making_molecules_explicit(self, raw_urea_times_two_top_fix):
        raw = deepcopy(raw_urea_times_two_top_fix)
        top = Topology(raw)
        assert len(top.atoms) == 16
        assert len(top.bonds) == 14
        assert len(top.pairs) == 0
        assert len(top.angles) == 0
        assert len(top.proper_dihedrals) == 16
        assert len(top.improper_dihedrals) == 6
        assert top.moleculetypes["Reactive"].atoms == top.atoms
        for i in range(8):
            n1 = str(i + 1)
            n2 = str(i + 9)
            a1 = top.atoms[n1]
            a2 = top.atoms[n2]
            assert a1.atom == a2.atom
            assert a1.type == a2.type
            assert a1.mass == a2.mass
            assert int(a2.nr) == int(a1.nr) + 8
            assert int(a2.cgnr) == int(a1.cgnr) + 8

    def test_merging_molecules(self, raw_two_different_ureas_top_fix):
        raw = deepcopy(raw_two_different_ureas_top_fix)
        top = Topology(raw)
        assert len(top.atoms) == 16
        assert len(top.bonds) == 14
        assert len(top.pairs) == 0
        assert len(top.angles) == 0
        assert len(top.proper_dihedrals) == 16
        assert len(top.improper_dihedrals) == 6
        assert top.moleculetypes["Reactive"].atoms == top.atoms
        for i in range(8):
            n1 = str(i + 1)
            n2 = str(i + 9)
            a1 = top.atoms[n1]
            a2 = top.atoms[n2]
            assert a1.atom == a2.atom
            assert a1.type == a2.type
            assert a1.mass == a2.mass
            if i == 7:
                assert a1.residue == "TOTALLYNOTUREA"
                assert a2.residue == "URE"
            else:
                assert a1.residue == "URE"
                assert a2.residue == "URE"
            assert int(a2.nr) == int(a1.nr) + 8
            assert int(a2.cgnr) == int(a1.cgnr) + 8
            assert a1.bound_to_nrs == [str(int(x) - 8) for x in a2.bound_to_nrs]
            if i < 7:
                assert int(a1.resnr) == int(a2.resnr) - 2
            else:
                assert int(a1.resnr) == int(a2.resnr) - 1

    def test_merging_solvent(self, hexala_rawtop_fix):
        top = Topology(
            deepcopy(hexala_rawtop_fix),
            is_reactive_predicate_f=lambda mol: mol in ["Protein", "SOL"],
        )
        assert len(top.atoms) == 38013
        assert top.atoms["73"].resnr == "9"
        assert top.atoms["74"].resnr == "9"
        assert top.atoms["75"].resnr == "9"
        assert top.atoms["76"].resnr == "10"
        assert top.atoms["77"].resnr == "10"
        assert top.atoms["78"].resnr == "10"
        assert top.atoms["79"].resnr == "11"
        assert top.atoms["80"].resnr == "11"
        assert top.atoms["81"].resnr == "11"


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
    def test_reindex_no_change(self, hexala_top_fix: Topology):
        org_top: Topology = deepcopy(hexala_top_fix)
        update = hexala_top_fix.reindex_atomnrs()

        # test produced mapping
        assert len(update.keys()) == 12
        for mapping in update.values():
            for k, v in mapping.items():
                assert k == v

        # test topology
        assert org_top.atoms == hexala_top_fix.atoms
        assert org_top.bonds == hexala_top_fix.bonds
        assert org_top.angles == hexala_top_fix.angles
        assert org_top.proper_dihedrals == hexala_top_fix.proper_dihedrals
        assert org_top.improper_dihedrals == hexala_top_fix.improper_dihedrals

    @given(atomindex=st.integers(min_value=1, max_value=72))
    @settings(suppress_health_check=[HealthCheck.function_scoped_fixture], deadline=400)
    def test_del_atom_hexala(self, hexala_top_fix, atomindex):
        top: Topology = deepcopy(hexala_top_fix)

        atom = top.atoms[str(atomindex)]
        # bonds
        bound_nrs = deepcopy(atom.bound_to_nrs)
        bound_atms = [top.atoms[i] for i in bound_nrs]
        # angles
        angles_to_delete = []
        angles_to_update = []
        for key in top.angles.keys():
            if str(atomindex) in key:
                angles_to_delete.append(key)
            else:
                angles_to_update.append(key)
        # proper dihedrals
        pd_to_delete = []
        pd_to_update = []
        for key in top.proper_dihedrals.keys():
            if str(atomindex) in key:
                pd_to_delete.append(key)
            else:
                pd_to_update.append(key)
        # improper dihedrals
        id_to_delete = []
        id_to_update = []
        for key in top.improper_dihedrals.keys():
            if str(atomindex) in key:
                id_to_delete.append(key)
            else:
                id_to_update.append(key)

        update = top.del_atom(str(atomindex), parameterize=False)
        rev_update = {v: k for k, v in update.items()}

        for nr, atm in zip(bound_nrs, bound_atms):
            if int(nr) > atomindex:
                assert int(atm.nr) == int(nr) - 1
            else:
                assert int(atm.nr) == int(nr)
            assert atm.is_radical

        assert atom not in top.atoms.values()
        assert len(atom.bound_to_nrs) == 0

        # angles
        for a_del in angles_to_delete:
            assert None in [update.get(a) for a in a_del]
        for a_up in angles_to_update:
            k1, k2, k3 = [update[a] for a in a_up]
            new = top.angles[k1, k2, k3]
            old = hexala_top_fix.angles[a_up]
            assert old.ai == rev_update.get(new.ai)
            assert old.aj == rev_update.get(new.aj)
            assert old.ak == rev_update.get(new.ak)

        # proper dihedrals
        for pd_del in pd_to_delete:
            assert None in tuple([update.get(a) for a in pd_del])
        for pd_up in pd_to_update:
            k1, k2, k3, k4 = [update[a] for a in pd_up]
            new = top.proper_dihedrals[k1, k2, k3, k4]
            old = hexala_top_fix.proper_dihedrals[pd_up]
            assert old.ai == rev_update.get(new.ai)
            assert old.aj == rev_update.get(new.aj)
            assert old.ak == rev_update.get(new.ak)
            assert old.al == rev_update.get(new.al)

            for d_key in old.dihedrals:
                old_d = old.dihedrals[d_key]
                new_d = new.dihedrals[d_key]
                assert new_d.ai == update.get(old_d.ai)
                assert new_d.aj == update.get(old_d.aj)
                assert new_d.ak == update.get(old_d.ak)
                assert new_d.al == update.get(old_d.al)

        # improper dihedrals
        for id_del in id_to_delete:
            assert None in tuple([update.get(a) for a in id_del])
        for id_up in id_to_update:
            k1, k2, k3, k4 = [update[a] for a in id_up]
            new = top.improper_dihedrals[k1, k2, k3, k4]
            old = hexala_top_fix.improper_dihedrals[id_up]
            assert old.ai == rev_update.get(new.ai)
            assert old.aj == rev_update.get(new.aj)
            assert old.ak == rev_update.get(new.ak)
            assert old.al == rev_update.get(new.al)

            for d_key in old.dihedrals:
                old_d = old.dihedrals[d_key]
                new_d = new.dihedrals[d_key]
                assert new_d.ai == update.get(old_d.ai)
                assert new_d.aj == update.get(old_d.aj)
                assert new_d.ak == update.get(old_d.ak)
                assert new_d.al == update.get(old_d.al)

    def test_break_bind_bond_hexala(self, hexala_top_fix):
        top = deepcopy(hexala_top_fix)
        og_top = deepcopy(top)

        bondindex = 24
        bond_key = list(top.bonds.keys())[bondindex]
        logging.info(f"bond_key: {bond_key}")
        assert top.bonds.get(bond_key) is not None

        top.break_bond(bond_key)
        assert top.bonds.get(bond_key) is None
        assert not top.validate_bond(top.atoms[bond_key[0]], top.atoms[bond_key[1]])
        with pytest.raises(ValueError):
            top.break_bond(bond_key)

        top.bind_bond(bond_key)
        assert top.bonds.get(bond_key) is not None
        assert top.validate_bond(top.atoms[bond_key[0]], top.atoms[bond_key[1]])
        with pytest.raises(ValueError):
            top.bind_bond(bond_key)

        assert top.bonds == og_top.bonds
        assert top.pairs == og_top.pairs
        assert top.angles == og_top.angles
        assert top.proper_dihedrals == og_top.proper_dihedrals
        assert top.improper_dihedrals == og_top.improper_dihedrals

    @given(bondindex=st.integers(min_value=0, max_value=70))
    @settings(
        suppress_health_check=[HealthCheck.function_scoped_fixture], deadline=1000
    )
    @pytest.mark.slow
    def test_break_bind_random_bond_hexala(self, hexala_top_fix, bondindex):
        top = deepcopy(hexala_top_fix)
        og_top = deepcopy(top)
        bond_key = list(top.bonds.keys())[bondindex]
        assert top.bonds.get(bond_key) is not None
        top.break_bond(bond_key)
        assert top.bonds.get(bond_key) is None
        top.bind_bond(bond_key)
        assert top.bonds.get(bond_key) is not None

        assert top.bonds == og_top.bonds
        assert top.pairs == og_top.pairs
        assert top.angles == og_top.angles
        assert top.proper_dihedrals == og_top.proper_dihedrals
        assert top.improper_dihedrals == og_top.improper_dihedrals

    def test_generate_topology_from_bound_to(self, hexala_top_fix):
        og_top = deepcopy(hexala_top_fix)
        newtop = deepcopy(hexala_top_fix)
        newtop.bonds.clear()
        newtop.pairs.clear()
        newtop.angles.clear()
        newtop.proper_dihedrals.clear()

        assert newtop.bonds == {}
        assert newtop.moleculetypes[REACTIVE_MOLECULEYPE].bonds == {}

        newtop._regenerate_topology_from_bound_to()

        assert newtop.bonds == og_top.bonds
        assert newtop.pairs == og_top.pairs
        assert newtop.angles == og_top.angles
        assert newtop.proper_dihedrals == og_top.proper_dihedrals

    @pytest.mark.slow
    @settings(max_examples=1, phases=[Phase.generate], deadline=600)
    @given(top_break=random_topology_and_break())
    def test_break_bind_bond_invertible(self, top_break):
        top, to_break = top_break
        og_top = deepcopy(top)
        assert top.bonds.get(to_break) is not None
        top.break_bond(to_break)
        assert top.bonds.get(to_break) is None
        top.bind_bond(to_break)
        assert top.bonds.get(to_break) is not None

        assert top.bonds == og_top.bonds
        assert top.pairs == og_top.pairs
        assert top.angles == og_top.angles
        assert top.proper_dihedrals == og_top.proper_dihedrals


class TestHexalaTopology:
    @pytest.fixture
    def top_break_29_35_fix(self, filedir) -> Topology:
        hexala_top = read_top(filedir / "hexala_break29-35.top")
        return Topology(hexala_top)

    @pytest.fixture
    def top_hat_34_to_29_fix(self, filedir) -> Topology:
        hexala_top = read_top(filedir / "hexala_hat34-29.top")
        return Topology(
            hexala_top,
        )

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

        atoms = get_reactive_section(topology, "atoms")
        bonds = get_reactive_section(topology, "bonds")
        pairs = get_reactive_section(topology, "pairs")
        angles = get_reactive_section(topology, "angles")
        dihedrals = get_reactive_section(topology, "dihedrals")
        assert dihedrals
        proper_dihedrals = [x for x in dihedrals if x[4] == "9"]
        improper_dihedrals = [x for x in dihedrals if x[4] == "4"]

        atoms_new = get_reactive_section(topology_new, "atoms")
        bonds_new = get_reactive_section(topology_new, "bonds")
        pairs_new = get_reactive_section(topology_new, "pairs")
        angles_new = get_reactive_section(topology_new, "angles")
        dihedrals_new = get_reactive_section(topology_new, "dihedrals")
        assert dihedrals_new
        proper_dihedrals_new = [x for x in dihedrals_new if x[4] == "9"]
        improper_dihedrals_new = [x for x in dihedrals_new if x[4] == "4"]

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
        proper_dihedraldiff = set([tuple(x[:4]) for x in proper_dihedrals]) - set(
            [tuple(x[:4]) for x in proper_dihedrals_new]
        )
        improper_dihedraldiff = set([tuple(x[:4]) for x in improper_dihedrals]) - set(
            [tuple(x[:4]) for x in improper_dihedrals_new]
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
        assert proper_dihedraldiff == set(
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
            ]
        )
        assert improper_dihedraldiff == set(
            [
                ("7", "9", "15", "17"),
                ("9", "17", "15", "16"),
            ]
        )

    def test_top_update_dict(self, raw_hexala_top_fix):
        raw = raw_hexala_top_fix
        raw_copy = deepcopy(raw)
        top = Topology(raw_copy)
        top._update_dict()
        assert top.top["dihedraltypes"]["content"] == raw["dihedraltypes"]["content"]
        assert (
            top.top["moleculetype_Reactive"]["subsections"]["dihedrals"]["content"]
            == raw["moleculetype_Protein"]["subsections"]["dihedrals"]["content"]
        )
        # "fix" the expected differences to test the rest
        assert raw["molecules"]["content"][0] == ["Protein", "1"]
        assert top.top["molecules"]["content"][0] == ["Reactive", "1"]
        raw["molecules"]["content"][0] = ["Reactive", "1"]

        assert raw["moleculetype_Protein"]["content"][0] == ["Protein", "3"]
        assert top.top["moleculetype_Reactive"]["content"][0] == ["Reactive", "3"]
        for section in ["atoms", "bonds", "angles", "dihedrals", "pairs"]:
            assert (
                top.top["moleculetype_Reactive"]["subsections"][section]["content"]
                == raw["moleculetype_Protein"]["subsections"][section]["content"]
            )

    def test_top_properties(self, hexala_top_fix):
        top = deepcopy(hexala_top_fix)

        # initial correct number of atoms and bonds
        assert len(top.atoms) == 72
        assert len(top.bonds) == 71

        # order is correct
        val = 0
        for atom_id in top.atoms.keys():
            assert int(atom_id) > val
            val = int(atom_id)

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
        protein = top.moleculetypes[REACTIVE_MOLECULEYPE]
        atomnr = "29"
        bonds = protein._get_atom_bonds(atomnr)
        angles = protein._get_atom_angles(atomnr)
        proper_dihedrals = protein._get_atom_proper_dihedrals(atomnr)
        improper_dihedrals = protein._get_atom_improper_dihedrals(atomnr, top.ff)

        assert len(bonds) == 4
        assert len(angles) == 13
        assert len(proper_dihedrals) == 25
        assert len(improper_dihedrals) == 3

        atomnr = "9"
        bonds = protein._get_atom_bonds(atomnr)
        angles = protein._get_atom_angles(atomnr)
        proper_dihedrals_center = protein._get_center_atom_dihedrals(atomnr)

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

    def test_hat_34_to_29_after_break(self, hexala_top_fix, top_hat_34_to_29_fix):
        """HAT of H at 34 to C at 29"""
        top = deepcopy(hexala_top_fix)
        top_ref = deepcopy(top_hat_34_to_29_fix)
        top.break_bond(("29", "35"))
        top.break_bond(("31", "34"))
        top.bind_bond(("34", "29"))

        # compare topologies
        assert len(top.bonds) == len(top_ref.bonds)
        assert len(top.pairs) == len(top_ref.pairs)
        assert len(top.angles) == len(top_ref.angles)
        assert len(top.proper_dihedrals) == len(top_ref.proper_dihedrals)
        # assert len(top.improper_dihedrals) == len(top_ref.improper_dihedrals)

        # inspect HAT hydrogen
        h = top.atoms["34"]
        assert h.is_radical == False
        # changes name to 'HX' because HA already exists
        # could be smart to change this behavior to get HA1 and HA2 in this case
        assert h.atom == "HX"
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

    def test_get_residue_by_bonding(self, hexala_top_fix):
        top = hexala_top_fix
        a = top.atoms["1"]
        res = get_residue_by_bonding(a, top.atoms)
        assert len(res) == 6
        for a in res.values():
            assert a.residue == "ACE"
        a = top.atoms["10"]
        res = get_residue_by_bonding(a, top.atoms)
        assert len(res) == 10
        for a in res.values():
            assert a.residue == "ALA"
        a = top.atoms["25"]
        res = get_residue_by_bonding(a, top.atoms)
        for a in res.values():
            assert a.residue == "ALA"


class TestPolymerFF:
    def test_sections_are_complete(self, filedir):
        path = filedir / "polymer/topol.top"
        raw_top = read_top(path)
        top = Topology(raw_top)
        assert len(top.ff.nonbond_params) == 6
        assert top.ff.nonbond_params[("B0", "B0")] == NonbondParamType(
            i="B0",
            j="B0",
            funct="1",
            c0="2.5",
            c1="2.5",
            id="B0---B0",
            id_sym="B0---B0",
        )
        assert top.ff.atomtypes["B1"] == AtomType(
            type="B1",
            id="B1",
            id_sym="B1",
            at_num="",
            mass="20000.0",
            charge="0.000",
            ptype="A",
            sigma="0.0",
            epsilon="0.0",
        )

    def test_nrexcl_from_the_correct_moleculetype(self, filedir):
        path = filedir / "polymer/topol.top"
        raw_top = read_top(path)
        top = Topology(raw_top)
        assert top.moleculetypes["Reactive"].nrexcl == "1"


class TestRadicalAla:
    @pytest.fixture
    def top_noprm_fix(self, filedir) -> Topology:
        hexala_top = read_top(filedir / "Ala_R_noprm.top")
        return Topology(hexala_top)

    @pytest.fixture
    def top_noprm_explicitR_fix(self, filedir) -> Topology:
        hexala_top = read_top(filedir / "Ala_R_noprm.top")
        return Topology(hexala_top, radicals="9")

    @pytest.fixture
    def top_prm_fix(self, filedir) -> Topology:
        hexala_top = read_top(filedir / "Ala_R_prm.top")
        return Topology(hexala_top)

    def test_is_radical(self, top_noprm_fix):
        assert top_noprm_fix.atoms["9"].is_radical == True
        assert top_noprm_fix.atoms["10"].is_radical == False

    def test_is_radical_explicit(self, top_noprm_explicitR_fix):
        assert top_noprm_explicitR_fix.atoms["9"].is_radical == True
        assert top_noprm_explicitR_fix.atoms["10"].is_radical == False


class TestChargeAssignment:
    @pytest.fixture
    def top_charged(self, filedir) -> Topology:
        top = read_top(filedir / "IMREE.top")
        return Topology(top)

    def test_break_neutral(self, top_charged: Topology):
        recipe_steps = [Break(atom_id_1="7", atom_id_2="9")]
        top_charged.break_bond((recipe_steps[0].atom_id_1, recipe_steps[0].atom_id_2))
        top_charged.update_partial_charges(recipe_steps)
        fragment1_nrs = ["7", "8"]
        fragment2_nrs = [str(i) for i in range(9, 26)]
        assert np.isclose(
            sum([float(top_charged.atoms[nr].charge) for nr in fragment1_nrs]), 0
        )
        assert np.isclose(
            sum([float(top_charged.atoms[nr].charge) for nr in fragment2_nrs]), 0
        )

    def test_break_negative(self, top_charged: Topology):
        recipe_steps = [Break(atom_id_1="69", atom_id_2="80")]
        top_charged.break_bond((recipe_steps[0].atom_id_1, recipe_steps[0].atom_id_2))
        top_charged.update_partial_charges(recipe_steps)
        fragment1_nrs = [str(i) for i in range(67, 80)]
        fragment2_nrs = ["80", "81"]
        assert np.isclose(
            sum([float(top_charged.atoms[nr].charge) for nr in fragment1_nrs]), -1
        )
        assert np.isclose(
            sum([float(top_charged.atoms[nr].charge) for nr in fragment2_nrs]), 0
        )

    def test_break_positive(self, top_charged: Topology):
        recipe_steps = [Break(atom_id_1="33", atom_id_2="36")]
        top_charged.break_bond((recipe_steps[0].atom_id_1, recipe_steps[0].atom_id_2))
        top_charged.update_partial_charges(recipe_steps)
        fragment1_nrs = [str(i) for i in range(26, 36)] + ["48", "49"]
        fragment2_nrs = [str(i) for i in range(36, 48)]
        assert np.isclose(
            sum([float(top_charged.atoms[nr].charge) for nr in fragment1_nrs]), 0
        )
        assert np.isclose(
            sum([float(top_charged.atoms[nr].charge) for nr in fragment2_nrs]), 1
        )

    def test_HAT_charge_assignment(self, top_charged: Topology):
        recipe_steps = [
            Break(atom_id_1="50", atom_id_2="51"),
            Bind(atom_id_1="52", atom_id_2="51"),
        ]
        top_charged.del_atom("53")
        top_charged.break_bond((recipe_steps[0].atom_id_1, recipe_steps[0].atom_id_2))
        top_charged.bind_bond((recipe_steps[1].atom_id_1, recipe_steps[1].atom_id_2))
        top_charged.update_partial_charges(recipe_steps)
        fragment1_nrs = [str(i) for i in range(50, 66)]
        assert np.isclose(
            sum([float(top_charged.atoms[nr].charge) for nr in fragment1_nrs]), 0
        )
