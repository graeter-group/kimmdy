from kimmdy.parsing import read_topol, read_plumed
from kimmdy import changemanager
import pytest
import os

from pathlib import Path
from copy import deepcopy

from kimmdy.reaction import Conversion, ConversionRecipe, ConversionType, ReactionOutcome, ReactionResult

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


def test_break_bond_plumed():
    plumeddat = read_plumed(Path('plumed.dat'))
    breakpair = ('9', '15')

    recipe = [Conversion(ConversionType.BREAK, breakpair)]

    changemanager.modify_plumed(
        recipe,
        Path("plumed.dat"),
        Path("plumed-mod.dat"),
        Path("distances.dat"),
    )

    newplumeddat = read_plumed(Path('plumed-mod.dat'))

    oldset = set(tuple(x["atoms"]) for x in plumeddat["distances"])
    newset = set(tuple(x["atoms"]) for x in newplumeddat["distances"])
    diffs = list(oldset - newset)
    assert len(diffs) == 1 and diffs[0] == breakpair


class TestTopologyMethods:
    input_f = Path(__file__).parent / "test_files/test_changemanager/hexala_out.top"
    topology = read_topol(input_f)
    topology = topol_split_dihedrals(topology)
    termdict = {
        "bonds": [["9", "10", "1"]],
        "pairs": [["5", "10", "1"], ["8", "10", "1"]],
        "angles": [["7", "9", "10", "1"], ["10", "9", "15", "1"]],
        "propers": [["5", "7", "9", "10", "9"]],
        "impropers": [["5", "9", "7", "8", "4"]],
    }

    def test_topol_remove(self):
        self.topology = changemanager.topol_remove_terms(self.topology, self.termdict)
        for key, val in self.termdict.items():
            for entry in val:
                assert entry not in self.topology[key]

    def test_topol_add(self):
        self.topology = changemanager.topol_add_terms(self.topology, self.termdict)

        for key, val in self.termdict.items():
            for entry in val:
                assert entry in self.topology[key]



class TestLocalGraphConstructMethods:
    input_f = Path(__file__).parent / "test_files/test_changemanager/hexala_out.top"
    topology = read_topol(input_f)
    topology = topol_split_dihedrals(topology)

    input_ff = Path(__file__).parent / "test_files/assets/amber99sb-star-ildnp.ff"
    heavy_idx = "9"

    test_graph = changemanager.LocalGraph(topology, heavy_idx, input_ff)

    def test_construct_graph(self):
        self.test_graph.construct_graph()
        assert len(self.test_graph.atoms) == 16
        assert len(self.test_graph.atoms_idx) == 16
        assert len(self.test_graph.bonds) == 15

    def test_order_lists(self):
        # need to run test_construct_graph beforehand
        self.test_graph.order_lists()
        val = 0
        for atom_idx in self.test_graph.atoms_idx:
            assert int(atom_idx) > val
            val = int(atom_idx)

        val = 0
        for bond in self.test_graph.bonds:
            assert int(bond[0]) >= val
            val = int(bond[0])
        return

    def test_update_atoms_list(self):
        self.test_graph.update_atoms_list()
        assert (
            len(self.test_graph.atoms_atomname) == 16
            and self.test_graph.atoms_atomname[0] == "CH3"
        )
        assert (
            len(self.test_graph.atoms_atomtype) == 16
            and self.test_graph.atoms_atomtype[0] == "CT"
        )
        assert (
            len(self.test_graph.atoms_resname) == 16
            and self.test_graph.atoms_resname[0] == "ACE"
        )

    def test_update_bound_to(self):
        self.test_graph.update_bound_to()
        for atom in self.test_graph.atoms:
            assert len(atom.bound_to) > 0
            for bonded_idx in atom.bound_to:
                bonded_pos = self.test_graph.atoms_idx.index(bonded_idx)
                assert any(
                    [x == atom.idx for x in self.test_graph.atoms[bonded_pos].bound_to]
                )


def test_build_PADs():
    input_f = Path(__file__).parent / "test_files/test_changemanager/hexala_out.top"
    topology = read_topol(input_f)
    topology = topol_split_dihedrals(topology)
    input_ff = Path(__file__).parent / "test_files/assets/amber99sb-star-ildnp.ff"

    new_graph = changemanager.LocalGraph(topology, "29", input_ff, None, 20)

    assert len(new_graph.atoms) == len(topology["atoms"]) - 9
    assert len(new_graph.bonds) == len(topology["bonds"]) - 1
    assert len(new_graph.pairs) == len(topology["pairs"]) - 1
    assert len(new_graph.angles) == len(topology["angles"]) - 1
    assert len(new_graph.proper_dihedrals) == len(topology["propers"]) - 1
    assert len(new_graph.improper_dihedrals) == len(topology["impropers"]) - 1


class TestLocalGraphAddRemoveMethods:
    input_f = Path(__file__).parent / "test_files/test_changemanager/hexala_out.top"
    topology = read_topol(input_f)
    topology = topol_split_dihedrals(topology)
    input_ff = Path(__file__).parent / "test_files/assets/amber99sb-star-ildnp.ff"
    test_graph = changemanager.LocalGraph(topology, "1", input_ff, None, 3)
    original_graph = deepcopy(test_graph)
    n_atoms = len(test_graph.get_atoms_property("idx"))
    n_bonds = len(test_graph.bonds)

    def test_add_atom(self):
        self.test_graph.add_atom(changemanager.Atom("10"))
        self.test_graph.add_atom(changemanager.Atom("11"))
        self.test_graph.add_atom(changemanager.Atom("12"))

        assert len(self.test_graph.atoms) == self.n_atoms + 3
        assert len(self.test_graph.atoms_idx) == self.n_atoms + 3

    def test_add_bond(self):
        self.test_graph.add_bond(["9", "10"])
        self.test_graph.add_bond(["10", "11"])
        self.test_graph.add_bond(["11", "12"])

        assert ["11", "12"] in self.test_graph.bonds
        assert len(self.test_graph.bonds) == self.n_bonds + 3

    def test_remove_terms(self):
        # TODO: rewrite
        pass


class TestLocalGraphFFMethods:
    input_ff = Path(__file__).parent / "test_files/assets/amber99sb-star-ildnp.ff"
    input_f = Path(__file__).parent / "test_files/test_changemanager/hexala_out.top"
    topology = read_topol(input_f)
    topology = topol_split_dihedrals(topology)
    test_graph = changemanager.LocalGraph(topology, "1", input_ff)

    full_graph = changemanager.LocalGraph(topology, "9", input_ff)
    full_graph.construct_graph()
    full_graph.order_lists()
    full_graph.update_atoms_list()
    full_graph.update_bound_to()
    full_graph.build_PADs()
    full_graph.order_lists()

    def test_ff_AA_todict(self):
        ffaminoacids = self.test_graph.ff_AA_todict()
        assert "ALA" in ffaminoacids.keys()
        assert all(
            x in ffaminoacids["GLY"]["atoms"]
            for x in [
                ["N", "N", "0.41570", "1"],
                ["H", "H", "0.27190", "2"],
                ["CA", "CT", "0.02520", "3"],
                ["HA1", "H1", "0.06980", "4"],
                ["HA2", "H1", "0.06980", "5"],
                ["C", "C", "0.59730", "6"],
                ["O", "O", "0.56790", "7"],
            ]
        )
        assert all(
            x in ffaminoacids["NME"]["bonds"]
            for x in [
                ["N", "H"],
                ["N", "CH3"],
                ["CH3", "HH31"],
                ["CH3", "HH32"],
                ["CH3", "HH33"],
                ["C", "N"],
            ]
        )

class TestLocalGraphParameterize:
    input_f = Path(__file__).parent / "test_files/test_changemanager/AlaCaR_out.top"
    input_ff = Path(__file__).parent / "test_files/assets/amber99sb-star-ildnp.ff"

    topology = read_topol(input_f)
    topology = topol_split_dihedrals(topology)
    full_graph = changemanager.LocalGraph(topology, "9", input_ff)
    full_graph.get_ff_sections()

    atom_terms = full_graph.get_terms_with_atom("9", center=True)
    atom_terms["pairs"].clear()

    def test_search_terms(self):
        termdict = {}
        termdict["bonds"] = self.full_graph.search_terms(
            self.full_graph.bonds, "9", False
        )
        termdict["propers"] = self.full_graph.search_terms(
            self.full_graph.proper_dihedrals, "9", False
        )
        termdict["centerpropers"] = self.full_graph.search_terms(
            self.full_graph.proper_dihedrals, "9", True
        )
        assert termdict["bonds"] == [["7", "9"], ["9", "10"], ["9", "14"]]
        assert termdict["propers"][0] == ["1", "5", "7", "9"]
        assert termdict["centerpropers"] == [
            ["5", "7", "9", "10"],
            ["5", "7", "9", "14"],
            ["8", "7", "9", "10"],
            ["8", "7", "9", "14"],
            ["7", "9", "10", "11"],
            ["7", "9", "10", "12"],
            ["7", "9", "10", "13"],
            ["14", "9", "10", "11"],
            ["14", "9", "10", "12"],
            ["14", "9", "10", "13"],
            ["7", "9", "14", "15"],
            ["7", "9", "14", "16"],
            ["10", "9", "14", "15"],
            ["10", "9", "14", "16"],
        ]

    def test_add_function_to_terms(self):
        termdict = {
            "bonds": [["b1"], ["b2"]],
            "pairs": [["p1"], ["p2"]],
            "angles": [["a1"], ["a2"]],
            "propers": [["p1"], ["p2"]],
            "impropers": [["i1"], ["i2"]],
        }
        termdict_funct = self.full_graph.add_function_to_terms(termdict)
        funct_dict = {
            "bonds": "1",
            "pairs": "1",
            "angles": "1",
            "propers": "9",
            "impropers": "4",
        }
        for key, val in termdict_funct.items():
            for entry in val:
                assert entry[-1] == funct_dict[key]

    def test_is_radical(self):
        assert self.full_graph.is_radical("9")
        assert not self.full_graph.is_radical("10")

    def test_parameterize_bonded_terms(self):
        terms_bond = self.full_graph.parameterize_bonded_terms(
            "bondtypes", self.full_graph.bonds
        )
        terms_angles = self.full_graph.parameterize_bonded_terms(
            "angletypes", self.full_graph.angles
        )
        assert all([len(x) == 5 for x in terms_bond])
        assert terms_angles[4] == ["5", "7", "9", "1", "121.900", "418.400"]  # C N CA

    def test_patch_bond(self):
        self.atom_terms["bonds"] = self.full_graph.parameterize_bonded_terms(
            "bondtypes", self.atom_terms["bonds"]
        )
        unpatched_terms = deepcopy(self.atom_terms["bonds"])
        newfrac = 0.98
        NCaR_offset = 0.006
        self.full_graph.patch_bond(
            self.atom_terms["bonds"], "9", newfrac, 0.955, NCaR_offset, CNR_offset=0.008
        )
        assert (
            pytest.approx(float(self.atom_terms["bonds"][0][3]), 0.002)
            == float(unpatched_terms[0][3]) * newfrac - NCaR_offset
        )
        assert (
            pytest.approx(float(self.atom_terms["bonds"][1][3]), 0.002)
            == float(unpatched_terms[1][3]) * newfrac
        )

    def test_patch_angle(self):
        self.atom_terms["angles"] = self.full_graph.parameterize_bonded_terms(
            "angletypes", self.atom_terms["angles"]
        )
        newtheteq = 117
        self.full_graph.patch_angle(self.atom_terms["angles"], "9", newtheteq)
        assert all([int(float(x[4])) == newtheteq for x in self.atom_terms["angles"]])

    def test_patch_dihedral_CA(self):
        # maybe also test whether dihedrals are well formed
        self.atom_terms["propers"] = [[*x, "9"] for x in self.atom_terms["propers"]]
        phivals = ["1.6279944", "21.068532", "1.447664"]
        psivals = ["6.556746", "20.284450", "0.297901"]
        self.atom_terms["propers"] = self.full_graph.patch_dihedral_CA(
            self.atom_terms["propers"], "9", phivals, psivals
        )
        assert self.atom_terms["propers"][0][6] == phivals[0]
        assert self.atom_terms["propers"][4][6] == psivals[1]

    def test_patch_improper(self):
        newphik = "43.93200"
        self.atom_terms["impropers"] = self.full_graph.patch_improper("9", newphik)
        assert self.atom_terms["impropers"] == [
            [
                "7",
                "10",
                "9",
                "14",
                "4",
                "180.0000000",
                str(newphik),
                "2",
                ";",
                "patched",
                "parameter",
            ]
        ]

    def test_parameterize_around_atom(self):
        atom_terms_rad = self.full_graph.parameterize_around_atom("9")
        atom_terms_nonrad = self.full_graph.parameterize_around_atom("10")
        for key, val in atom_terms_rad.items():
            if not key in ["pairs"]:
                assert any(term[-1].endswith("parameter") for term in val)

        for key, val in atom_terms_nonrad.items():
            if not key in ["pairs", "impropers"]:
                assert not all(term[-1].endswith("parameter") for term in val)
