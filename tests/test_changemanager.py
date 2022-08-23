from kimmdy.parsing import read_topol, read_plumed, topol_split_dihedrals
from kimmdy import changemanager

from pathlib import Path
from copy import deepcopy
import pytest


def test_break_bond_top():
    input_f = Path(__file__).parent / "test_files/test_changemanager/hexala_out.top"
    topology = read_topol(input_f)
    breakpair = (9, 15)
    breakpair = (str(breakpair[0]), str(breakpair[1]))

    topology_new = changemanager.break_bond_top(deepcopy(topology), breakpair)

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


def test_break_bond_plumed():
    input_f = Path(__file__).parent / "test_files/test_changemanager/plumed.dat"
    plumeddat = read_plumed(input_f)
    breakpair = (9, 15)
    breakpair = (str(breakpair[0]), str(breakpair[1]))

    newplumeddat = changemanager.break_bond_plumed(
        deepcopy(plumeddat), breakpair, "distances.dat"
    )

    oldset = set(tuple(x["atoms"]) for x in plumeddat["distances"])
    newset = set(tuple(x["atoms"]) for x in newplumeddat["distances"])
    diffs = list(oldset - newset)
    assert len(diffs) == 1 and diffs[0] == breakpair


def test_find_heavy():
    movepair = ["9", "10"]
    topology = {"bonds": []}
    for i in range(0, 20, 2):
        topology["bonds"].append([str(i), str(i + 1)])
    heavy_idx = changemanager.find_heavy(topology["bonds"], movepair[0])
    assert heavy_idx == "8"


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
        return

    def test_topol_change_at_an(self):
        atom = ["10", "HX", "H9", "ALA"]
        changemanager.topol_change_at_an(self.topology, atom)
        assert self.topology["atoms"][12][1] == "HX"
        assert self.topology["atoms"][12][4] == "H9"
        atom = ["20", "HC", "HB1", "ALA"]
        changemanager.topol_change_at_an(self.topology, atom)
        assert self.topology["atoms"][23][4] == "HB4"
        return


class TestLocalGraphConstructMethods:
    input_f = Path(__file__).parent / "test_files/test_changemanager/hexala_out.top"
    topology = read_topol(input_f)
    topology = topol_split_dihedrals(topology)

    input_ff = (
        Path(__file__).parent / "test_files/assets/amber99sb-star-ildnp.ff"
    )
    heavy_idx = "9"

    testGraph = changemanager.localGraph(topology, heavy_idx, input_ff)

    def test_construct_graph(self):
        self.testGraph.construct_graph()
        assert len(self.testGraph.AtomList) == 16
        assert len(self.testGraph.atoms_idx) == 16
        assert len(self.testGraph.BondList) == 15

    def test_order_lists(self):
        # need to run test_construct_graph beforehand
        self.testGraph.order_lists()
        val = 0
        for atom_idx in self.testGraph.atoms_idx:
            assert int(atom_idx) > val
            val = int(atom_idx)

        val = 0
        for bond in self.testGraph.BondList:
            assert int(bond[0]) >= val
            val = int(bond[0])
        return

    def test_update_atoms_list(self):
        self.testGraph.update_atoms_list()
        assert (
            len(self.testGraph.atoms_atomname) == 16
            and self.testGraph.atoms_atomname[0] == "CH3"
        )
        assert (
            len(self.testGraph.atoms_atomtype) == 16
            and self.testGraph.atoms_atomtype[0] == "CT"
        )
        assert (
            len(self.testGraph.atoms_resname) == 16
            and self.testGraph.atoms_resname[0] == "ACE"
        )

    def test_update_bound_to(self):
        self.testGraph.update_bound_to()
        for atom in self.testGraph.AtomList:
            assert len(atom.bound_to) > 0
            for bonded_idx in atom.bound_to:
                bonded_pos = self.testGraph.atoms_idx.index(bonded_idx)
                assert any(
                    [
                        x == atom.idx
                        for x in self.testGraph.AtomList[bonded_pos].bound_to
                    ]
                )


def test_build_PADs():
    input_f = Path(__file__).parent / "test_files/test_changemanager/hexala_out.top"
    topology = read_topol(input_f)
    topology = topol_split_dihedrals(topology)

    input_ff = (
        Path(__file__).parent / "test_files/assets/amber99sb-star-ildnp.ff"
    )
    heavy_idx = "9"

    newGraph = changemanager.localGraph(topology, "29", input_ff)
    newGraph.construct_graph(depth=20)
    newGraph.order_lists()
    newGraph.update_atoms_list()
    newGraph.update_bound_to()
    newGraph.build_PADs()

    assert len(newGraph.AtomList) == len(topology["atoms"]) - 9
    assert len(newGraph.BondList) == len(topology["bonds"]) - 1
    assert len(newGraph.PairList) == len(topology["pairs"]) - 1
    assert len(newGraph.AngleList) == len(topology["angles"]) - 1
    assert len(newGraph.ProperList) == len(topology["propers"]) - 1
    assert len(newGraph.ImproperList) == len(topology["impropers"]) - 1
    return


class TestLocalGraphAddRemoveMethods:
    testGraph = changemanager.localGraph(None, "1", None)

    def test_add_atom(self):
        self.testGraph.add_Atom(changemanager.Atom("2"))
        self.testGraph.add_Atom(changemanager.Atom("3"))
        self.testGraph.add_Atom(changemanager.Atom("4"))
        assert len(self.testGraph.AtomList) == 4
        assert len(self.testGraph.atoms_idx) == 4

    def test_add_bond(self):
        self.testGraph.add_Bond(["1", "2"])
        self.testGraph.add_Bond(["2", "3"])
        self.testGraph.add_Bond(["3", "4"])
        assert ["2", "3"] in self.testGraph.BondList
        assert len(self.testGraph.BondList) == 3

    def test_remove_terms(self):
        self.testGraph.topology = {"impropers": [["1", "3", "2", "4"]]}
        self.testGraph.order_lists()

        self.testGraph.update_bound_to()
        self.testGraph.build_PADs()

        rmvdir = {
            "bonds": [["1", "2"], ["2", "3"], ["3", "4"]],
            "pairs": [["1", "4"]],
            "angles": [["1", "2", "3"], ["2", "3", "4"]],
            "propers": [["1", "2", "3", "4"]],
            "impropers": [["1", "3", "2", "4"]],
        }
        self.testGraph.remove_terms(rmvdir)
        assert (
            self.testGraph.BondList
            == self.testGraph.PairList
            == self.testGraph.AngleList
            == self.testGraph.ProperList
            == self.testGraph.ImproperList
            == []
        )
        return


class TestLocalGraphFFMethods:
    input_ff = (
        Path(__file__).parent / "test_files/assets/amber99sb-star-ildnp.ff"
    )
    testGraph = changemanager.localGraph(None, "1", input_ff)

    input_f = Path(__file__).parent / "test_files/test_changemanager/hexala_out.top"
    topology = read_topol(input_f)
    topology = topol_split_dihedrals(topology)
    fullGraph = changemanager.localGraph(topology, "9", input_ff)
    fullGraph.construct_graph()
    fullGraph.order_lists()
    fullGraph.update_atoms_list()
    fullGraph.update_bound_to()
    fullGraph.build_PADs()
    fullGraph.order_lists()

    def test_ff_AA_todict(self):
        ffaminoacids = self.testGraph.ff_AA_todict()
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

    def test_get_H_FF_at_an(self):
        at_an1 = self.fullGraph.get_H_ff_at_an("9")
        assert at_an1[0] == "H1"
        assert at_an1[1] == "HA"
        at_an2 = self.fullGraph.get_H_ff_at_an("11")
        assert at_an2 == ("HC", "HB1")

    def test_compare_ff_impropers(self):
        self.fullGraph.ImproperList.append(
            ["5", "9", "7", "10", "4", "180.00", "4.60240", "2"]
        )
        adddict, rmvdict = self.fullGraph.compare_ff_impropers("9", "7")
        assert adddict["impropers"] == [["5", "9", "7", "8", "4"]]
        assert rmvdict["impropers"] == [
            ["5", "9", "7", "10", "4", "180.00", "4.60240", "2"]
        ]

    def test_get_ff_sections(self):
        self.testGraph.get_ff_sections()
        assert self.testGraph.ff["bondtypes"][1] == [
            "C",
            "C",
            "1",
            "0.1525",
            "259408.0",
            ";",
            "new99",
        ]
        assert self.testGraph.ff["angletypes"][-1] == [
            "CT",
            "CC",
            "CC",
            "1",
            "115.970",
            "541.070",
            ";",
            "AOK_parm.PYL",
        ]


class TestLocalGraphParameterize:
    input_f = Path(__file__).parent / "test_files/test_changemanager/AlaCaR_out.top"
    input_ff = (
        Path(__file__).parent / "test_files/assets/amber99sb-star-ildnp.ff"
    )

    topology = read_topol(input_f)
    topology = topol_split_dihedrals(topology)
    fullGraph = changemanager.localGraph(topology, "9", input_ff)
    fullGraph.construct_graph()
    fullGraph.order_lists()
    fullGraph.update_atoms_list()
    fullGraph.update_bound_to()
    fullGraph.build_PADs()
    fullGraph.order_lists()
    fullGraph.get_ff_sections()

    atom_terms = fullGraph.get_terms_with_atom("9", center=True)
    atom_terms["pairs"].clear()

    def test_search_terms(self):
        termdict = {}
        termdict["bonds"] = self.fullGraph.search_terms(
            self.fullGraph.BondList, "9", False
        )
        termdict["propers"] = self.fullGraph.search_terms(
            self.fullGraph.ProperList, "9", False
        )
        termdict["centerpropers"] = self.fullGraph.search_terms(
            self.fullGraph.ProperList, "9", True
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
        termdict_funct = self.fullGraph.add_function_to_terms(termdict)
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
        assert self.fullGraph.is_radical("9")
        assert not self.fullGraph.is_radical("10")

    def test_parameterize_bonded_terms(self):
        terms_bond = self.fullGraph.parameterize_bonded_terms(
            "bondtypes", self.fullGraph.BondList
        )
        terms_angles = self.fullGraph.parameterize_bonded_terms(
            "angletypes", self.fullGraph.AngleList
        )
        assert all([len(x) == 5 for x in terms_bond])
        assert terms_angles[4] == ["5", "7", "9", "1", "121.900", "418.400"]  # C N CA

    def test_patch_bond(self):
        self.atom_terms["bonds"] = self.fullGraph.parameterize_bonded_terms(
            "bondtypes", self.atom_terms["bonds"]
        )
        unpatched_terms = deepcopy(self.atom_terms["bonds"])
        newfrac = 0.98
        NCaR_offset = 0.006
        self.fullGraph.patch_bond(
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
        self.atom_terms["angles"] = self.fullGraph.parameterize_bonded_terms(
            "angletypes", self.atom_terms["angles"]
        )
        unpatched_terms = deepcopy(self.atom_terms["angles"])
        newtheteq = 117
        self.fullGraph.patch_angle(self.atom_terms["angles"], "9", newtheteq)
        assert all([int(float(x[4])) == newtheteq for x in self.atom_terms["angles"]])

    def test_patch_dihedral_CA(self):
        # maybe also test whether dihedrals are well formed
        self.atom_terms["propers"] = [[*x, "9"] for x in self.atom_terms["propers"]]
        unpatched_terms = deepcopy(self.atom_terms["propers"])
        phivals = ["1.6279944", "21.068532", "1.447664"]
        psivals = ["6.556746", "20.284450", "0.297901"]
        self.atom_terms["propers"] = self.fullGraph.patch_dihedral_CA(
            self.atom_terms["propers"], "9", phivals, psivals
        )
        assert self.atom_terms["propers"][0][6] == phivals[0]
        assert self.atom_terms["propers"][4][6] == psivals[1]

    def test_patch_improper(self):
        newphik = "43.93200"
        self.atom_terms["impropers"] = self.fullGraph.patch_improper("9", newphik)
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
        atom_terms_rad = self.fullGraph.parameterize_around_atom("9")
        atom_terms_nonrad = self.fullGraph.parameterize_around_atom("10")
        for key, val in atom_terms_rad.items():
            if not key in ["pairs"]:
                assert any(term[-1].endswith("parameter") for term in val)

        for key, val in atom_terms_nonrad.items():
            if not key in ["pairs", "impropers"]:
                assert not all(term[-1].endswith("parameter") for term in val)
