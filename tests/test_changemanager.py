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
