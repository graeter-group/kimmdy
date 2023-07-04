import pytest
import shutil
from pathlib import Path

from kimmdy.reaction import Break, Bind, Move, RecipeStep
from kimmdy.parsing import read_plumed
from kimmdy.changemanager import break_bond_plumed, modify_coords
from conftest import SlimFiles
import os


@pytest.fixture
def tmpdir(tmp_path) -> Path:
    dirname = "test_changemanager"
    try:
        filedir = Path(__file__).parent / "test_files" / dirname
        assetsdir = Path(__file__).parent / "test_files" / "assets"
    except NameError:
        filedir = Path("./tests/test_files") / dirname
        assetsdir = Path("./tests/test_files") / "assets"
    test_dir = tmp_path / dirname
    shutil.copytree(filedir, test_dir)
    Path(test_dir / "amber99sb-star-ildnp.ff").symlink_to(
        assetsdir / "amber99sb-star-ildnp.ff",
        target_is_directory=True,
    )
    return test_dir


def test_plumed_break(tmpdir):
    plumed = read_plumed(tmpdir / "plumed_nat.dat")
    plumed_break_ref = read_plumed(tmpdir / "plumed_break29-35.dat")

    breakpair = (29, 35)
    plumed_break = break_bond_plumed(plumed, breakpair, Path("distances.dat"))

    assert plumed_break["distances"] == plumed_break_ref["distances"]
    assert plumed_break["prints"] == plumed_break_ref["prints"]


def test_modify_coords_break(tmpdir):
    steps: list[RecipeStep] = [Break(28, 34)]
    files = SlimFiles(outputdir=tmpdir)
    files.input["trr"] = tmpdir / "pull.trr"
    files.input["tpr"] = tmpdir / "pull.tpr"
    files.input["top"] = tmpdir / "hexala_out.top"
    files.output["top"] = tmpdir / "topol_mod.top"
    run_parameter_growth = modify_coords(steps, files)
    assert run_parameter_growth
    assert files.input["top"] == files.outputdir / "top_merge.top"
    assert (files.outputdir / "top_merge.top").exists()


def test_modify_coords_move(tmpdir):
    steps = [Move(ix_to_move=0, new_coords=[[0.0, 0.0, 0.0], 100.0])]
    files = SlimFiles(outputdir=tmpdir)
    files.input["trr"] = tmpdir / "pull.trr"
    files.input["tpr"] = tmpdir / "pull.tpr"
    run_parameter_growth = modify_coords(steps, files)
    assert not run_parameter_growth
    assert files.output["trr"].exists()
    assert files.output["gro"].exists()
    # could check whether the coordinates were actually changed, probably using mda
    # could even randomize idx and coords

    # def test_parameterize_bonded_terms(self):
    #     terms_bond = self.full_graph.parameterize_bonded_terms(
    #         "bondtypes", self.full_graph.bonds
    #     )
    #     terms_angles = self.full_graph.parameterize_bonded_terms(
    #         "angletypes", self.full_graph.angles
    #     )
    #     assert all([len(x) == 5 for x in terms_bond])
    #     assert terms_angles[4] == ["5", "7", "9", "1", "121.900", "418.400"]  # C N CA

    # def test_patch_bond(self):
    #     self.atom_terms["bonds"] = self.full_graph.parameterize_bonded_terms(
    #         "bondtypes", self.atom_terms["bonds"]
    #     )
    #     unpatched_terms = deepcopy(self.atom_terms["bonds"])
    #     newfrac = 0.98
    #     NCaR_offset = 0.006
    #     self.full_graph.patch_bond(
    #         self.atom_terms["bonds"], "9", newfrac, 0.955, NCaR_offset, CNR_offset=0.008
    #     )
    #     assert (
    #         pytest.approx(float(self.atom_terms["bonds"][0][3]), 0.002)
    #         == float(unpatched_terms[0][3]) * newfrac - NCaR_offset
    #     )
    #     assert (
    #         pytest.approx(float(self.atom_terms["bonds"][1][3]), 0.002)
    #         == float(unpatched_terms[1][3]) * newfrac
    #     )

    # def test_patch_angle(self):
    #     self.atom_terms["angles"] = self.full_graph.parameterize_bonded_terms(
    #         "angletypes", self.atom_terms["angles"]
    #     )
    #     newtheteq = 117
    #     self.full_graph.patch_angle(self.atom_terms["angles"], "9", newtheteq)
    #     assert all([int(float(x[4])) == newtheteq for x in self.atom_terms["angles"]])

    # def test_patch_dihedral_CA(self):
    #     # maybe also test whether dihedrals are well formed
    #     self.atom_terms["propers"] = [[*x, "9"] for x in self.atom_terms["propers"]]
    #     phivals = ["1.6279944", "21.068532", "1.447664"]
    #     psivals = ["6.556746", "20.284450", "0.297901"]
    #     self.atom_terms["propers"] = self.full_graph.patch_dihedral_CA(
    #         self.atom_terms["propers"], "9", phivals, psivals
    #     )
    #     assert self.atom_terms["propers"][0][6] == phivals[0]
    #     assert self.atom_terms["propers"][4][6] == psivals[1]

    # def test_patch_improper(self):
    #     newphik = "43.93200"
    #     self.atom_terms["impropers"] = self.full_graph.patch_improper("9", newphik)
    #     assert self.atom_terms["impropers"] == [
    #         [
    #             "7",
    #             "10",
    #             "9",
    #             "14",
    #             "4",
    #             "180.0000000",
    #             str(newphik),
    #             "2",
    #             ";",
    #             "patched",
    #             "parameter",
    #         ]
    #     ]

    # def test_parameterize_around_atom(self):
    #     atom_terms_rad = self.full_graph.parameterize_around_atom("9")
    #     atom_terms_nonrad = self.full_graph.parameterize_around_atom("10")
    #     for key, val in atom_terms_rad.items():
    #         if not key in ["pairs"]:
    #             assert any(term[-1].endswith("parameter") for term in val)

    #     for key, val in atom_terms_nonrad.items():
    #         if not key in ["pairs", "impropers"]:
    #             assert not all(term[-1].endswith("parameter") for term in val)
