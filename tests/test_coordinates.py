import pytest
import shutil
from pathlib import Path
import os
import MDAnalysis as mda

from kimmdy.recipe import Break, Bind, Place, Relax, RecipeStep
from kimmdy.parsing import read_plumed, read_top, read_top
from kimmdy.plugins import BasicParameterizer
from conftest import SlimFiles
from kimmdy.tasks import TaskFiles
from kimmdy.topology.topology import Topology
from kimmdy.topology.atomic import Bond
from kimmdy.coordinates import (
    merge_top_slow_growth,
    get_explicit_or_type,
    place_atom,
    break_bond_plumed,
)


@pytest.fixture
def tmpdir(tmp_path) -> Path:
    dirname = "test_coordinates"
    try:
        filedir = Path(__file__).parent / "test_files" / dirname
        assetsdir = Path(__file__).parent / "test_files" / "assets"
    except NameError:
        filedir = Path("./tests/test_files") / dirname
        assetsdir = Path("./tests/test_files") / "assets"
    test_dir = tmp_path / dirname
    shutil.copytree(filedir, test_dir)

    return test_dir


## test coordinate changes


@pytest.fixture(scope="module")
def coordinates_files():
    filedir = Path(__file__).parent / "test_files" / "test_coordinates"
    ffdir = Path(__file__).parent / "test_files" / "assets" / "amber99sb-star-ildnp.ff"
    topA_path = filedir / "topol_stateA.top"
    topB_path = filedir / "topol_stateB.top"
    topFEP_path = filedir / "topol_FEP.top"
    if (filedir / "amber99sb-star-ildnp.ff").exists():
        (filedir / "amber99sb-star-ildnp.ff").unlink()
    (filedir / "amber99sb-star-ildnp.ff").symlink_to(
        ffdir,
        target_is_directory=True,
    )
    stateA = read_top(topA_path)
    stateB = read_top(topB_path)
    fep = read_top(topFEP_path)
    topA = Topology(stateA)
    topB = Topology(stateB)
    topFEP = Topology(fep)

    files = {
        "topA": topA,
        "topB": topB,
        "topFEP": topFEP,
        "topA_path": topA_path,
        "topB_path": topB_path,
        "topFEP_path": topFEP,
        "fep": fep,
    }
    yield files
    (filedir / "amber99sb-star-ildnp.ff").unlink()


def test_get_bondobj(coordinates_files):
    bond1_keys = ("17", "18")
    bond1obj = get_explicit_or_type(
        bond1_keys,
        coordinates_files["topA"].bonds[bond1_keys],
        coordinates_files["topA"].ff.bondtypes,
        coordinates_files["topA"],
    )

    bond2_keys = ("17", "19")
    bond2obj = get_explicit_or_type(
        bond2_keys,
        coordinates_files["topA"].bonds[bond2_keys],
        coordinates_files["topA"].ff.bondtypes,
        coordinates_files["topA"],
    )
    assert float(bond1obj.c0) == pytest.approx(0.10100) and float(
        bond1obj.c1
    ) == pytest.approx(363171.2)

    assert float(bond2obj.c0) == pytest.approx(0.13600) and float(
        bond2obj.c1
    ) == pytest.approx(282001.6)


def test_place_atom(tmpdir):
    files = TaskFiles(
        get_latest=lambda: f"DummyCallable",
    )
    files.input = {
        "tpr": tmpdir / "pull.tpr",
        "trr": tmpdir / "pull.trr",
    }
    files.outputdir = tmpdir

    step = Place(id_to_place="1", new_coords=(0, 0, 0))
    timespan = [(0.0, 100.0)]

    place_atom(files, step, timespan)

    assert files.output["trr"].exists()
    assert files.output["gro"].exists()

    u = mda.Universe(str(files.output["gro"]))
    coords = tuple(u.select_atoms(f"id {step.id_to_place}")[0].position)
    assert coords == step.new_coords
    print(coords)


## test plumed changes


def test_plumed_break(tmpdir):
    files = TaskFiles(
        get_latest=lambda: f"DummyCallable",
    )
    files.input = {
        "plumed": tmpdir / "plumed_nat.dat",
        "plumed_out": Path("distances.dat"),
    }
    files.outputdir = tmpdir
    recipe_steps = [Break(28, 34)]
    breakpair = (recipe_steps[0].atom_id_1, recipe_steps[0].atom_id_2)

    break_bond_plumed(files, breakpair, Path("plumed_mod.dat"))

    plumed_nat = read_plumed(files.input["plumed"])
    plumed_break_ref = read_plumed(tmpdir / "plumed_break29-35.dat")
    plumed_break_test = read_plumed(files.output["plumed"])

    assert plumed_break_test["distances"] == plumed_break_ref["distances"]
    assert plumed_break_test["prints"] == plumed_break_ref["prints"]
    assert len(plumed_break_test["distances"]) == len(plumed_nat["distances"]) - 1


## test topology changes
def test_merge_prm_top(coordinates_files):
    """this tests a topology merge for a HAT reaction from a Ca (nr 19) radical to a N (nr 26) radical"""
    topmerge = merge_top_slow_growth(
        coordinates_files["topA"], coordinates_files["topB"]
    )

    # write_top(
    #     topmerge.to_dict(),
    #     Path(
    #         "/hits/fast/mbm/hartmaec/kimmdys/kimmdy_main/tests/test_files/test_coordinates/topol_curr.top"
    #     ),
    # )

    assert topmerge.atoms == coordinates_files["topFEP"].atoms
    assert topmerge.bonds.keys() == coordinates_files["topFEP"].bonds.keys()
    assert topmerge.angles.keys() == coordinates_files["topFEP"].angles.keys()
    assert topmerge.pairs.keys() == coordinates_files["topFEP"].pairs.keys()
    assert (
        topmerge.proper_dihedrals.keys()
        == coordinates_files["topFEP"].proper_dihedrals.keys()
    )
    assert (
        topmerge.improper_dihedrals.keys()
        == coordinates_files["topFEP"].improper_dihedrals.keys()
    )

    assert topmerge.bonds[("19", "27")].funct == "3"
    assert topmerge.bonds[("26", "27")].funct == "3"
    assert topmerge.angles[("17", "19", "20")].c3 != None
    assert topmerge.proper_dihedrals[("15", "17", "19", "24")].dihedrals["3"].c5 == "3"
    assert topmerge.improper_dihedrals[("17", "20", "19", "24")].c5 == "2"
    # asster one dihedral merge improper/proper
