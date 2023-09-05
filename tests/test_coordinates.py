import pytest
import shutil
from pathlib import Path
import os
import MDAnalysis as mda

from kimmdy.recipe import Break, Bind, Place, Relax, RecipeStep
from kimmdy.parsing import read_plumed, read_top, read_top
from kimmdy.plugins import BasicParameterizer
from conftest import DummyFiles
from kimmdy.tasks import TaskFiles
from kimmdy.topology.topology import Topology
from kimmdy.topology.atomic import Bond
from kimmdy.coordinates import (
    merge_top_slow_growth,
    get_explicit_or_type,
    place_atom,
    break_bond_plumed,
)


## test coordinate changes
def test_get_bondobj(arranged_tmp_path):
    top_A = Topology(read_top(Path("topol_stateA.top")))

    bond1_keys = ("17", "18")
    bond1obj = get_explicit_or_type(
        bond1_keys,
        top_A.bonds[bond1_keys],
        top_A.ff.bondtypes,
        top_A.moleculetypes["Protein"],
    )

    bond2_keys = ("17", "19")
    bond2obj = get_explicit_or_type(
        bond2_keys,
        top_A.bonds[bond2_keys],
        top_A.ff.bondtypes,
        top_A.moleculetypes["Protein"],
    )
    assert float(bond1obj.c0) == pytest.approx(0.10100) and float(
        bond1obj.c1
    ) == pytest.approx(363171.2)

    assert float(bond2obj.c0) == pytest.approx(0.13600) and float(
        bond2obj.c1
    ) == pytest.approx(282001.6)


def test_place_atom(arranged_tmp_path):
    files = DummyFiles()
    files.input = {
        "tpr": arranged_tmp_path / "pull.tpr",
        "trr": arranged_tmp_path / "pull.trr",
    }
    files.outputdir = arranged_tmp_path

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


def test_plumed_break(arranged_tmp_path):
    files = TaskFiles(
        get_latest=lambda: f"DummyCallable",
    )
    files.input = {
        "plumed": arranged_tmp_path / "plumed_nat.dat",
        "plumed_out": Path("distances.dat"),
    }
    files.outputdir = arranged_tmp_path
    recipe_steps = [Break(28, 34)]
    breakpair = (recipe_steps[0].atom_id_1, recipe_steps[0].atom_id_2)

    break_bond_plumed(files, breakpair, Path("plumed_mod.dat"))

    plumed_nat = read_plumed(files.input["plumed"])
    plumed_break_ref = read_plumed(arranged_tmp_path / "plumed_break29-35.dat")
    plumed_break_test = read_plumed(files.output["plumed"])

    assert plumed_break_test["labeled_action"] == plumed_break_ref["labeled_action"]
    assert plumed_break_test["prints"] == plumed_break_ref["prints"]
    assert (
        len(plumed_break_test["labeled_action"])
        == len(plumed_nat["labeled_action"]) - 1
    )


## test topology changes
def test_merge_prm_top(arranged_tmp_path):
    """this tests a topology merge for a HAT reaction from a Ca (nr 19) radical to a N (nr 26) radical"""

    top_A = Topology(read_top(Path("topol_stateA.top")))
    top_B = Topology(read_top(Path("topol_stateB.top")))
    top_merge_ref = Topology(read_top(Path("topol_FEP.top")))

    top_merge = merge_top_slow_growth(top_A, top_B)

    # write_top(
    #     topmerge.to_dict(),
    #     Path(
    #         "/hits/fast/mbm/hartmaec/kimmdys/kimmdy_main/tests/test_files/test_coordinates/topol_curr.top"
    #     ),
    # )

    assert top_merge.atoms == top_merge_ref.atoms
    assert top_merge.bonds.keys() == top_merge_ref.bonds.keys()
    assert top_merge.angles.keys() == top_merge_ref.angles.keys()
    assert top_merge.pairs.keys() == top_merge_ref.pairs.keys()
    assert top_merge.proper_dihedrals.keys() == top_merge_ref.proper_dihedrals.keys()
    assert (
        top_merge.improper_dihedrals.keys() == top_merge_ref.improper_dihedrals.keys()
    )

    assert top_merge.bonds[("19", "27")].funct == "3"
    assert top_merge.bonds[("26", "27")].funct == "3"
    assert top_merge.angles[("17", "19", "20")].c3 != None
    assert top_merge.proper_dihedrals[("15", "17", "19", "24")].dihedrals["3"].c5 == "3"
    assert top_merge.improper_dihedrals[("17", "20", "19", "24")].c5 == "2"
    # assert one dihedral merge improper/proper
