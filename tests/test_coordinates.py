import pytest
from pathlib import Path
import re
import MDAnalysis as mda
import subprocess as sp
from kimmdy.recipe import Break, Place
from kimmdy.parsing import read_plumed, read_top
from kimmdy.constants import REACTIVE_MOLECULEYPE
from conftest import DummyFiles
from kimmdy.tasks import TaskFiles
from kimmdy.topology.atomic import Dihedral
from kimmdy.topology.topology import Topology
from kimmdy.coordinates import (
    merge_top_slow_growth,
    get_explicit_or_type,
    place_atom,
    break_bond_plumed,
)
from kimmdy.utils import truncate_sim_files


def test_get_explicit_MultipleDihedrals(arranged_tmp_path):
    top_a = Topology(read_top(Path("topol_stateA.top")))
    # has only one dihedral
    key1 = ("1", "5", "7", "8")
    assert top_a.proper_dihedrals[key1].dihedrals[""] == Dihedral(*key1, "9")

    # has multiple dihedrals
    key_mult = ("15", "17", "19", "24")
    assert top_a.proper_dihedrals[key_mult].dihedrals.keys() == {"1", "2", "3"}
    assert top_a.proper_dihedrals[key_mult].dihedrals["1"] == Dihedral(
        *key_mult, funct="9", c0="180.000000", c1="1.6279944", periodicity="1"
    )
    assert top_a.proper_dihedrals[key_mult].dihedrals["2"] == Dihedral(
        *key_mult, funct="9", c0="180.000000", c1="21.068532", periodicity="2"
    )
    assert top_a.proper_dihedrals[key_mult].dihedrals["3"] == Dihedral(
        *key_mult, funct="9", c0="180.000000", c1="1.447664", periodicity="3"
    )


def test_get_bondobj(arranged_tmp_path):
    top_a = Topology(read_top(Path("topol_stateA.top")))

    bond1_keys = ("17", "18")
    bond1obj = get_explicit_or_type(
        bond1_keys,
        top_a.bonds[bond1_keys],
        top_a.ff.bondtypes,
        top_a.moleculetypes[REACTIVE_MOLECULEYPE],
    )

    bond2_keys = ("17", "19")
    bond2obj = get_explicit_or_type(
        bond2_keys,
        top_a.bonds[bond2_keys],
        top_a.ff.bondtypes,
        top_a.moleculetypes[REACTIVE_MOLECULEYPE],
    )
    assert bond1obj is not None
    assert bond2obj is not None
    assert bond1obj.c0 is not None
    assert bond1obj.c1 is not None
    assert bond2obj.c0 is not None
    assert bond2obj.c1 is not None
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

    step = Place(id_to_place="1", new_coords=(1, 2, 3))

    place_atom(files, step, ttime=100)

    assert files.output["trr"].exists()
    assert files.output["gro"].exists()

    u = mda.Universe(str(files.output["gro"]))
    coords = tuple(u.select_atoms(f"id {step.id_to_place}")[0].position)
    assert coords == step.new_coords


def test_plumed_break(arranged_tmp_path):
    """
    test plumed changes
    """
    files = TaskFiles(
        get_latest=lambda: "DummyCallable",
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


def test_merge_prm_top(arranged_tmp_path):
    """this tests a topology merge for a HAT reaction from a Ca (nr 19) radical to a N (nr 26) radical"""

    top_A = Topology(read_top(Path("topol_stateA.top")))
    top_B = Topology(read_top(Path("topol_stateB.top")))
    top_merge_ref = Topology(read_top(Path("topol_FEP.top")))

    top_merge = merge_top_slow_growth(top_A, top_B)

    assert top_merge.atoms == top_merge_ref.atoms
    assert top_merge.bonds.keys() == top_merge_ref.bonds.keys()
    assert top_merge.angles.keys() == top_merge_ref.angles.keys()
    assert top_merge.pairs.keys() == top_merge_ref.pairs.keys()
    assert top_merge.proper_dihedrals.keys() == top_merge_ref.proper_dihedrals.keys()
    assert (
        top_merge.improper_dihedrals.keys() == top_merge_ref.improper_dihedrals.keys()
    )
    assert len(top_merge.exclusions) == 1
    assert top_merge.exclusions == top_merge_ref.exclusions

    assert top_merge.bonds[("19", "27")].funct == "3"
    assert top_merge.bonds[("26", "27")].funct == "3"
    assert top_merge.angles[("17", "19", "20")].c3 is not None
    assert top_merge.proper_dihedrals[("15", "17", "19", "24")].dihedrals["3"].c5 == "3"
    assert (
        top_merge.improper_dihedrals[("17", "20", "19", "24")].dihedrals["2"].c5 == "2"
    )
    # assert one dihedral merge improper/proper


@pytest.mark.require_gmx
def test_truncate_sim_files(arranged_tmp_path):
    files = DummyFiles()
    files.input = {
        "trr": arranged_tmp_path / "relax.trr",
        "xtc": arranged_tmp_path / "relax.xtc",
        "edr": arranged_tmp_path / "relax.edr",
        "gro": arranged_tmp_path / "relax.gro",
    }
    files.outputdir = arranged_tmp_path
    time = 5.2
    truncate_sim_files(files, time)

    for p in files.input.values():
        assert p.exists()
        assert p.with_name(p.name + ".tail").exists()

    p = sp.run(
        f"gmx -quiet -nocopyright check -f {files.input['trr']}",
        text=True,
        capture_output=True,
        shell=True,
    )
    # FOR SOME REASON gmx check writes in stderr instead of stdout
    m = re.search(r"Last frame.*time\s+(\d+\.\d+)", p.stderr)
    assert m, p.stderr
    last_time = m.group(1)
    assert last_time == "5.000"
