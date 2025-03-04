import copy
import logging
import numpy as np
from math import sqrt
import pytest
from pathlib import Path
import MDAnalysis as mda
from kimmdy.recipe import Break, Place
from kimmdy.parsing import read_plumed, read_top, write_top
from kimmdy.constants import DEFAULT_EDISSOC, REACTIVE_MOLECULEYPE
from conftest import DummyFiles
from kimmdy.tasks import TaskFiles
from kimmdy.topology.atomic import Bond, BondType, Dihedral, Pair
from kimmdy.topology.topology import Topology
from kimmdy.coordinates import (
    MoleculeTypeMerger,
    merge_top_slow_growth,
    get_explicit_or_type,
    place_atom,
    break_bond_plumed,
)


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
    assert isinstance(bond1obj, BondType)
    assert isinstance(bond2obj, Bond)
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
        get_latest=lambda: "DummyCallable",  # pyright: ignore
    )
    files.input = {
        "plumed": arranged_tmp_path / "plumed_nat.dat",
        "plumed_out": Path("distances.dat"),
    }
    files.outputdir = arranged_tmp_path
    recipe_steps = [Break(28, 34)]
    breakpair = (recipe_steps[0].atom_id_1, recipe_steps[0].atom_id_2)

    break_bond_plumed(files, breakpair, Path("plumed_mod.dat"))
    plumed_path = files.input["plumed"]
    assert plumed_path is not None
    plumed_nat = read_plumed(plumed_path)
    plumed_break_ref = read_plumed(arranged_tmp_path / "plumed_break29-35.dat")
    plumed_break_test = read_plumed(files.output["plumed"])

    assert plumed_break_test["labeled_action"] == plumed_break_ref["labeled_action"]
    assert plumed_break_test["prints"] == plumed_break_ref["prints"]
    assert (
        len(plumed_break_test["labeled_action"])
        == len(plumed_nat["labeled_action"]) - 1
    )


def test_morse_parameters(arranged_tmp_path):
    top_a = Topology(read_top(Path("topol_stateA.top")))
    top_b = Topology(read_top(Path("topol_stateB.top")))
    merger = MoleculeTypeMerger(
        top_a.moleculetypes[REACTIVE_MOLECULEYPE],
        top_b.moleculetypes[REACTIVE_MOLECULEYPE],
        ff=top_a.ff,
    )

    bond_key = ("19", "27")
    bond_a = merger.mol_a.bonds.get(bond_key)

    bond_b = merger.mol_b.bonds.get(bond_key)
    atomtypes_a: tuple[str, str] = tuple([merger.mol_a.atoms[atom_id].type for atom_id in bond_key])  # type: ignore
    atomtypes_b: tuple[str, str] = tuple(
        [
            (
                t
                if (t := getattr(merger.mol_b.atoms[id], "typeB", None))
                else merger.mol_b.atoms[id].type
            )
            for id in bond_key
        ]
    )  # type: ignore
    assert atomtypes_a == ("CT", "H")
    assert atomtypes_b == ("CT", "H1")

    depth_a, beta_a = merger._get_morse_parameters(atomtypes=atomtypes_a, bond=bond_a)
    depth_b, beta_b = merger._get_morse_parameters(atomtypes=atomtypes_b, bond=bond_b)
    type_key: tuple[str, str] = tuple(sorted(atomtypes_b))  # type: ignore
    d = DEFAULT_EDISSOC.get(type_key)
    assert d is None
    assert depth_a == 300
    assert depth_b == 300

    assert bond_b is not None
    k = bond_b.c1
    assert k is not None
    beta = sqrt(float(k) / (2 * depth_b))

    assert beta == pytest.approx(21.9393, rel=1e-5)
    assert beta_a == 19
    assert beta_b == pytest.approx(21.9393, rel=1e-5)


def test_lj_parameters(arranged_tmp_path):
    id1 = "20"
    id2 = "27"
    top_a = Topology(read_top(Path("topol_stateA.top")))
    top_b = Topology(read_top(Path("topol_stateB.top")))
    merger = MoleculeTypeMerger(
        top_a.moleculetypes[REACTIVE_MOLECULEYPE],
        top_b.moleculetypes[REACTIVE_MOLECULEYPE],
        ff=top_a.ff,
    )

    t1 = merger.mol_a.atoms[id1].type
    t2 = merger.mol_a.atoms[id2].type
    assert t1 == "CT"
    assert t2 == "H"
    type1 = merger.ff.atomtypes[t1]
    type2 = merger.ff.atomtypes[t2]
    sigma_a = 0.5 * (float(type1.sigma) + float(type2.sigma))
    epsilon_a = np.sqrt(float(type1.epsilon) * float(type2.epsilon))

    sigma, epsilon = merger._get_LJ_parameters(id1, id2, use_state_b=False)
    assert sigma == sigma_a
    assert epsilon == epsilon_a

    assert sigma == pytest.approx(0.22344, abs=1e-5)
    assert epsilon == pytest.approx(0.1734, abs=1e-5)

    id1 = "17"
    id2 = "27"
    sigma, epsilon = merger._get_LJ_parameters(id1, id2, use_state_b=False)
    assert sigma == pytest.approx(0.21595, abs=1e-5)
    assert epsilon == pytest.approx(0.21615, abs=1e-5)


def test_merge_hat_details(arranged_tmp_path, caplog):
    top_a = Topology(read_top(Path("topol_stateA.top")))
    top_b = Topology(read_top(Path("topol_stateB.top")))
    merger = MoleculeTypeMerger(
        top_a.moleculetypes[REACTIVE_MOLECULEYPE],
        top_b.moleculetypes[REACTIVE_MOLECULEYPE],
        ff=top_a.ff,
    )
    merger.merge()
    assert set([("19", "27")]) == merger.affected_interactions.bonds.added
    assert set([("26", "27")]) == merger.affected_interactions.bonds.removed
    assert (
        set(
            [
                ("24", "27"),
                ("17", "27"),
                ("20", "27"),
            ]
        )
        == merger.affected_interactions.angles.added
    )
    assert (
        set(
            [
                ("24", "27"),
                ("27", "28"),
            ]
        )
        == merger.affected_interactions.angles.removed
    )
    assert (
        set(
            [
                ("17", "27"),
                ("20", "27"),
            ]
        )
        == merger.affected_interactions.angles.added
        - merger.affected_interactions.angles.removed
    )
    assert (
        set(
            [
                ("24", "27"),
            ]
        )
        == merger.affected_interactions.angles.added
        & merger.affected_interactions.angles.removed
    )

    assert (
        set(
            [
                ("15", "27"),
                ("18", "27"),
                ("21", "27"),
                ("22", "27"),
                ("23", "27"),
                ("25", "27"),
                ("26", "27"),
            ]
        )
        == merger.affected_interactions.dihedrals.added
    )

    assert (
        set(
            [
                ("27", "29"),
                ("27", "30"),
                ("27", "34"),
                (
                    "19",
                    "27",
                ),  # this case is interesting, that's the dihedral that is removed, but the atoms are then in a bond
                (
                    "25",
                    "27",
                ),  # this is also interesting, because the dihedral is removed and added (as a different dihedral)
            ]
        )
        == merger.affected_interactions.dihedrals.removed
    )


def test_merge_hat_top(arranged_tmp_path, caplog):
    """Topology merge for a HAT reaction of H 27 from N 26 to C 19"""
    caplog.set_level(logging.INFO)
    top_a = Topology(read_top(Path("topol_stateA.top")))
    top_b = Topology(read_top(Path("topol_stateB.top")))
    top_merge_ref = Topology(read_top(Path("topol_FEP.top")))
    top_merge = merge_top_slow_growth(top_a=top_a, top_b=top_b)

    # for debugging
    # write_top(top_merge.to_dict(), Path("/tmp/kimmdtests_topol_merge.top"))

    assert top_merge.bonds[("19", "27")].funct == "3"
    assert top_merge.bonds[("26", "27")].funct == "3"
    assert top_merge.angles[("17", "19", "20")].c3 is not None

    assert top_merge.proper_dihedrals[("15", "17", "19", "24")].dihedrals["3"].c5 == "3"
    assert (
        top_merge.improper_dihedrals[("17", "20", "19", "24")].dihedrals["2"].c5 == "2"
    )

    assert top_merge.proper_dihedrals == top_merge_ref.proper_dihedrals
    # assert top_merge.improper_dihedrals == top_merge_ref.improper_dihedrals
    assert top_merge.atoms == top_merge_ref.atoms
    assert top_merge.bonds == top_merge_ref.bonds
    assert top_merge.angles == top_merge_ref.angles
    assert top_merge.reactive_molecule.pairs == top_merge_ref.reactive_molecule.pairs
    assert top_merge.pairs == top_merge_ref.pairs
    assert top_merge.exclusions == top_merge_ref.exclusions


def test_merge_small_hat_details(arranged_tmp_path, caplog):
    top = Topology(read_top(Path("AlaCaR_out.top")))
    top_a = copy.deepcopy(top)
    top_b = top

    top_b.break_bond(("7", "8"))
    top_b.bind_bond(("8", "9"))
    merger = MoleculeTypeMerger(
        top_a.moleculetypes[REACTIVE_MOLECULEYPE],
        top_b.moleculetypes[REACTIVE_MOLECULEYPE],
        ff=top_a.ff,
    )
    merger.merge()

    assert merger.affected_interactions.bonds.removed == {("7", "8")}
    assert merger.affected_interactions.bonds.added == {("8", "9")}
    assert merger.affected_interactions.angles.added == {
        ("8", "10"),
        ("7", "8"),
        ("8", "14"),
    }
    assert merger.affected_interactions.angles.removed == {("5", "8"), ("8", "9")}

    assert merger.affected_interactions.dihedrals.added == {
        ("8", "15"),
        ("8", "11"),
        ("8", "12"),
        ("8", "13"),
        ("8", "16"),
        ("5", "8"),
    }

    assert merger.affected_interactions.dihedrals.removed == {
        ("8", "10"),
        ("8", "14"),
        ("6", "8"),
        ("1", "8"),
    }

    added_bonds = (
        merger.affected_interactions.bonds.added
        - merger.affected_interactions.bonds.removed
    )
    removed_bonds = (
        merger.affected_interactions.bonds.removed
        - merger.affected_interactions.bonds.added
    )
    added_angles = (
        merger.affected_interactions.angles.added
        - merger.affected_interactions.angles.removed
    )
    removed_angles = (
        merger.affected_interactions.angles.removed
        - merger.affected_interactions.angles.added
    )
    swapping_angles = (
        merger.affected_interactions.angles.added
        & merger.affected_interactions.angles.removed
    )
    morphing_angles = merger.affected_interactions.angles.morphed
    added_dihedrals = (
        merger.affected_interactions.dihedrals.added
        - merger.affected_interactions.dihedrals.removed
    )
    removed_dihedrals = (
        merger.affected_interactions.dihedrals.removed
        - merger.affected_interactions.dihedrals.added
    )
    swapping_dihedrals = (
        merger.affected_interactions.dihedrals.added
        & merger.affected_interactions.dihedrals.removed
    )
    morphing_dihedrals = merger.affected_interactions.dihedrals.morphed

    assert len(added_bonds) == 1
    assert len(removed_bonds) == 1
    assert len(added_angles) == 3
    assert len(removed_angles) == 2
    assert len(swapping_angles) == 0
    assert len(added_dihedrals) == 6
    assert len(removed_dihedrals) == 4

    assert len(merger.helper_pairs) > 0

    # the atoms of the changing bond still stay in an angle
    # so they keep being excluded
    assert ("8", "9") not in merger.helper_pairs
    assert ("7", "8") not in merger.helper_pairs


def test_merge_small_hat_with_more_overlaps(arranged_tmp_path, caplog):
    top = Topology(read_top(Path("AlaCaR_out.top")))
    # move the radical to N16 before we start
    top.break_bond(("16", "17"))
    top.bind_bond(("9", "17"))

    top_a = copy.deepcopy(top)
    top_b = top

    # now the hat to investigate
    top_b.break_bond(("18", "21"))
    top_b.bind_bond(("16", "21"))

    merger = MoleculeTypeMerger(
        top_a.moleculetypes[REACTIVE_MOLECULEYPE],
        top_b.moleculetypes[REACTIVE_MOLECULEYPE],
        ff=top_a.ff,
    )
    merger.merge()

    assert merger.affected_interactions.bonds.removed == {("18", "21")}
    assert merger.affected_interactions.bonds.added == {("16", "21")}
    assert merger.affected_interactions.angles.added == {
        ("14", "21"),
        ("18", "21"),
    }
    assert merger.affected_interactions.angles.removed == {
        ("19", "21"),
        ("20", "21"),
        ("16", "21"),
    }

    assert merger.affected_interactions.dihedrals.added == {
        ("20", "21"),
        ("19", "21"),
        ("15", "21"),
        ("9", "21"),
    }

    assert merger.affected_interactions.dihedrals.removed == {
        ("14", "21"),
    }

    # 14 21 is interesting, because is get's removed as an angle
    # but added as a dihedral
    # so we expect a Pair transition from exlcusion (all 0.0000) to half-strengh 1-4 interactions
    # sigma and epsilon being for CR and H1
    merger.helper_pairs[("14", "21")] = Pair(
        "14", "21", "1", "0.00000", "0.00000", "0.22344", "0.17340"
    )
