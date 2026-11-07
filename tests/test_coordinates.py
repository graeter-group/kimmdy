from types import SimpleNamespace
import copy
import multiprocessing as mp
import logging
import numpy as np
from math import sqrt
import pytest
from pathlib import Path
import MDAnalysis as mda
from kimmdy.recipe import Bind, Break, Place
from kimmdy.parsing import read_top, write_top
from kimmdy.constants import REACTIVE_MOLECULEYPE
from conftest import DummyFiles
from kimmdy.topology.atomic import Angle, Atom, Bond, BondType, Dihedral, Pair
from kimmdy.topology.topology import Topology
from kimmdy.coordinates import (
    AffectedInteractions,
    MoleculeTypeMerger,
    merge_top_slow_growth,
    get_explicit_or_type,
    place_atom,
)


def test_merge_harmonic_to_urey_bradley_angle():
    merger = object.__new__(MoleculeTypeMerger)
    key = ("1", "2", "3")
    harmonic = Angle(*key, funct="1", c0="109.5", c1="300.0")
    urey_bradley = Angle(
        *key,
        funct="5",
        c0="112.0",
        c1="450.0",
        c2="0.245",
        c3="5000.0",
    )

    merged = merger._merge_urey_bradley_angles(key, harmonic, urey_bradley)

    assert merged.funct == "5"
    assert (merged.c0, merged.c1, merged.c2, merged.c3) == (
        "109.5",
        "300.0",
        "0.245",
        "0.00",
    )
    assert (merged.c4, merged.c5, merged.c6, merged.c7) == (
        "112.0",
        "450.0",
        "0.245",
        "5000.0",
    )


def test_implicit_interactions_ignore_generated_pairs():
    """A pair without numbers blocks the transferred-H retyping only when the force field does not generate pairs."""
    merger = object.__new__(MoleculeTypeMerger)
    merger.mol_a = SimpleNamespace(
        bonds={("1", "2"): Bond("1", "2", funct="1", c0="0.0960", c1="462750.0")},
        angles={},
        proper_dihedrals={},
        improper_dihedrals={},
        pairs={("2", "5"): Pair("2", "5", funct="1")},
    )

    # gen-pairs = yes: grompp builds the pair from the atom types, retyping is safe
    merger.ff = SimpleNamespace(defaults=[["1", "2", "yes", "0.5", "0.8333"]])
    assert merger._implicit_interactions("2") == []

    # gen-pairs = no (or missing): the pair needs a [ pairtypes ] entry, keep blocking
    merger.ff = SimpleNamespace(defaults=[["1", "2", "no", "0.5", "0.8333"]])
    assert merger._implicit_interactions("2") == ["pair 2-5"]
    merger.ff = SimpleNamespace(defaults=[["1", "2"]])
    assert merger._implicit_interactions("2") == ["pair 2-5"]

    # a bond without numbers still blocks, whatever gen-pairs says
    merger.ff = SimpleNamespace(defaults=[["1", "2", "yes", "0.5", "0.8333"]])
    merger.mol_a.bonds[("1", "2")] = Bond("1", "2", funct="1")
    assert merger._implicit_interactions("2") == ["bond 1-2"]


@pytest.mark.parametrize(
    "type_a, type_b, expected",
    [
        ("HO", "HC", "HC"),  # the B type has an LJ size: kept
        ("HO", "HO", "H1"),  # lands on an oxygen again: no LJ size in B
        ("HC", "HO", "H1"),  # moves onto an oxygen: no LJ size in B
    ],
)
def test_retype_transferred_h_uses_b_type(type_a, type_b, expected):
    """State B of the relax takes the B type of a moved H, H1 only where the B type has no LJ size."""
    merger = object.__new__(MoleculeTypeMerger)
    merger.retype_transferred_h = True
    merger.transferred_h_type = "H1"
    merger.ff = SimpleNamespace(
        atomtypes={
            "HO": SimpleNamespace(epsilon="0.00000e+00"),
            "HC": SimpleNamespace(epsilon="6.56888e-02"),
            "H1": SimpleNamespace(epsilon="6.56888e-02"),
        },
        defaults=[["1", "2", "yes", "0.5", "0.8333"]],
    )
    # H 2 moves from atom 1 to atom 3; merge_atoms gave it its B type
    h_a = Atom("2", type_a, "1", "RES", "H", "2", "0.4", "1.008")
    if type_b != type_a:
        h_a.typeB = type_b
    h_b = Atom("2", type_b, "1", "RES", "H", "2", "0.4", "1.008")
    merger.mol_a = SimpleNamespace(
        atoms={"2": h_a},
        bonds={
            ("1", "2"): Bond("1", "2", funct="1", c0="0.0960", c1="462750.0"),
            ("2", "3"): Bond("2", "3", funct="1", c0="0.1090", c1="284512.0"),
        },
        angles={},
        proper_dihedrals={},
        improper_dihedrals={},
        pairs={},
    )
    merger.mol_b = SimpleNamespace(
        atoms={"2": h_b},
        bonds={("2", "3"): Bond("2", "3", funct="1", c0="0.1090", c1="284512.0")},
    )
    merger.affected_interactions = AffectedInteractions()

    merger._retype_transferred_hydrogens()

    assert h_a.type == type_a  # state A keeps its type
    assert (h_a.typeB or h_a.type) == expected


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

    depth_a, beta_a = merger._get_morse_parameters(bond=bond_a)
    depth_b, beta_b = merger._get_morse_parameters(bond=bond_b)
    # always the default well depth
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
    top_a_ref = Topology(read_top(Path("topol_FEP.top")))
    merge_top_slow_growth(top_a=top_a, top_b=top_b)

    # for debugging
    # write_top(top_a.to_dict(), Path("/tmp/kimmdtests_topol_merge.top"))

    assert top_a.bonds[("19", "27")].funct == "3"
    assert top_a.bonds[("26", "27")].funct == "3"
    assert top_a.angles[("17", "19", "20")].c3 is not None

    assert top_a.proper_dihedrals[("15", "17", "19", "24")].dihedrals["3"].c5 == "3"
    assert top_a.improper_dihedrals[("17", "20", "19", "24")].dihedrals["2"].c5 == "2"

    assert top_a.proper_dihedrals == top_a_ref.proper_dihedrals
    # assert top_a.improper_dihedrals == top_a_ref.improper_dihedrals
    assert top_a.atoms == top_a_ref.atoms
    assert top_a.bonds == top_a_ref.bonds
    assert top_a.angles == top_a_ref.angles
    assert top_a.reactive_molecule.pairs == top_a_ref.reactive_molecule.pairs
    assert top_a.pairs == top_a_ref.pairs
    assert top_a.exclusions == top_a_ref.exclusions


@pytest.mark.parametrize("arranged_tmp_path", ["test_coordinates"], indirect=True)
def test_merge_uses_b_charges_in_both_states(arranged_tmp_path):
    """The relaxation topology carries the charges of B; only types and masses are interpolated."""
    top_a = Topology(read_top(Path("topol_stateA.top")))
    top_b = Topology(read_top(Path("topol_stateB.top")))
    merge_top_slow_growth(top_a=top_a, top_b=top_b)

    for nr, atom in top_a.atoms.items():
        assert atom.charge == top_b.atoms[nr].charge
        assert atom.chargeB in (None, atom.charge)
    assert top_a.atoms["26"].typeB is None  # only the charge differed
    assert top_a.atoms["27"].typeB == "H1"  # the moved hydrogen still changes type


@pytest.mark.parametrize("arranged_tmp_path", ["test_coordinates"], indirect=True)
def test_merge_morse_only_writes_b_types(arranged_tmp_path):
    """morse_only gives an atom its B type even if its charge did not change."""
    top_a = Topology(read_top(Path("topol_stateA.top")))
    top_b = Topology(read_top(Path("topol_stateB.top")))
    merge_top_slow_growth(top_a=top_a, top_b=top_b, morse_only=True)

    # the moved H 27 has the same charge in A and B, but goes from H to H1
    assert top_a.atoms["27"].type == "H"
    assert top_a.atoms["27"].typeB == "H1"


def _hat(top: Topology, donor: str, hydrogen: str, acceptor: str) -> None:
    """Move a hydrogen from donor to acceptor as a recipe segment does: bonds, then names and types."""
    top_before = copy.deepcopy(top)
    top.break_bond((donor, hydrogen))
    top.bind_bond((hydrogen, acceptor))
    top.update_names_LJ(
        top_before,
        [
            Break(atom_id_1=donor, atom_id_2=hydrogen),
            Bind(atom_id_1=hydrogen, atom_id_2=acceptor),
        ],
    )


def test_hat_keeps_impropers_of_acceptor_residue(arranged_tmp_path):
    """An H that moves into a residue with an atom of its name keeps that atom's impropers.

    H 17 (H of NME) moves to CA 9 of ALA, which has its own H 8. bind_bond
    builds the impropers of CA 9 while 17 is still named H; once
    update_names_LJ names it HA, the amide improper of ALA is (5 9 7 8) again.
    """
    top = Topology(read_top(Path("AlaCaR_out.top")))
    assert ("5", "9", "7", "8") in top.improper_dihedrals

    _hat(top, "16", "17", "9")

    assert top.atoms["17"].atom == "HA"
    assert ("5", "9", "7", "8") in top.improper_dihedrals
    assert ("5", "9", "7", "17") not in top.improper_dihedrals


def test_merge_small_hat_details(arranged_tmp_path, caplog):
    top = Topology(read_top(Path("AlaCaR_out.top")))
    top_a = copy.deepcopy(top)
    top_b = top

    _hat(top_b, "7", "8", "9")
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


def test_merge_restores_pairs_excluded_only_by_union_graph(arranged_tmp_path):
    top = Topology(read_top(Path("AlaCaR_out.top")))
    # Split the peptide so donor and acceptor are physically disconnected in
    # both endpoints, like a lipid donating H to a hydroxyl radical.
    top.break_bond(("14", "16"))
    top_a = copy.deepcopy(top)
    top_b = copy.deepcopy(top)

    top_b.break_bond(("7", "8"))
    top_b.bind_bond(("8", "16"))

    merger = MoleculeTypeMerger(
        top_a.moleculetypes[REACTIVE_MOLECULEYPE],
        top_b.moleculetypes[REACTIVE_MOLECULEYPE],
        ff=top_a.ff,
    )
    merger.merge()

    # These pairs are fully nonbonded in A and B, but the union graph contains
    # the artificial paths N7-H8-N16 and C5-N7-H8-N16.
    expected_restored_pairs = {("7", "16"), ("5", "16")}
    assert expected_restored_pairs <= merger.union_exclusion_helper_pairs
    for pair_key in expected_restored_pairs:
        pair = merger.mol_a.pairs[pair_key]
        assert pair.c0 is not None
        assert pair.c1 is not None
        assert pair.c2 is not None
        assert pair.c3 is not None
        assert pair_key in merger.mol_a.exclusions


def test_union_exclusion_reconciliation_keeps_connected_hat_behavior(
    arranged_tmp_path,
):
    top = Topology(read_top(Path("AlaCaR_out.top")))
    top_a = copy.deepcopy(top)
    top_b = copy.deepcopy(top)

    _hat(top_b, "7", "8", "9")

    merger = MoleculeTypeMerger(
        top_a.moleculetypes[REACTIVE_MOLECULEYPE],
        top_b.moleculetypes[REACTIVE_MOLECULEYPE],
        ff=top_a.ff,
    )
    merger.merge()

    assert merger.union_exclusion_helper_pairs == set()


def test_merge_small_hat_with_more_overlaps(arranged_tmp_path, caplog):
    top = Topology(read_top(Path("AlaCaR_out.top")))
    # move the radical to N16 before we start
    _hat(top, "16", "17", "9")

    top_a = copy.deepcopy(top)
    top_b = top

    # now the hat to investigate
    _hat(top_b, "18", "21", "16")

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


def test_merge_small_hat_harmonic_acceptor_h_term(arranged_tmp_path):
    """The optional relax helper adds one growing type-6 H-acceptor bond and nothing else."""

    def merged_dict(helper: bool) -> dict:
        """Move the H from atom 7 to atom 9, build the relax topology with or without the spring, and return it as text sections."""
        top = Topology(read_top(Path("AlaCaR_out.top")))
        top_a = copy.deepcopy(top)
        _hat(top, "7", "8", "9")
        merge_top_slow_growth(top_a, top, harmonic_acceptor_h_relax_term=helper)
        return top_a.to_dict()

    plain = merged_dict(False)
    helped = merged_dict(True)
    bonds_plain = plain["moleculetype_Reactive"]["subsections"]["bonds"]["content"]
    bonds_helped = helped["moleculetype_Reactive"]["subsections"]["bonds"]["content"]
    extra = [line for line in bonds_helped if line not in bonds_plain]
    assert len(bonds_helped) == len(bonds_plain) + 1
    assert len(extra) == 1
    ai, aj, funct, b0_a, k_a, b0_b, k_b = extra[0][:7]
    assert (ai, aj, funct) == ("8", "9", "6")
    assert float(k_a) == 0.0 and float(k_b) == 50000.0
    assert b0_a == b0_b
    for section in ("pairs", "exclusions", "angles", "atoms"):
        assert (
            plain["moleculetype_Reactive"]["subsections"][section]
            == helped["moleculetype_Reactive"]["subsections"][section]
        )


def _hat_states(case: str) -> tuple[Topology, Topology]:
    """Return the topologies before (A) and after (B) one HAT reaction of the test files.

    small: H 8 moves from atom 7 to atom 9 in AlaCaR_out.top.
    large: H 27 moves from N 26 to C 19 (topol_stateA.top -> topol_stateB.top).
    """
    if case == "small":
        top_b = Topology(read_top(Path("AlaCaR_out.top")))
        top_a = copy.deepcopy(top_b)
        _hat(top_b, "7", "8", "9")
        return top_a, top_b
    return Topology(read_top(Path("topol_stateA.top"))), Topology(read_top(Path("topol_stateB.top")))


@pytest.mark.skipif(
    "fork" not in mp.get_all_start_methods(),
    reason="parallel dihedral merging needs fork helper processes",
)
@pytest.mark.parametrize("case", ["small", "large"])
def test_parallel_dihedral_merge_matches_serial(arranged_tmp_path, case):
    """Merging the dihedrals on 4 helper processes gives the same relax topology as on one CPU.

    Compares the affected interactions (which steer the later pair handling),
    every merged proper and improper dihedral with all its values, the pairs
    and exclusions, and the written [ dihedrals ] lines including their order
    (both merges go through the dihedrals sorted by atom numbers).
    """

    def merged(parallel: bool):
        """Merge the reaction once, serially or in parallel; return the merger and the written dihedral lines."""
        top_a, top_b = _hat_states(case)
        merger = MoleculeTypeMerger(
            top_a.moleculetypes[REACTIVE_MOLECULEYPE],
            top_b.moleculetypes[REACTIVE_MOLECULEYPE],
            ff=top_a.ff,
            parallel_dihedral_merge=parallel,
            parallel_dihedral_workers=4,
        )
        merger.merge()

        # the written relax topology, from a second identical merge
        top_a, top_b = _hat_states(case)
        merge_top_slow_growth(
            top_a=top_a,
            top_b=top_b,
            parallel_dihedral_merge=parallel,
            parallel_dihedral_workers=4,
        )
        lines = top_a.to_dict()["moleculetype_Reactive"]["subsections"]["dihedrals"]["content"]
        return merger, [tuple(line) for line in lines]

    serial, serial_lines = merged(parallel=False)
    parallel, parallel_lines = merged(parallel=True)

    # the interactions the reaction touched
    assert parallel.affected_interactions == serial.affected_interactions
    assert serial.affected_interactions.dihedrals.added or serial.affected_interactions.dihedrals.removed

    # every merged dihedral with all its periodicities and values
    assert parallel.mol_a.proper_dihedrals == serial.mol_a.proper_dihedrals
    assert parallel.mol_a.improper_dihedrals == serial.mol_a.improper_dihedrals

    # pairs and exclusions, which follow from the affected interactions
    assert parallel.mol_a.pairs == serial.mol_a.pairs
    assert parallel.mol_a.exclusions == serial.mol_a.exclusions

    # the same written lines, in the same order
    assert parallel_lines == serial_lines
