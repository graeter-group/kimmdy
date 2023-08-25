from hypothesis import given, strategies as stace
from kimmdy.recipe import *
import csv
import pytest
from dataclasses import asdict


def test_combine_recipes():
    empty_step = BondOperation(1, 5)
    rp1a = Recipe([empty_step], rates=[1], timespans=[(0.0, 1.0)])
    rp1b = Recipe([empty_step], rates=[1], timespans=[(1.0, 2.0)])
    rp2 = Recipe([empty_step, empty_step], rates=[1], timespans=[(1.0, 3.0)])

    rp1a.combine_with(rp1b)
    assert rp1a.timespans == [(0.0, 1.0), (1.0, 2.0)]
    assert rp1a.rates == [1, 1]

    with pytest.raises(ValueError):
        rp1a.combine_with(rp2)


def test_compare_bond_operations():
    sp1 = BondOperation(1, 5)
    sp2 = BondOperation(2, 5)
    sp21 = BondOperation(atom_id_1=3, atom_id_2="6")
    sp22 = BondOperation(atom_ix_1=2, atom_id_2="6")

    assert sp1 == sp1
    assert sp1 != sp2
    assert sp2 == sp21
    assert sp2 == sp22


# TODO: make this a test of BondOperation
# def test_compare_place():
#     m1 = Place(1)
#     m11 = Place(ix_to_place=1)
#     m12 = Place(id_to_place=2)  #shouldnt this throw an error or something?
#     m13 = Place(id_to_place="2")

#     assert m1 == m11
#     assert m1 == m12
#     assert m1 == m13


# def test_move_class_initialization_mixed():
#     # Initialize with a mix of integers and strings
#     m1 = Place(ix_to_place=5, id_to_bind="4", ix_to_break=2)
#     m2 = Place(id_to_move="6", ix_to_bind=3, id_to_break="3")

#     # Instances should be equal
#     assert m1 == m2


# def test_move_class_initialization_negative():
#     # Initialize with non-matching integers and strings
#     m1 = Place(ix_to_place=5, ix_to_bind=3, ix_to_break=2)
#     m2 = Place(id_to_move="7", id_to_bind="4", id_to_break="3")

#     # Instances should not be equal
#     assert m1 != m2


# def test_move_class_initialization_partial():
#     # Initialize with only some values
#     m1 = Place(ix_to_place=5, ix_to_bind=3)
#     m2 = Place(id_to_move="6", id_to_bind="4")

#     # Instances should be equal
#     assert m1 == m2


# def test_move_class_initialization_partial_negative():
#     # Initialize with only some values, but different values
#     m1 = Place(ix_to_place=5, ix_to_bind=3)
#     m2 = Place(id_to_move="7", id_to_bind="4")

#     # Instances should not be equal
#     assert m1 != m2


# tests with random inputs
# @given(
#     ix_to_place=st.integers(min_value=0, max_value=1000),
#     ix_to_bind=st.integers(min_value=0, max_value=1000),
#     ix_to_break=st.integers(min_value=0, max_value=1000),
# )
# def test_move_class_initialization_integers(ix_to_place, ix_to_bind, ix_to_break):
#     # Initialize with random integers
#     m = Place(ix_to_place=ix_to_place, ix_to_bind=ix_to_bind, ix_to_break=ix_to_break)

#     # Check that properties match input
#     assert m.ix_to_place == ix_to_place
#     assert m.ix_to_bind == ix_to_bind
#     assert m.ix_to_break == ix_to_break
#     assert m.id_to_move == str(ix_to_place + 1)
#     assert m.id_to_bind == str(ix_to_bind + 1)
#     assert m.id_to_break == str(ix_to_break + 1)


# @given(
#     id_to_move=st.integers(min_value=1, max_value=1001).map(str),
#     id_to_bind=st.integers(min_value=1, max_value=1001).map(str),
#     id_to_break=st.integers(min_value=1, max_value=1001).map(str),
# )
# def test_move_class_initialization_strings(id_to_move, id_to_bind, id_to_break):
#     # Initialize with random strings
#     m = Place(id_to_move=id_to_move, id_to_bind=id_to_bind, id_to_break=id_to_break)

#     # Check that properties match input
#     assert m.id_to_move == id_to_move
#     assert m.id_to_bind == id_to_bind
#     assert m.id_to_break == id_to_break
#     assert m.ix_to_place == int(id_to_move) - 1
#     assert m.ix_to_bind == int(id_to_bind) - 1
#     assert m.ix_to_break == int(id_to_break) - 1


# from hypothesis import given, strategies as st


# @given(
#     ix_to_place=st.integers(min_value=0, max_value=1000),
#     ix_to_bind=st.integers(min_value=0, max_value=1000),
#     ix_to_break=st.integers(min_value=0, max_value=1000),
# )
# def test_move_class_initialization_integers(ix_to_place, ix_to_bind, ix_to_break):
#     # Initialize with random integers
#     m1 = Place(ix_to_place=ix_to_place, ix_to_bind=ix_to_bind, ix_to_break=ix_to_break)

#     # Initialize with corresponding strings
#     m2 = Place(
#         id_to_move=ix_to_place + 1,
#         id_to_bind=ix_to_bind + 1,
#         id_to_break=ix_to_break + 1,
#     )

#     # Compare instances
#     assert m1 == m2


# @given(
#     id_to_move=st.integers(min_value=1, max_value=1001).map(str),
#     id_to_bind=st.integers(min_value=1, max_value=1001).map(str),
#     id_to_break=st.integers(min_value=1, max_value=1001).map(str),
# )
# def test_move_class_initialization_strings(id_to_move, id_to_bind, id_to_break):
#     # Initialize with random strings
#     m1 = Place(id_to_move=id_to_move, id_to_bind=id_to_bind, id_to_break=id_to_break)

#     # Initialize with corresponding integers
#     m2 = Place(
#         ix_to_place=int(id_to_move) - 1,
#         ix_to_bind=int(id_to_bind) - 1,
#         ix_to_break=int(id_to_break) - 1,
#     )

#     # Compare instances
#     assert m1 == m2


@pytest.fixture
def recipe_collection():
    rps = [
        Recipe(
            [BondOperation(1, 5), BondOperation(2, 6), BondOperation(3, 7)],
            rates=[1],
            timespans=[(0.0, 1.0)],
        ),
        Recipe([BondOperation(1, 5)], rates=[1], timespans=[(0.0, 1.0)]),
        Recipe([BondOperation(1, 5)], rates=[2], timespans=[(1.0, 2.0)]),
        Recipe(
            [BondOperation(1, 5), BondOperation(1, 5)],
            rates=[1],
            timespans=[(2.0, 3.0)],
        ),
        Recipe([BondOperation(2, 6)], rates=[1], timespans=[(3.0, 4.0)]),
        Recipe(
            [BondOperation(1, 5), BondOperation(1, 5)],
            rates=[3],
            timespans=[(4.0, 5.0)],
        ),
        Recipe(
            [BondOperation(2, 6), BondOperation(3, 7)],
            rates=[1],
            timespans=[(4.0, 5.0)],
        ),
    ]
    return RecipeCollection(rps)


def test_aggregate_recipe_collection(recipe_collection):
    assert len(recipe_collection.recipes) == 7
    recipe_collection.aggregate_reactions()

    assert len(recipe_collection.recipes) == 5
    assert recipe_collection.recipes[1].timespans == [
        (0.0, 1.0),
        (1.0, 2.0),
    ]
    assert recipe_collection.recipes[2].timespans == [(2.0, 3.0), (4.0, 5.0)]


def test_recipe_collection_from_csv(tmp_path, recipe_collection):
    csv_path = tmp_path / "test_out.csv"
    recipe_collection.to_csv(csv_path)
    loaded = RecipeCollection.from_csv(csv_path)[0]
    assert loaded == recipe_collection


def test_recipe_collection_from_csv_picked(tmp_path, recipe_collection):
    csv_path = tmp_path / "test_out.csv"
    picked = recipe_collection.recipes[2]
    recipe_collection.to_csv(csv_path, picked_recipe=picked)
    loaded, loaded_pick = RecipeCollection.from_csv(csv_path)
    assert loaded == recipe_collection
    assert picked == loaded_pick


def test_recipe_collection_to_csv(tmp_path, recipe_collection):
    csv_path = tmp_path / "test_out.csv"
    recipe_collection.to_csv(csv_path)
    with open(csv_path, newline="") as f:
        reader = csv.DictReader(f)
        rows = [r for r in reader]
    for read_rp, org_rp in zip(rows, recipe_collection.recipes):
        org_rp_d = asdict(org_rp)
        keys_to_check = ["rates", "timespans"]
        for key in keys_to_check:
            org_val = str(org_rp_d[key])
            if org_val == "None":
                org_val = ""
            assert org_val == read_rp[key], f"{org_rp_d[key]} != {read_rp[key]}"


def test_recipe_collection_to_csv_picked(tmp_path, recipe_collection):
    csv_path = tmp_path / "test_out.csv"
    picked = recipe_collection.recipes[2]
    recipe_collection.to_csv(csv_path, picked_recipe=picked)
    with open(csv_path, newline="") as f:
        reader = csv.DictReader(f)
        rows = [r for r in reader]
    for read_rp, org_rp in zip(rows, recipe_collection.recipes):
        org_rp_d = asdict(org_rp)
        keys_to_check = ["rates", "timespans"]
        for key in keys_to_check:
            org_val = str(org_rp_d[key])
            if org_val == "None":
                org_val = ""
            assert org_val == read_rp[key], f"{org_rp_d[key]} != {read_rp[key]}"


def test_recipe_collection_to_dill(tmp_path, recipe_collection):
    csv_path = tmp_path / "test_out.dill"
    recipe_collection.to_dill(csv_path)
    loaded_rr = RecipeCollection.from_dill(csv_path)
    assert loaded_rr == recipe_collection
