import csv
from dataclasses import asdict

import pytest
from hypothesis import given
from hypothesis import strategies as st

from kimmdy import recipe


## Test RecipeSteps
# tests with random inputs
@given(
    ix_1=st.integers(min_value=0, max_value=1000),
    ix_2=st.integers(min_value=0, max_value=1000),
)
def test_bo_initialization_integers(ix_1, ix_2):
    # Initialize with random integers
    m = recipe.BondOperation(atom_ix_1=ix_1, atom_ix_2=ix_2)

    # Check that properties match input
    assert m.atom_ix_1 == ix_1
    assert m.atom_ix_2 == ix_2
    assert m.atom_id_1 == str(ix_1 + 1)
    assert m.atom_id_2 == str(ix_2 + 1)


@given(
    id_1=st.integers(min_value=1, max_value=1001).map(str),
    id_2=st.integers(min_value=1, max_value=1001).map(str),
)
def test_bo_initialization_strings(id_1, id_2):
    # Initialize with random strings
    m = recipe.BondOperation(atom_id_1=id_1, atom_id_2=id_2)

    # Check that properties match input
    assert m.atom_id_1 == id_1
    assert m.atom_id_2 == id_2
    assert m.atom_ix_1 == int(id_1) - 1
    assert m.atom_ix_2 == int(id_2) - 1


def test_bo_initialization_mixed():
    # Initialize with a mix of integers and strings
    m1 = recipe.BondOperation(atom_ix_1=5, atom_id_2="4")
    m2 = recipe.BondOperation(
        atom_id_1="6",
        atom_ix_2=3,
    )

    # Instances should be equal
    assert m1 == m2


@given(
    ix_1=st.integers(min_value=0, max_value=1000),
    ix_2=st.integers(min_value=0, max_value=1000),
)
def test_bo_initialization_separate(ix_1, ix_2):
    # Initialize with random integers
    m1 = recipe.BondOperation(atom_ix_1=ix_1, atom_ix_2=ix_2)

    # Initialize with corresponding strings
    m2 = recipe.BondOperation(
        atom_id_1=str(ix_1 + 1),
        atom_id_2=str(ix_2 + 1),
    )

    # Compare instances
    assert m1 == m2


def test_bo_initialization_unequal():
    # Initialize with non-matching integers and strings
    m1 = recipe.BondOperation(5, 3)
    m2 = recipe.BondOperation(1, 2)

    # Instances should not be equal
    assert m1 != m2


def test_bo_initialization_wrong_type():
    # Should raise an error because initialization is with the wrong type
    with pytest.raises(ValueError):
        recipe.BondOperation("1", "2")  # type: ignore
    with pytest.raises(ValueError):
        recipe.BondOperation(atom_id_1=0, atom_id_2=1)  # type: ignore


@given(
    ix_1=st.integers(min_value=0, max_value=1000),
    ix_2=st.integers(min_value=0, max_value=1000),
)
def test_bind_like_bo(ix_1, ix_2):
    # Initialize with random integers
    m = recipe.Bind(atom_ix_1=ix_1, atom_ix_2=ix_2)

    # Check that properties match input
    assert m.atom_ix_1 == ix_1
    assert m.atom_ix_2 == ix_2
    assert m.atom_id_1 == str(ix_1 + 1)
    assert m.atom_id_2 == str(ix_2 + 1)


@given(
    id_1=st.integers(min_value=1, max_value=1001).map(str),
    id_2=st.integers(min_value=1, max_value=1001).map(str),
)
def test_break_like_bo(id_1, id_2):
    # Initialize with random strings
    m = recipe.Break(atom_id_1=id_1, atom_id_2=id_2)

    # Check that properties match input
    assert m.atom_id_1 == id_1
    assert m.atom_id_2 == id_2
    assert m.atom_ix_1 == int(id_1) - 1
    assert m.atom_ix_2 == int(id_2) - 1


def test_place_initialization():
    m1 = recipe.Place(ix_to_place=1, new_coords=(0, 0, 0))
    m11 = recipe.Place(id_to_place="2", new_coords=(0, 0, 0))
    m2 = recipe.Place(ix_to_place=2, new_coords=(0, 0, 0))
    m3 = recipe.Place(ix_to_place=2, new_coords=(0, 1, 0))

    assert m1 == m11
    assert m1 != m2
    assert m2 != m3

    with pytest.raises(TypeError):
        recipe.Place(ix_to_place=1)  # type: ignore
    with pytest.raises(ValueError):
        recipe.Place(id_to_place=1, new_coords=(0, 0, 0))  # type: ignore


def test_relax_initialization():
    recipe.Relax()


def test_relax_comparison():
    r1 = recipe.Relax()
    r2 = recipe.Relax()
    assert r1 == r2


## Test Recipe
def test_combine_recipes():
    empty_step = recipe.BondOperation(1, 5)
    rp1a = recipe.Recipe([empty_step], rates=[1], timespans=[(0.0, 1.0)])
    rp1b = recipe.Recipe([empty_step], rates=[1], timespans=[(1.0, 2.0)])
    rp2 = recipe.Recipe([empty_step, empty_step], rates=[1], timespans=[(1.0, 3.0)])

    rp1a.combine_with(rp1b)
    assert rp1a.timespans == [(0.0, 1.0), (1.0, 2.0)]
    assert rp1a.rates == [1, 1]

    with pytest.raises(ValueError):
        rp1a.combine_with(rp2)


# ToDo: test other recipe methods


## Test RecipeCollection
@pytest.fixture
def recipe_collection():
    rps = [
        recipe.Recipe(
            [
                recipe.BondOperation(1, 5),
                recipe.BondOperation(2, 6),
                recipe.BondOperation(3, 7),
            ],
            rates=[1],
            timespans=[(0.0, 1.0)],
        ),
        recipe.Recipe([recipe.BondOperation(1, 5)], rates=[1], timespans=[(0.0, 1.0)]),
        recipe.Recipe([recipe.BondOperation(1, 5)], rates=[2], timespans=[(1.0, 2.0)]),
        recipe.Recipe(
            [recipe.BondOperation(1, 5), recipe.BondOperation(1, 5)],
            rates=[1],
            timespans=[(2.0, 3.0)],
        ),
        recipe.Recipe([recipe.BondOperation(2, 6)], rates=[1], timespans=[(3.0, 4.0)]),
        recipe.Recipe(
            [recipe.BondOperation(1, 5), recipe.BondOperation(1, 5)],
            rates=[3],
            timespans=[(4.0, 5.0)],
        ),
        recipe.Recipe(
            [recipe.BondOperation(2, 6), recipe.BondOperation(3, 7)],
            rates=[1],
            timespans=[(4.0, 5.0)],
        ),
    ]
    return recipe.RecipeCollection(rps)


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
    loaded = recipe.RecipeCollection.from_csv(csv_path)[0]
    assert loaded == recipe_collection


def test_recipe_collection_from_csv_picked(tmp_path, recipe_collection):
    csv_path = tmp_path / "test_out.csv"
    picked = recipe_collection.recipes[2]
    recipe_collection.to_csv(csv_path, picked_recipe=picked)
    loaded, loaded_pick = recipe.RecipeCollection.from_csv(csv_path)
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
    loaded_rr = recipe.RecipeCollection.from_dill(csv_path)
    assert loaded_rr == recipe_collection
