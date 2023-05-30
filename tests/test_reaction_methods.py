from kimmdy.reaction import *
import csv
import pytest
from dataclasses import asdict


def test_no_double_times():
    with pytest.raises(ValueError):
        Recipe([RecipeStep()], rates=[1.0, 2.0, 3.0], times=[1.0, 1.0, 2.0])


def test_combine_recipes():
    rp1a = Recipe([RecipeStep()], rates=[1], times=[1.0])
    rp1b = Recipe([RecipeStep()], rates=[1], times=[2.0])
    rp2 = Recipe([RecipeStep(), RecipeStep()], rates=[1], times=[3.0])

    rp1a.combine_with(rp1b)
    assert rp1a.times == [1, 2]
    assert rp1a.rates == [1, 1]

    with pytest.raises(ValueError):
        rp1a.combine_with(rp2)


@pytest.fixture
def recipe_collection():
    rps = [
        Recipe([RecipeStep(), RecipeStep(), RecipeStep()], rates=[1], times=[0.0]),
        Recipe([RecipeStep()], rates=[1], times=[1.0]),
        Recipe([RecipeStep()], rates=[1], times=[2.0]),
        Recipe([RecipeStep(), RecipeStep()], rates=[1], times=[3.0]),
        Recipe([RecipeStep()], rates=[1], times=[4.0]),
        Recipe([RecipeStep(), RecipeStep()], rates=[1], times=[5.0]),
    ]
    return RecipeCollection(rps)


def test_aggregate_recipe_collection(recipe_collection):
    recipe_collection.aggregate_reactions()

    assert len(recipe_collection.recipes) == 3
    assert recipe_collection.recipes[1].times == [1.0, 2.0, 4.0]
    assert recipe_collection.recipes[2].times == [3.0, 5.0]


def test_recipe_collection_to_csv(tmp_path, recipe_collection):
    csv_path = tmp_path / "test_out.csv"
    recipe_collection.to_csv(csv_path)
    with open(csv_path, newline="") as f:
        reader = csv.DictReader(f)
        rows = [r for r in reader]
    for read_rp, org_rp in zip(rows, recipe_collection.recipes):
        org_rp_d = asdict(org_rp)
        keys_to_check = ["rates", "times", "avg_rates", "avg_timespans"]
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
