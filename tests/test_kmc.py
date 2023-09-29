import pytest
from numpy.random import default_rng
from kimmdy.recipe import Recipe, RecipeCollection, Break, Bind
from kimmdy.kmc import rf_kmc, frm, KMCResult


@pytest.fixture
def recipe_collection():
    rps: list[Recipe] = [
        Recipe([Break(1, 2)], rates=[0.0], timespans=[(0.0, 1.0)]),
        Recipe([Bind(2, 3)], rates=[0.12, 0.0], timespans=[(0.0, 6.0), (6.0, 10.0)]),
        Recipe(
            [Break(3, 4), Bind(4, 5)],
            rates=[0.15, 0.06],
            timespans=[(2.0, 4.0), (4.0, 8.0)],
        ),
        Recipe([Break(4, 5), Bind(5, 6)], rates=[1.0], timespans=[(0.0, 0.0)]),
    ]
    return RecipeCollection(rps)


@pytest.fixture
def reference_KMC():
    return KMCResult(
        recipe=Recipe(
            [Bind(2, 3)], rates=[0.12, 0.0], timespans=[(0.0, 6.0), (6.0, 10.0)]
        ),
        time_delta=0.04032167624965666,
        reaction_probability=[0.0, 0.72, 0.54, 0.0],
    )


def compare_to_ref(result: KMCResult, reference: KMCResult):
    assert result.recipe == reference.recipe
    for i in range(len(reference.reaction_probability)):
        assert (
            abs(result.reaction_probability[i] - reference.reaction_probability[i])
            < 1e-9
        )


def test_rf_kmc_empty():
    KMC_dict = rf_kmc(RecipeCollection([]))
    assert KMC_dict.recipe == Recipe([], [], [])


def test_frm_empty():
    KMC_dict = frm(RecipeCollection([]))
    assert KMC_dict.recipe == Recipe([], [], [])


def test_rf_kmc_unlike_ref(reference_KMC):
    rng = default_rng(1)
    # first random numbers are array([0.51182162, 0.9504637])
    KMC_dict = rf_kmc(RecipeCollection([]), rng=rng)
    with pytest.raises(AssertionError):
        compare_to_ref(KMC_dict, reference_KMC)


def test_rf_kmc_calculation(recipe_collection, reference_KMC):
    rng = default_rng(1)
    # first random numbers are array([0.51182162, 0.9504637])
    KMC_dict = rf_kmc(recipe_collection, rng=rng)
    compare_to_ref(KMC_dict, reference_KMC)
    assert abs(KMC_dict.time_delta - reference_KMC.time_delta) < 1e-9


def test_frm_calculation(recipe_collection, reference_KMC):
    rng = default_rng(1)
    # first random numbers are array([0.51182162, 0.9504637 , 0.14415961, 0.94864945])
    # with the seed 0, the second reaction is picked, so that the time step for rfKMC and FRM are equal
    # because the same random value is used to determine the time step
    KMC_dict = frm(recipe_collection, rng=rng, MD_time=10)
    compare_to_ref(KMC_dict, reference_KMC)


def test_frm_no_event(recipe_collection):
    rng = default_rng(1)
    # first random numbers are array([0.51182162, 0.9504637 , 0.14415961, 0.94864945])
    new_recipes = recipe_collection.recipes[2:4]
    KMC_dict = frm(RecipeCollection(new_recipes), rng=rng, MD_time=None)
    assert KMC_dict.recipe == Recipe([], [], [])
