import pytest
import numpy as np
from numpy.random import default_rng
from kimmdy.recipe import Recipe, RecipeCollection, Break, Bind
from kimmdy.kmc import (
    KMCError,
    KMCReject,
    rf_kmc,
    frm,
    extrande,
    extrande_mod,
    KMCAccept,
    total_index_to_index_within_plugin,
)


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
def reference_KMC() -> KMCAccept:
    return KMCAccept(
        recipe=Recipe(
            [Bind(2, 3)], rates=[0.12, 0.0], timespans=[(0.0, 6.0), (6.0, 10.0)]
        ),
        time_delta=0.04032167624965666,
        reaction_probability=[0.0, 0.72, 0.54, 0.0],
        time_start=0,
        time_start_index=0,
    )


@pytest.fixture
def reference_extrande_KMC() -> KMCAccept:
    return KMCAccept(
        recipe=Recipe([Bind(2, 3)], rates=[0.12], timespans=[(0.0, 6.0)]),
        time_delta=0,
        time_start=3.85725338647224,
        time_start_index=0,
        reaction_probability=None,
    )


def compare_to_ref(result: KMCAccept, reference: KMCAccept):
    assert isinstance(result, KMCAccept)
    assert isinstance(reference, KMCAccept)
    assert result.recipe == reference.recipe
    assert reference.reaction_probability, "Reference reaction probability is None"
    assert reference.reaction_probability, "Reference reaction probability is None"
    assert result.reaction_probability, "Result reaction probability is None"
    for i in range(len(reference.reaction_probability)):
        assert (
            abs(result.reaction_probability[i] - reference.reaction_probability[i])
            < 1e-9
        )


def test_rf_kmc_empty():
    kmc = rf_kmc(RecipeCollection([]))
    assert isinstance(kmc, KMCError)


def test_frm_empty():
    kmc = frm(RecipeCollection([]))
    assert isinstance(kmc, KMCError)


def test_extrande_empty():
    kmc = extrande(RecipeCollection([]), 1.0)
    assert isinstance(kmc, KMCError)


def test_extrande_mod_empty():
    kmc = extrande_mod(RecipeCollection([]), 1.0)
    assert isinstance(kmc, KMCError)


def test_rf_kmc_unlike_ref(reference_KMC):
    rng = default_rng(1)
    # first random numbers are array([0.51182162, 0.9504637])
    kmc = rf_kmc(RecipeCollection([]), rng=rng)
    assert isinstance(kmc, KMCError)
    with pytest.raises(AssertionError):
        compare_to_ref(kmc, reference_KMC)


def test_extrande_calculation(recipe_collection, reference_extrande_KMC):
    rng = default_rng(1)
    # first random numbers are array([0.51182162, 0.9504637])
    kmc = extrande(recipe_collection, 1.0, rng=rng)
    assert isinstance(kmc, KMCAccept)
    assert kmc.recipe == reference_extrande_KMC.recipe
    assert abs(kmc.time_start - reference_extrande_KMC.time_start) < 1e-9
    assert abs(kmc.time_delta - reference_extrande_KMC.time_delta) < 1e-9


def test_extrande_mod_calculation(recipe_collection, reference_extrande_KMC):
    rng = default_rng(1)
    # first random numbers are array([0.51182162, 0.9504637])
    kmc = extrande_mod(recipe_collection, 1.0, rng=rng)
    assert isinstance(kmc, KMCAccept)
    assert kmc.recipe == reference_extrande_KMC.recipe
    assert kmc.time_start is not None
    assert (
        abs(kmc.time_start - 2.4632908726674225) < 1e-9
    )  # single result different, same distribution?
    assert abs(kmc.time_delta - reference_extrande_KMC.time_delta) < 1e-9


def test_rf_kmc_calculation(recipe_collection, reference_KMC):
    rng = default_rng(1)
    # first random numbers are array([0.51182162, 0.9504637])
    kmc = rf_kmc(recipe_collection, rng=rng)
    assert isinstance(kmc, KMCAccept)
    compare_to_ref(kmc, reference_KMC)
    assert abs(kmc.time_delta - reference_KMC.time_delta) < 1e-9


def test_frm_calculation(recipe_collection, reference_KMC):
    rng = default_rng(1)
    # first random numbers are array([0.51182162, 0.9504637 , 0.14415961, 0.94864945])
    # with the seed 0, the second reaction is picked, so that the time step for rfKMC and FRM are equal
    # because the same random value is used to determine the time step
    kmc = frm(recipe_collection, rng=rng, MD_time=10)
    assert isinstance(kmc, KMCAccept)
    compare_to_ref(kmc, reference_KMC)


def test_frm_no_event(recipe_collection):
    rng = default_rng(1)
    # first random numbers are array([0.51182162, 0.9504637 , 0.14415961, 0.94864945])
    new_recipes = recipe_collection.recipes[2:4]
    kmc = frm(RecipeCollection(new_recipes), rng=rng, MD_time=None)
    assert isinstance(kmc, KMCError)


def test_compare_extrande_extrande_mod(recipe_collection):
    rng = default_rng(1)
    extrande_list = [extrande(recipe_collection, 1.0, rng=rng) for _ in range(2000)]
    extrande_mod_list = [
        extrande_mod(recipe_collection, 1.0, rng=rng) for _ in range(2000)
    ]
    ext_ts = np.array([r.time_start for r in extrande_list if isinstance(r, KMCAccept)])
    extmod_ts = np.array(
        [r.time_start for r in extrande_mod_list if isinstance(r, KMCAccept)]
    )
    mask = np.nonzero(ext_ts != np.array(None))[0]
    mmask = np.nonzero(extmod_ts != np.array(None))[0]

    assert abs(ext_ts[mask].mean() - extmod_ts[mmask].mean()) < 0.1

    ext_rs = np.concatenate(
        [r.recipe.rates for r in extrande_list if isinstance(r, KMCAccept)]
    )
    extmod_rs = np.concatenate(
        [r.recipe.rates for r in extrande_mod_list if isinstance(r, KMCAccept)]
    )

    assert abs(ext_rs.mean() - extmod_rs.mean()) < 0.002


def test_total_index_to_index_within_plugin():
    ns = [3, 2, 4, 1]
    assert total_index_to_index_within_plugin(0, ns) == 0
    assert total_index_to_index_within_plugin(1, ns) == 1
    assert total_index_to_index_within_plugin(2, ns) == 2
    assert total_index_to_index_within_plugin(3, ns) == 0
    assert total_index_to_index_within_plugin(4, ns) == 1
    assert total_index_to_index_within_plugin(5, ns) == 0
