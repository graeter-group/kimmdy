from typing import Union

import logging
import numpy as np
from scipy.interpolate import interp1d
import scipy
from numpy.random import default_rng
from kimmdy.reaction import RecipeCollection, Recipe


def rfKMC(
    recipe_collection: RecipeCollection, rng: np.random.BitGenerator = default_rng()
) -> Recipe:
    """Rejection-Free Monte Carlo.
    takes RecipeCollection and choses a recipe based on the relative propensity of the events.

    Parameters
    ---------
    reaction_result: ReactionResult
        from which one will be choosen
    rng: np.random.BitGenerator
        to generate random numbers in the KMC step
    """
    # compare e.g. <https://en.wikipedia.org/wiki/Kinetic_Monte_Carlo#Rejection-free_KMC>

    # check for empty ReactionResult
    if len(recipe_collection.recipes) == 0:
        logging.warning("Empty ReactionResult; no reaction chosen")
        return Recipe(), 0, []

    # 0. Initialization
    constant_rates = []
    recipes_steps = []
    for recipe in recipe_collection.recipes:
        # 1.1 Calculate the propensity function for each reaction
        dt = [x[1] - x[0] for x in recipe.timespans]
        constant_rates.append(sum(map(lambda x, y: x * y, dt, recipe.rates)))
        recipes_steps.append(recipe.recipe_steps)
    # 2. Set the total propensity to the sum of individual propensities
    rates_cumulative = np.cumsum(constant_rates)
    total_rate = rates_cumulative[-1]
    # 3. Generate two independent uniform (0,1) random numbers u1,u2
    u = rng.random(2)
    logging.debug(
        f"Random values u: {u}, cumulative rates {rates_cumulative}, total rate {total_rate}"
    )

    # 4. Find the even to carry out, mu, using binary search (np.searchsorted)
    pos = np.searchsorted(rates_cumulative, u[0] * total_rate)
    chosen_recipe = recipes_steps[pos]
    logging.debug(f"Chosen Recipe: {chosen_recipe}")
    # 5. Calculate the time step associated with mu
    dt = 1 / total_rate * np.log(1 / u[1])

    return chosen_recipe, dt, constant_rates


def FRM(
    recipe_collection: RecipeCollection,
    rng: np.random.BitGenerator = default_rng(),
    MD_time: Union[float, None] = None,
) -> Recipe:
    """First Reaction Method variant of Kinetic Monte Carlo.
    takes RecipeCollection and choses a recipe based on which reaction would occur.

    Parameters
    ---------
    reaction_result: ReactionResult
        from which one will be choosen
    rng: np.random.BitGenerator
        to generate random numbers in the KMC step
    MD_time: Union[float, None] [ps]
        to compare conformational events with reaction events in the time domain
    """
    # compare e.g. <https://en.wikipedia.org/wiki/Kinetic_Monte_Carlo#Time-dependent_Algorithms>

    # check for empty ReactionResult
    if len(recipe_collection.recipes) == 0:
        logging.warning("Empty ReactionResult; no reaction chosen")
        return Recipe(), 0, []

    # 0. Initialization
    start_time = recipe_collection.recipes[0].times[0]
    end_time = recipe_collection.recipes[0].times[-1]
    resolution = 50  # [1/ps]
    samples = np.linspace(
        start_time, end_time, num=resolution * (end_time - start_time)
    )
    rates = np.empty((len(recipe_collection.recipes), len(samples)))
    for i, recipe in enumerate(recipe_collection.recipes):
        # 1.1 Calculate the propensity function for each reaction
        f = interp1d(recipe.times, recipe.rates, kind="linear", bounds_error=True)
        # does not deal with implicit 0-rates inbetween explicit rates
        # use of interpolation is to get uniform time steps?!
        rates[i] = trapezoid(samples, f(samples))
        recipes.append(recipe)
        # add option to convert avg_rates to rates??

    return chosen_recipe, dt, rates
