from typing import Union

import logging
import numpy as np
import scipy.integrate
from numpy.random import default_rng
from kimmdy.reaction import RecipeCollection, Recipe


def integrate(y, x):
    return scipy.integrate.trapezoid(y, x)


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
        return Recipe()

    # 0. Initialization
    rates = []
    recipes = []
    for recipe in recipe_collection.recipes:
        if recipe.avg_rates is not None:
            # 1.1 Calculate the propensity function for each reaction
            avg_weigth = [x[1] - x[0] for x in recipe.avg_timespans]
            overall_timespan = recipe.avg_timespans[-1][1] - recipe.avg_timespans[0][0]
            rates.append(
                sum(map(lambda x, y: x * y, avg_weigth, recipe.avg_rates))
                / overall_timespan
            )
        else:
            # 1.2 Calculate the propensity function for each reaction
            curr_rate = integrate(recipe.rates, recipe.times)
            rates.append(curr_rate if not np.isnan(curr_rate) else 0)
        recipes.append(recipe)
    # 2. Set the total propensity to the sum of individual propensities
    rates_cumulative = np.cumsum(rates)
    total_rate = rates_cumulative[-1]
    # 3. Generate two independent uniform (0,1) random numbers u1,u2
    u = rng.random(2)
    logging.debug(
        f"Random values u: {u}, cumulative rates {rates_cumulative}, total rate {total_rate}"
    )

    # 4. Find the even to carry out, mu, using binary search (np.searchsorted)
    pos = np.searchsorted(rates_cumulative, u[0] * total_rate)
    chosen_recipe = recipes[pos]
    logging.debug(f"Chosen Recipe: {chosen_recipe}")
    # 5. Calculate the time step associated with mu
    dt = 1 / total_rate * np.log(1 / u[1])

    return chosen_recipe, dt, rates


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
        return Recipe()
