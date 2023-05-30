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
    takes RecipeCollection and choses a recipe.

    Parameters
    ---------
    reaction_result: ReactionResult
        from which one will be choosen
    """
    # compare e.g. <https://en.wikipedia.org/wiki/Kinetic_Monte_Carlo#Rejection-free_KMC>

    # check for empty ReactionResult
    if len(recipe_collection.recipes) == 0:
        logging.warning("Empty ReactionResult; no reaction chosen")
        return Recipe()

    rates = []
    recipes = []
    for recipe in recipe_collection.recipes:
        if recipe.avg_rates is not None:
            avg_weigth = [x[1] - x[0] for x in recipe.avg_timespans]
            rates.append(sum(map(lambda x, y: x * y, avg_weigth, recipe.avg_rates)))
        else:
            curr_rate = integrate(recipe.rates, recipe.times)
            rates.append(curr_rate if not np.isnan(curr_rate) else 0)
        recipes.append(recipe)

    rates_cumulative = np.cumsum(rates)
    total_rate = rates_cumulative[-1]
    u = rng.random()  # u in [0.0,1.0)
    logging.debug(
        f"Random value u: {u}, cumulative rates {rates_cumulative}, total rate {total_rate}"
    )

    # Find the even to carry out, mu, using binary search (np.searchsorted)
    pos = np.searchsorted(rates_cumulative, u * total_rate)
    chosen_recipe = recipes[pos]
    logging.debug(f"Chosen Recipe: {chosen_recipe}")

    return chosen_recipe, rates
