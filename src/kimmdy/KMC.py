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
    for i, path in enumerate(recipe_collection.recipes):
        if path.avg_rates != None:
            avg_weigth = [x[1] - x[0] for x in path.avg_frames]
            rates.append(sum(map(lambda x, y: x * y, avg_weigth, path.avg_rates)))
        else:
            curr_rate = integrate(path.r_ts + 1e-30, path.ts)
            # nan case could be problematic
            rates.append(curr_rate if str(curr_rate) != "nan" else 0)
        recipes.append(path.recipe)

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
