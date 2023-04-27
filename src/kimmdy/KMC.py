import logging
import numpy as np
import scipy.integrate 
from numpy.random import default_rng
from kimmdy.reaction import ReactionResult,ConversionRecipe

def integrate(y,x):
    return scipy.integrate.trapezoid(y)

def rfKMC(
    reaction_result: ReactionResult, rng: np.random.BitGenerator = default_rng()
) -> ConversionRecipe:
    """Rejection-Free Monte Carlo.
    takes ReactionResults and choses a recipe.

    Parameters
    ---------
    reaction_result: ReactionResult
        from which one will be choosen
    """
    # compare e.g. <https://en.wikipedia.org/wiki/Kinetic_Monte_Carlo#Rejection-free_KMC>

    # check for empty ReactionResult
    if len(reaction_result.outcomes) == 0:
        logging.warning("Empty ReactionResult; no reaction chosen")
        return  ConversionRecipe()

    rates = []
    recipes = []
    for i,outcome in enumerate(reaction_result.outcomes):
        curr_rate = integrate(outcome.r_ts+1e-30,outcome.ts)
        # nan case could be problematic
        rates.append(curr_rate if str(curr_rate) != 'nan' else 0 )
        recipes.append(outcome.recipe)

    rates_cumulative = np.cumsum(rates)
    total_rate = rates_cumulative[-1]
    u = rng.random()  # u in [0.0,1.0)
    logging.debug(f"Random value u: {u}, cumulative rates {rates_cumulative}, total rate {total_rate}")

    # Find the even to carry out, mu, using binary search (np.searchsorted)
    pos = np.searchsorted(rates_cumulative,u*total_rate)
    chosen_recipe = recipes[pos]
    logging.debug(f"Chosen Recipe: {chosen_recipe}")

    return chosen_recipe, rates