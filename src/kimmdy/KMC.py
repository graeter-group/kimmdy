from typing import Union

import logging
import numpy as np
from scipy.interpolate import interp1d
import scipy
from numpy.random import default_rng
from kimmdy.reaction import RecipeCollection, Recipe

# In our system, the reaction rate r = (deterministic) reaction constant k
# = stochastic reaction constant c (from gillespie 1977)
# = propensity a (from Anderson 2007)
# because of the fundamental premise of chemical kinetics
# and because we have one reactant molecule


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
        # 1.1 Calculate the constant rate for each reaction
        dt = [x[1] - x[0] for x in recipe.timespans]
        constant_rates.append(sum(map(lambda x, y: x * y, dt, recipe.rates)))
        recipes_steps.append(recipe.recipe_steps)
    # 2. Set the total rate to the sum of individual rates
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
    time_step = 1 / total_rate * np.log(1 / u[1])

    return chosen_recipe, time_step, constant_rates


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
    recipes_steps = []
    tau = [MD_time] if MD_time else []
    # conformational change is an event that occurs during MD time
    # 1. Generate M independent uniform (0,1) random numbers
    u = rng.random(len(recipe_collection.recipes))
    # 2. Calculate the cumulative probabilities for each event
    for i, recipe in enumerate(recipe_collection.recipes):
        dt = [x[1] - x[0] for x in recipe.timespans]
        cumulative_probabilities = np.cumsum(map(lambda x, y: x * y, dt, recipe.rates))
        # 3. For each event k, find the time until the this reaction takes place tau_k,
        #    if it is during the time when the propensities are defined (simulation time)
        if p := np.log(1 / u[i]) <= cumulative_probabilities[-1]:
            # this means a reaction will take place during the defined time
            pos_time = np.searchsorted(cumulative_probabilities, p)
            tau.append(recipe.avg_timespans[pos_time][1])
            recipes_steps.append(recipe.recipe_steps)
    # 4. Find the event whose putative time is least
    pos_event = np.argmin(tau)
    chosen_recipe = recipes_steps[pos_event]
    time_step = tau[pos_event]

    return chosen_recipe, time_step
