from typing import Union

import logging
import numpy as np
from numpy.random import default_rng
from kimmdy.reaction import RecipeCollection, Recipe

# In our system, the reaction rate r = (deterministic) reaction constant k
# = stochastic reaction constant c (from gillespie 1977)
# = propensity a (from Anderson 2007)
# because of the fundamental premise of chemical kinetics
# and because we have one reactant molecule


def rfKMC(
    recipe_collection: RecipeCollection, rng: np.random.BitGenerator = default_rng()
) -> dict:
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
        return {"chosen_recipe": RecipeCollection([])}

    # 0. Initialization
    reaction_probability = []
    recipes = []
    for recipe in recipe_collection.recipes:
        # 1.1 Calculate the probability for each reaction
        dt = [x[1] - x[0] for x in recipe.timespans]
        reaction_probability.append(sum(map(lambda x, y: x * y, dt, recipe.rates)))
        recipes.append(recipe)
    # 2. Set the total rate to the sum of individual rates
    probability_cumulative = np.cumsum(reaction_probability)
    probability_sum = probability_cumulative[-1]
    # 3. Generate two independent uniform (0,1) random numbers u1,u2
    u = rng.random(2)
    logging.debug(
        f"Random values u: {u}, cumulative probability {probability_cumulative}, probability sum {probability_sum}"
    )

    # 4. Find the even to carry out, mu, using binary search (np.searchsorted)
    pos = np.searchsorted(probability_cumulative, u[0] * probability_sum)
    recipe = recipe_collection.recipes[pos]
    logging.info(f"Chosen Recipe: {recipe}")
    # 5. Calculate the time step associated with mu
    time_step = 1 / probability_sum * np.log(1 / u[1])

    return {
        "recipe_steps": recipe.recipe_steps,
        "time_step": time_step,
        "reaction_probability": reaction_probability,
    }


def FRM(
    recipe_collection: RecipeCollection,
    rng: np.random.BitGenerator = default_rng(),
    MD_time: Union[float, None] = None,
) -> dict:
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
        return {
            "recipe_steps": None,
            "time_step": None,
            "reaction_probability": None,
        }

    # 0. Initialization
    reaction_probability = []
    recipes = ["MD"] if MD_time else []
    tau = [MD_time] if MD_time else []
    # conformational change is an event that occurs during MD time
    # 1. Generate M independent uniform (0,1) random numbers
    u = rng.random(len(recipe_collection.recipes))
    # 2. Calculate the cumulative probabilities for each event
    for i, recipe in enumerate(recipe_collection.recipes):
        dt = [x[1] - x[0] for x in recipe.timespans]
        cumulative_probability = np.cumsum(
            list(map(lambda x, y: x * y, dt, recipe.rates))
        )
        reaction_probability.append(cumulative_probability[-1])
        # 3. For each event k, find the time until the this reaction takes place tau_k,
        #    if it is during the time when the propensities are defined (simulation time)
        if (p := np.log(1 / u[i])) <= cumulative_probability[-1]:
            # this means a reaction will take place during the defined time
            pos_time = np.searchsorted(cumulative_probability, p)
            tau.append(recipe.timespans[pos_time][1])
            recipes.append(recipe)
    # 4. Find the event whose putative time is least
    try:
        pos_event = np.argmin(tau)
        chosen_recipe = recipes[pos_event]
        time_step = tau[pos_event]
    except ValueError:
        logging.warning(
            f"FRM recipe selection did not work, probably tau: {tau} is empty."
        )
        return {"chosen_recipe": RecipeCollection([])}

    return {
        "recipe_steps": chosen_recipe.recipe_steps,
        "time_step": time_step,
        "reaction_probability": reaction_probability,
    }
