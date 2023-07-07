"""
Kinetic Monte Carlo (KMC) classes and functions.

In our system, the reaction rate r = (deterministic) reaction constant k
= stochastic reaction constant c (from gillespie 1977)
= propensity a (from Anderson 2007)
because of the fundamental premise of chemical kinetics
and because we have one reactant molecule
"""
from typing import Union
import logging
import numpy as np
from dataclasses import dataclass
from numpy.random import default_rng
from kimmdy.reaction import RecipeCollection, RecipeStep



@dataclass
class KMCResult:
    """The result of a KMC step. Similar to a Recipe but for the concrete realization of a reaction.

    Attributes
    ----------
    recipe_steps :
        Single sequence of RecipeSteps to build product
    reaction_probability :
        Integral of reaction propensity with respect to time
    time_step :
        Time step during which the reaction occurs
    """
    recipe_steps: Union[list[RecipeStep], None] = None
    reaction_probability: Union[list[float], None] = None
    time_step: Union[float, None] = None


def rf_kmc(
    recipe_collection: RecipeCollection, rng: np.random.BitGenerator = default_rng()
) -> KMCResult:
    """Rejection-Free Monte Carlo.
    takes RecipeCollection and choses a recipe based on the relative propensity of the events.


    Compare e.g. <https://en.wikipedia.org/wiki/Kinetic_Monte_Carlo#Rejection-free_KMC>

    Parameters
    ---------
    recipe_collection :
        from which one will be choosen
    rng :
        function to generate random numbers in the KMC step
    """

    # check for empty ReactionResult
    if len(recipe_collection.recipes) == 0:
        logging.warning("Empty ReactionResult; no reaction chosen")
        return KMCResult()

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

    return KMCResult(
        recipe_steps=recipe.recipe_steps,
        reaction_probability=reaction_probability,
        time_step=time_step,
    )


def frm(
    recipe_collection: RecipeCollection,
    rng: np.random.BitGenerator = default_rng(),
    MD_time: Union[float, None] = None,
) -> dict:
    """First Reaction Method variant of Kinetic Monte Carlo.
    takes RecipeCollection and choses a recipe based on which reaction would occur.

    Compare e.g. <https://en.wikipedia.org/wiki/Kinetic_Monte_Carlo#Time-dependent_Algorithms>

    Parameters
    ---------
    recipe_collection :
        from which one will be choosen
    rng :
        to generate random numbers in the KMC step
    MD_time :
        time [ps] to compare conformational events with reaction events in the time domain
    """

    # check for empty ReactionResult
    if len(recipe_collection.recipes) == 0:
        logging.warning("Empty ReactionResult; no reaction chosen")
        return KMCResult()

    # 0. Initialization
    reaction_probability = []
    recipes_steps = [["MD"]] if MD_time else []
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
            recipes_steps.append(recipe.recipe_steps)
    # 4. Find the event whose putative time is least
    try:
        pos_event = np.argmin(tau)
        chosen_steps = recipes_steps[pos_event]
        time_step = tau[pos_event]
    except ValueError:
        logging.warning(
            f"FRM recipe selection did not work, probably tau: {tau} is empty."
        )
        return KMCResult()

    return KMCResult(
        recipe_steps=chosen_steps,
        reaction_probability=reaction_probability,
        time_step=time_step,
    )
