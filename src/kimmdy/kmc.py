"""
Kinetic Monte Carlo (KMC) classes and functions.

In our system, the reaction rate r = (deterministic) reaction constant k
= stochastic reaction constant c (from gillespie 1977)
= propensity a (from Anderson 2007)
because of the fundamental premise of chemical kinetics
and because we have one reactant molecule
"""
from typing import Optional
import logging
import numpy as np
from dataclasses import dataclass, field
from numpy.random import default_rng
from kimmdy.recipe import RecipeCollection, Recipe

logger = logging.getLogger(__name__)


@dataclass
class KMCResult:
    """The result of a KMC step. Similar to a Recipe but for the concrete realization of a reaction.

    Attributes
    ----------
    recipe
        Single sequence of RecipeSteps to build product
    reaction_probability
        Integral of reaction propensity with respect to time
    time_delta
        MC time jump during which the reaction occurs [ps]
    time_start
        Time, from which the reaction starts. The reaction changes the
        geometry/topology of this timestep and continues from there. [ps]
    """

    recipe: Recipe = field(default_factory=lambda: Recipe([], [], []))
    reaction_probability: Optional[list[float]] = None
    time_delta: Optional[float] = None
    time_start: Optional[float] = None


def rf_kmc(
    recipe_collection: RecipeCollection, rng: np.random.BitGenerator = default_rng()
) -> KMCResult:
    """Rejection-Free Monte Carlo.
    Takes RecipeCollection and choses a recipe based on the relative propensity of the events.
    The 'start' time of the reaction is the time of the highest rate of the accepted reaction.

    Compare e.g. [Wikipedia KMC - rejection free](https://en.wikipedia.org/wiki/Kinetic_Monte_Carlo#Rejection-free_KMC)

    Parameters
    ---------
    recipe_collection
        from which one will be choosen
    rng
        function to generate random numbers in the KMC step
    """

    # check for empty ReactionResult
    if len(recipe_collection.recipes) == 0:
        logger.warning("Empty ReactionResult; no reaction chosen")
        return KMCResult()

    # 0. Initialization
    reaction_probability = []

    # 1. Calculate the probability for each reaction
    for recipe in recipe_collection.recipes:
        dt = [x[1] - x[0] for x in recipe.timespans]
        reaction_probability.append(sum(np.multiply(dt, recipe.rates)))

    # 2. Set the total rate to the sum of individual rates
    probability_cumulative = np.cumsum(reaction_probability)
    probability_sum = probability_cumulative[-1]

    # 3. Generate two independent uniform (0,1) random numbers u1,u2
    u = rng.random(2)
    logger.debug(
        f"Random values u: {u}, cumulative probability {probability_cumulative}, probability sum {probability_sum}"
    )

    # 4. Find the even to carry out, mu, using binary search (np.searchsorted)
    pos = np.searchsorted(probability_cumulative, u[0] * probability_sum)
    recipe = recipe_collection.recipes[pos]
    reaction_time = recipe.timespans[np.argmax(recipe.rates)][1]
    logger.info(f"Chosen Recipe: {recipe} at time {reaction_time}")

    # 5. Calculate the time step associated with mu
    time_delta = np.log(1 / u[1]) / probability_sum

    return KMCResult(
        recipe=recipe,
        reaction_probability=reaction_probability,
        time_delta=time_delta,
        time_start=reaction_time,
    )


def frm(
    recipe_collection: RecipeCollection,
    rng: np.random.BitGenerator = default_rng(),
    MD_time: Optional[float] = None,
) -> KMCResult:
    """First Reaction Method variant of Kinetic Monte Carlo.
    takes RecipeCollection and choses a recipe based on which reaction would occur.

    Compare e.g. [Wikipedia KMC - time dependent](https://en.wikipedia.org/wiki/Kinetic_Monte_Carlo#Time-dependent_Algorithms)

    Parameters
    ---------
    recipe_collection
        from which one will be choosen
    rng
        to generate random numbers in the KMC step
    MD_time
        time [ps] to compare conformational events with reaction events in the time domain
    """

    # TODO:
    # Currently, this looks at the complete traj and says something happens in this
    # integral or not. It should split it into smaller junks, since the rates
    # fluctuate. Then, still a whole integral is accepted, but one integral does not
    # contain many different conformations.

    # check for empty ReactionResult
    if len(recipe_collection.recipes) == 0:
        logger.warning("Empty ReactionResult; no reaction chosen")
        return KMCResult()

    # 0. Initialization
    reaction_probability = []
    # empty recipe for continuing the MD
    recipes = [Recipe([], [], [])] if MD_time else []
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
        time_delta = tau[pos_event]
    except ValueError:
        logger.warning(
            f"FRM recipe selection did not work, probably tau: {tau} is empty."
        )
        return KMCResult()

    return KMCResult(
        recipe=chosen_recipe,
        reaction_probability=reaction_probability,
        time_delta=time_delta,
        time_start=0,
    )
