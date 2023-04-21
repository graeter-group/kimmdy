import logging
from numpy.random import default_rng
rng = default_rng()
from kimmdy.reaction import ReactionResult,ConversionRecipe

def rfKMC(
    reaction_results: list[ReactionResult],
) -> ConversionRecipe:
    """Rejection-Free Monte Carlo.
    takes a list of ReactionResults and choses a recipe.

    Parameters
    ---------
    reaction_results: list[ReactionResults]
        from which one will be choosen
    """
    # compare e.g. <https://en.wikipedia.org/wiki/Kinetic_Monte_Carlo#Rejection-free_KMC>

    # flatten the list of rates form the reaction results
    rates = []
    recipes = []
    for reaction_result in reaction_results:
        for outcome in reaction_result:
            rates.append(outcome.rate)
            recipes.append(outcome.recipe)

    total_rate = sum(rates)
    t = rng.random()  # t in [0.0,1.0)
    logging.debug(f"Random value t: {t}, rates {rates}, total rate {total_rate}")
    rate_running_sum = 0

    # if nothing is choosen, return an empty ConversionRecipe
    result = ConversionRecipe()
    for i, rate in enumerate(rates):
        rate_running_sum += rate
        if (t * total_rate) <= rate_running_sum:
            result = recipes[i]
            break
    logging.debug(f"Result: {result}")

    return result