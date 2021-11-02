import logging
from kimmdy.reaction import Reaction, ConversionRecipe

class Homolysis(Reaction):
    """Homolytic bond breaking leading to 2 radicals"""
    # TODO
    def get_rates(self):
        logging.info("Generating rates for reaction: homolysis")
        rates = [0, 1, 42]
        return rates

    def get_recipe(self):
        logging.info("Generating convertion recipe for reaction: homolysis")
        recipe = ConversionRecipe("break", [(1, 2)])
        return recipe


