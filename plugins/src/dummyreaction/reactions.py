from kimmdy.reaction import RecipeCollection, ReactionPlugin
import logging


class dummy_reaction(ReactionPlugin):
    """Dummy reaction, returns empty RecipeCollection"""

    def get_recipe_collection(self, files) -> RecipeCollection:
        logging.info("Starting dummy reaction, will do nothing")
        return RecipeCollection([])
