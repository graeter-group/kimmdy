from kimmdy.recipe import RecipeCollection
from kimmdy.plugins import ReactionPlugin


class DummyReaction(ReactionPlugin):
    """Dummy reaction, returns empty RecipeCollection"""

    def get_recipe_collection(self, files) -> RecipeCollection:
        files.logger.info("Starting dummy reaction, will do nothing")
        return RecipeCollection([])
