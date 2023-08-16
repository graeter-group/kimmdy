from kimmdy.reaction import RecipeCollection, ReactionPlugin


class DummyReaction(ReactionPlugin):
    """Dummy reaction, returns empty RecipeCollection"""

    def get_recipe_collection(self, files) -> RecipeCollection:
        files.logger.info("Starting dummy reaction, will do nothing")
        return RecipeCollection([])
