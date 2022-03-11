from kimmdy.reaction import Reaction, ReactionResult, ConversionRecipe
import logging


class dummy_reaction(Reaction):
    """Dummy reaction, does not change the topology"""

    def get_reaction_result(self, files) -> ReactionResult:
        logging.info("Starting dummy reaction, will do nothing")
        return ReactionResult()

    @property
    def type_scheme(self):
        """Dict of types of possible entries in config.
        Used to read and check the input config.
        To not use this feature return empty dict
        """
        return dict()
