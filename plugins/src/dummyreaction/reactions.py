from kimmdy.reaction import Reaction, ReactionOutcome, ReactionResult, ConversionRecipe
import logging


class dummy_reaction(Reaction):
    """Dummy reaction, does not change the topology"""

    def get_reaction_result(self, files) -> ReactionResult:
        logging.info("Starting dummy reaction, will do nothing")
        return [ReactionOutcome([], 0)]
