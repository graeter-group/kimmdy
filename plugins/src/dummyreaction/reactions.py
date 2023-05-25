from kimmdy.reaction import ReactionResults, ReactionPlugin
import logging


class dummy_reaction(ReactionPlugin):
    """Dummy reaction, returns empty ReactionResults"""

    def get_reaction_result(self, files) -> ReactionResults:
        logging.info("Starting dummy reaction, will do nothing")
        return ReactionResults([])
