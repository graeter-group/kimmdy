"""Plugin base classes and basic instances thereof.
"""

from __future__ import annotations
from typing import TYPE_CHECKING
from abc import ABC, abstractmethod
import logging

if TYPE_CHECKING:
    from kimmdy.runmanager import RunManager
    from kimmdy.config import Config
    from kimmdy.recipe import RecipeCollection
    from kimmdy.tasks import TaskFiles
    from kimmdy.topology.topology import Topology


logger = logging.getLogger(__name__)


class ReactionPlugin(ABC):
    """Reaction base class

    Parameters
    ----------
    name : str
        Name of the reaction
    runmng : Runmanager
        RunManager instance
    """

    def __init__(self, name: str, runmng: RunManager):
        self.name = name
        self.runmng = runmng
        # sub config, settings of this specific reaction:
        self.config: Config = self.runmng.config.reactions.__getattribute__(self.name)

        logger.debug(f"Reaction {self.name} instatiated.")

    @abstractmethod
    def get_recipe_collection(self, files: TaskFiles) -> RecipeCollection:
        """Get a RecipeCollection as a result of the reaction.

        This is run as a [](`~kimmdy.tasks.Task`) in the RunManager.
        How the RecipeCollection is built is up to the reaction.
        It has access to the current state of the system via the
        runmanager `self.runmng` and the files.

        Parameters
        ----------
        files :
            TaskFiles instance
        """
        pass


class Parameterizer(ABC):
    type_scheme = dict()

    @abstractmethod
    def parameterize_topology(self, current_topology: Topology) -> Topology:
        pass


class BasicParameterizer(Parameterizer):
    """reconstruct base force field state"""

    def parameterize_topology(self, current_topology: Topology) -> None:
        """Do nothing,
        all necessary actions should already have happened in bind_bond and break_bond of Topology
        """
        pass
