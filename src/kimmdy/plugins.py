"""Plugin base classes and basic instances thereof.

Also discovers and loads KIMMDY plugins.
"""

from __future__ import annotations

import logging
import sys
from abc import ABC, abstractmethod
from typing import TYPE_CHECKING, Optional

if TYPE_CHECKING:
    from kimmdy.config import Config
    from kimmdy.recipe import RecipeCollection
    from kimmdy.runmanager import RunManager
    from kimmdy.tasks import TaskFiles
    from kimmdy.topology.topology import Topology


reaction_plugins: dict[str, type[ReactionPlugin]] = {}
broken_reaction_plugins: dict[str, Exception] = {}
parameterization_plugins: dict[str, type[Parameterizer]] = {}
broken_parameterization_plugins: dict[str, Exception] = {}


def discover_plugins():
    if sys.version_info > (3, 10):
        from importlib_metadata import entry_points

        discovered_reaction_plugins = entry_points(group="kimmdy.reaction_plugins")
        discovered_parameterization_plugins = entry_points(
            group="kimmdy.parameterization_plugins"
        )
    else:
        from importlib.metadata import entry_points

        discovered_reaction_plugins = entry_points()["kimmdy.reaction_plugins"]
        discovered_parameterization_plugins = entry_points()[
            "kimmdy.parameterization_plugins"
        ]

    for _ep in discovered_reaction_plugins:
        try:
            reaction_plugins[_ep.name] = _ep.load()
        except Exception as _e:
            broken_reaction_plugins[_ep.name] = _e

    for _ep in discovered_parameterization_plugins:
        try:
            parameterization_plugins[_ep.name] = _ep.load()
        except Exception as _e:
            broken_parameterization_plugins[_ep.name] = _e


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

        logger.debug(f"Reaction {self.name} instantiated.")

    def __repr__(self):
        return self.name

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
    def __init__(self, **kwargs):
        self.type_scheme = dict()

    @abstractmethod
    def parameterize_topology(
        self, current_topology: Topology, focus_nrs: Optional[set[str]]
    ) -> Topology:
        pass


class BasicParameterizer(Parameterizer):
    """reconstruct base force field state"""

    def parameterize_topology(
        self, current_topology: Topology, focus_nrs: Optional[set[str]] = None
    ) -> Topology:
        """Do nothing,
        all necessary actions should already have happened in bind_bond and break_bond of Topology
        """
        return current_topology
