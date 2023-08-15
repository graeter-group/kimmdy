from kimmdy.topology.topology import Topology
from abc import ABC, abstractmethod


class Parameterizer(ABC):
    type_scheme = dict()

    @abstractmethod
    def parameterize_topology(self, current_topology: Topology) -> Topology:
        pass


class BasicParameterizer(Parameterizer):
    """reconstruct base force field state"""

    def parameterize_topology(self, current_topology: Topology) -> None:
        # all necessary actions should already have happened in bind_bond and break_bond of Topology
        pass
