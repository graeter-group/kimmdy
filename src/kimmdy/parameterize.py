from kimmdy.topology.topology import Topology
from abc import ABC, abstractmethod


class Parameterizer(ABC):
    type_scheme = dict()

    @abstractmethod
    def parameterize_topology(self, current_topology: Topology) -> Topology:
        pass


class DefaultParameterizer(Parameterizer):
    """default for parameter patch: reconstruct base force field state"""

    def parameterize_topology(self, current_topology):
        raise NotImplementedError
