from kimmdy.topology.topology import Topology
from kimmdy.parameterize import Parameterizer


class GrappaInterface(Parameterizer):
    def parameterize_topology(self, current_topology: Topology) -> Topology:
        return current_topology
