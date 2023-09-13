"""Example useage of installed plugins.
"""
import sys
from importlib.metadata import entry_points


discovered_plugins = entry_points(group="kimmdy.reaction_plugins")
discovered_plugins, "dummyreaction" in discovered_plugins.names

print(discovered_plugins)
reaction = discovered_plugins["dummyreaction"].load()
print(reaction)
