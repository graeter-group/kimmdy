"""
KIMMDY __init__ file. 

Discovers and loads KIMMDY plugins.
"""


from __future__ import annotations
import sys

from kimmdy.reaction import ReactionPlugin
from kimmdy.parameterize import Parameterizer

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

reaction_plugins: dict[str, ReactionPlugin | Exception] = {}
for _ep in discovered_reaction_plugins:
    try:
        reaction_plugins[_ep.name] = _ep.load()
    except Exception as _e:
        reaction_plugins[_ep.name] = _e

parameterization_plugins: dict[str, Parameterizer | Exception] = {}
for _ep in discovered_parameterization_plugins:
    try:
        parameterization_plugins[_ep.name] = _ep.load()
    except Exception as _e:
        parameterization_plugins[_ep.name] = _e
