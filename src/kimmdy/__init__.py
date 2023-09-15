"""
KIMMDY __init__ file. 

Discovers and loads KIMMDY plugins.
"""


from __future__ import annotations
from typing import TYPE_CHECKING
import sys

if TYPE_CHECKING:
    from kimmdy.plugins import Parameterizer, ReactionPlugin

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

reaction_plugins: dict[str, type[ReactionPlugin]] = {}
broken_reaction_plugins: dict[str, Exception] = {}
for _ep in discovered_reaction_plugins:
    try:
        reaction_plugins[_ep.name] = _ep.load()
    except Exception as _e:
        broken_reaction_plugins[_ep.name] = _e

parameterization_plugins: dict[str, type[Parameterizer]] = {}
broken_parameterization_plugins: dict[str, Exception] = {}
for _ep in discovered_parameterization_plugins:
    try:
        parameterization_plugins[_ep.name] = _ep.load()
    except Exception as _e:
        broken_parameterization_plugins[_ep.name] = _e
