from __future__ import annotations
import sys

from kimmdy.reaction import ReactionPlugin

if sys.version_info > (3, 10):
    from importlib_metadata import entry_points

    discovered_plugins = entry_points(group="kimmdy.plugins")
else:
    from importlib.metadata import entry_points

    discovered_plugins = entry_points()["kimmdy.plugins"]

plugins: dict[str, ReactionPlugin | Exception] = {}
for _ep in discovered_plugins:
    try:
        plugins[_ep.name] = _ep.load()
    except Exception as _e:
        plugins[_ep.name] = _e
