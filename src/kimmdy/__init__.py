from __future__ import annotations
import sys

from kimmdy.reaction import Reaction

if sys.version_info > (3, 10):
    from importlib_metadata import entry_points
else:
    from importlib.metadata import entry_points

discovered_plugins = entry_points()["kimmdy.plugins"]
plugins: dict[str, Reaction | Exception] = {}
for _ep in discovered_plugins:
    try:
        plugins[_ep.name] = _ep.load()
    except Exception as _e:
        plugins[_ep.name] = _e
