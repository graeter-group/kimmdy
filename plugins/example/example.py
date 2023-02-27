# %%
import sys

if sys.version_info < (3, 10):
    from importlib_metadata import entry_points
else:
    from importlib.metadata import entry_points

# %%
discovered_plugins = entry_points(group="kimmdy.plugins")
discovered_plugins, "dummyreaction" in discovered_plugins.names
# %%
print(discovered_plugins)
reaction = discovered_plugins["dummyreaction"].load()
r = reaction()
type(r)

# %%
