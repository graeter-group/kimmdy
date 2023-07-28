"""
Handle the schema for the config file.
To  be used by the config module to validate the config file and set defaults
for the Config object.

Reserved keywords:
    - pytype
    - default
    - description
    - type
"""
import json
import importlib.resources as pkg_resources
import sys
from typing import Optional

# needed for eval of type_scheme from schema
# don't remove even if lsp says it's unused
import kimmdy
import pathlib
from pathlib import Path


class Sequence(list):
    """A sequence of tasks."""

    def __init__(self, tasks: list):
        list.__init__(self)
        for task in tasks:
            if isinstance(task, dict):
                for _ in range(task["mult"]):
                    assert isinstance(
                        task["tasks"], list
                    ), "Grouped tasks must be a list!"
                    self.extend(task["tasks"])
            else:
                self.append(task)

    def __repr__(self):
        return f"Sequence({list.__repr__(self)})"


def load_kimmdy_schema() -> dict:
    """Return the schema for the config file"""
    path = pkg_resources.files(kimmdy) / "kimmdy-yaml-schema.json"
    with path.open("rt") as f:
        schema = json.load(f)
    return schema


def load_plugin_schemas() -> dict:
    """Return the schemas for the plugins"""
    if sys.version_info > (3, 10):
        from importlib_metadata import entry_points

        discovered_plugins = entry_points(group="kimmdy.plugins")
    else:
        from importlib.metadata import entry_points

        discovered_plugins = entry_points()["kimmdy.plugins"]

    schemas = {}
    for entry_point in discovered_plugins:
        # get entry point of plugin
        plugin = entry_point.load()
        # get main module from that plugin
        plugin = plugin.__module__.split(".")[0]
        if plugin == "kimmdy":
            continue
        path = pkg_resources.files(plugin) / "kimmdy-yaml-schema.json"
        name = plugin.split(".")[-1]
        with path.open("rt") as f:
            schemas[name] = json.load(f)
    return schemas


def convert_schema_to_dict(dictionary: dict) -> dict:
    """Convert a dictionary from a raw json schema to a nested dictionary
    where each leaf entry is a dictionary with the "pytype" and "default".
    """
    result = {}
    properties = dictionary.get("properties")
    patternProperties = dictionary.get("patternProperties")
    if properties is None and patternProperties is not None:
        properties = patternProperties
    if properties is None:
        return result
    if patternProperties is not None:
        properties.update(patternProperties)
    for key, value in properties.items():
        if not isinstance(value, dict):
            continue
        result[key] = {}
        json_type = value.get("type")
        if json_type == "object":
            result[key] = convert_schema_to_dict(value)

        pytype = value.get("pytype")
        default = value.get("default")
        description = value.get("description")
        if pytype is not None:
            result[key]["pytype"] = eval(pytype)
        if default is not None:
            result[key]["default"] = default
        if description is not None:
            result[key]["description"] = description

    return result


def get_combined_scheme() -> dict:
    """Return the schema for the config file"""
    schema = load_kimmdy_schema()
    schemas = load_plugin_schemas()
    kimmdy_dict = convert_schema_to_dict(schema)
    plugin_dicts = {k: convert_schema_to_dict(schema) for k, schema in schemas.items()}
    for k, v in plugin_dicts.items():
        kimmdy_dict["reactions"].update({k: v})

    return kimmdy_dict


def prune(d: dict) -> dict:
    """Remove empty dicts from a nested dict"""
    if not isinstance(d, dict):
        return d
    return {
        k: v
        for k, v in ((k, prune(v)) for k, v in d.items())
        if v is not None and v != {}
    }


def flatten_scheme(scheme, section="") -> list:
    """recursively get properties and their desicripions from the scheme"""
    ls = []
    for key, value in scheme.items():
        if not isinstance(value, dict):
            continue
        if section:
            key = f"{section}.{key}"

        description = value.get("description", "")
        pytype = value.get("pytype")
        if pytype is not None:
            pytype = pytype.__name__
        else:
            pytype = ""
        default = value.get("default", "")

        ls.append((key, description, pytype, default))

        for k in value.keys():
            if k not in ["pytype", "default", "description", "type"]:
                k_esc = k
                if k == ".*":
                    k_esc = "\\*"
                s = f"{key}.{k_esc}"
                ls.extend(flatten_scheme(value[k], section=s))

    return ls


def generate_markdown_table(scheme, append=False):
    table = []
    if not append:
        table.append("| Option | Description | Type | Default |")
        table.append("| --- | --- | --- | --- | --- |")

    for key, pytype, description, default in scheme:
        if pytype == "":
            key = f"**{key}**"
        row = f"| {key} | {pytype} | {description} | {default} |"
        table.append(row)

    return "\n".join(table)
