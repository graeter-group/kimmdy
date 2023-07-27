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



def convert_schema_to_dict(
    dictionary: dict) -> dict:
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
        json_type = value.get('type')
        if json_type == 'object':
            result[key] = convert_schema_to_dict(value)
            continue

        pytype = value.get("pytype")
        default = value.get("default")
        if pytype is not None:
            result[key]["pytype"] = eval(pytype)
        if default is not None:
            result[key]["default"] = default

    return result


def config_schema() -> dict:
    """Return the schema for the config file"""
    path = pkg_resources.files(kimmdy) / "kimmdy-yaml-schema.json"
    with path.open("rt") as f:
        schema = json.load(f)
    return schema


def prune(d: dict) -> dict:
    """Remove empty dicts from a nested dict"""
    if not isinstance(d, dict):
        return d
    return {
        k: v
        for k, v in ((k, prune(v)) for k, v in d.items())
        if v is not None and v != {}
    }

def create_schemes() -> tuple[dict, dict]:
    """Create a type scheme from the schema"""

    type_scheme = {}
    default_scheme = {}
    schemas = get_all_plugin_schemas()
    for schema in schemas:
        plugin_type_scheme = convert_schema_to_type_dict(schema)
        plugin_type_scheme = prune(plugin_type_scheme)
        plugin_default_scheme = convert_schema_to_type_dict(
            schema, field="default", eval_field=False
        )
        type_scheme.update(plugin_type_scheme)
        default_scheme.update(plugin_default_scheme)

    schema = config_schema()
    main_type_scheme = convert_schema_to_type_dict(schema)
    main_type_scheme = prune(main_type_scheme)
    type_scheme.update(main_type_scheme)

    main_default_scheme = convert_schema_to_type_dict(
            schema, field="default", eval_field=False
        )
    main_default_scheme = prune(main_default_scheme)

    default_scheme.update(main_default_scheme)

    return type_scheme, default_scheme


def get_properties(schema, section=""):
    """recursively get propteries and their desicripions from the schema"""
    properties = []
    property_section = schema.get("properties")

    if property_section is None:
        patterns = schema.get("patternProperties")
        if patterns is not None:
            patterns = patterns.get(".*")
            section = f"{section}.\\*"
            property_section = patterns.get("properties")

    if not property_section:
        return properties
    for key, value in property_section.items():
        if section:
            key = f"{section}.{key}"
        description = value.get("description", "")
        pytype = value.get("pytype", "")
        properties.append((key, pytype, description))
        if value.get("type", "") == "object":
            properties.extend(get_properties(value, key))
    return properties


def get_plugin_schema(plugin) -> Optional[dict]:
    """Return the schema for the config file for a kimmdy reaction plugin."""
    path = pkg_resources.files(plugin) / "kimmdy-yaml-schema.json"
    if not path.is_file():
        return None
    schema = {}
    name = plugin.split(".")[-1]
    schema['properties'] = {}
    schema['properties']['reactions'] = {}
    schema['properties']['reactions']['properties'] = {}
    with path.open("rt") as f:
        schema['properties']['reactions']['properties'][name] = json.load(f)
    return schema


def get_all_plugin_schemas():
    """loop over all discovered plugins and print their schema"""
    if sys.version_info > (3, 10):
        from importlib_metadata import entry_points
        discovered_plugins = entry_points(group="kimmdy.plugins")
    else:
        from importlib.metadata import entry_points
        discovered_plugins = entry_points()["kimmdy.plugins"]

    schemas = []
    for entry_point in discovered_plugins:
        # get entry point of plugin
        plugin = entry_point.load()
        # get main module from that plugin
        plugin = plugin.__module__.split(".")[0]
        if plugin == "kimmdy":
            continue
        schema = get_plugin_schema(plugin)
        if schema is None:
            continue
        schemas.append(schema)
    return schemas


def get_overall_schema():
    """get the schema for the overall config file"""
    schema = config_schema()
    schemas = get_all_plugin_schemas()
    for schema in schemas:
        # schema.update(schema)
        pass
    return schema

def generate_markdown_table(schema, section="", append=False):
    properties = get_properties(schema, section)
    table = []
    if not append:
        table.append("| Option | Type | (_default_) Description |")
        table.append("| --- | --- | --- | --- |")

    for key, pytype, description in properties:
        if pytype == "":
            key = f"**{key}**"
        row = f"| {key} | {pytype} | {description} |"
        table.append(row)

    return "\n".join(table)


# type_scheme, default_scheme = create_schemes()

