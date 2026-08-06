import json
from typing_extensions import deprecated
from kimmdy import schema


def test_load_schema(arranged_tmp_path):
    schema.get_combined_scheme()


def test_flatten_scheme(arranged_tmp_path):
    with open(arranged_tmp_path / "kimmdy-yaml-schema.json", "r") as f:
        loaded = json.load(f)
    schema_dict = schema.convert_schema_to_dict(loaded)
    flat = schema.flatten_scheme(schema_dict)
    assert isinstance(flat, list)
    assert all([isinstance(v, dict) for v in flat])
    for d in flat:
        assert set(d.keys()) == set(
            [
                "default",
                "deprecated",
                "desc",
                "enum",
                "key",
                "type",
            ]
        )


def test_product_template_schema_definition():
    loaded = schema.load_kimmdy_schema()
    topology_properties = loaded["properties"]["topology"]["properties"]

    sources = topology_properties["product_template_sources"]
    assert sources["type"] == "array"
    assert sources["pytype"] == "Sequence"
    assert sources["default"] == ["rtp"]
    assert set(sources["items"]["enum"]) == {"itp", "rtp"}
    assert sources["uniqueItems"] is True

    itp_glob = topology_properties["product_template_itp_glob"]
    assert itp_glob["type"] == "string"
    assert itp_glob["pytype"] == "str"
    assert itp_glob["default"] == "*.itp"

    transitions = topology_properties["product_template_transitions"]
    assert transitions["type"] == "object"
    assert transitions["additionalProperties"] == {"type": "string"}
    assert transitions["default"] == {}
