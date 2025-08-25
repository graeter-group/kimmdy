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
