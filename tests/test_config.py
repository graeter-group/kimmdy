from pathlib import Path

import pytest
import yaml

from kimmdy.config import Config
from kimmdy.parsing import read_top
from kimmdy.runmanager import RunManager, get_existing_files
from kimmdy.topology.ff import FF
from kimmdy.topology.topology import Topology


@pytest.mark.require_gmx
def test_parse_config1_casting(arranged_tmp_path):
    config = Config(Path("config1.yml"))

    assert config.dryrun is True
    assert isinstance(config.dryrun, bool)
    assert config.max_tasks == 10

    assert isinstance(config.cwd, Path)
    assert isinstance(config.out, Path)

    assert isinstance(config.mds, Config)
    assert isinstance(config.mds.equilibrium1, Config)
    assert isinstance(config.plumed, Path)
    assert isinstance(config.mds.pull1.use_plumed, bool)
    assert isinstance(config.mds.pull1.mdp, Path)


def test_non_existent_sections_with_defaults_in_subsections_are_created(
    arranged_tmp_path,
):
    config = Config(Path("config1.yml"))

    assert config.log.file.name == "kimmdy.log"
    assert config.log.level == "INFO"


def test_no_sections_are_created_for_not_mentioned_reactions(arranged_tmp_path):
    config_1 = Config(Path("config1.yml"))
    config_2 = Config(Path("config2.yml"))

    assert config_1.reactions.homolysis
    assert config_1.edissoc.name == "edissoc.dat"
    assert len(config_1.reactions.__dict__) == 1

    assert len(config_2.reactions.__dict__) == 2
    assert config_2.reactions.homolysis
    assert config_2.edissoc.name == "edissoc.dat"
    assert config_2.reactions.hat_naive
    assert config_2.reactions.hat_naive.polling_rate == 1
    assert config_2.reactions.hat_naive.h_cutoff == 4


def test_subsections_with_defaults_are_kept(arranged_tmp_path):
    config = Config(Path("config2.yml"))

    assert config.log.file.name == "kimmdy.log"
    assert config.log.level == "debug"
    assert config.topology.reactive
    assert getattr(config.topology.reactive, "nrexcl", None) is None


def test_out_is_generated_from_name_if_not_set(arranged_tmp_path):
    config_1 = Config(Path("config1.yml"))
    config_2 = Config(Path("config2.yml"))

    assert config_1.name == "kimmdy"
    assert config_1.out.name == "test_config_1"

    assert config_2.name == "config2"
    assert config_2.out.name == "config2"


def test_general_settings_for_mds_are_set(arranged_tmp_path):
    config = Config(Path("config1.yml"))

    assert len(config.mds.__dict__) == 3
    # explicitly set:
    assert config.mds.equilibrium1.mdp.name == "pullf1500_equil.mdp"
    assert config.mds.pull1.use_plumed == True
    # '.*' defaults:
    assert config.mds.equilibrium1.use_plumed == False


def test_complains_if_plumed_used_but_not_set(arranged_tmp_path):
    with open(Path("config1.yml"), "r") as f:
        raw = yaml.safe_load(f)
    del raw["plumed"]

    with pytest.raises(
        AssertionError,
        match="Plumed requested in md section, but not defined at config root",
    ):
        _ = Config(opts=raw)


def test_parse_reaction_only(arranged_tmp_path):
    with open(Path("config1.yml"), "r") as f:
        raw = yaml.safe_load(f)
    raw["tpr"] = "pull.tpr"
    raw["trr"] = "pull.trr"
    del raw["mds"]
    del raw["changer"]["coordinates"]["md"]
    raw["sequence"] = ["homolysis"]

    config = Config(opts=raw)

    assert isinstance(config.tpr, Path)
    assert isinstance(config.trr, Path)
    assert isinstance(config.plumed, Path)
    assert config.sequence == ["homolysis"]


def test_parse_missing_mdp_file(arranged_tmp_path):
    with open(Path("config1.yml"), "r") as f:
        raw = yaml.safe_load(f)
    raw["mds"]["relax"]["mdp"] = "nonexisting.mdp"
    with pytest.raises(LookupError, match="File not found:"):
        Config(opts=raw)


def test_parse_missing_required_mdp_section(arranged_tmp_path):
    with open(Path("config1.yml"), "r") as f:
        raw = yaml.safe_load(f)
    del raw["mds"]["relax"]["mdp"]
    with pytest.raises(
        AssertionError, match="MD instance defined but contains no mdp file."
    ):
        Config(opts=raw)


def test_parse_sequence_missing_entry(arranged_tmp_path):
    with open(Path("config1.yml"), "r") as f:
        raw = yaml.safe_load(f)
    raw["sequence"][1]["tasks"][0] = "nonexistent_entry"
    with pytest.raises(
        AssertionError,
        match="Task nonexistent_entry listed in sequence, but not defined!",
    ):
        Config(opts=raw)


def test_parse_sequence_missing(arranged_tmp_path):
    with open(Path("config1.yml"), "r") as f:
        raw = yaml.safe_load(f)
    del raw["sequence"]
    with pytest.raises(AssertionError, match="No sequence defined in config!"):
        Config(opts=raw)


def test_parse_coordinates_bad_reference(arranged_tmp_path):
    with open(Path("config1.yml"), "r") as f:
        raw = yaml.safe_load(f)
    raw["changer"]["coordinates"]["md"] = "relax_nonexistent.mdp"
    with pytest.raises(
        AssertionError, match="Relax MD relax_nonexistent.mdp not in MD section!"
    ):
        Config(opts=raw)


def test_get_existing_files(arranged_tmp_path):
    config = Config(Path("config1.yml"))
    file_d = get_existing_files(config)
    assert set(file_d.keys()) == set(
        [
            "kimmdy.log",
            "top",
            "gro",
            "ndx",
            "plumed",
            "plumed_out",
            "pullf1500_equil.mdp",
            "pullf1500.mdp",
            "broken_equil_f1000.mdp",
            "edissoc.dat",
            "yml",
        ]
    )


def test_explicit_residuetypes(arranged_tmp_path):
    config = Config(Path("config3.yml"))
    ff = FF(top=read_top(config.top), residuetypes_path=config.residuetypes)

    assert len(ff.residuetypes.keys()) == 121
    assert ff.residuetypes["ALA"]


def test_explicit_radicals(arranged_tmp_path):
    config = Config(Path("config1.yml"))
    config.radicals = "1 2 3 4"
    top = Topology(read_top(config.top), radicals=config.radicals)

    assert len(top.radicals) == 4
    for radical in config.radicals.split(sep=" "):
        assert top.radicals.get(radical)


def test_config_retains_additional_kwargs(
    arranged_tmp_path,
):
    config = Config(Path("config4.yml"))
    assert config.changer.topology.parameterization_kwargs.test_arg == "hello"


def test_config_works_for_no_slow_growth_relax(
    arranged_tmp_path,
):
    config = Config(Path("conf5.yml"))
    runmgr = RunManager(config)
    assert runmgr


def test_config_parses_product_template_settings(arranged_tmp_path):
    with open(Path("config1.yml"), "r") as f:
        raw = yaml.safe_load(f)
    raw["topology"] = {
        "product_template_sources": ["itp", "rtp"],
        "product_template_itp_glob": "DOPA*.itp",
        "product_template_transitions": {"OHR": "OHRH"},
    }

    config = Config(opts=raw)

    assert list(config.topology.product_template_sources) == ["itp", "rtp"]
    assert config.topology.product_template_itp_glob == "DOPA*.itp"
    assert vars(config.topology.product_template_transitions) == {
        "OHR": "OHRH",
    }


@pytest.fixture
def schema_startup_config(tmp_path, monkeypatch):
    """Exercise Config construction without GROMACS or installed plugins."""
    from copy import deepcopy
    from kimmdy import config as config_module
    from kimmdy.schema import convert_schema_to_dict, load_kimmdy_schema

    scheme = convert_schema_to_dict(load_kimmdy_schema())
    scheme["reactions"]["test_rates"] = {
        "rate_equation": {
            "pytype": str,
            "enum": ["arrhenius", "eyring"],
            "default": "eyring",
        },
        "temperature": {"pytype": float, "exclusiveMinimum": 0, "default": 300},
        "coefficient": {"pytype": float, "minimum": 0, "default": 1},
        "unit": {"pytype": str, "enum": ["eV", "kcal/mol"], "default": "eV"},
        "required": ["temperature", "coefficient"],
    }
    monkeypatch.setattr(config_module, "reaction_plugins", {"test_rates": object()})
    monkeypatch.setattr(config_module, "check_file_exists", lambda **kwargs: None)
    monkeypatch.setattr(config_module, "check_gmx_version", lambda config: None)
    monkeypatch.setattr(config_module, "configure_logger", lambda config: None)
    base = {
        "cwd": tmp_path,
        "ff": tmp_path,
        "mds": {"sample": {"mdp": "sample.mdp"}},
        "changer": {"coordinates": {"md": "sample", "slow_growth": "full"}},
        "reactions": {"test_rates": {}},
        "sequence": ["sample", "reactions"],
    }

    def build(rate_options=None, overrides=None):
        opts = deepcopy(base)
        if rate_options is not None:
            opts["reactions"]["test_rates"] = rate_options
        opts.update(overrides or {})
        return Config(opts=opts, scheme=deepcopy(scheme))

    return build, scheme


@pytest.mark.parametrize(
    "options, message",
    [
        ({"rate_equation": "unknown"}, "test_rates.rate_equation must be one of"),
        ({"temperature": 0}, "test_rates.temperature must be > 0"),
        ({"temperature": -1}, "test_rates.temperature must be > 0"),
        ({"coefficient": -1}, "test_rates.coefficient must be >= 0"),
        ({"temperature": float("nan")}, "test_rates.temperature must be finite"),
        ({"temperature": float("inf")}, "test_rates.temperature must be finite"),
        ({"coefficient": float("-inf")}, "test_rates.coefficient must be finite"),
    ],
)
def test_schema_rejects_invalid_settings_during_initialization(
    schema_startup_config, options, message
):
    build, _ = schema_startup_config
    with pytest.raises(ValueError, match=message):
        build(options)


def test_schema_validates_defaults_and_allows_zero_coefficient(schema_startup_config):
    build, _ = schema_startup_config
    config = build({"coefficient": 0})
    assert config.reactions.test_rates.rate_equation == "eyring"
    assert config.reactions.test_rates.temperature == 300.0
    assert config.reactions.test_rates.coefficient == 0.0


@pytest.mark.parametrize("level", ["INFO", "debug"])
@pytest.mark.parametrize(
    "slow_growth, expected",
    [
        ("full", True),
        ("yes", True),
        (True, True),
        ("no", False),
        (False, False),
        ("", False),
        ("morse_only", "morse_only"),
    ],
)
def test_schema_preserves_enum_casing_and_slow_growth_aliases(
    schema_startup_config, level, slow_growth, expected
):
    build, _ = schema_startup_config
    config = build(
        {"unit": "eV", "temperature": 310, "coefficient": 1},
        overrides={
            "log": {"level": level},
            "changer": {"coordinates": {"md": "sample", "slow_growth": slow_growth}},
        },
    )
    assert config.log.level == level.lower()
    assert config.reactions.test_rates.unit == "ev"
    assert config.changer.coordinates.slow_growth == expected


def test_schema_does_not_require_unconfigured_reaction(schema_startup_config):
    build, scheme = schema_startup_config
    scheme["reactions"]["test_rates"]["required"].append("missing")
    config = build(overrides={"reactions": {}, "sequence": ["sample"]})
    assert config.reactions.get_attributes() == []


def test_schema_reports_required_setting_without_default(schema_startup_config):
    build, scheme = schema_startup_config
    scheme["reactions"]["test_rates"]["coefficient"].pop("default")
    with pytest.raises(ValueError, match="required.*test_rates.coefficient"):
        build()


def test_schema_keeps_deprecated_required_options_optional(schema_startup_config):
    build, scheme = schema_startup_config
    scheme["reactions"]["test_rates"]["old_option"] = {"deprecated": True}
    scheme["reactions"]["test_rates"]["required"].append("old_option")
    config = build()
    assert config.reactions.test_rates.rate_equation == "eyring"


def test_schema_keeps_value_of_typed_deprecated_option(
    schema_startup_config, monkeypatch
):
    from kimmdy.constants import CONFIG_LOGS

    monkeypatch.setitem(CONFIG_LOGS, "warnings", [])
    build, scheme = schema_startup_config
    scheme["reactions"]["test_rates"]["old_option"] = {
        "deprecated": True,
        "pytype": int,
        "description": "Old option.",
    }
    config = build({"old_option": "3"})
    assert config.reactions.test_rates.old_option == 3
    assert any(
        "config.reactions.test_rates.old_option" in m for m in CONFIG_LOGS["warnings"]
    )


def test_schema_keeps_value_of_untyped_deprecated_option(
    schema_startup_config, monkeypatch
):
    from kimmdy.constants import CONFIG_LOGS

    monkeypatch.setitem(CONFIG_LOGS, "warnings", [])
    build, scheme = schema_startup_config
    scheme["reactions"]["test_rates"]["old_option"] = {"deprecated": True}
    config = build({"old_option": "Some Value"})
    assert config.reactions.test_rates.old_option == "Some Value"
    assert any(
        "config.reactions.test_rates.old_option" in m for m in CONFIG_LOGS["warnings"]
    )


def test_plumed_options_are_deprecated_but_kept(schema_startup_config, monkeypatch):
    from kimmdy.constants import CONFIG_LOGS

    monkeypatch.setitem(CONFIG_LOGS, "warnings", [])
    build, _ = schema_startup_config
    config = build(
        overrides={
            "plumed": "plumed.dat",
            "mds": {"sample": {"mdp": "sample.mdp", "use_plumed": True}},
        }
    )
    assert isinstance(config.plumed, Path)
    assert config.mds.sample.use_plumed is True
    warnings = "\n".join(CONFIG_LOGS["warnings"])
    assert "config.plumed" in warnings
    assert "config.mds.sample.use_plumed" in warnings


def test_schema_checks_numeric_constraints_in_wildcard_sections(schema_startup_config):
    build, scheme = schema_startup_config
    scheme["mds"][".*"]["stride"] = {"pytype": int, "minimum": 1}
    with pytest.raises(ValueError, match="mds.sample.stride must be >= 1"):
        build(overrides={"mds": {"sample": {"mdp": "sample.mdp", "stride": 0}}})
