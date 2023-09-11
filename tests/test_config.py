from kimmdy.config import Config
from kimmdy.runmanager import get_existing_files
from pathlib import Path
import pytest


@pytest.mark.require_gmx
def test_parse_config1_casting(arranged_tmp_path):
    input_f = Path("config1.yml")
    config = Config(input_f)

    assert config.dryrun is True
    assert isinstance(config.dryrun, bool)
    assert config.max_tasks == 10

    assert isinstance(config.cwd, Path)
    assert isinstance(config.out, Path)

    assert isinstance(config.mds, Config)
    assert isinstance(config.mds.equilibrium1, Config)
    assert isinstance(config.plumed, Path)
    assert isinstance(config.mds.pull1.use_plumed, bool)


def test_non_existent_sections_with_defaults_in_subsections_are_created(
    arranged_tmp_path,
):
    path = Path("config1.yml")
    config = Config(path)
    assert config.log.file.name == "kimmdy.log"
    assert config.log.level == "INFO"


def test_no_sections_are_created_for_not_mentioned_reactions(arranged_tmp_path):
    path = Path("config1.yml")
    config = Config(path)
    assert config.reactions.homolysis
    assert config.reactions.homolysis.edis.name == "edissoc.dat"
    assert len(config.reactions.__dict__) == 1

    path = Path("config7.yml")
    config = Config(path)
    assert len(config.reactions.__dict__) == 2
    assert config.reactions.homolysis
    assert config.reactions.homolysis.edis.name == "edissoc.dat"
    assert config.reactions.hat_naive
    assert config.reactions.hat_naive.polling_rate == 1
    assert config.reactions.hat_naive.h_cutoff == 4


def test_subsections_with_defaults_are_kept(arranged_tmp_path):
    path = Path("config7.yml")
    config = Config(path)
    assert config.log.file.name == "kimmdy.log"
    assert config.log.level == "DEBUG"


def test_out_is_generated_from_name_if_not_set(arranged_tmp_path):
    path = Path("config1.yml")
    config = Config(path)
    assert config.name == "kimmdy"
    assert config.out.name == "test_config_1"

    path = Path("config7.yml")
    config = Config(path)
    assert config.name == "config7"
    assert config.out.name == "config7"


def test_general_settings_for_mds_are_set(arranged_tmp_path):
    path = Path("config1.yml")
    config = Config(path)
    assert len(config.mds.__dict__) == 3
    # explicitly set:
    assert config.mds.equilibrium1.mdp.name == "pullf1500_equil.mdp"
    assert config.mds.pull1.use_plumed == True
    # '.*' defaults:
    assert config.mds.equilibrium1.use_plumed == False


def test_complains_if_plumed_used_but_not_set(arranged_tmp_path):
    path = Path("config8.yml")
    with pytest.raises(
        AssertionError,
        match="Plumed requested in md section, but not defined at config root",
    ):
        _ = Config(path)


def test_parse_config2_start_with_reaction(arranged_tmp_path):
    input_f = Path("config2.yml")
    config = Config(input_f)

    assert isinstance(config.tpr, Path)
    assert isinstance(config.trr, Path)
    assert isinstance(config.plumed, Path)
    assert config.sequence == ["homolysis"]


def test_parse_config3_missing_mdp_file(arranged_tmp_path):
    input_f = Path("config3.yml")

    with pytest.raises(LookupError):
        Config(input_f)


def test_parse_config4_sequence_missing_entry(arranged_tmp_path):
    input_f = Path("config4.yml")

    with pytest.raises(AssertionError):
        Config(input_f)


def test_parse_config5_sequence_missing_entry_no_mds(arranged_tmp_path):
    input_f = Path("config5.yml")

    with pytest.raises(AssertionError):
        Config(input_f)


def test_parse_config6_changer_bad_reference(arranged_tmp_path):
    input_f = Path("config6.yml")

    with pytest.raises(AssertionError):
        Config(input_f)


def test_get_existing_files(arranged_tmp_path):
    input_f = Path("config1.yml")
    config = Config(input_f)
    file_d = get_existing_files(config)
    assert set(file_d.keys()) == set(
        [
            "",
            "top",
            "gro",
            "ndx",
            "ff",
            "plumed",
            "plumed_out",
            "pullf1500_equil.mdp",
            "pullf1500.mdp",
            "broken_equil_f1000.mdp",
            "edissoc.dat",
            "ffbonded.itp",
        ]
    )
