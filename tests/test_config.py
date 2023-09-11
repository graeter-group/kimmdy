from kimmdy.config import Config
from kimmdy.runmanager import get_existing_files
from pathlib import Path
import pytest


@pytest.mark.require_gmx
def test_parse_config1_casting(arranged_tmp_path):
    input_f = Path("config1.yml")

    assert input_f.exists(), "Input file not found"

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


def test_parse_config2_start_with_reaction(arranged_tmp_path):
    input_f = Path("config2.yml")

    assert input_f.exists(), "Input file not found"

    config = Config(input_f)

    assert isinstance(config.tpr, Path)
    assert isinstance(config.trr, Path)
    assert isinstance(config.plumed, Path)
    assert config.sequence == ["homolysis"]


def test_parse_config3_missing_mdp_file(arranged_tmp_path):
    input_f = Path("config3.yml")

    assert input_f.exists(), "Input file not found"

    with pytest.raises(LookupError):
        Config(input_f)


def test_parse_config4_sequence_missing_entry(arranged_tmp_path):
    input_f = Path("config4.yml")

    assert input_f.exists(), "Input file not found"

    with pytest.raises(AssertionError):
        Config(input_f)


def test_parse_config5_sequence_missing_entry_no_mds(arranged_tmp_path):
    input_f = Path("config5.yml")

    assert input_f.exists(), "Input file not found"

    with pytest.raises(AssertionError):
        Config(input_f)


def test_parse_config6_changer_bad_reference(arranged_tmp_path):
    input_f = Path("config6.yml")

    assert input_f.exists(), "Input file not found"

    with pytest.raises(AssertionError):
        Config(input_f)


def test_get_existing_files(arranged_tmp_path):
    input_f = Path("config1.yml")

    assert input_f.exists(), "Input file not found"

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
