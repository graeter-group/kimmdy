import py
from kimmdy.config import Config
from pathlib import Path
import pytest

def test_parse_config1_casting():
    c1 = Path(__file__).parent / 'test_files/config_test/config1.yml'
    assert c1.exists(), "Input file not found"

    config = Config(c1)

    assert config.dryrun is True
    assert isinstance(config.dryrun, bool)
    
    assert config.experiment == "Experiment 1"

    assert config.run == 1
    
    assert isinstance(config.cwd, Path)

    assert isinstance(config.plumed, Config)
    assert isinstance(config.plumed.dat, Path)

def test_parse_config2_missing_tpr_in_nvt():
    input_f = Path(__file__).parent / 'test_files/config_test/config2.yml'
    assert input_f.exists(), "Input file not found"
    
    with pytest.raises(AssertionError):
        config = Config(input_f)


def test_parse_config3_missing_npt_mdp_file():
    input_f = Path(__file__).parent / 'test_files/config_test/config3.yml'
    assert input_f.exists(), "Input file not found"

    with pytest.raises(LookupError):
        config = Config(input_f)


if __name__ == "__main__":
    test_parse_config1_casting()
    test_parse_config2_missing_tpr_in_nvt()
    test_parse_config3_missing_npt_mdp_file()
    