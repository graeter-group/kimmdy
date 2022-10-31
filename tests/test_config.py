from kimmdy.config import Config
from pathlib import Path
import pytest
import os


def test_parse_config1_casting():
    try:
        input_f = Path(__file__).parent / "test_files/config_test/config1.yml"
        os.chdir(input_f.parent)
        assert input_f.exists(), "Input file not found"

        config = Config(input_f)

        assert config.dryrun is True
        assert isinstance(config.dryrun, bool)
        assert config.iterations == 10

        assert isinstance(config.cwd, Path)
        assert isinstance(config.out, Path)

        assert isinstance(config.mds, Config)
        assert isinstance(config.mds.equilibrium1, Config)
        assert isinstance(config.mds.pull1.plumed, Config)
        assert isinstance(config.mds.pull1.plumed.dat, Path)
    finally:
        for d in input_f.parent.glob("test_config_1*"):
            # [f.unlink() for f in d.iterdir()]
            d.rmdir()


def test_parse_config2_missing_dat_in_plumed():
    try:
        input_f = Path(__file__).parent / "test_files/config_test/config2.yml"
        os.chdir(input_f.parent)
        assert input_f.exists(), "Input file not found"

        with pytest.raises(ValueError):
            config = Config(input_f)
            if not hasattr(config.mds.pull.plumed, "dat"):
                raise ValueError
    finally:
        for d in input_f.parent.glob("test_config_2*"):
            [f.unlink() for f in d.iterdir()]
            d.rmdir()


def test_parse_config3_missing_mdp_file():
    try:
        input_f = Path(__file__).parent / "test_files/config_test/config3.yml"
        os.chdir(input_f.parent)
        assert input_f.exists(), "Input file not found"

        with pytest.raises(LookupError):
            Config(input_f)
    finally:
        for d in input_f.parent.glob("test_config_3*"):
            [f.unlink() for f in d.iterdir()]
            d.rmdir()

def test_parse_config4_sequence_missing_entry():
    try:
        input_f = Path(__file__).parent / "test_files/config_test/config4.yml"
        os.chdir(input_f.parent)
        assert input_f.exists(), "Input file not found"

        with pytest.raises(AssertionError):
            Config(input_f)
    finally:
        for d in input_f.parent.glob("test_config_4*"):
            [f.unlink() for f in d.iterdir()]
            d.rmdir()

def test_parse_config5_sequence_missing_entry_no_mds():
    try:
        input_f = Path(__file__).parent / "test_files/config_test/config5.yml"
        os.chdir(input_f.parent)
        assert input_f.exists(), "Input file not found"

        with pytest.raises(AssertionError):
            Config(input_f)
    finally:
        for d in input_f.parent.glob("test_config_5*"):
            [f.unlink() for f in d.iterdir()]
            d.rmdir()

def test_parse_config6_changer_bad_reference():
    try:
        input_f = Path(__file__).parent / "test_files/config_test/config6.yml"
        os.chdir(input_f.parent)
        assert input_f.exists(), "Input file not found"

        with pytest.raises(AssertionError):
            Config(input_f)
    finally:
        for d in input_f.parent.glob("test_config_6*"):
            [f.unlink() for f in d.iterdir()]
            d.rmdir()

if __name__ == "__main__":
    test_parse_config1_casting()
    test_parse_config2_missing_dat_in_plumed()
    test_parse_config3_missing_mdp_file()
    test_parse_config4_sequence_missing_entry()
    test_parse_config5_sequence_missing_entry_no_mds()
    test_parse_config6_changer_bad_reference()
