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

        assert config.iterations == 3

        assert isinstance(config.cwd, Path)
        assert isinstance(config.out, Path)

        assert isinstance(config.plumed, Config)
        assert isinstance(config.plumed.dat, str)
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
            if not hasattr(config.plumed,'dat'):
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


if __name__ == "__main__":
    test_parse_config1_casting()
    test_parse_config2_missing_dat_in_plumed()
    test_parse_config3_missing_mdp_file()
