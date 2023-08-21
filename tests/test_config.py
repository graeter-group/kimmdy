from kimmdy.config import Config
from kimmdy.runmanager import get_existing_files
from pathlib import Path
import pytest
import os


def test_parse_config1_casting():
    input_f = Path(__file__).parent / "test_files/test_config/config1.yml"
    try:
        os.chdir(input_f.parent)
        assert input_f.exists(), "Input file not found"

        config = Config(input_f)

        assert config.dryrun is True
        assert isinstance(config.dryrun, bool)
        assert config.max_tasks == 10

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
    input_f = Path(__file__).parent / "test_files/test_config/config2.yml"
    try:
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
    input_f = Path(__file__).parent / "test_files/test_config/config3.yml"
    try:
        os.chdir(input_f.parent)
        assert input_f.exists(), "Input file not found"

        with pytest.raises(LookupError):
            Config(input_f)
    finally:
        for d in input_f.parent.glob("test_config_3*"):
            [f.unlink() for f in d.iterdir()]
            d.rmdir()


def test_parse_config4_sequence_missing_entry():
    input_f = Path(__file__).parent / "test_files/test_config/config4.yml"
    try:
        os.chdir(input_f.parent)
        assert input_f.exists(), "Input file not found"

        with pytest.raises(AssertionError):
            Config(input_f)
    finally:
        for d in input_f.parent.glob("test_config_4*"):
            [f.unlink() for f in d.iterdir()]
            d.rmdir()


def test_parse_config5_sequence_missing_entry_no_mds():
    input_f = Path(__file__).parent / "test_files/test_config/config5.yml"
    try:
        os.chdir(input_f.parent)
        assert input_f.exists(), "Input file not found"

        with pytest.raises(AssertionError):
            Config(input_f)
    finally:
        for d in input_f.parent.glob("test_config_5*"):
            [f.unlink() for f in d.iterdir()]
            d.rmdir()


def test_parse_config6_changer_bad_reference():
    input_f = Path(__file__).parent / "test_files/test_config/config6.yml"
    try:
        os.chdir(input_f.parent)
        assert input_f.exists(), "Input file not found"

        with pytest.raises(AssertionError):
            Config(input_f)
    finally:
        for d in input_f.parent.glob("test_config_6*"):
            [f.unlink() for f in d.iterdir()]
            d.rmdir()


def test_parse_config7_start_with_reaction():
    input_f = Path(__file__).parent / "test_files/test_config/config7.yml"
    try:
        os.chdir(input_f.parent)
        assert input_f.exists(), "Input file not found"

        config = Config(input_f)
        
        assert isinstance(config.tpr, Path)
        assert isinstance(config.trr, Path)
        assert isinstance(config.plumed, Config)
        assert isinstance(config.plumed.dat, Path)
        assert config.sequence == ["homolysis"]


    finally:
        for d in input_f.parent.glob("test_config_7*"):
            [f.unlink() for f in d.iterdir()]
            d.rmdir()


def test_get_existing_files():
    input_f = Path(__file__).parent / "test_files/test_config/config1.yml"
    try:
        os.chdir(input_f.parent)
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
                "plumed.dat",
                "pullf1500_equil.mdp",
                "pullf1500.mdp",
                "broken_equil_f1000.mdp",
                "edissoc.dat",
                "ffbonded.itp",
            ]
        )

    finally:
        for d in input_f.parent.glob("test_config_1*"):
            [f.unlink() for f in d.iterdir()]
            d.rmdir()


if __name__ == "__main__":
    test_parse_config1_casting()
    test_parse_config2_missing_dat_in_plumed()
    test_parse_config3_missing_mdp_file()
    test_parse_config4_sequence_missing_entry()
    test_parse_config5_sequence_missing_entry_no_mds()
    test_parse_config6_changer_bad_reference()
