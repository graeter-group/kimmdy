import pytest
import shutil
import os
from pathlib import Path
from dataclasses import dataclass, field

from kimmdy.runmanager import RunManager
from kimmdy.topology.topology import Topology
from kimmdy.config import Config
from kimmdy.tasks import TaskFiles
from kimmdy.parsing import read_top, read_json
from typing import Callable

# name of this file is fixed for pytest to recognize the fixture without importing


@dataclass
class SlimFiles(TaskFiles):
    input: dict[str, Path] = field(default_factory=dict)
    output: dict[str, Path] = field(default_factory=dict)
    outputdir: Path = Path()
    get_latest: Callable = lambda x: x


@pytest.fixture
def generic_rmgr(tmp_path):
    shutil.copytree(
        Path(__file__).parent / "test_files/test_homolysis",
        tmp_path,
        dirs_exist_ok=True,
    )
    os.chdir(tmp_path)
    Path(tmp_path / "amber99sb-star-ildnp.ff").symlink_to(
        Path(__file__).parent / "test_files/assets/amber99sb-star-ildnp.ff",
        target_is_directory=True,
    )
    return RunManager(Config(tmp_path / "kimmdy.yml"))


@pytest.fixture(scope="function")
def generic_topology():
    filedir = Path(__file__).parent / "test_files/test_integration/hat_naive"
    top_path = filedir / "Ala_out.top"
    top_dict = read_top(top_path)
    top = Topology(top_dict)
    return top


@pytest.fixture
def generic_parameter_input():
    return read_json(
        Path(__file__).parent / "test_files/assets/GrAPPa_input_alanine.json"
    )


@pytest.fixture
def generic_parameter_output():
    return read_json(
        Path(__file__).parent / "test_files/assets/GrAPPa_output_alanine.json"
    )
