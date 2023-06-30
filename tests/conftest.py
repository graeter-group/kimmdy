import pytest
import shutil
import os
from pathlib import Path
from dataclasses import dataclass, field

from kimmdy.runmanager import RunManager
from kimmdy.config import Config
from kimmdy.tasks import TaskFiles

# name of this file is fixed for pytest to recognize the fixture without importing


@dataclass
class SlimFiles:
    input: dict[str, Path] = field(default_factory=dict)
    output: dict[str, Path] = field(default_factory=dict)
    outputdir: Path = Path()


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
