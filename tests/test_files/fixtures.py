import pytest
import shutil
import os
from pathlib import Path

from kimmdy.runmanager import RunManager
from kimmdy.config import Config
from kimmdy.tasks import TaskFiles


@pytest.fixture
def generic_rmgr(tmp_path):
    shutil.copytree(Path(__file__).parent / "test_files/test_homolysis", tmp_path)
    os.chdir(tmp_path)
    Path(tmp_path / "amber99sb-star-ildnp.ff").symlink_to(
        Path(__file__).parent / "test_files/assets/amber99sb-star-ildnp.ff",
        target_is_directory=True,
    )
    rmgr = RunManager(Config(tmp_path / "kimmdy.yml"))
    files = TaskFiles(rmgr)
    files.input["top"] = Path("topol.top")
    for f_path in ["plumed.dat", "distances.dat", "edissoc.dat", "ffbonded.itp"]:
        files.input[f_path] = tmp_path / f_path
