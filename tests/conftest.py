"""File for pytest configuration in python and fixture definition.

The name 'conftest.py' is recognized by pytest to execute it before tests.
"""

import pytest
import shutil
import os
from pathlib import Path
from dataclasses import dataclass
from typing import Callable
from kimmdy.plugins import discover_plugins

from kimmdy.tasks import TaskFiles
from kimmdy.utils import get_gmx_dir


## create pytest mark decorators ##
#  register marker
def pytest_configure(config):
    # register an additional mark
    config.addinivalue_line(
        "markers", "require_gmx: mark test to run if gmx is executable"
    )
    config.addinivalue_line(
        "markers", "require_grappa: mark test to run if grappa is available"
    )


# look for mark and define mark action
def pytest_runtest_setup(item):
    require_gmx = [mark for mark in item.iter_markers(name="require_gmx")]
    if require_gmx:
        if not get_gmx_dir():
            pytest.skip(
                f"{item.name} skipped. Command 'gmx' not found, can't test gmx dir parsing."
            )

    require_grappa = [mark for mark in item.iter_markers(name="require_grappa")]
    if require_grappa:
        try:
            import grappa
        except ModuleNotFoundError:
            pytest.skip(f"{item.name} skipped. Can't import module grappa.")


## fixtures for setup and teardown ##
@pytest.fixture
def arranged_tmp_path(tmp_path: Path, request: pytest.FixtureRequest):
    """Arrange temporary directory for tests.

    With files for the test and a symlink to forcefield.
    """

    discover_plugins()

    # if fixture was parameterized, use this for directory with input files
    if hasattr(request, "param"):
        file_dir = Path(__file__).parent / "test_files" / request.param
    # else use stem of requesting file to find directory with input files
    else:
        file_dir = Path(__file__).parent / "test_files" / request.path.stem
    # arrange tmp_path
    shutil.copytree(file_dir, tmp_path, dirs_exist_ok=True)
    assetsdir = Path(__file__).parent / "test_files" / "assets"
    Path(tmp_path / "amber99sb-star-ildnp.ff").symlink_to(
        assetsdir / "amber99sb-star-ildnp.ff",
        target_is_directory=True,
    )
    # change cwd to tmp_path
    os.chdir(tmp_path.resolve())
    return tmp_path


## dummy classes ##
@dataclass
class DummyFiles(TaskFiles):
    get_latest: Callable = lambda: f"DummyCallable"


## general object fixtures ##
