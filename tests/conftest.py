"""File for pytest configuration in python and fixture definition.

The name 'conftest.py' is recognized by pytest to execute it before tests.
"""

import json
import os
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

import pytest
from kimmdy.plugins import discover_plugins

from kimmdy.tasks import TaskFiles
from kimmdy.utils import get_gmx_dir


def pytest_addoption(parser):
    parser.addoption(
        "--gromacs-alias",
        default="gmx",
        help="GROMACS executable used by integration tests (default: gmx)",
    )
    parser.addoption(
        "--mdrun-prefix",
        default="",
        help="Command prefix used for integration-test mdrun calls",
    )


@pytest.fixture(scope="session")
def gromacs_launch(pytestconfig):
    """Return explicit GROMACS integration-test launch overrides."""
    return (
        pytestconfig.getoption("gromacs_alias"),
        pytestconfig.getoption("mdrun_prefix"),
    )


def _configure_integration_yaml(config_path, gromacs_launch):
    """Patch a copied integration YAML without changing its source fixture."""
    gromacs_alias, mdrun_prefix = gromacs_launch
    lines = [
        line
        for line in config_path.read_text().splitlines()
        if not line.startswith("mdrun_prefix:")
    ]
    launch_lines = [f"gromacs_alias: {json.dumps(gromacs_alias)}"]
    if mdrun_prefix:
        launch_lines.append(f"mdrun_prefix: {json.dumps(mdrun_prefix)}")

    alias_index = next(
        (i for i, line in enumerate(lines) if line.startswith("gromacs_alias:")),
        None,
    )
    if alias_index is None:
        lines = launch_lines + lines
    else:
        lines[alias_index : alias_index + 1] = launch_lines

    config_path.write_text("\n".join(lines) + "\n")


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
        gromacs_alias = item.config.getoption("gromacs_alias")
        if not get_gmx_dir(gromacs_alias):
            pytest.skip(
                f"{item.name} skipped. Command {gromacs_alias!r} not found, "
                "can't test GROMACS directory parsing."
            )

    require_grappa = [mark for mark in item.iter_markers(name="require_grappa")]
    if require_grappa:
        try:
            import grappa
        except ModuleNotFoundError:
            pytest.skip(f"{item.name} skipped. Can't import module grappa.")


## fixtures for setup and teardown ##
@pytest.fixture
def arranged_tmp_path(
    tmp_path: Path, request: pytest.FixtureRequest, gromacs_launch
):
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
    if hasattr(request, "param") and str(request.param).startswith(
        "test_integration/"
    ):
        for config_path in tmp_path.glob("*.yml"):
            _configure_integration_yaml(config_path, gromacs_launch)
    assetsdir = Path(__file__).parent / "test_files" / "assets"
    if not (tmp_path / "amber99sb-star-ildnp.ff").exists():
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
