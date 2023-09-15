"""
Tests for default_reactions reaction plugin.

Is skipped if KIMMDY was installed without the plugin.
"""
import pytest
from kimmdy import reaction_plugins
from kimmdy.runmanager import RunManager

pytest.importorskip("homolysis")
pytest.importorskip("hat_naive")
pytest.importorskip("dummyreaction")

from dataclasses import dataclass
from typing import Callable
import os
import shutil
from pathlib import Path
import numpy as np
from homolysis.reaction import Homolysis
from kimmdy.config import Config
from kimmdy.recipe import Break
from kimmdy.topology.topology import Topology
from kimmdy.parsing import (
    read_plumed,
    read_top,
    read_distances_dat,
    read_edissoc,
)
from kimmdy.utils import (
    get_atomnrs_from_plumedid,
    get_atomtypes_from_atomnrs,
    get_bondprm_from_atomtypes,
    morse_transition_rate,
)
from kimmdy.tasks import TaskFiles


@dataclass
class DummyRunmanager(RunManager):
    top: Topology
    config: Config


@dataclass
class DummyFiles(TaskFiles):
    get_latest: Callable = lambda: f"DummyCallable"


@pytest.fixture
def homolysis_files(tmp_path: Path):
    file_dir = Path(__file__).parent / "test_default_reactions"
    shutil.copytree(file_dir, tmp_path, dirs_exist_ok=True)
    os.chdir(tmp_path)
    top = Topology(read_top(Path("topol.top")))
    plumed = read_plumed(Path("plumed.dat"))
    distances = read_distances_dat(Path("distances.dat"))
    ffbonded = read_top(Path("ffbonded.itp"))
    edissoc = read_edissoc(Path("edissoc.dat"))

    assetsdir = Path(__file__).parents[2] / "tests" / "test_files" / "assets"
    Path(tmp_path / "amber99sb-star-ildnp.ff").symlink_to(
        assetsdir / "amber99sb-star-ildnp.ff",
        target_is_directory=True,
    )

    files = {
        "top": top,
        "plumed": plumed,
        "distances": distances,
        "ffbonded": ffbonded,
        "edissoc": edissoc,
    }

    return files


## test homolysis
def test_get_atomnrs(homolysis_files):
    atomnrs = get_atomnrs_from_plumedid("d1", homolysis_files["plumed"])
    assert atomnrs == ["7", "9"]


def test_get_atomntypes(homolysis_files):
    atomnrs = ["7", "9"]
    atomtypes = get_atomtypes_from_atomnrs(atomnrs, homolysis_files["top"])
    assert atomtypes == ["N", "CT"]


def test_lookup_bondprm(homolysis_files):
    b0, kb, E_dis = get_bondprm_from_atomtypes(
        ["CT", "C"], homolysis_files["ffbonded"], homolysis_files["edissoc"]
    )
    assert abs(b0 - 0.15220) < 1e-9
    assert abs(kb - 265265.6) < 1e-9
    assert abs(E_dis - 348) < 1e-9


def test_fail_lookup_bondprm(homolysis_files):
    with pytest.raises(
        KeyError, match="Did not find dissociation energy for atomtypes"
    ):
        b0, kb, e_dis = get_bondprm_from_atomtypes(
            ["X", "Z"],
            homolysis_files["ffbonded"],
            homolysis_files["edissoc"],
        )


def test_morse_transition_rate(homolysis_files):
    b0, kb, es_dis = get_bondprm_from_atomtypes(
        ["CT", "C"], homolysis_files["ffbonded"], homolysis_files["edissoc"]
    )

    rs_ref = list(np.linspace(0.9, 1.3, 8) * b0)
    ks, fs = morse_transition_rate(rs_ref, b0, es_dis, kb)

    ks_ref = np.asarray(
        [
            0.00000000e00,
            0.00000000e00,
            6.34487872e-42,
            1.01004484e-11,
            2.67853959e-03,
            2.07275944e-01,
            2.88000000e-01,
            2.88000000e-01,
        ]
    )
    fs_ref = np.asarray(
        [
            -6327.85756373,
            -2095.88999762,
            541.22525767,
            2101.46367022,
            2944.48244374,
            3318.63831687,
            3396.91825041,
            3396.91825041,
        ]
    )

    assert all(np.isclose(ks, ks_ref))
    assert all(np.isclose(fs, fs_ref))


def test_get_recipe_collection(homolysis_files):
    config = Config(Path("kimmdy.yml"))
    rmgr = DummyRunmanager(homolysis_files["top"], config)

    files = DummyFiles()
    files.input["plumed"] = Path("plumed.dat")
    files.input["plumed_out"] = Path("distances.dat")
    files.input["edis"] = Path("edissoc.dat")
    files.input["itp"] = Path("ffbonded.itp")
    r = Homolysis(name="homolysis", runmng=rmgr)
    rc = r.get_recipe_collection(files)

    plumed = read_plumed(files.input["plumed"])
    assert len(rc.recipes) == len(plumed["labeled_action"])
    for recipe in rc.recipes:
        assert len(recipe.recipe_steps) == 2
        assert type(recipe.recipe_steps[0]) == type(Break(1, 2))
        assert len(recipe.rates) == 1
        assert type(recipe.rates[0]) in [float, np.float32, np.float64]
        assert len(recipe.timespans) == 1
        assert len(recipe.timespans[0]) == 2
        for time in recipe.timespans[0]:
            assert type(time) in [float, np.float32, np.float64]


## test hat_naive
# TODO:

## test dummyreaction
# TODO:
