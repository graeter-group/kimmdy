"""
Tests for homolysis reaction plugin.

Assumes KIMMDY was installed with the plugin.
E.g. via pip install -r requirements.txt.
"""
import pytest
from pathlib import Path
import numpy as np
from homolysis.reaction import Homolysis
from kimmdy.recipe import Break
from kimmdy.parsing import (
    read_plumed,
    read_top,
    read_distances_dat,
    read_edissoc,
    read_rtp,
)
from kimmdy.utils import (
    get_atominfo_from_plumedid,
    get_bondprm_from_atomtypes,
    morse_transition_rate,
)
from kimmdy.tasks import TaskFiles


@pytest.fixture
def homolysis_files():
    filedir = Path(__file__).parent / "test_files/test_homolysis"
    top = read_top(filedir / "topol.top")
    plumed = read_plumed(filedir / "plumed.dat")
    distances = read_distances_dat(filedir / "distances.dat")
    ffbonded = read_top(filedir / "ffbonded.itp")
    edissoc = read_edissoc(filedir / "edissoc.dat")

    files = {
        "top": top,
        "plumed": plumed,
        "distances": distances,
        "ffbonded": ffbonded,
        "edissoc": edissoc,
    }
    return files


def test_lookup_atominfo(homolysis_files):
    atomtypes, atomid = get_atominfo_from_plumedid(
        "d1", homolysis_files["plumed"], homolysis_files["top"]
    )
    assert atomtypes == frozenset(["N", "CT"])
    assert atomid == ["7", "9"]


def test_lookup_bondprm(homolysis_files):
    b0, kb, E_dis = get_bondprm_from_atomtypes(
        frozenset(["CT", "C"]), homolysis_files["ffbonded"], homolysis_files["edissoc"]
    )
    assert abs(b0 - 0.15220) < 1e-9
    assert abs(kb - 265265.6) < 1e-9
    assert abs(E_dis - 348) < 1e-9


def test_fail_lookup_bondprm(homolysis_files):
    with pytest.raises(KeyError):
        b0, kb, e_dis = get_bondprm_from_atomtypes(
            frozenset(["X", "Z"]),
            homolysis_files["ffbonded"],
            homolysis_files["edissoc"],
        )


def test_morse_transition_rate(homolysis_files):
    b0, kb, es_dis = get_bondprm_from_atomtypes(
        frozenset(["CT", "C"]), homolysis_files["ffbonded"], homolysis_files["edissoc"]
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


def test_get_recipe_collection(generic_rmgr):
    # curr_path = Path().cwd()

    files = TaskFiles(generic_rmgr)
    files.input["top"] = Path("topol.top")
    files.input["plumed"] = Path("plumed.dat")
    files.input["plumed_out"] = Path("distances.dat")
    files.input["edis"] = Path("edissoc.dat")
    files.input["itp"] = Path("ffbonded.itp")
    r = Homolysis(name="homolysis", runmng=generic_rmgr)
    rc = r.get_recipe_collection(files)

    plumed = read_plumed(files.input["plumed"])
    assert len(rc.recipes) == len(plumed["distances"])
    for recipe in rc.recipes:
        assert len(recipe.recipe_steps) == 1
        assert type(recipe.recipe_steps[0]) == type(Break(1, 2))
        assert len(recipe.rates) == 1
        assert type(recipe.rates[0]) in [float, np.float32, np.float64]
        assert len(recipe.timespans) == 1
        assert len(recipe.timespans[0]) == 2
        for time in recipe.timespans[0]:
            assert type(time) in [float, np.float32, np.float64]
