from kimmdy.parsing import read_top, write_top
from kimmdy.topology.topology import Topology
from kimmdy.topology.atomic import Bond
from kimmdy.coordinates import merge_top_parameter_growth, get_atomicobj
from conftest import SlimFiles


import pytest
from pathlib import Path


@pytest.fixture(scope="module")
def coordinates_files():
    filedir = Path(__file__).parent / "test_files" / "test_coordinates"
    ffdir = Path(__file__).parent / "test_files" / "assets" / "amber99sb-star-ildnp.ff"
    topA_path = filedir / "topol_stateA.top"
    topB_path = filedir / "topol_stateB.top"
    topFEP_path = filedir / "topol_FEP.top"
    if (filedir / "amber99sb-star-ildnp.ff").exists():
        (filedir / "amber99sb-star-ildnp.ff").unlink()
    (filedir / "amber99sb-star-ildnp.ff").symlink_to(
        ffdir,
        target_is_directory=True,
    )
    stateA = read_top(topA_path)
    stateB = read_top(topB_path)
    fep = read_top(topFEP_path)
    topA = Topology(stateA)
    topB = Topology(stateB)
    topFEP = Topology(fep)

    files = {
        "topA": topA,
        "topB": topB,
        "topFEP": topFEP,
        "topA_path": topA_path,
        "topB_path": topB_path,
        "topFEP_path": topFEP,
    }
    yield files
    (filedir / "amber99sb-star-ildnp.ff").unlink()


def test_get_bondobj(coordinates_files):
    bond1_keys = ["17", "18"]
    bond1obj = get_atomicobj(bond1_keys, Bond, coordinates_files["topA"])

    bond2_keys = ["17", "19"]
    bond2obj = get_atomicobj(bond2_keys, Bond, coordinates_files["topA"])
    assert float(bond1obj.c0) == pytest.approx(0.10100) and float(
        bond1obj.c1
    ) == pytest.approx(363171.2)

    assert float(bond2obj.c0) == pytest.approx(0.13600) and float(
        bond2obj.c1
    ) == pytest.approx(282001.6)


def test_merge_prm_top(coordinates_files):
    topmerge = merge_top_parameter_growth(
        coordinates_files["topA"], coordinates_files["topB"]
    )
    # topmerge_path = coordinates_files["topB_path"].parent / "top_merge.top"
    # topdict = topmerge.to_dict()
    # write_top(topdict, topmerge_path)
    # topFEP does not work as a reference, the file must be changed for this test to work
    # try:
    assert topmerge.atoms == coordinates_files["topFEP"].atoms
    assert topmerge.bonds == coordinates_files["topFEP"].bonds
    assert topmerge.angles == coordinates_files["topFEP"].angles
    assert topmerge.pairs == coordinates_files["topFEP"].pairs
    # finally:
    #     topmerge_path.unlink()
