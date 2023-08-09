from kimmdy.parsing import read_top, write_top
from kimmdy.topology.topology import Topology
from kimmdy.topology.atomic import Bond
from kimmdy.coordinates import merge_top_parameter_growth, get_explicit_or_type
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
        "fep": fep,
    }
    yield files
    (filedir / "amber99sb-star-ildnp.ff").unlink()


def test_get_bondobj(coordinates_files):
    bond1_keys = ("17", "18")
    bond1obj = get_explicit_or_type(
        bond1_keys,
        coordinates_files["topA"].bonds[bond1_keys],
        coordinates_files["topA"].ff.bondtypes,
        coordinates_files["topA"],
    )

    bond2_keys = ("17", "19")
    bond2obj = get_explicit_or_type(
        bond2_keys,
        coordinates_files["topA"].bonds[bond2_keys],
        coordinates_files["topA"].ff.bondtypes,
        coordinates_files["topA"],
    )
    assert float(bond1obj.c0) == pytest.approx(0.10100) and float(
        bond1obj.c1
    ) == pytest.approx(363171.2)

    assert float(bond2obj.c0) == pytest.approx(0.13600) and float(
        bond2obj.c1
    ) == pytest.approx(282001.6)


def test_merge_prm_top(coordinates_files):
    """this tests a topology merge for a HAT reaction from a Ca (nr 19) radical to a N (nr 26) radical"""
    topmerge = merge_top_parameter_growth(
        coordinates_files["topA"], coordinates_files["topB"]
    )

    # write_top(
    #     topmerge.to_dict(),
    #     Path(
    #         "/hits/fast/mbm/hartmaec/kimmdys/kimmdy_main/tests/test_files/test_coordinates/topol_curr.top"
    #     ),
    # )

    assert topmerge.atoms == coordinates_files["topFEP"].atoms
    assert topmerge.bonds.keys() == coordinates_files["topFEP"].bonds.keys()
    assert topmerge.angles.keys() == coordinates_files["topFEP"].angles.keys()
    assert topmerge.pairs.keys() == coordinates_files["topFEP"].pairs.keys()
    assert (
        topmerge.proper_dihedrals.keys()
        == coordinates_files["topFEP"].proper_dihedrals.keys()
    )
    assert (
        topmerge.improper_dihedrals.keys()
        == coordinates_files["topFEP"].improper_dihedrals.keys()
    )

    assert topmerge.bonds[("19", "27")].funct == "3"
    assert topmerge.bonds[("26", "27")].funct == "3"
    assert topmerge.angles[("17", "19", "20")].c3 != None
    assert topmerge.proper_dihedrals[("15", "17", "19", "24")].dihedrals["3"].c5 == "3"
    assert topmerge.improper_dihedrals[("17", "20", "19", "24")].c5 == "2"
    # asster one dihedral merge improper/proper
