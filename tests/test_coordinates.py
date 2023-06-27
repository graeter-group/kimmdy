from kimmdy.parsing import read_top, write_top
from kimmdy.tasks import TaskFiles
from kimmdy.topology.topology import Topology, TopologyDict
from kimmdy.topology.atomic import Bond
from kimmdy.coordinates import merge_top_prmgrowth, get_bondobj

import pytest
from pathlib import Path


@pytest.fixture(scope="module")
def coordinates_files():
    filedir = Path(__file__).parent / "test_files" / "test_coordinates"
    ffdir = Path(__file__).parent / "test_files" / "assets" / "amber99sb-star-ildnp.ff"
    topA_path = filedir / "topol_stateA.top"
    topB_path = filedir / "topol_stateB.top"
    stateA = read_top(topA_path)
    stateB = read_top(topB_path)
    fep = read_top(filedir / "topol_FEP.top")
    topA = Topology(stateA, ffdir)
    topB = Topology(stateB, ffdir)
    topFEP = Topology(fep, ffdir)

    files = {
        "topA": topA,
        "topB": topB,
        "topFEP": topFEP,
        "topA_path": topA_path,
        "topB_path": topB_path,
    }
    return files


def test_get_bondobj(coordinates_files):
    bond1_keys = ["17", "18"]
    bond1 = Bond.from_top_line("17       18       1   ".split())
    bond1obj = get_bondobj(bond1_keys, bond1, coordinates_files["topA"])

    bond2_keys = ["17", "19"]
    bond2 = Bond.from_top_line("17       19       1        0.13600  282001.6".split())
    bond2obj = get_bondobj(bond2_keys, bond2, coordinates_files["topA"])
    assert bond1obj.c0 == "0.10100" and bond1obj.c1 == "363171.2"
    assert bond2obj.c0 == "0.13600" and bond2obj.c1 == "282001.6"


def test_merge_prm_top(generic_rmgr, coordinates_files):
    # needs a proper runmgr
    files = TaskFiles(generic_rmgr)
    files.input["top"] = coordinates_files["topA_path"]
    files.output["top"] = coordinates_files["topB_path"]
    topmerge = merge_top_prmgrowth(files)
    topdict = topmerge.to_dict()
    write_top(
        topdict, coordinates_files["topB_path"].parent / "top_merge.top"
    )
    # topFEP does not work as a reference, the file must be changed for this test to work
    assert topmerge == coordinates_files["topFEP"]