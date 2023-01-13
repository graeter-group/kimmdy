#%%
#%autoreload
import os
import string
from hypothesis import given, strategies as st
from kimmdy import parsing
from pathlib import Path
from copy import deepcopy


def set_dir():
    try:
        test_dir = Path(__file__).parent / "test_files/test_parsing"
    except NameError:
        test_dir = Path("./tests/test_files/test_parsing")
    os.chdir(test_dir)


set_dir()

#%%
#### Example file urea.gro ####
# from <https://manual.gromacs.org/documentation/current/reference-manual/topologies/topology-file-formats.html>
# should parse
def test_parser_doesnt_crash_on_example():
    set_dir()
    urea_path = Path("urea.gro")
    top = parsing.read_topol(urea_path)
    assert isinstance(top, dict)


#%%
#### Parsing it's own output should return the same top on urea.gro ####
def test_doubleparse_urea():
    set_dir()
    urea_path = Path("urea.gro")
    top = parsing.read_topol(urea_path)
    p = Path("pytest_urea.top")
    parsing.write_topol(top, p)
    top2 = parsing.read_topol(p)
    p2 = Path("pytest_urea2.top")
    parsing.write_topol(top2, p2)
    top3 = parsing.read_topol(p2)
    assert top2 == top3


#%%
#### Parsing should be invertible ####
allowed_text = st.text(
    string.ascii_letters + string.digits + "!\"$%&'()*+,-./:<=>?@\\^_`{|}~", min_size=1
)


@given(
    d=st.dictionaries(
        allowed_text,
        st.lists(st.lists(allowed_text, min_size=1), min_size=1),
        min_size=1,
    )
)
def test_parser_invertible(d):
    p = Path("tmp/pytest_topol.top")
    p.parent.mkdir(exist_ok=True)
    parsing.write_topol(d, p)
    d2 = parsing.read_topol(p)
    assert d == d2


#%%
def test_parse_xml_ff():
    set_dir()
    ff_path = Path("amber99sb_trunc.xml")
    d = parsing.read_xml_ff(ff_path)
    refdict = {
        "HEAD": {},
        "AtomTypes": {
            "N": {"element": "N", "mass": 14.00672},
            "H": {"element": "H", "mass": 1.007947},
        },
        "HarmonicBondForce": {
            "H_N": {"length": 0.101, "k": 363171.2},
            "C_*_*_O": {"periodicity1": 2.0, "phase1": 3.14159265359, "k1": 43.932},
        },
        "NonbondedForce": {
            "N": {"charge": -0.4157, "sigma": 0.324999852378, "epsilon": 0.71128},
            "H": {"charge": 0.2719, "sigma": 0.106907846177, "epsilon": 0.0656888},
        },
    }
    assert d == refdict
