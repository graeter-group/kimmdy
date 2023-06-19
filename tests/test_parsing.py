# %%
# %autoreload
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


# %%
#### Example file urea.gro ####
# from <https://manual.gromacs.org/documentation/current/reference-manual/topologies/topology-file-formats.html>
# should parse
def test_parser_doesnt_crash_on_example():
    set_dir()
    urea_path = Path("urea.gro")
    top = parsing.read_top(urea_path)
    assert isinstance(top, dict)


# %%
def test_doubleparse_urea():
    """Parsing it's own output should return the same top on urea.gro"""
    set_dir()
    urea_path = Path("urea.gro")
    top = parsing.read_top(urea_path)
    p = Path("pytest_urea.top")
    parsing.write_top(top, p)
    top2 = parsing.read_top(p)
    p2 = Path("pytest_urea2.top")
    parsing.write_top(top2, p2)
    top3 = parsing.read_top(p2)
    assert top2 == top3


def test_parsing_includes_as_blocks():
    set_dir()
    urea_path = Path("urea.gro")
    top = parsing.read_top(urea_path)
    assert top["includes"] is not None


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
    parsing.write_top(d, p)
    d2 = parsing.read_top(p)
    assert d == d2


# %%
def test_parse_xml_ff():
    set_dir()
    ff_path = Path("amber99sb_trunc.xml")
    xml = parsing.read_xml_ff(ff_path)

    atomtypes = xml.find("AtomTypes")
    atomtypes.findall("Type")
    assert atomtypes.findall("Type")[0].attrib == {
        "class": "N",
        "element": "N",
        "mass": "14.00672",
    }
    assert atomtypes.findall("Type")[1].attrib == {
        "class": "H",
        "element": "H",
        "mass": "1.007947",
    }
