# %%
# %autoreload
import os
import re
import string
from hypothesis import given, strategies as st
from kimmdy import parsing
from pathlib import Path
from copy import deepcopy
import pytest


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


#### Parsing should be invertible ####
allowed_text = st.text(
    string.ascii_letters + string.digits + "!\"$%&'()*+,-./:<=>?@\\^_`{|}~", min_size=1
)


@given(
    # a list of lists that correspond to a sections of a top file
    sections=st.lists(
        st.lists(
            allowed_text,
            min_size=2,
            max_size=5,
        ),
        min_size=1,
        max_size=5,
    )
)
def test_parser_invertible(sections):
    # flatten list of lists of strings to list of strings with subsection headers
    # use first element of each section as header
    for s in sections:
        header = s[0]
        header = re.sub(r"\d", "x", header)
        s[0] = f"[ {header} ]"
    ls = [l for s in sections for l in s]
    print(ls)
    p = Path("tmp/pytest_topol.top")
    p2 = Path("tmp/pytest_topol2.top")
    p.parent.mkdir(exist_ok=True)
    with open(p, "w") as f:
        f.write("\n".join(ls))
    top = parsing.read_top(p)
    parsing.write_top(top, p2)
    top2 = parsing.read_top(p2)
    assert top == top2


@given(ls=st.lists(allowed_text))
def test_parser_fails_without_sections(ls):
    p = Path("tmp/pytest_topol.top")
    p.parent.mkdir(exist_ok=True)
    with open(p, "w") as f:
        f.writelines(ls)
    with pytest.raises(ValueError):
        parsing.read_top(p)
    assert True


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
