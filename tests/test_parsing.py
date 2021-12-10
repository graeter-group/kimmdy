from hypothesis import given, strategies as st
import string
from kimmdy import parsing
from pathlib import Path

allowed_text = st.text(string.ascii_letters +
                       string.digits +
                       '!"$%&\'()*+,-./:<=>?@\\^_`{|}~', min_size=1)

@given(d=st.dictionaries(allowed_text, st.lists(allowed_text, min_size=1)))
def test_parser_invertible(d):
    p = Path("/tmp/pytest_topol.top")
    parsing.write_topol(d, p)
    d2 = parsing.read_topol(p)
    assert d == d2

