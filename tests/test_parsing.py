import re
import string
from pathlib import Path

import pytest
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

from kimmdy import parsing
from kimmdy.constants import AA3
from kimmdy.utils import get_gmx_dir


## test topology parser
def test_parser_doesnt_crash_on_example(arranged_tmp_path, caplog):
    """Example file urea.top
    from <https://manual.gromacs.org/documentation/current/reference-manual/topologies/topology-file-formats.html>
    """
    top_dict = parsing.read_top(Path("urea.top"))
    assert isinstance(top_dict, dict)


def test_doubleparse_urea(arranged_tmp_path):
    """Parsing it's own output should return the same top on urea.top"""
    top_dict = parsing.read_top(Path("urea.top"))

    top2_path = Path("pytest_urea.top")
    parsing.write_top(top_dict, top2_path)
    top2_dict = parsing.read_top(top2_path)

    top3_path = Path("pytest_urea2.top")
    parsing.write_top(top2_dict, top3_path)
    top3_dict = parsing.read_top(top3_path)

    assert top2_dict == top3_dict


@pytest.mark.skipif(
    not get_gmx_dir("gmx"),
    reason="Command 'gmx' not found, can't test gmx dir parsing.",
)
def test_ff_includes_with_gmxdir(arranged_tmp_path):
    top_dict = parsing.read_top(Path("urea.top"))
    gmx_dir = get_gmx_dir("gmx")
    assert gmx_dir is not None, "gmx dir not found"
    tip3_dict = parsing.read_top(gmx_dir / "top" / "amber99.ff" / "tip3p.itp")

    assert top_dict["atomtypes"]
    assert top_dict["bondtypes"]
    assert top_dict["moleculetype_SOL"] == tip3_dict["moleculetype_SOL"]


def test_ff_includes_with_ff_in_cwd(arranged_tmp_path):
    top_dict = parsing.read_top(Path("hexala.top"))
    ions_dict = parsing.read_top(Path("amber99sb-star-ildnp.ff/ions.itp"))
    assert top_dict["atomtypes"]
    assert top_dict["bondtypes"]
    for ion in [
        "IB+",
        "CA",
        "CL",
        "NA",
        "MG",
        "K",
        "RB",
        "CS",
        "LI",
        "ZN",
    ]:
        assert top_dict[f"moleculetype_{ion}"] == ions_dict[f"moleculetype_{ion}"]


# test whether topology parsing is invertible
allowed_text = st.text(
    string.ascii_letters + string.digits + "!\"$%&'()+,-./:<=>?@\\^_`{|}~", min_size=1
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
@settings(suppress_health_check=[HealthCheck.function_scoped_fixture], deadline=500)
def test_parser_invertible(sections, arranged_tmp_path):
    # flatten list of lists of strings to list of strings with subsection headers
    # use first element of each section as header
    for s in sections:
        header = s[0]
        header = re.sub(r"\d", "x", header)
        s[0] = f"[ {header} ]"
    ls = [l for s in sections for l in s]

    top_path = Path("topol.top")
    top2_path = Path("topol2.top")
    with open(top_path, "w") as f:
        f.write("\n".join(ls))
    top_dict = parsing.read_top(top_path)
    parsing.write_top(top_dict, top2_path)
    top2_dict = parsing.read_top(top2_path)
    assert top_dict == top2_dict


@given(ls=st.lists(allowed_text))
@settings(suppress_health_check=[HealthCheck.function_scoped_fixture], deadline=500)
def test_parser_fails_without_sections(ls, arranged_tmp_path):
    p = Path("random_content.top")
    with open(p, "w") as f:
        f.writelines(ls)
    with pytest.raises(ValueError):
        parsing.read_top(p)


## test ff file parsing
def test_parse_aminoacids_read_top():
    aminoacids_path = (
        Path(__file__).parent
        / "test_files"
        / "assets"
        / "amber99sb-star-ildnp.ff"
        / "aminoacids.rtp"
    )
    aminoacids_dict = parsing.read_top(aminoacids_path, use_gmx_dir=False)
    for aminoacid in AA3:
        assert (
            entry := aminoacids_dict.get(aminoacid)
        ), f"Aminoacid {aminoacid} not in {aminoacids_path.name}"
        ref_subsections = ["atoms", "bonds", "impropers"]
        subsections = list(entry["subsections"].keys())

        assert all(
            x in subsections for x in ref_subsections
        ), f"Aminoacid {aminoacid} does not have the subsections {ref_subsections} but {subsections}"
        assert all(len(x) == 4 for x in entry["subsections"]["atoms"]["content"])
        assert all(len(x) == 2 for x in entry["subsections"]["bonds"]["content"])
        assert all(
            len(x) in [4, 7] for x in entry["subsections"]["impropers"]["content"]
        )


def test_parse_ffbonded_read_top():
    ffbonded_path = (
        Path(__file__).parent
        / "test_files"
        / "assets"
        / "amber99sb-star-ildnp.ff"
        / "ffbonded.itp"
    )
    ffbonded_dict = parsing.read_top(ffbonded_path)
    ref_sections = ["bondtypes", "angletypes", "dihedraltypes"]
    ffbonded_sections = list(ffbonded_dict.keys())

    assert all(
        x in ffbonded_sections for x in ref_sections
    ), f"Sections {ref_sections} should be in ffbonded sections: {ffbonded_sections}"
    assert all(
        len(x) == 5 for x in ffbonded_dict["bondtypes"]["content"]
    ), "Unexpected number of elements in bondtypes"
    assert all(
        len(x) == 6 for x in ffbonded_dict["angletypes"]["content"]
    ), "Unexpected number of elements in angletypes"
    assert all(
        len(x) == 8 for x in ffbonded_dict["dihedraltypes"]["content"]
    ), "Unexpected number of elements in dihedraltypes"


## test plumed parsing
def test_plumed_read(arranged_tmp_path):
    plumed_dict = parsing.read_plumed(Path("plumed.dat"))

    assert set(plumed_dict.keys()) == set(["labeled_action", "prints", "other"])
    assert isinstance(plumed_dict["labeled_action"], dict)
    assert len(plumed_dict["labeled_action"]) == 12
    assert isinstance(plumed_dict["prints"], list)
    assert plumed_dict["prints"][0]["FILE"] == Path("distances.dat")


def test_plumed_write_identity(arranged_tmp_path):
    plumed_dict = parsing.read_plumed(Path("plumed.dat"))
    plumed_mod_path = Path("plumed_mod.dat")
    parsing.write_plumed(plumed_dict, plumed_mod_path)
    plumed_mod_dict = parsing.read_plumed(plumed_mod_path)

    assert plumed_mod_dict == plumed_dict


## test json parsing

## test misc file parsing


def test_edissoc_read(arranged_tmp_path):
    edissoc_dict = parsing.read_edissoc(Path("edissoc.dat"))
    assert len(edissoc_dict.keys()) == 14
    for k, v in edissoc_dict.items():
        for kv, vv in v.items():
            assert isinstance(kv, frozenset)
            assert isinstance(vv, float)


def test_marker_file_parsing(tmp_path: Path):
    parsing.write_time_marker(tmp_path / "marker_file1", "event1")
    parsing.write_time_marker(tmp_path / "marker_file2", "event1")
    parsing.write_time_marker(tmp_path / "marker_file2", "event2")
    parsing.write_time_marker(tmp_path / "marker_file2", "event3")
    parsing.write_time_marker(tmp_path / "marker_file2", "event2")

    es1, ts1 = parsing.read_time_marker(tmp_path / "marker_file1")
    es2, ts2 = parsing.read_time_marker(tmp_path / "marker_file2")

    assert len(es1) == len(ts1) == 1
    assert len(es2) == len(ts2) == 4

    assert "event1" in es1
    assert "event0" not in es1

    assert "event0" not in es2
    assert "event1" in es2
    assert "event2" in es2
    assert "event3" in es2

    assert (ts2[1] - ts2[0]).total_seconds() > 0
