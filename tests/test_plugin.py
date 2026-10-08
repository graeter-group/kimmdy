from kimmdy import plugin_utils
from kimmdy.parsing import read_top
from kimmdy.topology.topology import Topology
from pathlib import Path
import pytest

# backbone bonds of hexalanine_homolysis as in its former plumed.dat
# (backbone minus "C-N H* O*", the kimmdy-create-plumed default)
HEXALA_BACKBONE_BONDS = {
    ("7", "9"),
    ("9", "15"),
    ("17", "19"),
    ("19", "25"),
    ("27", "29"),
    ("29", "35"),
    ("37", "39"),
    ("39", "45"),
    ("47", "49"),
    ("49", "55"),
    ("57", "59"),
    ("59", "65"),
}


def test_bondstate_csv_io(arranged_tmp_path):
    bondstat = plugin_utils.bondstats_from_csv(".kimmdy.bondstats")
    assert isinstance(bondstat, dict)
    assert all(
        [
            k in d.keys()
            for k in ["b0", "delta_d", "harmonic_f", "mean_d", "mean_f", "plumed_id"]
            for d in bondstat.values()
        ]
    )

    out = arranged_tmp_path / "bondstats_out"
    plugin_utils.bondstats_to_csv(bondstat, out)

    with open(out, "r") as f:
        written_lines = f.readlines()
    with open(".kimmdy.bondstats", "r") as f:
        original_lines = f.readlines()
    assert all([w == r for w, r in zip(written_lines, original_lines)])


def test_bondstats_csv_without_harmonic_force(tmp_path):
    path = tmp_path / "bondstats_7col"
    path.write_text(
        "ai,aj,plumed_id,mean_d,mean_f,delta_d,b0,harmonic_f\n"
        "7,9,d0,0.146791,422.989854,0.001891,0.144900\n"
    )
    bondstats = plugin_utils.bondstats_from_csv(path)
    assert bondstats[("7", "9")]["harmonic_f"] is None
    assert bondstats[("7", "9")]["mean_d"] == pytest.approx(0.146791)

    out = tmp_path / "bondstats_out"
    plugin_utils.bondstats_to_csv(bondstats, out)
    assert out.read_text() == path.read_text().rstrip("\n")


@pytest.mark.parametrize(
    "arranged_tmp_path", ["test_integration/hexalanine_homolysis"], indirect=True
)
def test_select_bonds_matches_create_plumed_defaults(arranged_tmp_path):
    top = Topology(read_top(Path("hexala_out.top")))
    backbone = {
        nr for nr, atom in top.atoms.items() if atom.atom in {"N", "CA", "C", "O"}
    }
    bonds = plugin_utils.select_bonds(
        top=top, atom_nrs=backbone, exclude=["C-N", "H*", "O*"]
    )
    assert set(bonds) == HEXALA_BACKBONE_BONDS


@pytest.mark.parametrize(
    "arranged_tmp_path", ["test_integration/hexalanine_homolysis"], indirect=True
)
def test_select_bonds_exclusions(arranged_tmp_path):
    top = Topology(read_top(Path("hexala_out.top")))
    assert len(plugin_utils.select_bonds(top=top)) == len(top.bonds)

    bonds = plugin_utils.select_bonds(top=top, exclude=["H*", "O"])
    names = [{top.atoms[a].atom, top.atoms[b].atom} for a, b in bonds]
    assert names
    assert not any(n.startswith("H") or n == "O" for pair in names for n in pair)
    assert {"C", "N"} in names

    bonds = plugin_utils.select_bonds(top=top, exclude=["C-N"])
    names = [{top.atoms[a].atom, top.atoms[b].atom} for a, b in bonds]
    assert {"C", "N"} not in names
