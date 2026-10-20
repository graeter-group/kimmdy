from kimmdy import plugin_utils
import pytest


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


def test_edissoc_layered_fallback():
    edissoc = {
        "general": {("C", "N"): 500.0, ("CA", "CB"): 350.0},
        "ALA": {("C", "N"): 510.0},
    }
    # residue section first
    assert plugin_utils.get_edissoc_from_atomnames(["N", "C"], edissoc, "ALA") == 510.0
    # residue section exists, but pair only in general
    assert (
        plugin_utils.get_edissoc_from_atomnames(["CB", "CA"], edissoc, "ALA") == 350.0
    )
    # residue section missing
    assert (
        plugin_utils.get_edissoc_from_atomnames(["CA", "CB"], edissoc, "GLY") == 350.0
    )
    # pair nowhere
    assert plugin_utils.get_edissoc_from_atomnames(["C", "O"], edissoc, "ALA") == 400.0
    # residue section missing and no general section
    assert (
        plugin_utils.get_edissoc_from_atomnames(["C", "N"], {"ALA": {}}, "GLY") == 400.0
    )
