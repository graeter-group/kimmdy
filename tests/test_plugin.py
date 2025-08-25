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
