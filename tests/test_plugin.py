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
