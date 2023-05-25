from kimmdy.reaction import *
import csv
import pytest
from dataclasses import asdict


def test_no_double_frames():
    with pytest.raises(ValueError):
        ReactionPath([Conversion()], rates=[1.0, 2.0, 3.0], frames=[1, 1, 2])


def test_combine_reaction_paths():
    rp1a = ReactionPath([Conversion()], rates=[1], frames=[1])
    rp1b = ReactionPath([Conversion()], rates=[1], frames=[2])
    rp2 = ReactionPath([Conversion(), Conversion()], rates=[1], frames=[3])

    rp1a.combine_with(rp1b)
    assert rp1a.frames == [1, 2]
    assert rp1a.rates == [1, 1]

    with pytest.raises(ValueError):
        rp1a.combine_with(rp2)


@pytest.fixture
def reaction_result():
    rps = [
        ReactionPath([Conversion(), Conversion(), Conversion()], rates=[1], frames=[0]),
        ReactionPath([Conversion()], rates=[1], frames=[1]),
        ReactionPath([Conversion()], rates=[1], frames=[2]),
        ReactionPath([Conversion(), Conversion()], rates=[1], frames=[3]),
        ReactionPath([Conversion()], rates=[1], frames=[4]),
        ReactionPath([Conversion(), Conversion()], rates=[1], frames=[5]),
    ]
    return ReactionResults(rps)
    

def test_aggregate_reaction_result(reaction_result):
    reaction_result.aggregate_reactions()

    assert len(reaction_result.reaction_paths) == 3
    assert reaction_result.reaction_paths[1].frames == [1, 2, 4]
    assert reaction_result.reaction_paths[2].frames == [3, 5]


def test_reaction_result_to_csv(tmp_path, reaction_result):
    csv_path = tmp_path / "test_out.csv"
    reaction_result.to_csv(csv_path)
    with open(csv_path, newline="") as f:
        reader = csv.DictReader(f)
        rows = [r for r in reader]
    for read_rp, org_rp in zip(rows, reaction_result.reaction_paths):
        org_rp_d = asdict(org_rp)
        keys_to_check = ["rates", "frames", "avg_rates", "avg_frames"]
        for key in keys_to_check:
            org_val = str(org_rp_d[key])
            if org_val == "None":
                org_val = ""
            assert org_val == read_rp[key], f"{org_rp_d[key]} != {read_rp[key]}"


def test_reaction_result_to_dill(tmp_path, reaction_result):
    csv_path = tmp_path / "test_out.dill"
    reaction_result.to_dill(csv_path)
    loaded_rr = ReactionResults.from_dill(csv_path)
    assert loaded_rr == reaction_result



