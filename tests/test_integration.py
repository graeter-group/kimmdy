import os
import subprocess as sp
from pathlib import Path

import pytest

from kimmdy.cmd import kimmdy_run
from kimmdy.constants import MARK_DONE, MARK_FINISHED
from kimmdy.parsing import read_top, write_top
from kimmdy.plugins import parameterization_plugins
from kimmdy.topology.topology import Topology
from kimmdy.utils import get_task_directories


def read_last_line(file):
    with open(file, "rb") as f:
        try:  # catch OSError in case of a one line file
            f.seek(-2, os.SEEK_END)
            while f.read(1) != b"\n":
                f.seek(-2, os.SEEK_CUR)
        except OSError:
            f.seek(0)
        return f.readline().decode()


@pytest.mark.parametrize(
    "arranged_tmp_path", (["test_integration/emptyrun"]), indirect=True
)
def test_integration_emptyrun(arranged_tmp_path):
    Path("emptyrun.txt").touch()
    with pytest.raises(ValueError):
        kimmdy_run()
    assert len(list(Path.cwd().glob("emptyrun_001/*"))) == 3
    assert not (arranged_tmp_path / "minimal" / MARK_FINISHED).exists()


@pytest.mark.parametrize(
    "arranged_tmp_path", (["test_integration/minimal_input_files"]), indirect=True
)
def test_integration_valid_input_files(arranged_tmp_path):
    kimmdy_run()
    assert "Finished running last task" in read_last_line(Path("minimal.kimmdy.log"))
    assert (arranged_tmp_path / "minimal" / MARK_FINISHED).exists()
    assert len(list(Path.cwd().glob("minimal/*"))) == 4


@pytest.mark.parametrize(
    "arranged_tmp_path", (["test_integration/minimal_input_files"]), indirect=True
)
def test_grompp_with_kimmdy_topology(arranged_tmp_path):
    raw_top = read_top(Path("minimal.top"))
    top = Topology(raw_top)
    top_dict = top.to_dict()
    write_top(top.to_dict(), Path("output.top"))
    assert sp.run(
        [
            "gmx",
            "grompp",
            "-f",
            "minimal.mdp",
            "-c",
            "minimal.gro",
            "-p",
            "output.top",
            "-o",
            "minial.tpr",
        ],
        check=True,
    )


@pytest.mark.require_grappa
@pytest.mark.parametrize(
    "arranged_tmp_path", (["test_integration/triplehelix_pull"]), indirect=True
)
def test_grappa_partial_parameterization(arranged_tmp_path):
    parametrizer = parameterization_plugins["grappa"]()

    # parameterize everything
    top = Topology(read_top(Path("topol.top")), parametrizer=parametrizer)
    top.needs_parameterization = True
    top.update_parameters()

    # parameterize around delete
    top_del = Topology(read_top(Path("topol.top")), parametrizer=parametrizer)
    top_del.needs_parameterization = True
    top_del.update_parameters()
    top_del.break_bond(("1", "5"))
    top_del.needs_parameterization = True
    top_del.parameterization_focus_ids = set(["1", "5"])
    top_del.update_parameters()

    # check values are changed around homolysis atoms
    assert top.bonds[("5", "6")] != top_del.bonds[("5", "6")]
    assert top.angles[("6", "5", "7")] != top_del.angles[("6", "5", "7")]
    assert (
        top.proper_dihedrals[("13", "11", "14", "15")]
        != top_del.proper_dihedrals[("13", "11", "14", "15")]
    )
    assert (
        top.improper_dihedrals[("8", "5", "7", "9")]
        != top_del.improper_dihedrals[("8", "5", "7", "9")]
    )

    # check values are the same further away
    assert top.bonds[("200", "203")] == top_del.bonds[("200", "203")]
    assert top.angles[("228", "230", "232")] == top_del.angles[("228", "230", "232")]
    assert (
        top.proper_dihedrals[("1714", "1713", "1715", "1717")]
        == top_del.proper_dihedrals[("1714", "1713", "1715", "1717")]
    )
    assert (
        top.improper_dihedrals[("349", "357", "355", "356")]
        == top_del.improper_dihedrals[("349", "357", "355", "356")]
    )


@pytest.mark.parametrize(
    "arranged_tmp_path",
    (["test_integration/hexalanine_single_reaction"]),
    indirect=True,
)
def test_integration_single_reaction(arranged_tmp_path):
    kimmdy_run(input=Path("kimmdy.yml"))
    assert "Finished running last task" in read_last_line(
        Path("single_reaction_000.kimmdy.log")
    )
    assert len(list(Path.cwd().glob("single_reaction_000/*"))) == 8


@pytest.mark.parametrize(
    "arranged_tmp_path",
    (["test_integration/hexalanine_single_reaction"]),
    indirect=True,
)
def test_integration_just_reactions(arranged_tmp_path):
    kimmdy_run(input=Path("alternative_kimmdy.yml"))
    assert "Finished running last task" in read_last_line(
        Path("single_reaction_000.kimmdy.log")
    )
    assert len(list(Path.cwd().glob("single_reaction_000/*"))) == 8


@pytest.mark.slow
@pytest.mark.parametrize(
    "arranged_tmp_path", (["test_integration/alanine_hat_naive"]), indirect=True
)
def test_integration_hat_naive_reaction(arranged_tmp_path):
    kimmdy_run()
    assert "Finished running last task" in read_last_line(
        Path("alanine_hat_000.kimmdy.log")
    )
    assert len(list(Path.cwd().glob("alanine_hat_000/*"))) == 16


@pytest.mark.slow
@pytest.mark.parametrize(
    "arranged_tmp_path", (["test_integration/hexalanine_homolysis"]), indirect=True
)
def test_integration_homolysis_reaction(arranged_tmp_path):
    kimmdy_run()
    assert "Finished running last task" in read_last_line(
        Path("hexalanine_homolysis_000.kimmdy.log")
    )
    assert len(list(Path.cwd().glob("hexalanine_homolysis_000/*"))) == 13


@pytest.mark.slow
@pytest.mark.parametrize(
    "arranged_tmp_path", (["test_integration/triplehelix_pull"]), indirect=True
)
def test_integration_pull(arranged_tmp_path):
    kimmdy_run()
    assert "Finished running last task" in read_last_line(Path("kimmdy_001.kimmdy.log"))
    assert len(list(Path.cwd().glob("kimmdy_001/*"))) == 12


@pytest.mark.require_grappa
@pytest.mark.slow
@pytest.mark.parametrize(
    "arranged_tmp_path",
    (["test_integration/charged_peptide_homolysis_hat_naive"]),
    indirect=True,
)
def test_integration_whole_run(arranged_tmp_path):
    kimmdy_run()
    assert "Finished running last task" in read_last_line(Path("kimmdy_001.kimmdy.log"))
    assert len(list(Path.cwd().glob("kimmdy_001/*"))) == 25


@pytest.mark.slow
@pytest.mark.parametrize(
    "arranged_tmp_path", (["test_integration/alanine_hat_naive"]), indirect=True
)
def test_integration_restart(arranged_tmp_path):
    run_dir = Path("alanine_hat_000")
    kimmdy_run(input=Path("kimmdy_restart.yml"))
    n_files_original = len(list(run_dir.glob("*")))

    # restart already finished run
    kimmdy_run(input=Path("kimmdy_restart.yml"))
    assert "already finished" in read_last_line(Path("alanine_hat_000.kimmdy.log"))

    # try restart from stopped md
    task_dirs = get_task_directories(run_dir)
    (task_dirs[-1] / MARK_DONE).unlink()
    (arranged_tmp_path / run_dir / MARK_FINISHED).unlink()
    kimmdy_run(input=Path("kimmdy_restart.yml"))
    n_files_continue_md = len(list(run_dir.glob("*")))

    assert "Finished running last task" in read_last_line(
        Path("alanine_hat_000.kimmdy.log")
    )
    assert n_files_original == n_files_continue_md == 17

    # try restart from finished md
    task_dirs = get_task_directories(run_dir)
    (task_dirs[-4] / MARK_DONE).unlink()
    (arranged_tmp_path / run_dir / MARK_FINISHED).unlink()
    kimmdy_run(input=Path("kimmdy_restart.yml"))
    n_files_restart = len(list(run_dir.glob("*")))

    assert "Finished running last task" in read_last_line(
        Path("alanine_hat_000.kimmdy.log")
    )
    assert n_files_original == n_files_restart == 17


@pytest.mark.slow
@pytest.mark.parametrize(
    "arranged_tmp_path", (["test_integration/alanine_hat_naive"]), indirect=True
)
def test_integration_file_usage(arranged_tmp_path):
    """Do a kimmdy run
    and verify that at each task the correct files are used
    specifically when writing out gro and trr files when
    applying a reaction.
    """
    kimmdy_run()
    histfile = Path("alanine_hat_000/kimmdy.history").read_text().split("\n\n")
    header, blocks = histfile[0], histfile[1:]
    tasks = {}
    for block in blocks:
        lines = block.split("\n")
        name = lines[0].removeprefix("Task: ")
        tasks[name] = {"input": {}, "output": {}}
        section = "input"
        for l in lines[1:]:
            if l.startswith("Input:"):
                section = "input"
            elif l.startswith("Output:"):
                section = "output"
            elif l.startswith(" "):
                tasks[name][section][l.split(": ")[0].strip()] = l.split(": ")[
                    1
                ].strip()

    assert tasks["5_apply_recipe"]["output"]["gro"] == "6_relax/relax.gro"
    assert tasks["5_apply_recipe"]["output"]["trr"] == "6_relax/relax.trr"
    assert tasks["5_apply_recipe"]["output"]["xtc"] == "6_relax/relax.xtc"
    assert tasks["5_apply_recipe"]["output"]["top"] == "5_apply_recipe/Ala_out.top"

    assert tasks["6_relax"]["input"]["top"] == "5_apply_recipe/Ala_out_relax.top"
    assert tasks["6_relax"]["input"]["gro"] == "2_equilibrium/.kimmdy_reaction.gro"
    assert tasks["6_relax"]["input"]["trr"] == "2_equilibrium/equilibrium.trr"
    assert tasks["6_relax"]["output"]["trr"] == "6_relax/relax.trr"
    assert tasks["6_relax"]["output"]["xtc"] == "6_relax/relax.xtc"

    assert tasks["7_equilibrium"]["input"]["gro"] == "6_relax/relax.gro"
    assert tasks["7_equilibrium"]["input"]["trr"] == "6_relax/relax.trr"
    assert tasks["7_equilibrium"]["output"]["trr"] == "7_equilibrium/equilibrium.trr"
    assert tasks["7_equilibrium"]["output"]["cpt"] == "7_equilibrium/equilibrium.cpt"
    assert tasks["7_equilibrium"]["output"]["xtc"] == "7_equilibrium/equilibrium.xtc"
