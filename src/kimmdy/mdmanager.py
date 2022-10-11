from __future__ import annotations
import logging
from kimmdy.tasks import TaskFiles
from kimmdy.utils import run_shell_cmd


def dummy_step(files):
    run_shell_cmd("pwd>./pwd.pwd", files.outputdir)
    return files


def write_mdp():
    pass

def md(files: TaskFiles, addgrompp: str, addmdrun: str):
    """General MD simulation
    """
    # Theory: Should cpt be a default input?
    # TODO: enable md on cluster
    # TODO: Could be built to construct commands in a more general way

    outputdir = files.outputdir

    top = files.input["top"]
    gro = files.input["gro"]
    mdp = files.input["mdp"]
    idx = files.input["idx"]

    tpr = outputdir / "equil.tpr"           ## correct name?!?!?!
    outgro = outputdir / "equil.gro"
    outtrr = outputdir / "equil.trr"

    files.output = {"tpr": tpr, "gro": outgro, "trr": outtrr}

    gromppcmd = f"gmx grompp -p {top} -c {gro} -f {mdp} -n {idx} -o {tpr} -maxwarn 5" + addgrompp
    mdruncmd = f"gmx mdrun -v -s {tpr} -c {outgro} -o {outtrr}" + addmdrun

    run_shell_cmd(gromppcmd, outputdir)
    run_shell_cmd(mdruncmd, outputdir)

    return files

# TODO:  add continue function if simulation stops?!
