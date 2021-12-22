from __future__ import annotations
import logging
from kimmdy.tasks import TaskFiles
from kimmdy.utils import run_shell_cmd


def dummy_step(files):
    run_shell_cmd("pwd>./pwd.pwd", files.outputdir)

    return files


def write_mdp():
    pass


def minimzation(files: TaskFiles) -> TaskFiles:
    outputdir = files.outputdir
    if not outputdir:
        m = f"No output directory given for the current task"
        logging.error(m)
        raise FileNotFoundError(m)

    top = files.input["top"]
    gro = files.input["gro"]
    mdp = files.input["mdp"]
    outgro = outputdir / "minimized.gro"
    outgro = outputdir / "minimized.gro"
    tpr = outputdir / "minimization.tpr"

    files.output = {"tpr": tpr, "outgro": outgro}

    run_shell_cmd(f"gmx grompp -p {top} -c {gro} -r {gro} -f {mdp} -o {tpr}", outputdir)
    run_shell_cmd(f"gmx mdrun -s {tpr} -c {outgro}", outputdir)

    return files


def equilibration(files: TaskFiles) -> TaskFiles:
    outputdir = files.outputdir
    if not outputdir:
        m = f"No output directory given for the current task"
        logging.error(m)
        raise FileNotFoundError(m)

    gro = files.input["gro"]
    top = files.input["top"]
    mdp = files.input["mdp"]
    tpr = outputdir / "equil.tpr"
    outgro = outputdir / "equil.gro"

    files.output = {"tpr": tpr, "gro": outgro}

    run_shell_cmd(f"gmx grompp -p {top} -c {gro} -r {gro} -f {mdp} -o {tpr}", outputdir)
    run_shell_cmd(f"gmx mdrun -v -s {tpr} -c {outgro}", outputdir)

    return files


def equilibrium(files: TaskFiles):
    """equilibrium before pulling md"""
    outputdir = files.outputdir
    if not outputdir:
        m = f"No output directory given for the current task"
        logging.error(m)
        raise FileNotFoundError(m)

    top = files.input["top"]
    gro = files.input["gro"]
    mdp = files.input["mdp"]
    idx = files.input["idx"]
    tpr = outputdir / "equil.tpr"
    outgro = outputdir / "equil.gro"
    outtrr = outputdir / "equil.trr"

    files.output = {"tpr": tpr, "gro": outgro, "trr": outtrr}

    run_shell_cmd(
        f"gmx grompp -p {top} -c {gro} -f {mdp} -n {idx} -o {tpr} -maxwarn 5",
        outputdir,
    )
    run_shell_cmd(f"gmx mdrun -v -s {tpr} -c {outgro} -o {outtrr}", outputdir)

    return TaskFiles


def production(files: TaskFiles) -> TaskFiles:
    """normal pulling md"""
    outputdir = files.outputdir
    if not outputdir:
        m = f"No output directory given for the current task"
        logging.error(m)
        raise FileNotFoundError(m)

    top = files.input["top"]
    gro = files.input["gro"]
    mdp = files.input["mdp"]
    idx = files.input["idx"]
    cpt = files.input["cpt"]
    plumed_dat = files.input["plumed.dat"]
    tpr = outputdir / "prod.tpr"
    outgro = outputdir / "prod.gro"
    outtrr = outputdir / "prod.trr"

    files.output = {"tpr": tpr, "gro": outgro, "trr": outtrr}

    run_shell_cmd(
        f" gmx grompp -p {top} -c {gro} -f {mdp} -n {idx} -t {cpt} -o {tpr} -maxwarn 5",
        outputdir,
    )
    run_shell_cmd(
        f"gmx mdrun -v -s {tpr} -c {outgro} -plumed {plumed_dat} -o {outtrr}",
        outputdir,
    )

    return files


def relaxation(files: TaskFiles) -> TaskFiles:
    """equil after break -->> called from changer"""
    outputdir = files.outputdir
    if not outputdir:
        m = f"No output directory given for the current task"
        logging.error(m)
        raise FileNotFoundError(m)

    top = files.input["top"]
    gro = files.input["gro"]
    mdp = files.input["mdp"]
    idx = files.input["idx"]
    cpt = files.input["cpt"]
    tpr = outputdir / "relax.tpr"
    outgro = outputdir / "relax.gro"
    outtrr = outputdir / "relax.trr"

    files.output = {"tpr": tpr, "gro": outgro, "trr": outtrr}

    run_shell_cmd(
        f"gmx grompp -p {top} -c {gro} -f {mdp} -n {idx} -t {cpt} -o {tpr} -maxwarn 5",
        outputdir,
    )
    run_shell_cmd(f"gmx mdrun -v -s {tpr} -c {outgro} -o {outtrr}", outputdir)

    return files
