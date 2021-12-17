from __future__ import annotations
import logging
from kimmdy.runmanager import TaskFiles
from kimmdy.utils import run_shell_cmd


def dummy_step(files):
    run_shell_cmd("pwd>./pwd.pwd", files.outputdir)


def write_mdp():
    pass


def minimzation(files: TaskFiles):
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
    # TODO not sure if it is necessary to return the files here.
    # will they already have changed in the calling scope due to mutability?
    return files


def equilibration(files):
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


def equilibrium(out_dir, top, gro, mdp, idx):
    """equilibrium before pulling md"""
    tpr = out_dir / "equil.tpr"
    outgro = out_dir / "equil.gro"
    outtrr = out_dir / "equil.trr"
    run_shell_cmd(
        f"gmx grompp -p {top} -c {gro} -f {mdp} -n {idx} -o {tpr} -maxwarn 5",
        out_dir,
    )
    run_shell_cmd(f"gmx mdrun -v -s {tpr} -c {outgro} -o {outtrr}", out_dir)


def production(out_dir, top, gro, mdp, idx, cpt, plumed_dat):
    """normal pulling md"""
    tpr = out_dir / "prod.tpr"
    outgro = out_dir / "prod.gro"
    outtrr = out_dir / "prod.trr"
    run_shell_cmd(
        f" gmx grompp -p {top} -c {gro} -f {mdp} -n {idx} -t {cpt} -o {tpr} -maxwarn 5",
        out_dir,
    )
    run_shell_cmd(
        f"gmx mdrun -v -s {tpr} -c {outgro} -plumed {plumed_dat} -o {outtrr}",
        out_dir,
    )


def relaxation(out_dir, top, gro, mdp, idx, cpt):
    """equil after break -->> called from changer"""
    tpr = out_dir / "relax.tpr"
    outgro = out_dir / "relax.gro"
    outtrr = out_dir / "relax.trr"
    run_shell_cmd(
        f"gmx grompp -p {top} -c {gro} -f {mdp} -n {idx} -t {cpt} -o {tpr} -maxwarn 5",
        out_dir,
    )
    run_shell_cmd(f"gmx mdrun -v -s {tpr} -c {outgro} -o {outtrr}", out_dir)
