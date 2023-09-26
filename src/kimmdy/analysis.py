"""Analysis tools for KIMMDY runs.
For command line usage, run `kimmdy-analysis -h`.
"""
from typing import Union
from pathlib import Path
import MDAnalysis as mda
import subprocess as sp
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn.objects as so
import argparse
from seaborn import axes_style
import pandas as pd

from kimmdy.utils import run_shell_cmd
from kimmdy.parsing import read_json, write_json
from kimmdy.recipe import RecipeCollection


def get_analysis_dir(dir: Path) -> Path:
    """Get analysis directory for a KIMMDY run.

    Creates the directory if it does not exist.

    Parameters
    ----------
    dir
        Directory of KIMMDY run

    Returns
    -------
        Path to analysis directory
    """
    out = dir / "analysis"
    out.mkdir(exist_ok=True)
    return out


def get_step_directories(dir: Path, steps: Union[list[str], str] = "all") -> list[Path]:
    """
    create list of subdirectories that match the steps.
    If steps is "all", all subdirectories are returned.

    Parameters
    ----------
    dir
        Directory to search for subdirectories
    steps
        List of steps e.g. ["equilibrium", "production"]. Or a string "all" to return all subdirectories
    """
    directories = sorted(
        [p for p in dir.glob("*_*/") if p.is_dir()],
        key=lambda p: int(p.name.split("_")[0]),
    )
    if steps == "all":
        matching_directories = directories
    else:
        matching_directories = list(
            filter(lambda d: d.name.split("_")[1] in steps, directories)
        )

    if not matching_directories:
        raise ValueError(
            f"Could not find directories {steps} in {dir}. Thus, no trajectories can be concatenated"
        )

    return matching_directories


def concat_traj(dir: str, steps: Union[list[str], str], open_vmd: bool = False):
    """Find and concatenate trajectories (.xtc files) from a KIMMDY run into one trajectory.
    The concatenated trajectory is centered and pbc corrected.

    Parameters
    ----------
    dir
        Directory to search for subdirectories
    steps
        List of steps e.g. ["equilibrium", "production"]. Or a string "all" to return all subdirectories
    open_vmd
        Open concatenated trajectory in VMD
    """
    run_dir = Path(dir).expanduser().resolve()
    analysis_dir = get_analysis_dir(run_dir)

    directories = get_step_directories(run_dir, steps)

    out_xtc = analysis_dir / "concat.xtc"

    ## gather trajectories
    trajectories = []
    tprs = []
    gros = []
    for d in directories:
        trajectories.extend(d.glob("*.xtc"))
        tprs.extend(d.glob("*.tpr"))
        gros.extend(d.glob("*.gro"))

    assert (
        len(trajectories) > 0
    ), f"No trrs found to concatenate in {run_dir} with subdirectory names {steps}"

    trajectories = [str(t) for t in trajectories]

    ## write concatenated trajectory
    tmp_xtc = str(out_xtc.with_name("tmp.xtc"))
    run_shell_cmd(
        f"gmx trjcat -f {' '.join(trajectories)} -o {tmp_xtc} -cat",
        cwd=run_dir,
    )
    run_shell_cmd(
        f"echo '1 0' | gmx trjconv -f {tmp_xtc} -s {tprs[0]} -o {str(out_xtc)} -center -pbc mol",
        cwd=run_dir,
    )
    run_shell_cmd(f"rm {tmp_xtc}", cwd=run_dir)
    if open_vmd:
        gro = str(gros[0])
        run_shell_cmd(f"cp {gro} {str(analysis_dir)}", cwd=run_dir)
        run_shell_cmd(f"vmd {gro} {str(out_xtc)}", cwd=run_dir)


def plot_energy(
    dir: str, steps: Union[list[str], str], terms: list[str], open_plot: bool = False
):
    """Plot GROMACS energy for a KIMMDY run.

    Parameters
    ----------
    dir
        Directory to search for subdirectories
    steps
        List of steps e.g. ["equilibrium", "production"]. Or a string "all" to return all subdirectories.
        Default is "all".
    terms
        Terms from gmx energy that will be plotted. Uses 'Potential' by default.
    open_plot :
        Open plot in default system viewer.
    """
    run_dir = Path(dir).expanduser().resolve()
    xvg_entries = ["time"] + terms
    terms_str = "\n".join(terms)

    subdirs_matched = get_step_directories(run_dir, steps)

    analysis_dir = get_analysis_dir(run_dir)
    xvgs_dir = analysis_dir / "energy_xvgs"
    xvgs_dir.mkdir(exist_ok=True)

    ## gather energy files
    edrs = []
    for d in subdirs_matched:
        edrs.extend(d.glob("*.edr"))
    assert (
        len(edrs) > 0
    ), f"No GROMACS energy files in {run_dir} with subdirectory names {steps}"

    energy = {}
    for k in xvg_entries:
        energy[k] = []

    energy["step"] = []
    energy["step_ix"] = []

    time_offset = 0
    for i, edr in enumerate(edrs):
        ## write energy .xvg file
        xvg = str(xvgs_dir / edr.parents[0].with_suffix(".xvg").name)
        step_name = edr.parents[0].name.split("_")[1]

        run_shell_cmd(
            f"echo '{terms_str} \n\n' | gmx energy -f {str(edr)} -o {xvg}",
            cwd=run_dir,
        )

        ## read energy .xvg file
        with open(xvg, "r") as f:
            for line in f:
                if line[0] not in ["@", "#"]:
                    energy["step"].append(step_name)
                    energy["step_ix"].append(i)
                    for k, v in zip(xvg_entries, line.split()):
                        if k == "time":
                            energy[k].append(float(v) + time_offset)
                        else:
                            energy[k].append(float(v))

        time_offset = energy["time"][-1]

    df = pd.DataFrame(energy).melt(
        id_vars=["time", "step", "step_ix"], value_vars=terms
    )
    step_names = df[df["variable"] == terms[0]]
    step_names["variable"] = "Step"
    step_names["value"] = step_names["step_ix"]
    df = pd.concat([df, step_names], ignore_index=True)

    p = (
        so.Plot(df, x="time", y="value")
        .add(so.Line())
        .facet(row="variable")
        .share(x=True, y=False)
        .theme({**axes_style("white")})
        .label(x="Time [ps]", y="Energy [kJ/mol]")
    )
    p.plot(pyplot=True)

    step_names = step_names.groupby(["step", "step_ix"]).first().reset_index()
    for t, v, s in zip(step_names["time"], step_names["value"], step_names["step"]):
        plt.axvline(x=t, color="black", linestyle="--")
        plt.text(x=t, y=v + 0.5, s=s, fontsize=6)

    ax = plt.gca()
    steps_y_axis = [c for c in ax.get_children() if isinstance(c, mpl.axis.YAxis)][0]
    steps_y_axis.set_visible(False)
    output_path = str(run_dir / "analysis" / "energy.png")
    plt.savefig(output_path, dpi=300)

    if open_plot:
        sp.call(("xdg-open", output_path))


def radical_population(
    dir: str,
    steps: Union[list[str], str] = "all",
    select_atoms: str = "protein",
    open_plot: bool = False,
    open_vmd: bool = False,
):
    """Plot population of radicals for a KIMMDY run.

    Parameters
    ----------
    dir
        KIMMDY run directory to be analysed.
    steps
        List of steps e.g. ["equilibrium", "production"]. Or a string "all" to return all subdirectories.
        Default is "all".
    select_atoms
        Atoms chosen for radical population analysis, default is protein (uses MDAnalysis selection syntax)
    open_plot
        Open plot in default system viewer.
    open_vmd
        Open a pdb in VMD with the radical occupation as B-factors.
    """
    run_dir = Path(dir).expanduser().resolve()
    analysis_dir = get_analysis_dir(run_dir)

    subdirs_matched = get_step_directories(run_dir, steps)
    radical_jsons = run_dir.glob("**/radicals.json")

    radical_info = {"time": [], "radicals": []}
    for radical_json in radical_jsons:
        data = read_json(radical_json)
        for k in radical_info.keys():
            radical_info[k].append(data[k])

    write_json(radical_info, analysis_dir / "radical_population.json")

    ## get atoms from gro file
    gro = None
    for subdir in subdirs_matched:
        gro = list(subdir.glob("*.gro"))
        if gro:
            break
    if not gro:
        raise ValueError(f"No .gro file found in {run_dir}")

    u = mda.Universe(str(gro[0]), format="gro")
    atoms = u.select_atoms(select_atoms)
    atoms_identifier = [
        "-".join(x)
        for x in list(
            zip(
                [str(resid) for resid in atoms.resids],
                [str(name) for name in atoms.names],
            )
        )
    ]
    atom_ids = atoms.ids

    ## plot fingerprint
    counts = {i: 0.0 for i in atom_ids}
    n_states = len(radical_info["time"])
    for state in range(n_states):
        for idx in radical_info["radicals"][state]:
            if int(idx) in counts.keys():
                counts[int(idx)] += 1 / n_states

    plt.bar(x=atoms_identifier, height=counts.values())
    plt.xlabel("Atom identifier")
    plt.ylabel("Fractional Radical Occupancy")
    plt.ylim(0, 1)
    plt.xticks(atoms_identifier, rotation=90, ha="right")
    plt.tight_layout()
    output_path = str(analysis_dir / "radical_population_fingerprint.png")
    plt.savefig(output_path, dpi=300)

    u.add_TopologyAttr("tempfactors")
    atoms = u.select_atoms(select_atoms)
    atoms.tempfactors = list(counts.values())
    protein = u.select_atoms("protein")
    pdb_output = f"{analysis_dir}/radical_population.pdb"
    protein.write(pdb_output)

    if open_plot:
        sp.call(("xdg-open", output_path))

    if open_vmd:
        run_shell_cmd(f"vmd {pdb_output}", cwd=analysis_dir)


def plot_rates(dir: str):
    """Plot rates of all possible reactions for each 'decide_recipe' step.

    Parameters
    ----------
    dir
        Directory of KIMMDY run
    """
    run_dir = Path(dir).expanduser().resolve()
    analysis_dir = get_analysis_dir(run_dir)

    for recipes in run_dir.glob("*decide_recipe/recipes.csv"):
        rc, picked_rp = RecipeCollection.from_csv(recipes)
        i = recipes.parent.name.split("_")[0]
        rc.plot(analysis_dir / f"{i}_reaction_rates.svg", highlight_r=picked_rp)


def plot_runtime(dir: str, log_str: str, datefmt: str):
    """Plot runtime of all tasks.

    Parameters
    ----------
    dir
        Directory of KIMMDY run
    logfile
        KIMMDY logfile
    datefmt
        Date format in the KIMMDY logfile
    """

    def time_from_logline(line: str, sep: int):
        return datetime.strptime(" ".join(line.split()[:sep]), datefmt)

    import re
    from datetime import datetime

    run_dir = Path(dir).expanduser().resolve()
    analysis_dir = get_analysis_dir(run_dir)

    log_path = run_dir / log_str
    datefmt_substrings = len(datefmt.split())

    start_str = " INFO: Start "
    runstart_str = " INFO: Start run\n"
    end_str = " INFO: Done with "

    with open(log_path, "r") as f:
        log = f.readlines()

    starttime = time_from_logline(log[0], datefmt_substrings)
    endtime = time_from_logline(log[-1], datefmt_substrings)
    walltime = endtime - starttime
    print(starttime, endtime)
    print(walltime)

    # regex to find start of task

    # regex to find end of task
    # match both
    #
    task_dict = {}
    tasks = []
    runtimes = []
    for line in log:
        if start_str in line:
            if runstart_str in line:
                continue
            print(line)
            linesplit = line.split()
            # task name comes after datefmt in logging format
            task_dict[linesplit[datefmt_substrings]] = time_from_logline(
                line, datefmt_substrings
            )
        if end_str in line:
            print(line)
            linesplit = line.split()
            runtime = (
                time_from_logline(line, datefmt_substrings)
                - task_dict[linesplit[datefmt_substrings]]
            )

            tasks.append(linesplit[datefmt_substrings])
            runtimes.append(runtime.total_seconds())

    print(tasks, runtimes)
    plt.bar(tasks, runtimes)
    plt.savefig(analysis_dir / "test_runtime.png")


def get_analysis_cmdline_args() -> argparse.Namespace:
    """Parse command line arguments.

    Returns
    -------
        Parsed command line arguments
    """
    parser = argparse.ArgumentParser(
        description="Welcome to the KIMMDY analysis module. Use this module to analyse existing KIMMDY runs.",
    )
    subparsers = parser.add_subparsers(metavar="module", dest="module")

    parser_trjcat = subparsers.add_parser(
        name="trjcat", help="Concatenate trajectories of a KIMMDY run"
    )
    parser_trjcat.add_argument(
        "dir", type=str, help="KIMMDY run directory to be analysed."
    )
    parser_trjcat.add_argument(
        "--steps",
        "-s",
        nargs="*",
        default="all",
        help=(
            "Apply analysis method to subdirectories with these names. Uses all subdirectories by default."
        ),
    )
    parser_trjcat.add_argument(
        "--open-vmd",
        action="store_true",
        help="Open VMD with the concatenated trajectory.",
    )

    parser_plot_energy = subparsers.add_parser(
        name="plot_energy", help="Plot GROMACS energy for a KIMMDY run"
    )
    parser_plot_energy.add_argument(
        "dir", type=str, help="KIMMDY run directory to be analysed."
    )
    parser_plot_energy.add_argument(
        "--steps",
        "-s",
        nargs="*",
        default="all",
        help=(
            "Apply analysis method to subdirectories with these names. Uses all subdirectories by default."
        ),
    )
    parser_plot_energy.add_argument(
        "--terms",
        "-t",
        nargs="*",
        default=["Potential"],
        help=(
            "Terms from gmx energy that will be plotted. Uses 'Potential' by default."
        ),
    )
    parser_plot_energy.add_argument(
        "--open-plot", action="store_true", help="Open plot in default system viewer."
    )

    parser_radical_population = subparsers.add_parser(
        name="radical_population",
        help="Plot population of radicals for one or multiple KIMMDY run(s)",
    )
    parser_radical_population.add_argument(
        "dir", type=str, help="KIMMDY run directory to be analysed."
    )
    parser_radical_population.add_argument(
        "--select_atoms",
        "-a",
        type=str,
        help="Atoms chosen for radical population analysis, default is protein (uses MDAnalysis selection syntax)",
        default="protein",
    )
    parser_radical_population.add_argument(
        "--steps",
        "-s",
        nargs="*",
        default="all",
        help=(
            "Apply analysis method to subdirectories with these names. Uses all subdirectories by default."
        ),
    )
    parser_radical_population.add_argument(
        "--open-plot", action="store_true", help="Open plot in default system viewer."
    )
    parser_radical_population.add_argument(
        "--open-vmd",
        action="store_true",
        help="Open VMD with the concatenated trajectory."
        "To view the radical occupancy per atom, add a representation with the beta factor as color.",
    )

    parser_plot_rates = subparsers.add_parser(
        name="plot_rates",
        help="Plot rates of all possible reactions after a MD run. Rates must have been saved!",
    )
    parser_plot_rates.add_argument(
        "dir", type=str, help="KIMMDY run directory to be analysed."
    )

    parser_runtime = subparsers.add_parser(
        name="runtime",
        help="Plot runtime of the tasks of a kimmdy run.",
    )
    parser_runtime.add_argument(
        "dir", type=str, help="KIMMDY run directory to be analysed."
    )
    parser_runtime.add_argument(
        "--logfile", type=str, default="kimmdy.log", help="KIMMDY logfile."
    )
    parser_runtime.add_argument(
        "--datefmt",
        type=str,
        default="%d-%m-%y %H:%M:%S",
        help="Date format in the KIMMDY logfile.",
    )

    return parser.parse_args()


def entry_point_analysis():
    """Analyse existing KIMMDY runs."""
    args = get_analysis_cmdline_args()

    if args.module == "trjcat":
        concat_traj(args.dir, args.steps, args.open_vmd)
    elif args.module == "plot_energy":
        plot_energy(args.dir, args.steps, args.terms, args.open_plot)
    elif args.module == "radical_population":
        radical_population(
            args.dir, args.steps, args.select_atoms, args.open_plot, args.open_vmd
        )
    elif args.module == "plot_rates":
        plot_rates(args.dir)
    elif args.module == "runtime":
        plot_runtime(args.dir, args.logfile, args.datefmt)
    else:
        print(
            "No analysis module specified. Use -h for help and a list of available modules."
        )
