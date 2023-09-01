from typing import Union
from pathlib import Path
from math import isclose
import MDAnalysis as mda
import subprocess as sp
import matplotlib.pyplot as plt
import seaborn.objects as so
import argparse
from seaborn import axes_style
import pandas as pd

from kimmdy.utils import run_shell_cmd
from kimmdy.parsing import read_json, write_json
from kimmdy.recipe import RecipeCollection


def get_step_directories(dir: Path, steps: Union[list[str], str] = "all") -> list[Path]:
    """
    create list of subdirectories that match the steps.
    If steps is "all", all subdirectories are returned.

    Parameters
    ----------
    dir :
        Directory to search for subdirectories
    steps :
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
    dir :
        Directory to search for subdirectories
    steps :
        List of steps e.g. ["equilibrium", "production"]. Or a string "all" to return all subdirectories
    open_vmd :
        Open concatenated trajectory in VMD
    """
    run_dir = Path(dir).expanduser().resolve()

    directories = get_step_directories(run_dir, steps)

    ## create output dir
    analysis_dir = run_dir / "analysis"
    analysis_dir.mkdir(exist_ok=True)
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


def plot_energy(dir: str, steps: Union[list[str], str], terms: list[str], open_plot: bool = False):
    run_dir = Path(dir).expanduser().resolve()
    xvg_entries = ["time"] + terms

    subdirs_matched = get_step_directories(run_dir, steps)

    ## create output dir
    (run_dir / "analysis").mkdir(exist_ok=True)
    xvgs_dir = run_dir / "analysis" / "energy_xvgs"
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

    ## write energy .xvg files
    for edr in edrs:
        print(edr.parents[0].name + ".xvg")
        xvg = str(xvgs_dir / edr.parents[0].with_suffix(".xvg").name)
        terms_str = "\n".join(terms)
        run_shell_cmd(
            f"echo '{terms_str} \n\n' | gmx energy -f {str(edr)} -o {xvg}",
            cwd=run_dir,
        )

        ## read energy .xvg files
        with open(xvg, "r") as f:
            for line in f:
                if line[0] not in ["@", "#"]:
                    for k, v in zip(xvg_entries, line.split()):
                        energy[k].append(float(v))

    ## plot energy
    snapshot = range(len(energy['time']))
    steps_starts = [i for i in snapshot if isclose(energy["time"][i], 0)]
    steps_names = [str(edr.parents[0].name).split("_")[1] for edr in edrs]

    df = pd.DataFrame(energy).melt(id_vars=["time"], value_vars=terms)
    steps = pd.DataFrame({"time": steps_starts, "name": steps_names})

    p = (so.Plot(df, x = "time", y = "value")
        .facet(row="variable").share(x=True, y=False)
        .add(so.Line())
        # .add(so.Path(color="black"), data=steps)
        .theme({**axes_style("white")})
        .label(x = "Time [ps]", y = "Energy [kJ/mol]")
    )
    plt.axvline(x=5, color="black", linestyle="--")

    output_path = str(run_dir / "analysis" / "energy.png")
    p.save(output_path, dpi=300)

    # open png file with default system viewer
    if open_plot:
        sp.call(('xdg-open', output_path))



def radical_population(rundir: str, select_atoms: str):
    # TODO: weigh radical population by time

    ## set up directory to store radical information
    radical_info = {"time": [], "radicals": []}

    for curr_dir in rundir[::-1]:
        run_dir = Path(curr_dir).expanduser().resolve()

        ## find .gro file
        subdirs_sorted = sorted(
            list(filter(lambda d: d.is_dir(), run_dir.glob("*_*/"))),
            key=lambda p: int(p.name.split("_")[0]),
        )
        for subdir in subdirs_sorted:
            gro = list(subdir.glob("*.gro"))
            if gro:
                break
        assert gro

        ## gather radical info
        radical_jsons = run_dir.glob("**/radicals.json")
        # print(list(radical_jsons))

        ## parse radical info
        for radical_json in radical_jsons:
            data = read_json(radical_json)
            for k in radical_info.keys():
                radical_info[k].append(data[k])

    ## create output dir (only goes to first mentioned run_dir)
    (run_dir / "analysis").mkdir(exist_ok=True)
    out = run_dir / "analysis"
    out = Path(out).expanduser()

    ## write gathered radical info
    write_json(radical_info, out / "radical_population.json")

    ## get info from gro file
    u = mda.Universe(str(gro[0]), format="gro")
    print(u)
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
    atoms_id = atoms.ids
    print(atoms_identifier)
    print(atoms_id)

    ## plot fingerprint
    counts = {i: 0.0 for i in atoms_id}
    n_states = len(radical_info["time"])
    for state in range(n_states):
        for idx in radical_info["radicals"][state]:
            if int(idx) in counts.keys():
                counts[int(idx)] += 1 / n_states
    print(counts)

    plt.bar(x=atoms_identifier, height=counts.values())
    plt.xlabel("Atom identifier")
    plt.ylabel("Fractional Radical Occupancy")
    plt.ylim(0, 1)
    plt.xticks(atoms_identifier, rotation=90, ha="right")
    plt.tight_layout()
    plt.savefig(str(run_dir / "analysis" / "radical_population_fingerprint"), dpi=300)

    u.add_TopologyAttr("tempfactors")
    atoms = u.select_atoms(select_atoms)
    print(list(counts.values()))
    print(atoms.tempfactors)
    atoms.tempfactors = list(counts.values())
    protein = u.select_atoms("protein")
    protein.write(str(out))


def plot_rates(rundir: list):
    for curr_dir in rundir:
        run_dir = Path(curr_dir).expanduser().resolve()

        ## create output dir (only goes to first mentioned run_dir)
        (run_dir / "analysis").mkdir(exist_ok=True)
        out = run_dir / "analysis"

        for recipes in run_dir.glob("*decide_recipe/recipes.csv"):
            rc, picked_rp = RecipeCollection.from_csv(recipes)
            i = recipes.parent.name.split("_")[0]
            rc.plot(out / f"{i}_reaction_rates.svg", highlight_r=picked_rp)

def get_analysis_cmdline_args() -> argparse.Namespace:
    """Parse command line arguments.

    Returns
    -------
    :
        Parsed command line arguments
    """
    parser = argparse.ArgumentParser(
        description="Welcome to the KIMMDY analysis module. Use this module to analyse existing KIMMDY runs.",
    )
    subparsers = parser.add_subparsers(required=True, metavar="module", dest="module")

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
    parser_trjcat.add_argument("--open-vmd", action="store_true", help="Open VMD with the concatenated trajectory.")

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
    parser_plot_energy.add_argument("--open-plot", action="store_true", help="Open plot in default system viewer.")

    parser_radical_population = subparsers.add_parser(
        name="radical_population",
        help="Plot population of radicals for one or multiple KIMMDY run(s)",
    )
    parser_radical_population.add_argument(
        "dir", nargs="+", help="KIMMDY run directory to be analysed. Can be multiple."
    )
    parser_radical_population.add_argument(
        "--select_atoms",
        "-a",
        type=str,
        help="Atoms chosen for radical population analysis, default is protein (uses MDAnalysis selection syntax)",
        default="protein",
    )

    parser_plot_rates = subparsers.add_parser(
        name="plot_rates",
        help="Plot rates of all possible reactions after a MD run. Rates must have been saved!",
    )
    parser_plot_rates.add_argument(
        "dir", nargs="+", help="KIMMDY run directory to be analysed. Can be multiple."
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
        radical_population(args.dir, args.radical_population)
    elif args.module == "plot_rates":
        plot_rates(args.dir)

