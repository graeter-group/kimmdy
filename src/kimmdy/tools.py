"""
Standalone tools that are complementary to KIMMDY.
"""

from pathlib import Path
import shutil
import argparse

from kimmdy.topology.topology import Topology
from kimmdy.parsing import read_top, write_top
from kimmdy.utils import run_gmx, run_shell_cmd
from kimmdy.plugins import parameterization_plugins
from kimmdy.plugins import discover_plugins


def build_examples(restore: str):
    """Build example directories for KIMMDY from integration tests.

    Parameters
    ----------
    restore
        If True, overwrite input files in existing example directories. If "hard", also delete output files.
    """
    print(f'build_examples with restore option "{restore}"')
    overwrite = True if restore else False

    example_directories = [
        "alanine_hat_naive",
        "hexalanine_homolysis",
        "charged_peptide_homolysis_hat_naive",
        "triplehelix_pull",
        "hexalanine_single_reaction",
    ]
    root_dir = Path(__file__).parents[2]
    testfiles_path = root_dir / "tests" / "test_files" / "test_integration"
    examples_path = root_dir / "example"
    assets_path = Path(".") / ".." / ".." / "tests" / "test_files" / "assets"

    for directory in example_directories:
        try:
            print("Building example", directory, "...", end="")
            src = testfiles_path / directory
            dest = examples_path / directory
            if restore == "hard":
                shutil.rmtree(dest, ignore_errors=True)
            print(src, dest)
            shutil.copytree(src, dest, dirs_exist_ok=overwrite)
            (dest / "amber99sb-star-ildnp.ff").unlink(missing_ok=True)
            (dest / "amber99sb-star-ildnp.ff").symlink_to(
                assets_path / "amber99sb-star-ildnp.ff",
                target_is_directory=True,
            )
            print("done building examples")
        except FileExistsError as e:
            raise FileExistsError(
                f"Could not build example directory {directory} because it already exists. Try the --restore option to still build it."
            ) from e


def get_build_example_cmdline_args() -> argparse.Namespace:
    """Parse command line arguments.

    Returns
    -------
    Namespace
        parsed command line arguments
    """
    parser = argparse.ArgumentParser(description="Build examples for KIMMDY.")
    parser.add_argument(
        "-r",
        "--restore",
        const=True,
        nargs="?",
        type=str,
        help="Overwrite input files in existing example directories, use keyword 'hard' to also delete output files.",
    )
    return parser.parse_args()


def entry_point_build_examples():
    """Build examples from the command line."""
    args = get_build_example_cmdline_args()
    build_examples(args.restore)


def remove_hydrogen(
    gro: str,
    top: str,
    nr: str,
    parameterize: bool,
    equilibrate: bool,
    gmx_mdrun_flags: str = "",
):
    """Remove one hydrogen from a gro and top file to create a radical.

    Parameters
    ----------
    gro
        Path to GROMACS gro file
    top
        Path to GROMACS top file
    nr
        Atom number as indicated in the GROMACS gro and top file
    parameterize
        Parameterize topology with grappa after removing hydrogen.
    equilibrate
        Do a minimization and equilibration with GROMACS. Uses mdp files from kimmdy assets.
    """
    gro_path = Path(gro)
    top_path = Path(top)

    topology = Topology(read_top(top_path))

    ## check for input validity
    assert (atom_type := topology.atoms[nr].type).startswith(
        "H"
    ), f"Wrong atom type {atom_type} with nr {nr} for remove hydrogen, should start with 'H'"

    ## deal with top file, order is important here
    # [heavy_nr] = topology.atoms[nr].bound_to_nrs

    # parameterize with grappa
    if parameterize:
        discover_plugins()

        if "grappa" in parameterization_plugins.keys():
            topology.parametrizer = parameterization_plugins["grappa"]()

        else:
            raise KeyError(
                "No grappa in parameterization plugins. Can't continue to parameterize molecule"
            )

    update_map = topology.del_atom(nr, parameterize=parameterize)

    ## write top file
    top_stem = top_path.stem
    top_outpath = top_path.with_stem(top_stem + f"_d{nr}")
    write_top(topology.to_dict(), top_outpath)

    ## deal with gro file
    with open(gro_path, "r") as f:
        gro_raw = f.readlines()

    if len(gro_raw[3].split()) > 6:
        print("gro file likely contains velocities. These will be discarded.")

    gro_stem = gro_path.stem
    gro_outpath = gro_path.with_stem(gro_stem + f"_d{nr}")
    with open(gro_outpath, "w") as f:
        f.write(gro_raw[0])
        f.write(f"   {str(int(gro_raw[1])-1)}\n")
        for line in gro_raw[2:-1]:
            linesplit = line.split()
            if val := update_map.get(linesplit[2]):
                linesplit[2] = val
            else:
                continue
            f.write("{:>8s}   {:>4s} {:>4s} {:>7s} {:>7s} {:>7s}\n".format(*linesplit))
            # format not exactly as defined by GROMACS but should be good enough
        f.write(gro_raw[-1])

    ## minimize and equilibrate system using GROMACS
    if equilibrate:
        cwd = top_outpath.parent
        run_shell_cmd(
            f"cp {str(Path(__file__).parents[2] / 'tests'/'test_files'/'assets'/'md'/'*')} .",
            cwd=cwd,
        )
        run_gmx(
            f"gmx editconf -f {gro_outpath.name} -o {gro_outpath.stem}_box.gro -c -d 1.0 -bt dodecahedron",
            cwd=cwd,
        )
        run_gmx(
            f"gmx solvate -cp {gro_outpath.stem}_box.gro -p {top_outpath.name} -o {gro_outpath.stem}_solv.gro",
            cwd=cwd,
        )
        run_gmx(
            f"gmx grompp -f ions.mdp -c {gro_outpath.stem}_solv.gro -p {top_outpath.name} -o {gro_outpath.stem}_genion.tpr",
            cwd=cwd,
        )
        run_gmx(
            f"echo 'SOL' | gmx genion -s {gro_outpath.stem}_genion.tpr -p {top_outpath.name} -o {gro_outpath.stem}_ion.gro -conc 0.15 -neutral",
            cwd=cwd,
        )
        run_gmx(
            f"gmx grompp -f minim.mdp -c {gro_outpath.stem}_ion.gro -p {top_outpath.name} -o {gro_outpath.stem}_min.tpr -maxwarn 2",
            cwd=cwd,
        )
        run_gmx(f"gmx mdrun -deffnm {gro_outpath.stem}_min" + gmx_mdrun_flags, cwd=cwd)
        run_gmx(
            f"gmx grompp -f nvt.mdp -c {gro_outpath.stem}_min.gro -p {top_outpath.name} -o {gro_outpath.stem}_nvt.tpr -maxwarn 2",
            cwd=cwd,
        )
        run_gmx(f"gmx mdrun -deffnm {gro_outpath.stem}_nvt" + gmx_mdrun_flags, cwd=cwd)
        run_gmx(
            f"gmx grompp -f npt.mdp -c {gro_outpath.stem}_nvt.gro -p {top_outpath.name} -o {gro_outpath.stem}_npt.tpr -maxwarn 2",
            cwd=cwd,
        )
        run_gmx(
            f"gmx mdrun -v -deffnm {gro_outpath.stem}_npt" + gmx_mdrun_flags, cwd=cwd
        )


def get_remove_hydrogen_cmdline_args() -> argparse.Namespace:
    """
    parse cmdline args for remove_hydrogen
    """
    parser = argparse.ArgumentParser(
        description="Welcome to the KIMMDY remove hydrogen module",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("gro", help="GROMACS gro file")
    parser.add_argument("top", help="GROMACS top file")
    parser.add_argument(
        "nr", help="Atom number as indicated in the GROMACS gro and top file"
    )
    parser.add_argument(
        "-p",
        "--parameterize",
        action="store_true",
        help="Parameterize topology with grappa after removing hydrogen.",
        default=False,
    )
    parser.add_argument(
        "-e",
        "--equilibrate",
        action="store_true",
        help="Do a minimization and equilibration with GROMACS. Uses mdp files from kimmdy assets.",
        default=False,
    )
    return parser.parse_args()


def entry_point_remove_hydrogen():
    """Remove hydrogen by atom nr in a gro and topology file"""
    args = get_remove_hydrogen_cmdline_args()

    remove_hydrogen(args.gro, args.top, args.nr, args.parameterize, args.equilibrate)


## dot graphs
def topology_to_edgelist(top: Topology):
    """Convert a topology to a list of edges for a dot graph."""
    ls = []
    for b in top.bonds:
        ai = b[0]
        aj = b[1]
        x = top.atoms[ai]
        y = top.atoms[aj]
        left = f'"{ai} {x.type}" '
        right = f'"{aj} {y.type}" '
        ls.append(f"{left} -- {right};")
        if x.is_radical:
            ls.append(f'{left} [color="red"]')
        if y.is_radical:
            ls.append(f'{right} [color="red"]')
    return ls


def edgelist_to_dot_graph(ls: list[str], overlap: str = "true"):
    """Convert a list of edges to a dot graph."""
    header = f"""
        graph G {{
          layout=neato
          overlap={overlap}
          node [shape="circle"]
    """
    tail = """
        }
    """
    body = "\n".join(ls)
    return header + body + tail


def top_to_graph(top: Topology, overlap: str = "true"):
    """Convert a topology to a dot graph.

    Can be used in notebooks to write a dot file and render with quarto.
    """
    return edgelist_to_dot_graph(topology_to_edgelist(top), overlap)


def write_top_as_dot(top: Topology, path: str, overlap: str = "true"):
    """Write a topology as a dot graph to a file.

    Can be used in notebooks to write a dot file and render with quarto.
    """
    with open(path, "w") as f:
        f.write(top_to_graph(top, overlap))
