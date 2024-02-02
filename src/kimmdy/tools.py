"""
Standalone tools that are complementary to KIMMDY.
"""

from pathlib import Path
import shutil
import argparse
from typing import Optional

from kimmdy.topology.topology import Topology
from kimmdy.parsing import read_top, write_top
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


def modify_top(
    top_path: str,
    out_path: str,
    parameterize: bool,
    removeH: Optional(list[int]),
    grofile: Optional(str),
    gmx_mdrun_flags: str = "",
):
    """Modify topology in various ways.

    Parameters
    ----------
    top_path
        Path to GROMACS top file
    out_path
        Output topology file name
    parameterize
        Parameterize topology with grappa after removing hydrogen
    removeH
        Remove one or more hydrogens by atom nrs in the top file
    grofile
        GROMACS gro file; keep compatible with top file
    """
    print(
        "Changing topology\n"
        f"top: {top_path}\n"
        f"output: {out_path}\n"
        f"parameterize: {parameterize}\n"
        f"remove hydrogen: {removeH}\n"
        f"optional gro: {grofile}\n"
    )

    top = Topology(read_top(Path(top_path)))

    # remove hydrogen
    if removeH:
        for nr in removeH:
            # check for input validity
            if (atom_type := top.atoms[str(nr)].type).startswith("H"):
                update_map = top.del_atom(str(nr), parameterize=parameterize)
            else:
                print(
                    f"Wrong atom type {atom_type} with nr {nr} for remove hydrogen, should start with 'H'."
                )
                continue

    # parameterize with grappa
    if parameterize:
        # load grappa
        discover_plugins()
        if "grappa" in parameterization_plugins.keys():
            top.parametrizer = parameterization_plugins["grappa"]()
        else:
            raise KeyError(
                "No grappa in parameterization plugins. Can't continue to parameterize molecule"
            )
        # require parameterization when writing topology to dict
        top.needs_parameterization = True

    # write top file
    write_top(top.to_dict(), Path(out_path))

    # deal with gro file
    if grofile:
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
                f.write(
                    "{:>8s}   {:>4s} {:>4s} {:>7s} {:>7s} {:>7s}\n".format(*linesplit)
                )
                # format not exactly as defined by GROMACS but should be good enough
            f.write(gro_raw[-1])


def get_modify_top_cmdline_args() -> argparse.Namespace:
    """
    parse cmdline args for modify_top
    """
    parser = argparse.ArgumentParser(
        description="Welcome to the KIMMDY modify top module",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("top", help="GROMACS top file")
    parser.add_argument("out", help="Output top file name")

    parser.add_argument(
        "-p",
        "--parameterize",
        action="store_true",
        help="Parameterize topology with grappa.",
        default=False,
    )
    parser.add_argument(
        "-r",
        "--removeH",
        help="Remove one or more hydrogens by atom nrs in the top file.",
        nargs="+",
        type=int,
    )
    parser.add_argument(
        "-c",
        "--grofile",
        help="If necessary, also apply actions on gro file to create a compatible gro/top file pair.",
        type=str,
    )
    return parser.parse_args()


def entry_point_modify_top():
    """Modify topology file in various ways"""
    args = get_modify_top_cmdline_args()

    modify_top(args.top, args.out, args.parameterize, args.removeH, args.grofile)


# dot graphs
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
