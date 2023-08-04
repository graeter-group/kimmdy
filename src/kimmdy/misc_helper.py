"""
Miscelaneous utilitiies that didn't fit anywhere else for now.
"""
from pathlib import Path
import shutil
import argparse
import subprocess
from typing import Union
from kimmdy.topology.topology import Topology


def topology_to_edgelist(top: Topology):
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
    return edgelist_to_dot_graph(topology_to_edgelist(top), overlap)


def concat_traj(
    run_dir: Union[Path, str], out: Union[Path, str] = None, run_types=None
):
    """Find and concatenate trajectories from KIMMDY runs.

    Parameters
    ----------
    run_dir : Path|str
        Directory containing directories of multiple tasks.
    out : Path
        File Path into the output trr will be written.
        Default `concat.trr` in given run_dir
    run_types : list, optional
        List of tasks to get trrs from. If a task is not in this list
        it will be skipped. By default None
    """
    run_dir = Path(run_dir).expanduser().resolve()
    if out is None:
        out = run_dir / "concat.trr"
    out = Path(out).expanduser()
    assert out.suffix == ".trr", "Output file should be a trr file."
    dirs = sorted(
        list(filter(lambda d: d.is_dir(), run_dir.iterdir())),
        key=lambda p: int(p.name.split("_")[0]),
    )

    if run_types is not None:
        dirs = list(filter(lambda d: d.name.split("_")[1] in run_types, dirs))

    trrs = []
    for d in dirs:
        trrs.extend(d.glob("*.trr"))
    trrs = list(filter(lambda p: "rotref" not in p.stem, trrs))
    trrs = [str(t) for t in trrs]
    newline = "\n\t"
    assert (
        len(trrs) > 0
    ), f"No trrs found to concatenate in \n{newline.join([str(t) for t in dirs])}"

    command = f"gmx trjcat -f {' '.join(trrs)} -o {str(out)} -cat".split(" ")

    subprocess.run(command, cwd=run_dir)


def _build_examples(args: argparse.Namespace):
    overwrite = True if args.restore else False

    rename = {
        "hat_naive": "alanine",
        "homolysis": "hexalanine",
        "whole_run": "charged_peptide",
    }  # ,'pull':'tripelhelix'}
    basedir = Path(__file__).parents[2]
    testpath = basedir / "tests" / "test_files" / "test_integration"
    assetpath = basedir / "tests" / "test_files" / "assets"
    examplepath = basedir / "example"

    for k, v in rename.items():
        try:
            srcpath = testpath / k
            dstpath = examplepath / v
            if args.restore == "hard":
                shutil.rmtree(dstpath, ignore_errors=True)
            shutil.copytree(srcpath, dstpath, dirs_exist_ok=overwrite)
            (dstpath / "amber99sb-star-ildnp.ff").unlink(missing_ok=True)
            (dstpath / "amber99sb-star-ildnp.ff").symlink_to(
                (assetpath / "amber99sb-star-ildnp.ff"),
                target_is_directory=True,
            )
        except FileExistsError as e:
            raise FileExistsError(
                f"Could not build example directory {v} because it already exists. Try the --restore option to still build it."
            ) from e
