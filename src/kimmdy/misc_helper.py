from pathlib import Path
import subprocess
from typing import Union
from kimmdy.topology.topology import Topology


def topology_to_edgelist(top: Topology):
    ls = []
    for b in sorted(top.bonds):
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


def concat_traj(run_dir: Union[Path, str], out: Union[Path, str], run_types=None):
    """Find and concatenate trajectories from KIMMDY runs.

    Parameters
    ----------
    run_dir : Path
        Directory containing directories of multiple tasks.
    out : Path
        File Path into the output trr will be written.
        Can be relative to run_dir
    run_types : list, optional
        List of tasks to get trrs from. If a task is not in this list
        it will be skipped. By default None
    """
    run_dir = Path(run_dir).expanduser().resolve()
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
    trrs = [str(t) for t in trrs]

    command = f"gmx trjcat -f {' '.join(trrs)} -o {str(out)} -cat".split(" ")

    subprocess.run(command, cwd=run_dir)
