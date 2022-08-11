from pathlib import Path
import subprocess
from typing import Union


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
