from pathlib import Path
import subprocess


def concat_traj(run_dir: Path, out: Path, run_types=None):
    """Find and concatenate trajectories from KIMMDY runs.

    Parameters
    ----------
    run_dir : Path
        Directory containing directories of multiple tasks.
    out : Path
        File Path into the output trr will be written.
    run_types : list, optional
        List of tasks to get trrs from. If a task is not in this list
        it will be skipped. By default None
    """
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

    assert out.suffix == ".trr", "Output file should be a trr file."
    command = f"gmx trjcat -f {' '.join(trrs)} -o {str(out)} -cat".split(" ")

    subprocess.run(command, cwd=run_dir)
