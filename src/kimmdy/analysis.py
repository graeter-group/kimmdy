from typing import Union
from pathlib import Path
import subprocess
import argparse

from kimmdy.utils import run_shell_cmd


def concat_traj(args: argparse.Namespace):
    """Find and concatenate trajectories (.xtc files) from KIMMDY runs."""
    run_dir = Path(args.dir).expanduser().resolve()
    steps: Union[list, str] = args.steps
    print(run_dir)
    print(steps)

    ## check if step argument is valid
    if not isinstance(steps, list):
        if not steps in ["all"]:
            raise ValueError(f"Steps argument {steps} can not be dealt with.")

    ## create list of subdirectories of run_dir that match the ones named in steps
    print(list(run_dir.glob("*_*/")))
    subdirs_sorted = sorted(
        list(filter(lambda d: d.is_dir(), run_dir.glob("*_*/"))),
        key=lambda p: int(p.name.split("_")[0]),
    )
    print(subdirs_sorted)
    if steps == "all":
        steps = list(set([x.name.split("_")[1] for x in subdirs_sorted]))
    subdirs_matched = list(
        filter(lambda d: d.name.split("_")[1] in steps, subdirs_sorted)
    )
    print(subdirs_matched)

    ## check if there are matched directories
    if not subdirs_matched:
        raise ValueError(
            f"Could not find directories {steps} in {run_dir}. Thus, no trajectories can be concatenated"
        )

    ## create output dir
    (run_dir / "analysis").mkdir(exist_ok=True)
    out = run_dir / "analysis" / "concat.xtc"
    out = Path(out).expanduser()

    ## gather trajectories
    trajectories = []
    tprs = []
    for d in subdirs_matched:
        trajectories.extend(d.glob("*.xtc"))
        tprs.extend(d.glob("*.tpr"))

    # trajectories = list(filter(lambda p: "rotref" not in p.stem, trajectories))
    trajectories = [str(t) for t in trajectories]
    assert (
        len(trajectories) > 0
    ), f"No trrs found to concatenate in {run_dir} with subdirectory names {steps}"

    ## write concatenated trajectory
    run_shell_cmd(
        f"gmx trjcat -f {' '.join(trajectories)} -o {str(out.with_name('tmp.xtc'))} -cat",
        cwd=run_dir,
    )
    run_shell_cmd("sleep 0.5s", cwd=run_dir)
    run_shell_cmd(
        f"echo '1 0' | gmx trjconv -f {str(out.with_name('tmp.xtc'))} -s {tprs[0]} -o {str(out)} -center -pbc mol",
        cwd=run_dir,
    )
