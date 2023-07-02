# %%
import argparse
import logging
from pathlib import Path
import dill
from kimmdy.config import Config
from kimmdy.misc_helper import concat_traj
from kimmdy.runmanager import RunManager
from kimmdy.utils import check_gmx_version, increment_logfile
import sys


if sys.version_info > (3, 10):
    from importlib_metadata import version
else:
    from importlib.metadata import version


def get_cmdline_args():
    """Parse command line arguments and configure logger.

    Returns
    -------
    Namespace
        parsed command line arguments
    """
    parser = argparse.ArgumentParser(description="Welcome to KIMMDY")
    parser.add_argument(
        "--version", action="version", version=f'KIMMDY {version("kimmdy")}'
    )
    parser.add_argument(
        "--input", "-i", type=str, help="kimmdy input file", default="kimmdy.yml"
    )
    parser.add_argument(
        "--loglevel",
        "-l",
        type=str,
        help="logging level (CRITICAL, ERROR, WARNING, INFO, DEBUG)",
        default="DEBUG",
    )
    parser.add_argument(
        "--logfile", "-f", type=str, help="logfile", default="kimmdy.log"
    )
    parser.add_argument("--checkpoint", "-c", type=str, help="checkpoint file")
    parser.add_argument(
        "--concat",
        type=Path,
        nargs="?",
        const=True,
        help=(
            "Concatenate trrs of this run"
            "Optionally, the run directory can be give"
            "Will save as concat.trr in current directory"
        ),
    )
    return parser.parse_args()


def configure_logging(args, color=False):
    """Configure logging.

    Configures the logging module with optional colorcodes
    for the terminal.
    """

    increment_logfile(Path(args.logfile))
    if color:
        logging.addLevelName(logging.INFO, "\033[35mINFO\033[00m")
        logging.addLevelName(logging.ERROR, "\033[31mERROR\033[00m")
        logging.addLevelName(logging.WARNING, "\033[33mWARN\033[00m")
        format = "\033[34m %(asctime)s\033[00m: %(levelname)s: %(message)s"
    else:
        format = "%(asctime)s: %(levelname)s: %(message)s"

    logging.basicConfig(
        level=getattr(logging, args.loglevel.upper()),
        handlers=[
            logging.FileHandler(args.logfile, encoding="utf-8", mode="w"),
            logging.StreamHandler(sys.stdout),
        ],
        format=format,
        datefmt="%d-%m-%Y %H:%M",
    )


def _run(args):
    configure_logging(args)

    logging.info("Welcome to KIMMDY")
    logging.info("KIMMDY is running with these command line options:")
    logging.info(args)

    if args.concat:
        logging.info("KIMMDY will concatenate trrs and exit.")

        run_dir = Path().cwd()
        if type(args.concat) != bool:
            run_dir = args.concat
        concat_traj(run_dir)
        exit()

    if args.checkpoint:
        logging.info("KIMMDY is starting from a checkpoint.")
        with open(args.checkpoint, "rb") as f:
            runmgr = dill.load(f)
            runmgr.from_checkpoint = True
    else:
        config = Config(args.input)
        logging.debug(config)
        runmgr = RunManager(config)
        logging.debug("Using system GROMACS:")
        logging.debug(check_gmx_version(config))

    runmgr.run()


def kimmdy_run(
    input: Path = Path("kimmdy.yml"),
    loglevel: str = "DEBUG",
    logfile: Path = Path("kimmdy.log"),
    checkpoint: str = "",
    concat: bool = False,
):
    """Run KIMMDY from python."""
    args = argparse.Namespace(
        input=input,
        loglevel=loglevel,
        logfile=logfile,
        checkpoint=checkpoint,
        concat=concat,
    )
    _run(args)
    logging.shutdown()


def kimmdy():
    """Run KIMMDY from the command line.

    The configuration is gathered from the input file,
    which is `kimmdy.yml` by default.
    """
    args = get_cmdline_args()
    _run(args)
    logging.shutdown()


# %%
if __name__ == "__main__":
    kimmdy_run()
