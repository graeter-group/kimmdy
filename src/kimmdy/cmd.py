import argparse
import logging
from pathlib import Path
from kimmdy.config import Config
from kimmdy.runmanager import RunManager
from kimmdy.utils import check_gmx_version, increment_logfile
import sys


def get_cmdline_args():
    """Parse command line arguments and configure logger.

    Returns
    -------
    Namespace
        parsed command line arguments
    """
    parser = argparse.ArgumentParser(description="Welcome to KIMMDY")
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
        format="\033[34m %(asctime)s\033[00m: %(levelname)s: %(message)s"
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

    config = Config(args.input)
    logging.debug(config)

    logging.debug("Using system GROMACS:")
    logging.debug(check_gmx_version(config))

    runmgr = RunManager(config)
    runmgr.run()


def kimmdy_run(
    input: Path = Path("kimmdy.yml"),
    loglevel: str = "DEBUG",
    logfile: Path = Path("kimmdy.log"),
):
    """Run KIMMDY from python."""
    args = argparse.Namespace(input=input, loglevel=loglevel, logfile=logfile)
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


#%%
if __name__ == "__main__":
    kimmdy_run()
