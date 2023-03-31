import logging
from typing import Optional, Union
from pathlib import Path
import dill
from kimmdy.config import Config
import sys
from kimmdy.runmanager import RunManager
from kimmdy.utils import check_gmx_version
from kimmdy.config import BaseConfig, get_config


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
    return parser.parse_args()


def configure_logging(config: BaseConfig):
    """Configure logging.

    Configures the logging module with optional colorcodes
    for the terminal.
    """

    config.logging.logfile = increment_logfile(Path(config.logging.logfile))

    if color:
        logging.addLevelName(logging.INFO, "\033[35mINFO\033[00m")
        logging.addLevelName(logging.ERROR, "\033[31mERROR\033[00m")
        logging.addLevelName(logging.WARNING, "\033[33mWARN\033[00m")
        format = "\033[34m %(asctime)s\033[00m: %(levelname)s: %(message)s"
    else:
        format = " %(asctime)s: %(levelname)s: %(message)s"

    logging.basicConfig(
        level=getattr(logging, config.logging.loglevel.upper()),
        handlers=[
            logging.FileHandler(config.logging.logfile, encoding="utf-8", mode="w"),
            logging.StreamHandler(sys.stdout),
        ],
        format=format,
        datefmt="%d-%m-%Y %H:%M",
    )


def _run(opts: Union[Config, dict] = {}):
    config = get_config(opts)
    configure_logging(config)

    logging.info("Welcome to KIMMDY")
    logging.info("KIMMDY is running with these options:")
    logging.debug(config)

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


def kimmdy_run(opts: Union[Config, dict] = {}):
    """Run KIMMDY from python."""
    _run(opts)
    logging.shutdown()


def kimmdy():
    """Run KIMMDY from the command line.

    The configuration is gathered from the input file,
    which is `kimmdy.yml` by default.
    """
    _run()
    logging.shutdown()


if __name__ == "__main__":
    kimmdy_run()
