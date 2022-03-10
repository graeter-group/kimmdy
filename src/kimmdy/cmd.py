import argparse
import logging
from kimmdy.config import Config
from kimmdy.runmanager import RunManager
from kimmdy.utils import check_gmx_version
import sys


def get_cmdline_args():
    """Parse command line arguments and configure logger"""
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
    return parser.parse_args("")


def configure_logging(args, color=True):
    if color:
        logging.addLevelName(logging.INFO, "\033[35mINFO\033[00m")
        logging.addLevelName(logging.ERROR, "\033[31mERROR\033[00m")
        logging.addLevelName(logging.WARNING, "\033[33mWARN\033[00m")
    logging.basicConfig(
        level=getattr(logging, args.loglevel.upper()),
        handlers=[
            logging.FileHandler(args.logfile, encoding="utf-8", mode="w"),
            logging.StreamHandler(sys.stdout),
        ],
        format="\033[34m %(asctime)s\033[00m: %(levelname)s: %(message)s",
        datefmt="%d-%m-%Y %H:%M",
    )


def kimmdy():
    """Run KIMMDY with a configuration generated form the specified input file."""
    args = get_cmdline_args()
    configure_logging(args)

    logging.info("Welcome to KIMMDY")

    logging.info("KIMMDY is running with these command line options:")
    logging.info(args)

    config = Config(args.input)

    logging.debug("Using system GROMACS:")
    logging.debug(check_gmx_version())

    runmgr = RunManager(config)
    runmgr.run()


if __name__ == "__main__":
    kimmdy()
