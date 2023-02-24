import logging
from typing import Optional, Union
from kimmdy.config import Config
import sys
from kimmdy.runmanager import RunManager
from kimmdy.utils import check_gmx_version
from kimmdy.config import BaseConfig, get_config


def configure_logging(config: BaseConfig):
    """Configure logging.

    Configures the logging module with optional colorcodes
    for the terminal.
    """

    if config.logging.logfile.exists():
        log_curr = config.logging.logfile
        while log_curr.exists():
            out_end = log_curr.name.strip("#")[-3:]
            if out_end.isdigit():
                log_curr = log_curr.with_name(
                    f"#{log_curr.name[:-3]}{int(out_end)+1:03}#"
                )
            else:
                log_curr = log_curr.with_name(f"#{log_curr.name}_001#")

        config.logging.logfile.rename(log_curr)

    if config.logging.color:
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

    logging.debug("Using system GROMACS:")
    logging.debug(check_gmx_version(config))

    runmgr = RunManager(config)
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
