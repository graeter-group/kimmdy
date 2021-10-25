import argparse
import logging
import sys

# TODO: The current trend seems to be less hierarchical module structures,
# so I made cmd a file instead of a folder.
# We can discuss this.


def get_args():
    """Parse command line arguments and configure logger"""
    parser = argparse.ArgumentParser(description="Welcome to KIMMDY")
    parser.add_argument(
        "--input", "-i", type=str, help="kimmdy input file", default="kimmdy.yml"
    )
    parser.add_argument(
        "--loglevel",
        "-l",
        type=str,
        help="logging level (DEBUG, INFO, WARNING, ERROR)",
        default="DEBUG",
    )
    parser.add_argument(
        "--logfile", "-f", type=str, help="logfile", default="kimmdy.log"
    )

    args = parser.parse_args()

    logging.basicConfig(
        encoding="utf-8",
        level=getattr(logging, args.loglevel.upper()),
        handlers=[logging.FileHandler(args.logfile), logging.StreamHandler(sys.stdout)],
    )

    return args


def kimmdy_run():
    """Run KIMMDY with a configuration generated form the specified input file."""
    args = get_args()
    logging.info("KIMMDY is running with options:")
    logging.info(args)


if __name__ == "__main__":
    kimmdy_run()

