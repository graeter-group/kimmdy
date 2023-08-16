"""
Functions for starting KIMMDY either from python or the command line.
Also initialized logging and configuration.
"""
import argparse
import logging
from os import chmod
from pathlib import Path
import textwrap
import dill
from kimmdy.config import Config
from kimmdy.misc_helper import _build_examples
from kimmdy.runmanager import RunManager
from kimmdy.utils import check_gmx_version, increment_logfile
import importlib.resources as pkg_resources
import sys
from glob import glob

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
        "--logfile", "-f", type=Path, help="logfile", default="kimmdy.log"
    )
    parser.add_argument(
        "--checkpoint", "-p", type=str, help="start KIMMDY from a checkpoint file"
    )
    parser.add_argument(
        "--from-latest-checkpoint",
        "-c",
        action="store_true",
        help="continue. Start KIMMDY from the latest checkpoint file",
    )
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

    # flag to show available plugins
    parser.add_argument(
        "--show-plugins", action="store_true", help=("List available plugins")
    )

    # flag to print path to yaml schema
    parser.add_argument(
        "--show-schema-path",
        action="store_true",
        help=(
            "Print path to yaml schema for use with yaml-language-server e.g. in VSCode and Neovim"
            "# yaml-language-server: $schema=/path/to/kimmdy-yaml-schema.json"
        ),
    )

    # flag to print an example jobscript for slurm hpc clusters
    parser.add_argument(
        "--generate-jobscript",
        action="store_true",
        help=(
            """
        Instead of running KIMMDY directly, generate at jobscript.sh for slurm HPC clusters.
        You can then run this jobscript with sbatch jobscript.sh
        """
        ),
    )

    return parser.parse_args()


def get_analysis_cmdline_args():
    """
    concat :
    Don't perform a full KIMMDY run but instead concatenate trajectories
    from a previous run.
    """
    parser = argparse.ArgumentParser(
        description="Welcome to the KIMMDY analysis module"
    )
    subparsers = parser.add_subparsers(required=True, metavar="module", dest="module")

    ## trjcat
    parser_trjcat = subparsers.add_parser(
        name="trjcat", help="Concatenate trajectories of a KIMMDY run"
    )
    parser_trjcat.add_argument(
        "dir", type=str, help="KIMMDY run directory to be analysed."
    )
    parser_trjcat.add_argument(
        "--steps",
        "-s",
        nargs="*",
        default="all",
        help=(
            "Apply analysis method to subdirectories with these names. Uses all subdirectories by default"
        ),
    )

    ## plot_energy
    parser_plot_energy = subparsers.add_parser(
        name="plot_energy", help="Plot GROMACS energy for a KIMMDY run"
    )
    parser_plot_energy.add_argument(
        "dir", type=str, help="KIMMDY run directory to be analysed."
    )
    parser_plot_energy.add_argument(
        "--steps",
        "-s",
        nargs="*",
        default="all",
        help=(
            "Apply analysis method to subdirectories with these names. Uses all subdirectories by default"
        ),
    )
    parser_plot_energy.add_argument(
        "--terms",
        "-t",
        nargs="*",
        default=["Potential"],
        help=(
            "Terms from gmx energy that will be plotted. Uses 'Potential' by default"
        ),
    )

    ## radical population
    parser_radical_population = subparsers.add_parser(
        name="radical_population",
        help="Plot population of radicals for one or multiple KIMMDY run(s)",
    )
    parser_radical_population.add_argument(
        "dir", nargs="+", help="KIMMDY run directory to be analysed. Can be multiple."
    )
    parser_radical_population.add_argument(
        "--select_atoms",
        "-a",
        type=str,
        help="Atoms chosen for radical population analysis, default is protein (uses MDAnalysis selection syntax)",
        default="protein",
    )

    # plot rates at each decision step
    parser_plot_rates = subparsers.add_parser(
        name="plot_rates",
        help="Plot rates of all possible reactions after a MD run. Rates must have been saved!",
    )
    parser_plot_rates.add_argument(
        "dir", nargs="+", help="KIMMDY run directory to be analysed. Can be multiple."
    )

    return parser.parse_args()


class longFormatter(logging.Formatter):
    def format(self, record):
        saved_name = record.name  # save and restore for other formatters if desired
        parts = saved_name.split(".")
        if len(parts) > 1:
            record.name = parts[0][0] + "." + ".".join(p[:10] for p in parts[1:])
        else:
            record.name = parts[0]
        result = super().format(record)
        record.name = saved_name
        return result


def configure_logger(log_path, log_level):
    """Configure logging.

    Configures the logging module with optional colorcodes
    for the terminal.

    Parameters
    ----------
    log_path : Path
        File path to log file
    log_level : str
        Log
    """

    increment_logfile(log_path)

    log_conf = {
        "version": 1,
        "formatters": {
            "short": {
                "format": "%(name)-15s %(levelname)s: %(message)s",
                "datefmt": "%H:%M",
            },
            "full": {
                "format": "%(asctime)s %(name)-17s %(levelname)s: %(message)s",
                "datefmt": "%d-%m-%y %H:%M",
            },
            "full_cut": {
                "()": longFormatter,
                "format": "%(asctime)s %(name)-12s %(levelname)s: %(message)s",
                "datefmt": "%d-%m-%y %H:%M",
            },
        },
        "handlers": {
            "cmd": {
                "class": "logging.StreamHandler",
                "formatter": "short",
            },
            "file": {
                "class": "logging.FileHandler",
                "formatter": "full_cut",
                "filename": log_path,
            },
            "null": {
                "class": "logging.NullHandler",
            },
        },
        "loggers": {
            "kimmdy": {
                "level": log_level.upper(),
                "handlers": ["cmd", "file"],
            },
        },
        # Mute others, e.g. tensorflow, matplotlib
        "root": {
            "level": "CRITICAL",
            "handlers": ["null"],
        },
    }
    logging.config.dictConfig(log_conf)


def _run(args: argparse.Namespace):
    """Run kimmdy.

    Parameters
    ----------
    args :
        Command line arguments.
    """
    configure_logger(args.logfile, args.loglevel)

    if args.show_plugins:
        from kimmdy import (
            discovered_reaction_plugins,
            discovered_parameterization_plugins,
        )

        print("Available reaction plugins:")
        for plugin in discovered_reaction_plugins:
            print(plugin)

        print("Available parameterization plugins:")
        for plugin in discovered_parameterization_plugins:
            print(plugin)
        exit()

    if args.show_schema_path:
        path = pkg_resources.files("kimmdy") / "kimmdy-yaml-schema.json"
        print(f"{path}")

        exit()

    logger = logging.getLogger("kimmdy")
    logger.info("Welcome to KIMMDY")
    logger.info("KIMMDY is running with these command line options:")
    logger.info(args)

    if args.generate_jobscript:
        config = Config(args.input)
        runmgr = RunManager(config)
        runmgr.write_one_checkoint()

        content = f"""
        #!/bin/env bash
        #SBATCH --job-name={config.name}
        #SBATCH --output=kimmdy-job.log
        #SBATCH --error=kimmdy-job.log
        #SBATCH --time=24:00:00
        #SBATCH -N 1
        #SBATCH --ntasks-per-node=40
        #SBATCH --mincpus=40
        #SBATCH --exclusive
        #SBATCH --cpus-per-task=1
        #SBATCH --gpus 1
        #SBATCH --mail-type=ALL
        # #SBATCH -p <your-partition>.p
        # #SBATCH --mail-user=<your-email


        # Setup up your environment here
        # modules.sh might load lmod modules, set environment variables, etc.
        source ./_modules.sh

        CYCLE=24

        START=$(date +"%s")

        kimmdy --input {args.input} -c

        END=$(date +"%s")

        LEN=$((END-START))
        HOURS=$((LEN/3600))

        echo "$LEN seconds ran"
        echo "$HOURS full hours ran"

        let "CYCLE--"
        if [ $HOURS -lt $CYCLE ]; then
          echo "last cycle was just $HOURS h long, KIMMDY is done."
          exit 3
        else
          echo "cycle resubmitting"
          sbatch -J {config.name} ./jobscript.sh
          exit 2
        fi
        """
        path = "jobscript.sh"
        with open(path, "w") as f:
            f.write(textwrap.dedent(content))

        chmod(path, 0o755)

        exit()

    if args.from_latest_checkpoint:
        config = Config(args.input, no_increment_output_dir=True)
        cpts = glob(f"{config.name}*kimmdy.cpt")
        if not cpts:
            logging.error("No checkpoints found.")
            exit()
        cpts.sort()
        args.checkpoint = cpts[-1]

    if args.checkpoint:
        logger.info("KIMMDY is starting from a checkpoint.")
        with open(args.checkpoint, "rb") as f:
            runmgr = dill.load(f)
            runmgr.from_checkpoint = True
    else:
        config = Config(args.input)
        logger.debug(config)
        runmgr = RunManager(config)
        logger.debug("Using system GROMACS:")
        logger.debug(check_gmx_version(config))

    runmgr.run()


def kimmdy_run(
    input: Path = Path("kimmdy.yml"),
    loglevel: str = "DEBUG",
    logfile: Path = Path("kimmdy.log"),
    checkpoint: str = "",
    from_latest_checkpoint: bool = False,
    concat: bool = False,
    show_plugins: bool = False,
    show_schema_path: bool = False,
    generate_jobscript: bool = False,
):
    """Run KIMMDY from python.

    Parameters
    ----------
    input :
        kimmdy input yml file.
    loglevel :
        Loglevel. One of ["INFO", "WARNING", "MESSAGE", "DEBUG"]
    logfile :
        File path of the logfile.
    checkpoint :
        File path if a kimmdy.cpt file to restart KIMMDY from a checkpoint.
    from_latest_checkpoint :
        Start KIMMDY from the latest checkpoint.
    concat :
        Don't perform a full KIMMDY run but instead concatenate trajectories
        from a previous run.
    show_plugins :
        Show available plugins and exit.
    show_schema_path :
        Print path to yaml schema for use with yaml-language-server e.g. in VSCode and Neovim
    generate_jobscript :
        Instead of running KIMMDY directly, generate at jobscript.sh for slurm HPC clusters
    """
    args = argparse.Namespace(
        input=input,
        loglevel=loglevel,
        logfile=logfile,
        checkpoint=checkpoint,
        from_latest_checkpoint=from_latest_checkpoint,
        concat=concat,
        show_plugins=show_plugins,
        show_schema_path=show_schema_path,
        generate_jobscript=generate_jobscript,
    )
    _run(args)
    logging.shutdown()


def get_build_example_args():
    """Parse command line arguments.

    Returns
    -------
    Namespace
        parsed command line arguments
    """
    parser = argparse.ArgumentParser(description="Build examples for KIMMDY.")
    parser.add_argument(
        "-r",
        "--restore",
        const=True,
        nargs="?",
        help="Overwrite input files in existing example directories, use keyword 'hard' to also delete output files.",
    )
    return parser.parse_args()


def build_examples():
    """Build examples from the command line."""
    args = get_build_example_args()
    _build_examples(args)
    pass


def analysis():
    """Analyse existing KIMMDY runs."""

    args = get_analysis_cmdline_args()
    print(args)
    from kimmdy.analysis import concat_traj, plot_energy, radical_population, plot_rates

    if args.module == "trjcat":
        concat_traj(args)
    elif args.module == "plot_energy":
        plot_energy(args)
    elif args.module == "radical_population":
        radical_population(args)
    elif args.module == "plot_rates":
        plot_rates(args)


def kimmdy():
    """Run KIMMDY from the command line.

    The configuration is gathered from the input file,
    which is `kimmdy.yml` by default.
    """
    args = get_cmdline_args()
    _run(args)
    logging.shutdown()


if __name__ == "__main__":
    kimmdy_run()
