from dataclasses import dataclass
from enum import Enum
from omegaconf.errors import ValidationError, MissingMandatoryValue
from omegaconf import OmegaConf
from pathlib import Path
import logging
import sys


def check_file_exists(p: str):
    path = Path(p)
    if not path.exists():
        m = "File not found: " + p
        logging.error(m)
        raise LookupError(m)


class Sequence(list):
    def __init__(self, tasks: list):
        list.__init__(self)
        for task in tasks:
            if isinstance(task, dict):
                for _ in range(task["mult"]):
                    assert isinstance(
                        task["tasks"], list
                    ), "Grouped tasks must be a list!"
                    self.extend(task["tasks"])
            else:
                self.append(task)


type_scheme = {
    "experiment": str,
    ""
    "dryrun": bool,
    "iterations": int,
    "out": Path,
    "ff": Path,
    "top": Path,
    "gro": Path,
    "idx": Path,
    "plumed": {"dat": Path, "distances": Path},
    "minimization": {"mdp": Path, "tpr": Path},
    "equilibration": {
        "nvt": {"mdp": Path, "tpr": Path},
        "npt": {"mdp": Path, "tpr": Path},
    },
    "equilibrium": {"mdp": Path},
    "prod": {"mdp": Path},
    "changer": {"coordinates": {"md": {"mdp": Path}}},
    "reactions": {},
    "sequence": Sequence,
}

# classes for static code analysis
class PlumedConfig:
    dat: str
    distances: str


@dataclass
class MinimizationConfig:
    mdp: str
    tpr: str


@dataclass
class NvtConfig:
    mdp: str
    tpr: str


@dataclass
class NptConfig:
    mdp: str
    tpr: str


@dataclass
class EquilibrationConfig:
    nvt: NvtConfig
    npt: NptConfig


@dataclass
class MdConfig:
    mdp: str


@dataclass
class CoordinatesConfig:
    md: MdConfig


@dataclass
class ChangerConfig:
    coordinates: CoordinatesConfig


@dataclass
class HomolysisConfig:
    edis: str
    bonds: str


@dataclass
class ReactionsConfig:
    homolysis: HomolysisConfig


@dataclass
class ProdConfig:
    mdp: str


@dataclass
class LoggingConf:
    logfile: str = "kimmdy.log"
    loglevel: str = "DEBUG"
    color: bool = True


@dataclass
class BaseConfig:
    run: int
    experiment: str
    name: str
    dryrun: bool
    iterations: int
    out: str
    ff: str
    top: str
    gro: str
    idx: str
    plumed: PlumedConfig
    minimization: MinimizationConfig
    equilibration: EquilibrationConfig
    equilibrium: MdConfig
    changer: ChangerConfig
    reactions: ReactionsConfig
    prod: ProdConfig
    logging: LoggingConf = LoggingConf()


def get_config() -> BaseConfig:
    """Get the configuration object.
    The BaseConfig is merged with command line arguments and the configuration read from the yaml file.
    All settings read from the input file are accessible through nested attributes.
    """
    args = sys.argv[1:]  # first argument is the program name
    arglist = [a.replace("-", "") for a in args]

    if "h" in arglist or "help" in arglist:
        print(
            """
        Welcome to KIMMDY!
        Usage:
        -  `kimmdy help` for this help text
        -  `kimmdy` to run kimmdy looking for a kimmdy.yml configuration file
        -  `kimmdy <path>` to use a different configuration file
        -  `kimmdy <...>` with configuration keyword=value pairs
            to overwrite the configuration from the command line.
            Nested arguments in kimmdy.yaml can be accessed with `.`.
            e.g. `kimmdy logging.loglevel=DEBUG`.
        """
        )
        exit()

    try:
        base_conf = OmegaConf.structured(BaseConfig)
        input = "kimmdy.yaml"
        if arglist and (".yml" in arglist[0] or ".yaml" in arglist[0]):
            input = arglist[0]
            arglist = arglist[1:]
        yaml_conf = OmegaConf.load(input)
        cli_conf = OmegaConf.from_cli(arglist)
        conf = OmegaConf.merge(base_conf, yaml_conf, cli_conf)
    except MissingMandatoryValue as e:
        logging.error("Missing configuration:")
        logging.error(e)
        raise MissingMandatoryValue(e)
    except ValidationError as e:
        logging.error("Error parsing configuration:")
        logging.error(e)
        raise ValidationError(e)

    # The OmegaConf config object will have like an instance of our BaseConfig,
    # but the type checker does not know that.
    # By telling it that the return type of this function
    # is BaseConfig, we get static type checking in the
    # rest of the program, but have to ignore the 2 following warnings.
    validate_config(conf)
    return conf


def validate_config(conf: BaseConfig):
    for f in [
        conf.ff,
        conf.gro,
        conf.idx,
        conf.top,
        conf.equilibrium.mdp,
        conf.prod.mdp,
    ]:
        check_file_exists(f)
