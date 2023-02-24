from dataclasses import dataclass, field
from typing import Optional, Union
from omegaconf.errors import ValidationError, MissingMandatoryValue
from omegaconf import OmegaConf
import omegaconf as omg
from pathlib import Path
from kimmdy import plugins
import logging
import sys
from kimmdy.reaction import Reaction


def check_file_exists(p: Path):
    if not p.exists():
        m = "File not found: " + str(p.resolve())
        logging.error(m)
        raise LookupError(m)


@dataclass
class SingleTaskConfig:
    tasks: str
    mult: int = 1


@dataclass
class SequenceConfig:
    tasks: list[SingleTaskConfig] = field(default_factory=list)


@dataclass
class TaskConfig:
    tasks: Union[SingleTaskConfig,SequenceConfig]


@dataclass
class PlumedConfig:
    dat: str
    distances: str


class Sequence(list):
    """A sequence of tasks."""

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
    logfile: Path = Path("kimmdy.log")
    loglevel: str = "DEBUG"
    color: bool = False


@dataclass
class BaseConfig():
    run: int
    experiment: str
    name: str
    dryrun: bool
    iterations: int
    out: Path
    ff: Path
    ffpatch: Optional[Path]
    top: Path
    gro: Path
    idx: Path
    minimization: MinimizationConfig
    equilibration: EquilibrationConfig
    equilibrium: MdConfig
    changer: ChangerConfig
    reactions: ReactionsConfig
    prod: ProdConfig
    plumed: Optional[PlumedConfig]
    sequence: list
    logging: LoggingConf = LoggingConf()


@dataclass
class Config(BaseConfig):
    """Internal representation of the configuration generated
    from the input file, which enables validation before running
    and computationally expensive operations.
    All settings read from the input file are accessible through nested attributes.
    """


def get_config(opts: Union[Config, dict] = {}) -> Config:
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
            Nested arguments in kimmdy.yml can be accessed with `.`.
            e.g. `kimmdy logging.loglevel=DEBUG`.
        """
        )
        exit()

    plugin_configs = []
    if len(plugins) > 0:
        logging.info("Loading Plugins")
        for plg_name, plugin in plugins.items():
            logging.debug(f"Loading {plg_name}")
            if isinstance(plugin, Exception):
                logging.warn(
                    f"Plugin {plg_name} could not be loaded!\n{plugin}\n"
                )
            # this next type error can be ignored
            print(plugin)
            if issubclass(plugin, Reaction):
                c = OmegaConf.create(plugin.type_scheme)
                if c:
                    plugin_configs.append(c)

    try:
        base_conf = OmegaConf.structured(BaseConfig)
        input = "kimmdy.yml"
        if arglist and (".yml" in arglist[0] or ".yml" in arglist[0]):
            input = arglist[0]
            arglist = arglist[1:]
        yaml_conf = OmegaConf.load(input)
        cli_conf = OmegaConf.from_cli(arglist)
        conf = OmegaConf.merge(base_conf, *plugin_configs, yaml_conf, cli_conf, opts)
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
        check_file_exists(Path(f))

        # TODO:
        # # Validate sequence
        # if isinstance(attr, Sequence):
        #     for task in attr:
        #         assert hasattr(
        #             self, task
        #         ), f"Task {task} listed in sequence, but not defined!"
        #
        # # Validate reaction plugins
        # if attr_name == "reactions":
        #     for r in attr.get_attributes():
        #         assert r in (ks := list(plugins.keys())), (
        #             f"Error: Reaction plugin {r} not found!\n"
        #             + f"Available plugins: {ks}"
        #         )
