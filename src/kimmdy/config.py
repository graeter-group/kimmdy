import yaml
import logging
from pathlib import Path
from dataclasses import dataclass
from kimmdy import plugins
from kimmdy.reaction import Reaction


def check_file_exists(p: Path):
    if not p.exists():
        m = "File not found: " + str(p.resolve())
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
    "run": int,
    "dryrun": bool,
    "iterations": int,
    "out": Path,
    "ff": Path,
    "top": Path,
    "gro": Path,
    "idx": Path,
    "plumed": bool,
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
# class PlumedConfig:
#     dat: Path
#     distances: Path


class MinimizationConfig:
    mdp: Path
    tpr: Path


class NvtConfig:
    mdp: Path
    tpr: Path


class NptConfig:
    mdp: Path
    tpr: Path


class EquilibrationConfig:
    nvt: NvtConfig
    npt: NptConfig


class MdConfig:
    mdp: Path


class CoordinatesConfig:
    md: MdConfig


class ChangerConfig:
    coordinates: CoordinatesConfig


class HomolysisConfig:
    edis: Path
    bonds: Path


class ReactionsConfig:
    homolysis: HomolysisConfig


class ProdConfig:
    mdp: Path


class SequenceConfig(list):
    tasks: list


@dataclass
class Config:
    """Internal representation of the configuration generated
    from the input file, which enables validation before running
    and computationally expensive operations.
    All settings read from the input file are accessible through nested attributes.

    Parameters
    ----------
    input_file : Path
        Path to the config yaml file.
    recursive_dict : dict
        For internal use only, used in reading settings in recursively.
    type_scheme : dict
        dict containing types for casting and validating settings.

    """

    run: int
    experiment: str
    name: str
    dryrun: bool
    iterations: int
    out: Path
    ff: Path
    top: Path
    gro: Path
    idx: Path
    plumed: bool
    minimization: MinimizationConfig
    equilibration: EquilibrationConfig
    equilibrium: MdConfig
    changer: ChangerConfig
    reactions: ReactionsConfig
    prod: ProdConfig
    sequence: SequenceConfig

    def __init__(
        self, input_file: Path = None, recursive_dict=None, type_scheme=type_scheme
    ):
        if input_file is None and recursive_dict is None:
            m = "Error: No input file was provided!"
            logging.error(m)
            raise ValueError(m)

        if input_file is not None and not isinstance(input_file, Path):
            logging.debug(
                "Config input file was not type pathlib.Path, attempting conversion.."
            )
            input_file = Path(input_file)

        self.type_scheme = type_scheme
        if self.type_scheme is None:
            self.type_scheme = {}

        # top level before initialization here
        if input_file is not None:
            with open(input_file, "r") as f:
                raw = yaml.safe_load(f)
            self.raw = raw
            if self.raw is None:
                m = "Error: Could not read input file"
                logging.error(m)
                raise ValueError(m)
            recursive_dict = raw

            if len(plugins) > 0:
                logging.info("Loading Plugins:")
                for plg_name, plugin in plugins.items():
                    if isinstance(plugin, Exception):
                        logging.warn(
                            f"Plugin {plg_name} could not be loaded!\n{plugin}\n"
                        )
                    if issubclass(plugin, Reaction):
                        self.type_scheme["reactions"].update(plugin().type_scheme)
                        print(self.type_scheme["reactions"])

        # building config recursively
        if recursive_dict is not None:
            for name, val in recursive_dict.items():
                if isinstance(val, dict):
                    val = Config(
                        recursive_dict=val, type_scheme=self.type_scheme.get(name)
                    )
                logging.debug(f"Set attribute: {name}, {val}")
                self.__setattr__(name, val)

        # top level after initialization here
        if input_file is not None:
            self.cwd = (
                Path(cwd)
                if (cwd := self.raw.get("cwd"))
                else input_file.parent.resolve()
            )
            self.out = (
                Path(out) if (out := self.raw.get("out")) else self.cwd / self.name
            )
            # make sure self.out is empty
            while self.out.exists():
                logging.debug(f"Output dir {self.out} exists, incrementing name")
                out_end = self.out.name[-3:]
                if out_end.isdigit():
                    self.out = self.out.with_name(
                        f"{self.out.name[:-3]}{int(out_end)+1:03}"
                    )
                else:
                    self.out = self.out.with_name(self.out.name + "_001")
            self.out.mkdir()
            logging.info(f"Created output dir {self.out}")

            if not hasattr(self, "ff"):
                ffs = list(self.cwd.glob("*.ff"))
                assert len(ffs) == 1, "Wrong count of forcefields"
                assert ffs[0].is_dir(), "Forcefield should be a directory!"
                self.ff = ffs[0].resolve()

            self._cast_types()
            self._validate()

    def get_attributes(self):
        """Get a list of all attributes"""
        repr = self.__dict__.copy()
        repr.pop("type_scheme")
        return list(repr.keys())

    def __repr__(self):
        repr = self.__dict__.copy()
        repr.pop("type_scheme")
        return str(repr)

    def attr(self, attribute):
        """Get the value of a specific attribute.
        Alias for self.__getattribute__
        """
        return self.__getattribute__(attribute)

    def _cast_types(self):
        """Casts types defined in `type_scheme` to raw attributes."""
        attr_names = filter(lambda s: s[0] != "_", self.__dir__())
        for attr_name in attr_names:
            to_type = self.type_scheme.get(attr_name)
            attr = self.__getattribute__(attr_name)

            if to_type is not None:
                if attr is None:
                    raise ValueError(
                        f"ERROR in inputfile: Missing settings for {attr_name}"
                    )

                # nested:
                if isinstance(to_type, dict):
                    attr._cast_types()
                    continue
                # wrap single element if it should be a list
                elif to_type is list:
                    self.__setattr__(attr_name, [attr])
                try:
                    self.__setattr__(attr_name, to_type(attr))
                except ValueError as e:
                    m = (
                        f"Error: Attribute {attr_name} with value {attr} "
                        + f"and type {type(attr)} cound not be converted to {to_type}!"
                    )
                    logging.error(m)
                    raise ValueError(e)

            else:
                logging.debug(
                    f"{to_type} conversion found for attribute {attr_name} and not executed."
                )

    def _validate(self):
        """Validates attributes read from config file."""
        try:
            attr_names = filter(lambda s: s[0] != "_", self.__dir__())
            for attr_name in attr_names:
                logging.debug(f"validating: {attr_name}")
                attr = self.attr(attr_name)
                if isinstance(attr, Config):
                    attr._validate()

                # Check files from scheme
                if isinstance(attr, Path):
                    self.__setattr__(attr_name, attr.resolve())
                    if not str(attr) in [
                        "distances.dat"
                    ]:  # distances.dat wouldn't exist prior to the run
                        logging.debug(attr)
                        check_file_exists(attr)

                # Check config for consistency
                # if attr_name == "plumed":
                #     for necessary_f in ["dat", "distances"]:
                #         assert (
                #             necessary_f in attr.__dir__()
                #         ), f"{necessary_f} for {attr_name} is missing in config!"

                # Validate sequence
                if isinstance(attr, Sequence):
                    for task in attr:
                        assert hasattr(
                            self, task
                        ), f"Task {task} listed in sequence, but not defined!"

                # Validate reaction plugins
                if attr_name == "reactions":
                    for r in attr.get_attributes():
                        assert r in (ks := list(plugins.keys())), (
                            f"Error: Reaction plugin {r} not found!\n"
                            + f"Available plugins: {ks}"
                        )

        except AssertionError as e:
            logging.error(f"Validating input failed!\n{e}")
            raise ValueError(e)
