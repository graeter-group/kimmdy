import yaml
import logging
from pathlib import Path
from dataclasses import dataclass
from pprint import pprint


def check_file_exists(p: Path):
    if not p.exists():
        m = "File not found: " + str(p)
        logging.error(m)
        raise LookupError(m)


type_scheme = {
    "dryrun": bool,
    "experiment": str,
    "run": int,
    "cwd": Path,
    "top": Path,
    "gro": Path,
    "plumed": {"dat": Path, "distances": Path},
    "minimization": {"mdp": Path, "tpr": Path},
    "equilibration": {
        "nvt": {"mdp": Path, "tpr": Path},
        "npt": {"mdp": Path, "tpr": Path},
    },
    "reactions": list,
}


class Config:
    """
    Internal representation of the configuration generated
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

    cwd: Path

    def __init__(
        self, input_file: Path = None, recursive_dict=None, type_scheme=type_scheme
    ): 
        if input_file is None and recursive_dict is None:
            m = "Error: No input file was provided!"
            logging.error(m)
            raise ValueError(m)
        
        if input_file is not None and not isinstance(input_file, Path):
            logging.warn("Warning: Config input file was not type pathlib.Path, attemptin conversion..")
            Path(input_file)

        self.type_scheme = type_scheme
        if self.type_scheme is None:
            self.type_scheme = {}

        if input_file is not None:
            with open(input_file, "r") as f:
                raw = yaml.safe_load(f)
                self.raw = raw
                if self.raw is None:
                    m = "Error: Could not read input file"
                    logging.error(m)
                    raise ValueError(m)
                recursive_dict = raw

        if recursive_dict is not None:
            for name, val in recursive_dict.items():
                if isinstance(val, dict):
                    val = Config(
                        recursive_dict=val, type_scheme=self.type_scheme.get(name)
                    )
                logging.debug(f"Set attribute: {name}, {val}")
                self.__setattr__(name, val)

        if input_file is not None:
            Config.cwd = Path(cwd) if (cwd := raw.get("cwd")) else input_file.parent
            self._cast_types()
            self._validate()

    def __repr__(self):
        repr = self.__dict__.copy()
        repr.pop("type_scheme")
        return str(repr)

    def attr(self, attribute):
        """Alias for self.__getattribute__"""
        return self.__getattribute__(attribute)

    def _cast_types(self):
        """Casts types defined in `type_scheme` to raw attributes."""
        attr_names = filter(lambda s: s[0] != "_", self.__dir__())
        for attr_name in attr_names:
            to_type = self.type_scheme.get(attr_name)
            attr = self.__getattribute__(attr_name)

            if to_type is not None:
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
                logging.info(
                    f"{to_type} conversion found for attribute {attr_name} and not executed."
                )

    def _validate(self):
        """Validates attributes read from config file."""
        attr_names = filter(lambda s: s[0] != "_", self.__dir__())
        for attr_name in attr_names:
            logging.debug(f"validating: {attr_name}")
            attr = self.__getattribute__(attr_name)
            if isinstance(attr, Config):
                attr._validate()

            # Check files from scheme
            if isinstance(attr, Path):
                if not attr.is_absolute():
                    attr = Config.cwd / attr
                    self.__setattr__(attr_name, attr)
                check_file_exists(attr)

            # Check config for consistency
            if attr_name in ["nvt", "npt"]:
                for necessary_f in ["mdp", "tpr"]:
                    assert (
                        necessary_f in attr.__dir__()
                    ), f"{necessary_f} for {attr_name} is missing in config!"

            # Checks
            if attr_name == "reactions":
                logging.info(f"DUMMY VALIDATION: There are {len(attr)} reactions!")
