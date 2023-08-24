"""
Read and validate kimmdy.yml configuration files
and package into a parsed format for internal use.
"""
from __future__ import annotations
from typing import Any, Optional
import yaml
import logging
from pathlib import Path
from kimmdy import reaction_plugins
from kimmdy.schema import Sequence, get_combined_scheme
from kimmdy.utils import get_gmx_dir

logger = logging.getLogger(__name__)

GMX_BUILTIN_FF_DIR = get_gmx_dir() / "top"
"""Path to gromacs data directory with the built-in forcefields."""


def check_file_exists(p: Path):
    if not p.exists():
        m = "File not found: " + str(p.resolve())
        logger.error(m)
        raise LookupError(m)


class Config:
    """Internal representation of the configuration generated
    from the input file, which enables validation before running
    and computationally expensive operations.

    Parameters
    ----------
    input_file :
        Path to the config yaml file.
    recursive_dict :
        For internal use only, used in reading settings in recursively.
    scheme :
        dict containing types and defaults for casting and validating settings.
    section :
        current section e.g. to determine the level of recursion in nested configs
        e.g. "config", "config.mds" or "config.reactions.homolysis"
    """

    # override get and set attributes to satisy type checker and
    # acknowledge that we don't actually statically type-check the attributes
    def __getattribute__(self, name) -> Any:
        return object.__getattribute__(self, name)

    def __setattr__(self, name, value: Any):
        object.__setattr__(self, name, value)

    def __init__(
        self,
        input_file: Path | None = None,
        recursive_dict: dict | None = None,
        scheme: dict | None = None,
        section: str = "config",
        no_increment_output_dir: bool = False,
    ):
        # failure case: no input file and no values from dictionary
        if input_file is None and recursive_dict is None:
            m = "No input file was provided for Config"
            logger.error(m)
            raise ValueError(m)

        # initial scheme
        if scheme is None:
            scheme = get_combined_scheme()

        # read initial input file
        if input_file is not None:
            with open(input_file, "r") as f:
                raw = yaml.safe_load(f)
            if raw is None or not isinstance(raw, dict):
                m = "Could not read input file"
                logger.error(m)
                raise ValueError(m)
            recursive_dict = raw

        assert recursive_dict is not None
        assert scheme is not None
        for k, v in recursive_dict.items():
            if isinstance(v, dict):
                # recursive case
                general_subscheme = scheme.get(".*")
                subscheme = scheme.get(k)
                if subscheme is None:
                    subscheme = general_subscheme
                if general_subscheme is not None:
                    general_subscheme.update(subscheme)
                    subscheme = general_subscheme
                assert subscheme is not None, (k, v, scheme)
                subsection = f"{section}.{k}"
                subconfig = Config(
                    recursive_dict=v, scheme=subscheme, section=subsection
                )
                self.__setattr__(k, subconfig)
            else:
                # base case for recursion

                # merge ".*" and specific schemes
                global_opts = scheme.get(".*")
                opts = scheme.get(k)
                if opts is None and global_opts is None:
                    m = f"Unknown option {section}.{k} found in config file."
                    logger.error(m)
                    raise ValueError(m)
                if opts is None:
                    opts = {}
                if global_opts is not None:
                    opts.update(global_opts)
                pytype = opts.get("pytype")
                if pytype is None:
                    m = f"No type found for {k}"
                    logger.error(m)
                    raise ValueError(m)

                # cast to type
                v = pytype(v)
                self.__setattr__(k, v)

        # set defaults for attributes
        for k, v in scheme.items():
            if "pytype" not in v:
                # recursive case, skip.
                # handled at base case
                continue
            if not hasattr(self, k):
                # get default if not set in yaml
                default = v.get("default")
                pytype = v.get("pytype")
                if default is None:
                    m = f"Option required but no default found for {section}.{k}"
                    logger.debug(m)
                    continue
                if pytype is None:
                    m = f"No type found for default value of {section}.{k}: {default}"
                    logger.error(m)
                    raise ValueError(m)
                default = pytype(default)

                self.__setattr__(k, default)

        # validate on initial construction
        if section == "config":
            self._validate(no_increment_output_dir=no_increment_output_dir)

    def _validate(self, section: str = "config", no_increment_output_dir: bool = False):
        """Validates config."""
        logger.info(f"Validating Config")
        logger.info(f"Validating Config, section {section}")

        # globals / interconnected
        if section == "config":
            # forcfield
            ffdir = self.ff
            if ffdir == Path("*.ff"):
                ffs = list(self.cwd.glob("*.ff"))
                if len(ffs) > 1:
                    logger.warn(
                        f"Found {len(ffs)} forcefields in cwd, using first one: {ffs[0]}"
                    )
                assert ffs[0].is_dir(), "Forcefield should be a directory!"
                ffdir = ffs[0].resolve()
            elif not ffdir.exists():
                print(ffdir)
                ffdir = GMX_BUILTIN_FF_DIR / ffdir
                if not ffdir.exists():
                    m = f"Could not find forcefield {ffdir} in cwd or gromacs data directory"
                    logger.error(m)
                    raise AssertionError(m)
            self.ff = ffdir

            # Validate changer reference
            if hasattr(self, "changer"):
                if hasattr(self.changer, "coordinates"):
                    if "md" in self.changer.coordinates.__dict__.keys():
                        assert (
                            self.changer.coordinates.md in self.mds.__dict__.keys()
                        ), f"Relax MD {self.changer.coordinates.md} not in MD section!"

            # Validate reaction plugins
            if hasattr(self, "reactions"):
                for reaction_name, reaction_config in self.reactions.__dict__.items():
                    assert reaction_name in (ks := list(reaction_plugins.keys())), (
                        f"Error: Reaction plugin {reaction_name} not found!\n"
                        + f"Available plugins: {ks}"
                    )
                    if isinstance(reaction_plugins[reaction_name], Exception):
                        raise reaction_plugins[reaction_name]

            # Validate sequence
            assert hasattr(self, "sequence"), "No sequence defined!"
            for task in self.sequence:
                if not hasattr(self, task):
                    if hasattr(self, "mds"):
                        if hasattr(self.mds, task):
                            continue
                    if hasattr(self.reactions, task):
                        continue
                    raise AssertionError(
                        f"Task {task} listed in sequence, but not defined!"
                    )

            # Validate dat file only once defined
            if hasattr(self, "plumed"):
                if hasattr(self, "md"):
                    if hasattr(self.md, "plumed"):
                        raise AssertionError(
                            "plumed dat file defined multiple times. When "
                            "running md and loading existing measurements, "
                            "only define it once in the md section."
                        )

            if not hasattr(self, "out"):
                self.out = self.cwd / self.name

            # make sure self.out is empty
            while not no_increment_output_dir and self.out.exists():
                logger.debug(f"Output dir {self.out} exists, incrementing name")
                out_end = self.out.name[-3:]
                if out_end.isdigit():
                    self.out = self.out.with_name(
                        f"{self.out.name[:-3]}{int(out_end)+1:03}"
                    )
                else:
                    self.out = self.out.with_name(self.out.name + "_001")
            if not no_increment_output_dir:
                self.out.mkdir()
                logger.info(f"Created output dir {self.out}")

        # individual attributes, recursively
        for name, attr in self.__dict__.items():
            if type(attr) is Config:
                attr._validate(section=f"{section}.{name}")
                continue

            # Check files from scheme
            elif isinstance(attr, Path):
                path = attr
                path = path.resolve()
                self.__setattr__(name, path)
                if not path.is_dir():
                    check_file_exists(path)

    def attr(self, attribute):
        """Get the value of a specific attribute.
        Alias for self.__getattribute__
        """
        return self.__getattribute__(attribute)

    def get_attributes(self):
        """Get a list of all attributes"""
        return list(self.__dict__.keys())

    def __repr__(self):
        repr = self.__dict__
        return str(repr)
