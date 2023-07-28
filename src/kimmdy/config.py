"""
Read and validate kimmdy.yml configuration files
and package into a parsed format for internal use.
"""
from __future__ import annotations
from typing import Any, Optional
import yaml
import logging
from pathlib import Path, PosixPath
from kimmdy import plugins
from kimmdy.schema import Sequence, get_combined_scheme
from kimmdy.utils import get_gmx_dir


GMX_BUILTIN_FF_DIR = get_gmx_dir() / "top"
"""Path to gromacs data directory with the built-in forcefields."""


def check_file_exists(p: Path):
    if not p.exists():
        m = "File not found: " + str(p.resolve())
        logging.error(m)
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
    ):
        # failure case: no input file and no values from dictionary
        if input_file is None and recursive_dict is None:
            m = "No input file was provided for Config"
            logging.error(m)
            raise ValueError(m)

        if scheme is None:
            scheme = get_combined_scheme()

        if input_file is not None:
            with open(input_file, "r") as f:
                raw = yaml.safe_load(f)
            if raw is None or not isinstance(raw, dict):
                m = "Could not read input file"
                logging.error(m)
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
                assert subscheme is not None
                subconfig = Config(recursive_dict=v, scheme=subscheme)
                self.__setattr__(k, subconfig)
            else:
                # base case for recursion

                # merge ".*" and specific schemes
                global_opts = scheme.get(".*")
                opts = scheme.get(k)
                if opts is None and global_opts is None:
                    m = f"No scheme found for {k}"
                    logging.error(m)
                    raise ValueError(m)
                if opts is None:
                    opts = {}
                if global_opts is not None:
                    opts.update(global_opts)
                pytype = opts.get("pytype")
                if pytype is None:
                    m = f"No type found for {k}"
                    logging.error(m)
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
                    m = f"No default found for {k}"
                    logging.debug(m)
                    continue
                if pytype is None:
                    m = f"No type found for default value {k}: {default}"
                    logging.error(m)
                    raise ValueError(m)
                default = pytype(default)

                self.__setattr__(k, default)

    def validate(self, section: str = "config"):
        """Validates config."""
        logging.info(f"Validating Config")
        logging.info(f"Validating Config, section {section}")

        # globals / interconnected
        if section == "config":
            # forcfield
            ffdir = self.ff
            if ffdir == Path("*.ff"):
                ffs = list(self.cwd.glob("*.ff"))
                if len(ffs) > 1:
                    logging.warn(
                        f"Found {len(ffs)} forcefields in cwd, using first one: {ffs[0]}"
                    )
                assert ffs[0].is_dir(), "Forcefield should be a directory!"
                ffdir = ffs[0].resolve()
            elif not ffdir.exists():
                print(ffdir)
                ffdir = GMX_BUILTIN_FF_DIR / ffdir
                if not ffdir.exists():
                    m = f"Could not find forcefield {ffdir} in cwd or gromacs data directory"
                    logging.error(m)
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
                    assert reaction_name in (ks := list(plugins.keys())), (
                        f"Error: Reaction plugin {reaction_name} not found!\n"
                        + f"Available plugins: {ks}"
                    )

            # Validate sequence
            assert hasattr(self, "sequence"), "No sequence defined!"
            for task in self.sequence:
                if not hasattr(self, task):
                    if hasattr(self, "mds"):
                        if hasattr(self.mds, task):
                            continue
                    raise AssertionError(
                        f"Task {task} listed in sequence, but not defined!"
                    )

            if not hasattr(self, "out"):
                self.out = self.cwd / self.name

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

        # individual attributes, recursively
        for name, attr in self.__dict__.items():
            if type(attr) is Config:
                attr.validate(section=f"{section}.{name}")
                continue

            # Check files from scheme
            elif isinstance(attr, Path):
                path = attr
                # distances.dat wouldn't exist prior to the run
                if not str(attr) in ["distances.dat"] and not path.is_dir():
                    check_file_exists(path)

    def __repr__(self):
        repr = self.__dict__
        return str(repr)
