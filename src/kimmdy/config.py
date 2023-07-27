"""
Read and validate kimmdy.yml configuration files
and package into a parsed format for internal use.
"""
from __future__ import annotations
from typing import Any, Optional
import yaml
import logging
from pathlib import Path
from kimmdy import plugins
from kimmdy.schema import Sequence, get_combined_scheme


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
        scheme: dict | None = None
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
                # caste type and validate attribute
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

                # cast type
                v = pytype(v)
                self.__setattr__(k, v)

        # set defaults for attributes
        for k, v in scheme.items():
            if "pytype" not in v:
                # recursive case, skip
                continue
            if not hasattr(self, k):
                default = v.get("default")
                if default is None:
                    m = f"No default found for {k}"
                    logging.warning(m)
                    continue
                self.__setattr__(k, default)


    def validate(self):
        """Validates attributes read from config file."""
        logging.info(f"Validating Config")

        for attr in self.__dict__.items():
            print(attr)
            print(type(attr))
            if isinstance(attr, Config):
                attr.validate()

            # Check files from scheme
            if isinstance(attr, Path):
                path = attr.resolve()
                # TODO: do we resolve the path here for self as well?
                # distances.dat wouldn't exist prior to the run
                if not str(attr) in ["distances.dat"]:
                    check_file_exists(path)

            # Validate sequence
            if isinstance(attr, Sequence):
                for task in attr:
                    if not hasattr(self, task):
                        if hasattr(self, "mds"):
                            if hasattr(self.mds, task):
                                continue
                        raise AssertionError(
                            f"Task {task} listed in sequence, but not defined!"
                        )

            # Validate changer reference
            if attr == "raw":
                if hasattr(self, "changer"):
                    if hasattr(self.changer, "coordinates"):
                        if "md" in self.changer.coordinates.get_attributes():
                            assert (
                                self.changer.coordinates.md
                                in self.mds.get_attributes()
                            ), f"Relax MD {self.changer.coordinates.md} not in MD section!"

            # Validate reaction plugins
            if attr == "reactions":
                for reaction_name, reaction_config in attr.__dict__.items():
                    assert reaction_name in (ks := list(plugins.keys())), (
                        f"Error: Reaction plugin {reaction_name} not found!\n"
                        + f"Available plugins: {ks}"
                    )

        # if input_file is not None:
        #     self._set_defaults()
        #     if cwd := self.raw.get("cwd"):
        #         cwd = Path(cwd)
        #     else:
        #         cwd = input_file.parent.resolve()
        #     self.cwd = cwd
        #     if out := self.raw.get("out"):
        #         out = Path(out)
        #     else:
        #         out = self.cwd / self.name
        #     self.out = out
        #     # make sure self.out is empty
        #     while self.out.exists():
        #         logging.debug(f"Output dir {self.out} exists, incrementing name")
        #         out_end = self.out.name[-3:]
        #         if out_end.isdigit():
        #             self.out = self.out.with_name(
        #                 f"{self.out.name[:-3]}{int(out_end)+1:03}"
        #             )
        #         else:
        #             self.out = self.out.with_name(self.out.name + "_001")
        #     self.out.mkdir()
        #     logging.info(f"Created output dir {self.out}")

            # TODO: move to defaults and then valdiation
            # if not hasattr(self, "ff"):
            #     ffs = list(self.cwd.glob("*.ff"))
            #     if len(ffs) < 1:
            #         # TODO: it can work with a ff in the gromacs data dir
            #         # need to re-add the `ff` option but change a bit
            #         # to unify with read_top
            #         raise FileNotFoundError(
            #             "No forcefield found in cwd, please provide one!"
            #         )
            #     if len(ffs) > 1:
            #         logging.warn(
            #             f"Found {len(ffs)} forcefields in cwd, using first one: {ffs[0]}"
            #         )
            #     assert ffs[0].is_dir(), "Forcefield should be a directory!"
            #     self.ff = ffs[0].resolve()

            # TODO: why is this commented out?
            # assert hasattr(self,'mds'), "MD section not defined in config file!"
            # for attribute in self.mds.get_attributes():
            #     self.mds.attr(attribute).mdp = Path(self.mds.attr(attribute).mdp)


    def __repr__(self):
        repr = self.__dict__
        return str(repr)

