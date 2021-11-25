import yaml
import logging
from pathlib import Path
from dataclasses import dataclass


def check_file_exists(p: Path):
    if not p.exists():
        m = "File not found: " + str(p)
        logging.error(m)
        raise LookupError(m)


@dataclass
class Config:
    """
    Internal representation of the configuration generated
    from the input file, which enables validation before running
    and computationally expensive operations.
    """

    raw: dict
    dryrun: bool
    experiment: str
    cwd: Path
    top: Path
    gro: Path
    idx: Path
    plumed: dict
    #minimization: dict
    #equilibration: dict
    equilibrium: dict
    prod: dict
    changer: dict
    reactions: dict 

    def __init__(self, input_file: Path):
        with open(input_file, "r") as f:
            raw = yaml.safe_load(f)
            self.raw = raw
            if self.raw is None:
                raise ValueError("Could not read input file")

        self.dryrun = raw.get("dryrun")
        self.experiment = raw.get("experiment")
        self.iterations = raw.get("iterations")
        self.cwd = Path(cwd) if (cwd := raw.get("cwd")) else Path.cwd()
        self.top = Path(raw.get("top"))
        self.gro = Path(raw.get("gro"))
        self.idx = Path(raw.get("idx"))
        self.plumed = raw.get("plumed")
        self.minimization = raw.get("minimization")
        self.equilibration = raw.get("equilibration")
        self.equilibrium = raw.get("equilibrium")
        self.prod = raw.get("prod")
        self.changer = raw.get("changer")
        self.reactions = raw.get("reactions")

        

        self._validate()

    def _validate(self):
        for p in [self.top, self.gro]:
            check_file_exists(Path(self.cwd,p))
