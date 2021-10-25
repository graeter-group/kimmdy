import yaml
import logging
from pathlib import Path
from dataclasses import dataclass


@dataclass
class Config:
    raw: dict
    experiment: str
    cwd: Path
    reactions: list[str]

    def __init__(self, input_file: Path):
        with open(input_file, 'r') as f:
            try:
                raw = yaml.safe_load(f)
                self.raw = raw
            except yaml.YAMLError as e:
                print(e)

        self.experiment = raw['experiment']
        self.cwd = Path(cwd) if (cwd := raw.get('cwd')) else Path.cwd()
        # TODO: validate reactions have a matching module and are loaded
        self.reactions = raw.get('reactions')
