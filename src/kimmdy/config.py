import yaml
import logging
from pathlib import Path
from dataclasses import dataclass


@dataclass
class Config:
    raw: dict
    cwd: Path = Path.cwd()

    def __init__(self, input_file: Path):
        with open(input_file, "r") as f:
            try:
                self.aw = yaml.safe_load(f)
            except yaml.YAMLError as e:
                print(e)

