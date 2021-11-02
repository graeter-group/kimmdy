import logging
import queue
from pathlib import Path
from enum import Enum, auto
from dataclasses import dataclass
from typing import Callable
from kimmdy.config import Config
from kimmdy.reaction import ConversionType, ConversionRecipe


class State(Enum):
    # TODO: find out which states we can have and what to call them
    IDLE = auto()
    MD = auto()
    REACTION = auto()


def default_decision_strategy(rates):
    return "some reaction"

@dataclass
class Task: 
    f: Callable
    kwargs: dict

@dataclass
class RunManager:
    config: Config
    tasks: queue.Queue[Task]
    iteration: int
    trj: Path
    top: Path
    measurements: Path
    rates: list
    state: State

    def __init__(self, input_file: Path):
        self.config = Config(input_file)
        self.tasks = queue.Queue()
        self.state = State.IDLE
        self.iteration = 0

    def run(self):
        logging.info("Start run")

        self.tasks.put(Task(
            self._run_md_minim,
            {''}
        ))

        self.tasks.put(Task(
            self._run_md_eq,
            {}
        ))

    def __next__(self):
        t = self.tasks.get()
        if t:
            t.f(**t.kwargs)

    def _run_md_minim(self):
        logging.info("Start minimization md")
        topfile = self.config.top
        grofile = self.config.gro
        outgro = 'min' + str(n) + '.gro'
        mdpfile = self.config.minization.mdp
        tprfile = "em" + str(n) + ".tpr"
        auto.do_energy_minimisation(grofile, topfile, mdpfile, tprfile, outgro)

    def _run_md_eq(self):
        logging.info("Start equilibration md")

    def _run_md_prod(self):
        logging.info("Start production md")

    def _run_recipe(self, ReactionRecipe):
        logging.info("Start production md")

    def _query_reactions(self):
        logging.info("Query reactions")

    def _decide_reaction(self, decision_strategy=default_decision_strategy):
        logging.info("Decide on a reaction")
        winner = decision_strategy(self.rates)
        return winner


