import logging
import queue
from pathlib import Path
from enum import Enum, auto
from dataclasses import dataclass
from typing import Callable
from kimmdy.utils import run_shell_cmd
from kimmdy.config import Config
from kimmdy.reaction import ConversionType, ConversionRecipe


class State(Enum):
    IDLE = auto()
    MD = auto()
    REACTION = auto()
    DONE = auto()


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
    state: State
    trj: Path
    top: Path
    measurements: Path
    rates: list
    # md: MDManager

    def __init__(self, input_file: Path):
        self.config = Config(input_file)
        self.tasks = queue.Queue()
        self.iteration = 0
        self.state = State.IDLE
        self.trj = self.config.cwd / Path("prod_" + str(self.iteration) + ".trj")
        self.top = self.config.top
        self.measurements = Path("measurements")
        self.rates = []
        self.md = MDManager()

    def run(self):
        logging.info("Start run")

        self.tasks.put(Task(self._run_md_minim, {}))
        self.tasks.put(Task(self._run_md_eq, {"ensemble": "nvt"}))
        self.tasks.put(Task(self._run_md_eq, {"ensemble": "npt"}))

        while self.state is not State.DONE:
            next(self)

    def __iter__(self):
        return self

    def __next__(self):
        if self.tasks.empty():
            self.state = State.DONE
            return
        t = self.tasks.get()
        t.f(**t.kwargs)

    def _run_md_minim(self):
        logging.info("Start minimization md")
        self.state = State.MD
        self.md.minimzation(self)
        self.state = State.IDLE
        logging.info("Done minimizing")

    def _run_md_eq(self, ensemble):
        logging.info("Start equilibration md")
        self.state = State.MD
        self.md.equilibration(self, ensemble)
        self.state = State.IDLE
        logging.info("Done equilibrating")

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


class MDManager:
    def __init__(self):
        pass

    def write_mdp(self):
        pass

    def minimzation(self, runmgr: RunManager):
        topfile = runmgr.config.top
        grofile = runmgr.config.gro
        outgro = "min" + str(runmgr.iteration) + ".gro"
        mdpfile = runmgr.config.minimization["mdp"]
        tprfile = "em" + str(runmgr.iteration) + ".tpr"

        if runmgr.config.dryrun:
            logging.info("Pretending to run minimization")
            return

        run_shell_cmd(
            f"gmx grompp -p {topfile} -c {grofile} -r {grofile} -f {mdpfile} -o {tprfile}"
        )
        run_shell_cmd(f"gmx mdrun -s {tprfile} -c {outgro}")

    def equilibration(self, runmgr: RunManager, ensemble: str):
        topfile = runmgr.config.top
        grofile = runmgr.config.gro
        outgro = ensemble + '_' + str(runmgr.iteration) + ".gro"
        mdpfile = runmgr.config.equilibration[ensemble]["mdp"]
        tprfile = ensemble + '_' + str(runmgr.iteration) + ".tpr"

        if runmgr.config.dryrun:
            logging.info("Pretending to run equilibration")
            return

        run_shell_cmd(
            f"gmx grompp -p {topfile} -c {grofile} -r {grofile} -f {mdpfile} -o {tprfile}"
        )
        run_shell_cmd(f"gmx mdrun -v -s {tprfile} -c {outgro}")

    def production(self, runmgr: RunManager):
        logging.info("GROMACS go brrrrrr!")
