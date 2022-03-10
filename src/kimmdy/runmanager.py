from __future__ import annotations
import logging
import queue
from enum import Enum, auto
from typing import Callable
from kimmdy.config import Config
from kimmdy.reactions.homolysis import Homolysis
from kimmdy.reaction import ReactionResult, ConversionRecipe
import kimmdy.mdmanager as md
import kimmdy.changemanager as changer
from kimmdy.tasks import Task, TaskFiles, TaskMapping
from pprint import pformat
import random
from kimmdy import plugins

# file types of which there will be multiple files per type
AMBIGUOUS_SUFFS = [".dat", ".xvg", ".log", ".trr"]


def default_decision_strategy(
    reaction_results: list[ReactionResult],
) -> ConversionRecipe:
    """Rejection-Free Monte Carlo.
    takes a list of ReactionResults and choses a recipe.

    Parameters
    ---------
    reaction_reults: list[ReactionResults]
        from which one will be choosen
    """
    # compare e.g. <https://en.wikipedia.org/wiki/Kinetic_Monte_Carlo#Rejection-free_KMC>

    # flatten the list of rates form the reaction results
    rates = []
    recipes = []
    for reaction in reaction_results:
        for rate in reaction.rates:
            rates.append(rate)
        for recipe in reaction.recipes:
            recipes.append(recipe)

    total_rate = sum(rates)
    random.seed()
    t = random.random()  # t in [0.0,1.0)
    rate_running_sum = 0

    # if nothing is choosen, return an empty ConversionRecipe
    result = ConversionRecipe()
    for i in range(len(rates)):
        rate_running_sum += rates[i]
        if (t * total_rate) <= rate_running_sum:
            result = recipes[i]

    return result


class State(Enum):
    """State of the system.
    one of IDLE, MD, REACTION, DONE.
    """

    IDLE = auto()
    MD = auto()
    REACTION = auto()
    DONE = auto()


class RunManager:
    """The RunManager, a central piece.
    Manages the queue of tasks, communicates with the
    rest of the program and keeps track of global state.
    """

    reaction_results: list[ReactionResult]

    def __init__(self, config: Config):
        self.config = config
        self.tasks: queue.Queue[Task] = queue.Queue()  # tasks from config
        self.crr_tasks: queue.Queue[Task] = queue.Queue()  # current tasks
        self.iteration = 0
        self.iterations = self.config.iterations
        self.state = State.IDLE
        self.reaction_results = []
        # self.measurements = Path("measurements")
        # self.structure = self.config.gro
        # self.top = self.config.top
        # self.trj = self.config.cwd / ("prod_" + str(self.iteration) + ".trj")
        # self.plumeddat = self.config.plumed.dat
        # self.plumeddist = self.config.plumed.distances
        self.filehist: list[TaskFiles] = []
        self.filehist.append(
            TaskFiles(
                input={
                    "top": self.config.top,
                    "gro": self.config.gro,
                    "idx": self.config.idx,
                    "plumed.dat": self.config.plumed.dat,
                    "distances.dat": self.config.plumed.distances,
                }
            )
        )

        self.task_mapping: TaskMapping = {
            "equilibrium": self._run_md_equil,
            "prod": self._run_md_prod,
            "minimization": self._run_md_minim,
            "relax": self._run_md_relax,
            "reactions": self._query_reactions,
        }

        logging.debug("Configuration from input file:")
        logging.debug(pformat(self.config.__dict__))

    def get_latest(self, suffix: str):
        """Returns path to latest file of given type.
        For .dat files (in general ambiguous extensions) use full file name.
        Errors if file is nof found.
        """
        # TODO we shouldn't have to step backwards through the complete
        # history of IO actions. We should just be able to maintain a dictionary
        # of the latest file for each type.
        if suffix in AMBIGUOUS_SUFFS:
            logging.warn(f"{suffix} ambiguous! Please specify full file name!")
        for step in reversed(self.filehist):
            for path in list(step.input.values()) + list(step.output.values()):
                if suffix in str(path):
                    return path
        else:
            m = f"File {suffix} requested but not found!"
            logging.error(m)
            raise FileNotFoundError(m)

    def run(self):
        logging.info("Start run")
        logging.info("Building task list")

        for task in self.config.sequence:
            logging.debug(f"Put Task: {self.task_mapping[task]}")
            self.tasks.put(Task(self.task_mapping[task]))

        while not (self.state is State.DONE or self.iteration >= self.iterations):
            next(self)

        logging.info(
            f"Stop running tasks, state: {self.state}, iteration:{self.iteration}, max:{self.iterations}"
        )
        logging.info(f"History:\n{pformat(self.filehist)}")

    def __iter__(self):
        return self

    def __next__(self):
        if self.tasks.empty() and self.crr_tasks.empty():
            self.state = State.DONE
            return
        if not self.crr_tasks.empty():
            task = self.crr_tasks.get()
        else:
            task = self.tasks.get()
            self.iteration += 1
        if self.config.dryrun:
            logging.info(f"Pretending to run: {task.name} with args: {task.kwargs}")
            return
        files = task()

        if files.outputdir:
            # list files written by the task
            for path in files.outputdir.iterdir():
                suffix = path.suffix
                if suffix in AMBIGUOUS_SUFFS:
                    suffix = path.name
                files.output[suffix] = path

        # TODO Ohhh, I might be trying to implement the same thing twice
        # Need to clarify what part is responsible of keeping track of
        # the files a task needs and creates!
        self.filehist.append(files)

    def _dummy(self):
        logging.info("Start dummy task")
        files = TaskFiles()
        files.outputdir = self.config.out / f"dummy_{self.iteration}"
        files.outputdir.mkdir()
        files.input = {
            "top": self.get_latest("top"),
            "gro": self.get_latest("gro"),
        }
        md.dummy_step(files)
        return files

    def _run_md_equil(self) -> TaskFiles:
        logging.info("Start equilibration MD")
        self.state = State.MD
        files = TaskFiles()
        outputdir = self.config.out / f"equil_{self.iteration}"
        outputdir.mkdir()
        (outputdir / self.config.ff.name).symlink_to(self.config.ff)
        files.outputdir = outputdir
        files.input = {
            "top": self.get_latest("top"),
            "gro": self.get_latest("gro"),
            "mdp": self.config.equilibrium.mdp,
            "idx": self.config.idx,
        }
        md.equilibrium(files)
        logging.info("Done equilibrating")
        return files

    def _run_md_minim(self) -> TaskFiles:
        logging.info("Start minimization md")
        self.state = State.MD
        files = TaskFiles()
        outputdir = self.config.out / f"min_{self.iteration}"
        outputdir.mkdir()
        (outputdir / self.config.ff.name).symlink_to(self.config.ff)
        files.outputdir = outputdir
        files.input = {
            "top": self.get_latest("top"),
            "gro": self.get_latest("gro"),
            "mdp": self.config.minimization.mdp,
        }
        # perform step
        files = md.minimzation(files)
        logging.info("Done minimizing")
        return files

    def _run_md_eq(self) -> TaskFiles:
        # TODO combine w/ other equilibration?
        logging.info("Start _run_md_eq MD")
        self.state = State.MD
        files = TaskFiles()
        outputdir = self.config.out / f"equil_{self.iteration}"
        outputdir.mkdir()
        (outputdir / self.config.ff.name).symlink_to(self.config.ff)
        files.outputdir = outputdir
        files.input = {
            "top": self.get_latest("top"),
            "mdp": self.config.equilibrium.mdp,
            "gro": self.get_latest("gro"),
        }
        files = md.equilibration(files)
        logging.info("Done equilibrating")
        return files

    def _run_md_prod(self) -> TaskFiles:
        logging.info("Start production MD")
        self.state = State.MD
        files = TaskFiles()
        outputdir = self.config.out / f"prod_{self.iteration}"
        outputdir.mkdir()
        (outputdir / self.config.ff.name).symlink_to(self.config.ff)
        files.outputdir = outputdir
        files.input = {
            "top": self.get_latest("top"),
            "gro": self.get_latest("gro"),
            "mdp": self.config.prod.mdp,
            "idx": self.config.idx,
            "cpt": self.get_latest("cpt"),
            "plumed.dat": self.get_latest("plumed.dat"),
        }
        files = md.production(files)
        logging.info("Done minimizing")
        return files

    def _run_md_relax(self) -> TaskFiles:
        logging.info("Start relaxation MD")
        self.state = State.MD
        files = TaskFiles()
        outputdir = self.config.out / f"prod_{self.iteration}"
        outputdir.mkdir()
        (outputdir / self.config.ff.name).symlink_to(self.config.ff)
        files.outputdir = outputdir
        files.input = {
            "top": self.get_latest("top"),
            "gro": self.get_latest("gro"),
            "mdp": self.config.changer.coordinates.md.mdp,
            "idx": self.config.idx,
            "cpt": self.get_latest("cpt"),
        }
        files = md.relaxation(files)
        logging.info("Done simulating")
        return files

    def _query_reactions(self):
        logging.info("Query reactions")
        self.state = State.REACTION
        files = TaskFiles()

        reactions = self.config.reactions.get_attributes()

        for react_name in reactions:
            reaction = plugins.get(react_name)
            if reaction is None:
                logging.warning(
                    f"Reaction {react_name} could not be executed! Plugin not found!"
                )
                continue

            # TODO: Make this general for all reactions.
            # Maybe with a dict keeping all the newest files.
            files.input = {
                "plumed.dat": self.get_latest("plumed.dat"),
                "distances.dat": self.get_latest("distances.dat"),
                "top": self.get_latest("top"),
                "ffbonded.itp": self.config.reactions.homolysis.bonds,
                "edissoc.dat": self.config.reactions.homolysis.edis,
            }
            self.reaction_results.append(reaction().get_reaction_result(files))

        logging.info("Reaction done")
        return files

    def _decide_reaction(
        self,
        decision_strategy: Callable[
            [list[ReactionResult]], ConversionRecipe
        ] = default_decision_strategy,
    ):
        logging.info("Decide on a reaction")
        self.chosen_recipe = decision_strategy(self.reaction_results)
        logging.info("Chosen recipe is:")
        logging.info(self.chosen_recipe)
        self.crr_tasks.put(Task(self._run_recipe))
        return None, None

    def _run_recipe(self) -> TaskFiles:
        logging.info(f"Start Recipe in step {self.iteration}")
        logging.info(f"Breakpair: {self.chosen_recipe.atom_idx}")

        files = TaskFiles()
        outputdir = self.config.out / f"recipe_{self.iteration}"
        outputdir.mkdir()
        (outputdir / self.config.ff.name).symlink_to(self.config.ff)
        files.outputdir = outputdir

        files.input = {
            "top": self.get_latest("top"),
            "plumed.dat": self.get_latest("plumed.dat"),
        }

        newtop = outputdir / "topol_mod.top"
        newplumeddat = outputdir / "plumed.dat"
        newplumeddist = outputdir / "distances.dat"
        files.output = {
            "top": newtop,
            "plumed.dat": newplumeddat,
            "dist.dat": newplumeddist,
        }

        changer.modify_top(self.chosen_recipe, files.input["top"], newtop)
        logging.info(f"Wrote new topology to {newtop.parts[-3:]}")
        changer.modify_plumed(
            self.chosen_recipe, files.input["plumed.dat"], newplumeddat, newplumeddist
        )
        logging.info(f"Wrote new plumedfile to {newplumeddat.parts[-3:]}")
        logging.info(f"Looking for md in {self.config.changer.coordinates.__dict__}")
        # TODO clean this up, maybe make function for this in config
        if hasattr(self.config, "changer"):
            if hasattr(self.config.changer, "coordinates"):
                if hasattr(self.config.changer.coordinates, "md"):
                    self.crr_tasks.put(Task(self._run_md_relax))
        return files
