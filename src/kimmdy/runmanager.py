from __future__ import annotations
from dataclasses import field
import logging
from pathlib import Path
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
AMBIGUOUS_SUFFS = ["dat", "xvg", "log", "trr"]


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

    def __init__(self, config: Config):
        self.config = config
        self.tasks: queue.Queue[Task] = queue.Queue()  # tasks from config
        self.crr_tasks: queue.Queue[Task] = queue.Queue()  # current tasks
        self.iteration = 0
        self.iterations = self.config.iterations
        self.state = State.IDLE
        self.reaction_results: list[ReactionResult] = []
        self.latest_files: dict[str, Path] = {
            "top": self.config.top,
            "gro": self.config.gro,
            "idx": self.config.idx,
            "plumed.dat": self.config.plumed.dat,
            "distances.dat": self.config.plumed.distances,
        }

        self.filehist: list[dict[str, TaskFiles]] = [
            {"setup": TaskFiles(input=self.latest_files)}
        ]

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
        Errors if file is not found.
        """
        logging.debug("Getting latest suffix: " + suffix)
        try:
            path = self.latest_files[suffix]
            logging.debug("Found: " + str(path))
            return path
        except Exception:
            m = f"File {suffix} requested but not found!"
            logging.error(m)
            raise FileNotFoundError(m)

    def run(self):
        logging.info("Start run")
        logging.info("Build task list")

        for task in self.config.sequence:
            logging.debug(f"Put Task: {self.task_mapping[task]}")
            self.tasks.put(Task(self.task_mapping[task]))

        while not (self.state is State.DONE or self.iteration >= self.iterations):
            next(self)

        logging.info(
            f"Stop running tasks, state: {self.state}, iteration:{self.iteration}, max:{self.iterations}"
        )
        logging.info("History:")
        for x in self.filehist:
            for taskname, taskfiles in x.items():
                logging.info(
                    f"""
                Task: {taskname} with output directory: {taskfiles.outputdir}
                Task: {taskname}, input:\n{pformat(taskfiles.input)}
                Task: {taskname}, output:\n{pformat(taskfiles.output)}
                """
                )

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
            logging.info(f"Pretend to run: {task.name} with args: {task.kwargs}")
            return
        logging.debug("Start task: " + pformat(task))
        files = task()
        self._discover_output_files(task.name, files)

    def _discover_output_files(self, taskname, files: TaskFiles):
        """Discover further files written by a task.

        and add those files to the `files` as well as
        the file history and latest files.
        """
        # discover other files written by the task
        for path in files.outputdir.iterdir():
            suffix = path.suffix[1:]
            if suffix in AMBIGUOUS_SUFFS:
                suffix = path.name
            files.output[suffix] = files.outputdir / path

        logging.debug("Update latest files with: ")
        logging.debug(pformat(files.output))
        self.latest_files.update(files.output)
        logging.debug("Append to file history")
        self.filehist.append({taskname: files})

    def _create_task_directory(self, prefix: str) -> TaskFiles:
        files = TaskFiles()
        files.outputdir = self.config.out / f"{prefix}_{self.iteration}"
        files.outputdir.mkdir()
        (files.outputdir / self.config.ff.name).symlink_to(self.config.ff)
        return files

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
        files = self._create_task_directory("equilibration")
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
        logging.info("Setup _run_md_minim")
        self.state = State.MD
        files = self._create_task_directory("minimization")
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
        logging.info("Setup _run_md_eq MD")
        self.state = State.MD
        files = self._create_task_directory("equilibrium")
        files.input = {
            "top": self.get_latest("top"),
            "mdp": self.config.equilibrium.mdp,
            "gro": self.get_latest("gro"),
        }
        files = md.equilibration(files)
        logging.info("Done")
        return files

    def _run_md_prod(self) -> TaskFiles:
        logging.info("Setup _run_md_prod")
        self.state = State.MD
        files = self._create_task_directory("production")
        files.input = {
            "top": self.get_latest("top"),
            "gro": self.get_latest("gro"),
            "mdp": self.config.prod.mdp,
            "idx": self.config.idx,
            "cpt": self.get_latest("cpt"),
            "plumed.dat": self.get_latest("plumed.dat"),
        }
        files = md.production(files)
        logging.info("Done with production MD")
        return files

    def _run_md_relax(self) -> TaskFiles:
        logging.info("Start _run_md_relax")
        self.state = State.MD
        files = self._create_task_directory("relaxation")
        files.input = {
            "top": self.get_latest("top"),
            "gro": self.get_latest("gro"),
            "mdp": self.config.changer.coordinates.md.mdp,
            "idx": self.config.idx,
            "cpt": self.get_latest("cpt"),
        }
        files = md.relaxation(files)
        logging.info("Done with relaxation MD")
        return files

    def _query_reactions(self):
        logging.info("Query reactions")
        self.state = State.REACTION
        files = TaskFiles()

        reactions = self.config.reactions.get_attributes()

        for react_name in reactions:
            reaction = plugins[react_name]

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
        files = self._create_task_directory("recipe")
        files.input = {
            "top": self.get_latest("top"),
            "plumed.dat": self.get_latest("plumed.dat"),
        }

        files.output = {
            "top": files.outputdir / "topol_mod.top",
            "plumed.dat": files.outputdir / "plumed.dat",
            "distances.dat": files.outputdir / "distances.dat",
        }

        changer.modify_top(self.chosen_recipe, files.input["top"], files.output["top"])
        logging.info(f'Wrote new topology to {files.output["top"].parts[-3:]}')
        changer.modify_plumed(
            self.chosen_recipe,
            files.input["plumed.dat"],
            files.output["plumed.dat"],
            files.output["plumed.dist"],
        )
        logging.info(f'Wrote new plumedfile to {files.output["plumed.dat"].parts[-3:]}')
        logging.info(f"Looking for md in {self.config.changer.coordinates.__dict__}")
        # TODO clean this up, maybe make function for this in config
        if hasattr(self.config, "changer"):
            if hasattr(self.config.changer, "coordinates"):
                if hasattr(self.config.changer.coordinates, "md"):
                    self.crr_tasks.put(Task(self._run_md_relax))
        return files
