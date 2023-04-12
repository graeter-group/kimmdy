from __future__ import annotations
import logging
from pathlib import Path
import dill
import queue
from enum import Enum, auto
from typing import Callable
from kimmdy import config
from kimmdy.config import Config
from kimmdy.utils import increment_logfile
from kimmdy.parsing import read_topol
from kimmdy.reaction import ConversionType, Reaction, ReactionResult, ConversionRecipe
import kimmdy.mdmanager as md
import kimmdy.changemanager as changer
from kimmdy.tasks import Task, TaskFiles, TaskMapping
from kimmdy.utils import run_shell_cmd
from pprint import pformat
import random
from kimmdy import plugins
from kimmdy.topology.topology import Topology

# file types of which there will be multiple files per type
AMBIGUOUS_SUFFS = ["dat", "xvg", "log"]
# are there cases where we have multiple trr files?


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
    for reaction_result in reaction_results:
        for outcome in reaction_result:
            rates.append(outcome.rate)
            recipes.append(outcome.recipe)

    total_rate = sum(rates)
    random.seed()
    t = random.random()  # t in [0.0,1.0)
    logging.debug(f"Random value t: {t}, rates {rates}, total rate {total_rate}")
    rate_running_sum = 0

    # if nothing is choosen, return an empty ConversionRecipe
    result = ConversionRecipe()
    for i, rate in enumerate(rates):
        rate_running_sum += rate
        if (t * total_rate) <= rate_running_sum:
            result = recipes[i]
            break
    logging.debug(f"Result: {result}")

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
        self.from_checkpoint = False
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
            "trr": "",
            "edr": "",
        }
        self.histfile = increment_logfile(Path(f"{self.config.out}_history.log"))
        self.cptfile = increment_logfile(Path(f"{self.config.out}_kimmdy.cpt"))
        try:
            _ = self.config.ffpatch
        except AttributeError:
            self.config.ffpatch = None
        self.top = Topology(
            read_topol(self.config.top), self.config.ff, self.config.ffpatch
        )

        self.filehist: list[dict[str, TaskFiles]] = [
            {"setup": TaskFiles(runmng=self, input=self.latest_files)}
        ]

        self.task_mapping: TaskMapping = {
            "md": self._run_md,
            "reactions": [
                self._query_reactions,
                self._decide_reaction,
                self._run_recipe,
                self._relaxation,
            ],
        }

        # Instantiate reactions
        self.reactions = []
        react_names = self.config.reactions.get_attributes()
        # logging.info("Instantiating Reactions:", *react_names)
        for react_name in react_names:
            r = plugins[react_name]
            reaction = r(react_name, self)
            self.reactions.append(reaction)

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

        if not self.from_checkpoint:
            # allows for mapping one config entry to multiple tasks
            for entry in self.config.sequence:
                if entry in self.config.mds.get_attributes():
                    task = self.task_mapping["md"]
                    logging.info(f"Put Task: {task}")
                    self.tasks.put(Task(task, kwargs={"instance": entry}))
                else:
                    for task in self.task_mapping[entry]:
                        logging.info(f"Put Task: {task}")
                        self.tasks.put(Task(task))

        while not (self.state is State.DONE or self.iteration >= self.iterations):
            logging.info("Write checkpoint before next task")
            with open(self.cptfile, "wb") as f:
                dill.dump(self, f)
            next(self)

        logging.info(
            f"Stop running tasks, state: {self.state}, iteration:{self.iteration}, max:{self.iterations}"
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
        logging.info("Starting task: " + pformat(task))
        files = task()
        self._discover_output_files(task.name, files)

    def _discover_output_files(self, taskname, files: TaskFiles):
        """Discover further files written by a task.

        and add those files to the `files` as well as
        the file history and latest files.
        """
        # discover other files written by the task
        if hasattr(files, "outputdir"):
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

            logging.info(f"Current task files:")
            m = f"""
            Task: {taskname} with output directory: {files.outputdir}
            Task: {taskname}, input:\n{pformat(files.input)}
            Task: {taskname}, output:\n{pformat(files.output)}
            """
            logging.info(m)
            with open(self.histfile, "a") as f:
                f.write(m)

    def _create_task_directory(self, postfix: str) -> TaskFiles:
        """Creates TaskFiles object, output directory and symlinks ff."""
        files = TaskFiles(self)
        files.outputdir = self.config.out / f"{self.iteration}_{postfix}"
        files.outputdir.mkdir(exist_ok=self.from_checkpoint)
        if not (files.outputdir / self.config.ff.name).exists():
            (files.outputdir / self.config.ff.name).symlink_to(self.config.ff)
        return files

    def _dummy(self):
        logging.info("Start dummy task")
        files = TaskFiles(self)
        files.outputdir = self.config.out / f"{self.iteration}_dummy"
        files.outputdir.mkdir()
        run_shell_cmd("pwd>./pwd.pwd", files.outputdir)
        return files

    def _run_md(self, instance) -> TaskFiles:
        """General MD simulation"""
        logging.info(f"Start MD {instance}")
        self.state = State.MD

        files = self._create_task_directory(f"{instance}")
        md_config = self.config.mds.attr(instance)
        gmx_alias = self.config.gromacs_alias

        logging.warning(self.latest_files)
        top = files.input["top"]
        gro = files.input["gro"]
        trr = files.input["trr"]
        edr = files.input["edr"]
        mdp = md_config.mdp
        idx = files.input["idx"]

        outputdir = files.outputdir
        # make maxh and ntomp accessible?
        maxh = 24
        ntomp = 2

        grompp_cmd = f"{gmx_alias} grompp -p {top} -c {gro} -f {mdp} -n {idx} -o {instance}.tpr -maxwarn 5"
        # only appends these lines if there are trr and edr files
        if trr and edr:
            grompp_cmd += f" -t {trr} -e {edr}"
        mdrun_cmd = f"{gmx_alias} mdrun -s {instance}.tpr -cpi {instance}.cpt -x {instance}.xtc -o {instance}.trr -cpo {instance}.cpt -c {instance}.gro -g {instance}.log -e {instance}.edr -px {instance}_pullx.xvg -pf {instance}_pullf.xvg -ro {instance}-rotation.xvg -ra {instance}-rotangles.log -rs {instance}-rotslabs.log -rt {instance}-rottorque.log -maxh {maxh} -dlb yes -ntomp {ntomp}"
        # like this, the previous checkpoint file would not be used, -t and -e options from grompp
        # replace the checkpoint file if gen_vel = no in the mdp file

        if "plumed" in md_config.get_attributes():
            mdrun_cmd += f" -plumed {md_config.plumed.dat}"
            files.output = {"plumed.dat": md_config.plumed.dat}
            # add plumed.dat to output to indicate it as current plumed.dat file

        run_shell_cmd(grompp_cmd, outputdir)
        run_shell_cmd(mdrun_cmd, outputdir)

        logging.info("Done with MD {instance}")
        return files

    def _query_reactions(self):
        logging.info("Query reactions")
        self.state = State.REACTION
        # empty list for every new round of queries
        self.reaction_results: list[ReactionResult] = []

        for reaction in self.reactions:
            # TODO: refactor into Task
            files = self._create_task_directory(reaction.name)

            self.reaction_results.append(reaction.get_reaction_result(files))

        logging.info("Reaction done")
        return files  # necessary?

    def _decide_reaction(
        self,
        decision_strategy: Callable[
            [list[ReactionResult]], ConversionRecipe
        ] = default_decision_strategy,
    ):
        logging.info("Decide on a reaction")
        logging.debug(f"Available reaction results: {self.reaction_results}")
        self.chosen_recipe = decision_strategy(self.reaction_results)
        logging.info("Chosen recipe is:")
        logging.info(self.chosen_recipe)
        return

    def _run_recipe(self) -> TaskFiles:
        logging.info(f"Start Recipe in step {self.iteration}")
        logging.info(f"Recipe: {self.chosen_recipe}")

        files = self._create_task_directory("recipe")

        files.output = {"top": files.outputdir / "topol_mod.top"}

        logging.debug(f"Chose recipe: {self.chosen_recipe}")
        changer.modify_top(
            self.chosen_recipe,
            files.input["top"],
            files.output["top"],
            files.input["ff"],
            self.config.ffpatch,
            self.top,
        )
        logging.info(f'Wrote new topology to {files.output["top"].parts[-3:]}')

        if "plumed.dat" in self.latest_files:
            files.output["plumed.dat"] = files.outputdir / "plumed_mod.dat"
            changer.modify_plumed(
                self.chosen_recipe,
                files.input["plumed.dat"],
                files.output["plumed.dat"],
                files.input["distances.dat"],
            )
            logging.info(
                f'Wrote new plumedfile to {files.output["plumed.dat"].parts[-3:]}'
            )
        logging.info("Reaction done")
        return files

    def _relaxation(self) -> TaskFiles:
        # this task generates no files but is only there to generate subtasks
        logging.info(f"Start Relaxation in step {self.iteration}")
        logging.info(f"Type of relaxation: {self.config.changer.coordinates}")

        if hasattr(self.config.changer.coordinates, "md"):
            self.crr_tasks.put(
                Task(
                    self._run_md,
                    kwargs={"instance": self.config.changer.coordinates.md},
                )
            )
        # TODO add atom placement method from Kai as relaxation option
        logging.info(f"Relaxation done!")
