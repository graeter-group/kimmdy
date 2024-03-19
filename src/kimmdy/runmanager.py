"""
The Runmanager is the main entry point of the program.

It manages the queue of tasks, communicates with the
rest of the program and keeps track of global state.
"""

from __future__ import annotations

import logging
import queue
import time
from copy import copy, deepcopy
from datetime import timedelta
from enum import Enum, auto
from functools import partial
from pathlib import Path
from pprint import pformat
from typing import Optional
import re

import dill

from kimmdy.config import Config
from kimmdy.coordinates import break_bond_plumed, merge_top_slow_growth, place_atom
from kimmdy.kmc import KMCResult, extrande, extrande_mod, frm, rf_kmc
from kimmdy.parsing import read_top, write_json, write_top
from kimmdy.plugins import (
    BasicParameterizer,
    ReactionPlugin,
    parameterization_plugins,
    reaction_plugins,
)
from kimmdy.recipe import Bind, Break, CustomTopMod, Place, RecipeCollection, Relax
from kimmdy.tasks import Task, TaskFiles, get_plumed_out
from kimmdy.topology.topology import Topology
from kimmdy.topology.utils import get_is_reactive_predicate_f
from kimmdy.utils import run_gmx, truncate_sim_files

logger = logging.getLogger(__name__)

# file types of which there will be multiple files per type
AMBIGUOUS_SUFFS = ["dat", "xvg", "log", "itp", "mdp"]
# file strings which to ignore
IGNORE_SUBSTR = ["_prev.cpt", r"step\d+[bc]\.pdb"]
# are there cases where we have multiple trr files?


class State(Enum):
    """State of the system.
    one of IDLE, MD, REACTION, SETUP, DONE.
    """

    SETUP = auto()
    IDLE = auto()
    MD = auto()
    REACTION = auto()
    DONE = auto()


def get_existing_files(config: Config) -> dict:
    """Initialize latest_files with every existing file defined in config"""
    file_d = {}
    attr_names = config.get_attributes()
    for attr_name in attr_names:
        attr = getattr(config, attr_name)
        if isinstance(attr, Path):
            if attr.is_dir() or not attr.exists():
                continue
            # suffix without dot
            key = attr.suffix[1:]

            # AMBIGUOUS_SUFFS -> key whole name
            if key in AMBIGUOUS_SUFFS:
                key = attr.name

            if attr_name == "plumed":
                key = "plumed"
                file_d["plumed_out"] = attr.parent / get_plumed_out(attr)

            file_d[key] = attr
        elif isinstance(attr, Config):
            file_d.update(get_existing_files(attr))
    return file_d


class RunManager:
    """The Runmanager is the main entry point of the program.

    Manages the queue of tasks, communicates with the
    rest of the program and keeps track of global state.

    Attributes
    ----------
    config
        The configuration object.
    from_checkpoint
        Whether the runmanager was initialized from a checkpoint.
    tasks
        Tasks from config.
    crr_tasks
        Current tasks.
    iteration
        Current iteration.
    state
        Current state of the system.
    recipe_collection
        Collection of recipes.
    latest_files
        Dictionary of latest files.
    histfile
        Path to history file.
    cptfile
        Path to checkpoint file.
    top
        Topology object.
    filehist
        List of dictionaries of TaskFiles.
    task_mapping
        Mapping of task names to runmanager methods.
    reaction_plugins
        List of initialized reaction plugins used in the sequence.
    """

    def __init__(self, config: Config):
        self.config: Config = config
        self.from_checkpoint: bool = False
        self.tasks: queue.Queue[Task] = queue.Queue()  # tasks from config
        self.crr_tasks: queue.Queue[Task] = queue.Queue()  # current tasks
        self.iteration: int = -1  # start at -1 to have iteration 0 be the initial setup
        self.state: State = State.IDLE
        self.recipe_collection: RecipeCollection = RecipeCollection([])
        self.kmcresult: Optional[KMCResult] = None
        self.time: float = 0.0  # [ps]
        self.latest_files: dict[str, Path] = get_existing_files(config)
        logger.debug(f"Initialized latest files:\n{pformat(self.latest_files)}")
        self.histfile: Path = self.config.out / "kimmdy.history"
        self.cptfile: Path = self.config.out / "kimmdy.cpt"
        self.kmc_algorithm: str

        try:
            if self.config.changer.topology.parameterization == "basic":
                parameterizer = BasicParameterizer()
            else:
                parameterizer = parameterization_plugins[
                    self.config.changer.topology.parameterization
                ]()
        except KeyError as e:
            raise KeyError(
                f"The parameterization tool chosen in the configuration file: "
                f"'{self.config.changer.topology.parameterization}' can not be found in "
                f"the parameterization plugins: {list(parameterization_plugins.keys())}"
            ) from e

        self.top = Topology(
            top=read_top(self.config.top, self.config.ff),
            parametrizer=parameterizer,
            is_reactive_predicate_f=get_is_reactive_predicate_f(
                self.config.topology.reactive
            ),
            radicals=getattr(self.config, "radicals", None),
            residuetypes_path=getattr(self.config, "residuetypes", None),
        )
        self.filehist: list[dict[str, TaskFiles]] = [
            {"setup": TaskFiles(self.get_latest)}
        ]

        # Initialize reaction plugins used in the sequence
        self.reaction_plugins: list[ReactionPlugin] = []
        for name in self.config.reactions.get_attributes():
            logger.debug(f"Initializing reaction: {name}")
            Plugin = reaction_plugins[name]
            reaction_plugin = Plugin(name, self)
            self.reaction_plugins.append(reaction_plugin)

        self.kmc_mapping = {
            "extrande": extrande,
            "rfkmc": rf_kmc,
            "frm": frm,
            "extrande_mod": extrande_mod,
        }

        self.task_mapping = {
            "md": {"f": self._run_md, "kwargs": {}, "out": None},
            "reactions": [
                {"f": self._place_reaction_tasks, "kwargs": {}, "out": None},
                {
                    "f": self._decide_recipe,
                    "kwargs": {},
                    "out": "decide_recipe",
                },
                {"f": self._apply_recipe, "kwargs": {}, "out": "apply_recipe"},
            ],
        }
        """Mapping of task names to functions and their keyword arguments."""

    def run(self):
        logger.info("Start run")
        self.start_time = time.time()

        if self.from_checkpoint:
            logger.info("KIMMDY is starting from a checkpoint.")
        else:
            self._setup_tasks()

        while (
            self.state is not State.DONE
            and (self.iteration <= self.config.max_tasks or self.config.max_tasks == 0)
            and (
                (time.time() - self.start_time) / 3600 < self.config.max_hours
                or self.config.max_hours == 0
            )
        ):
            if self.config.write_checkpoint:
                logger.info("Writing checkpoint before next task")
                with open(self.cptfile, "wb") as f:
                    dill.dump(self, f)
            next(self)

        logger.info(
            f"Finished running tasks, state: {self.state} after "
            f"{timedelta(seconds=(time.time() - self.start_time))}"
        )

    def _setup_tasks(self):
        """Populates the tasks queue.
        Allows for mapping one sequence entry in the config to multiple tasks
        """
        logger.info("Building task list")
        # setup task
        task = Task(
            self,
            f=self._setup,
            kwargs={},
            out="setup",
        )
        self.tasks.put(task)
        # configured sequence
        for step in self.config.sequence:
            if step in self.config.mds.get_attributes():
                # entry is a type of MD
                md = self.task_mapping["md"]
                kwargs: dict = copy(md["kwargs"])
                kwargs.update({"instance": step})
                task = Task(
                    self,
                    f=md["f"],
                    kwargs=kwargs,
                    out=step,
                )
                self.tasks.put(task)

            elif step in self.config.reactions.get_attributes():
                # entry is a single reaction
                task_list = copy(self.task_mapping["reactions"])
                # 0 is place_reaction_tasks
                task_list[0] = copy(task_list[0])
                task_list[0]["kwargs"] = {"selected": step}

                for task_kwargs in task_list:
                    self.tasks.put(Task(self, **task_kwargs))
            elif step == "reactions":
                # check all reactions
                for task_kwargs in self.task_mapping["reactions"]:
                    self.tasks.put(Task(self, **task_kwargs))
            else:
                m = f"Unknown task encountered in the sequence: {step}"
                logger.error(m)
                raise ValueError(m)
        logger.info(f"Task list build:\n{pformat(list(self.tasks.queue), indent=8)}")

    def write_one_checkpoint(self):
        """Just write the first checkpoint and then exit

        Used to generate a starting point for jobscripts on hpc clusters
        that can easily self-submit after a timelimit was exceeded.
        """
        logger.info("Initial setup for first checkpoint")
        self._setup_tasks()
        with open(self.cptfile, "wb") as f:
            dill.dump(self, f)
        logger.info(f"Wrote checkpointfile to: {self.cptfile}, ")

    def get_latest(self, suffix: str):
        """Returns path to latest file of given type.

        For .dat files (in general ambiguous extensions) use full file name.
        Errors if file is not found.
        """
        logger.debug("Getting latest suffix: " + suffix)
        try:
            path = self.latest_files[suffix]
            logger.debug("Found: " + str(path))
            return path
        except Exception:
            m = f"File {suffix} requested but not found!"
            logger.error(m)
            raise FileNotFoundError(m)

    def __iter__(self):
        return self

    def __next__(self):
        current_time = time.time()
        if self.tasks.empty() and self.crr_tasks.empty():
            self.state = State.DONE
            return
        if not self.crr_tasks.empty():
            task = self.crr_tasks.get()
        else:
            task = self.tasks.get()
        if self.config.dryrun:
            logger.info(f"Pretend to run: {task.name} with args: {task.kwargs}")
            return
        logger.info(f"Starting task: {task.name} with args: {task.kwargs}")
        files = task()
        logger.info(
            f"Finished task: {task.name} after "
            f"{timedelta(seconds=(time.time() - current_time))}"
        )
        if files is not None:
            self._discover_output_files(task.name, files)

    def _discover_output_files(
        self, taskname: str, files: TaskFiles
    ) -> Optional[TaskFiles]:
        """Discover further files written by a task.

        and add those files to the `files` as well as
        the file history and latest files.
        """
        # discover other files written by the task
        if hasattr(files, "outputdir"):
            # check whether double suffs are properly defined in files by the task
            discovered_files = [
                p
                for p in files.outputdir.iterdir()
                if not any(re.search(s, p.name) for s in IGNORE_SUBSTR)
            ]
            suffs = [p.suffix[1:] for p in discovered_files]
            counts = [suffs.count(s) for s in suffs]
            for suff, c in zip(suffs, counts):
                if c != 1 and suff not in AMBIGUOUS_SUFFS:
                    if files.output.get(suff) is None:
                        e = (
                            "ERROR: Task produced multiple files with same suffix but "
                            "did not define with which to continue!\n"
                            f"Task {taskname}, Suffix {suff} found {c} times"
                        )
                        logger.error(e)
                        raise RuntimeError(e)

            # discover output files
            for path in files.outputdir.iterdir():
                suffix = path.suffix[1:]
                if suffix in AMBIGUOUS_SUFFS:
                    suffix = path.name
                # don't overwrite manually added keys in files.output
                if files.output.get(suffix) is not None:
                    continue
                files.output[suffix] = files.outputdir / path

            # remove double entries
            if "plumed" in files.input.keys():
                if plumed := files.output.get("plumed"):
                    files.output.pop(plumed.name)
                if plumed_out := files.output.get("plumed_out"):
                    files.output.pop(plumed_out.name)

            logger.debug(f"Update latest files with:\n{pformat(files.output)}")
            self.latest_files.update(files.output)
            self.filehist.append({taskname: files})

            m = f"""
            Task: {taskname} with output directory: {files.outputdir}
            Task: {taskname}, input:\n{pformat(files.input)}
            Task: {taskname}, output:\n{pformat(files.output)}
            """
            with open(self.histfile, "a") as f:
                f.write(m)

            return files
        else:
            logger.debug("No output directory found for task: " + taskname)
            return None

    def _setup(self, files: TaskFiles) -> TaskFiles:
        """A setup task to collect files processed by kimmdy such as the topology"""
        logger = files.logger
        logger.info("Start setup task")
        self.state = State.SETUP
        logger.info("Writing initial topology after parsing")

        if self.config.parameterize_at_setup:
            focus_nrs = set(self.top.atoms.keys())
            self.top.needs_parameterization = True
            self.top.update_parameters(focus_nrs)

        write_top(self.top.to_dict(), files.outputdir / self.config.top.name)
        files.output["top"] = files.outputdir / self.config.top.name
        logger.info("Done with setup")
        return files

    def _run_md(self, instance: str, files: TaskFiles) -> TaskFiles:
        """General MD simulation"""
        logger = files.logger
        logger.info(f"Start MD {instance}")
        self.state = State.MD

        md_config = self.config.mds.attr(instance)
        gmx_alias = self.config.gromacs_alias
        gmx_mdrun_flags = self.config.gmx_mdrun_flags

        top = files.input["top"]
        gro = files.input["gro"]
        files.input["mdp"] = md_config.mdp
        mdp = files.input["mdp"]
        ndx = files.input["ndx"]

        outputdir = files.outputdir

        grompp_cmd = (
            f"{gmx_alias} grompp -p {top} -c {gro} "
            f"-f {mdp} -n {ndx} -o {instance}.tpr -maxwarn 5"
        )

        # optional files for grompp:
        if self.latest_files.get("trr") is not None:
            trr = files.input["trr"]
            grompp_cmd += f" -t {trr}"
        ## disable use of edr for now
        # if self.latest_files.get("edr") is not None:
        #     edr = files.input["edr"]
        #     grompp_cmd += f" -e {edr}"

        mdrun_cmd = (
            f"{gmx_alias} mdrun -s {instance}.tpr -cpi {instance}.cpt "
            f"-x {instance}.xtc -o {instance}.trr -cpo {instance}.cpt "
            f"-c {instance}.gro -g {instance}.log -e {instance}.edr "
            f"-px {instance}_pullx.xvg -pf {instance}_pullf.xvg "
            f"-ro {instance}-rotation.xvg -ra {instance}-rotangles.log "
            f"-rs {instance}-rotslabs.log -rt {instance}-rottorque.log "
            f"{gmx_mdrun_flags}  "
        )

        if getattr(md_config, "use_plumed"):
            mdrun_cmd += f" -plumed {files.input['plumed']}"

            plumed_out = files.outputdir / get_plumed_out(files.input["plumed"])
            files.output["plumed_out"] = plumed_out

        # specify trr to prevent rotref trr getting set as standard trr
        files.output["trr"] = files.outputdir / f"{instance}.trr"
        logger.debug(f"grompp cmd: {grompp_cmd}")
        logger.debug(f"mdrun cmd: {mdrun_cmd}")
        run_gmx(grompp_cmd, outputdir)
        run_gmx(mdrun_cmd, outputdir)
        logger.info(f"Done with MD {instance}")
        return files

    def _place_reaction_tasks(self, selected: Optional[str] = None) -> None:
        logger.info("Start query reactions")
        self.state = State.REACTION
        # empty list for every new round of queries
        self.recipe_collection: RecipeCollection = RecipeCollection([])

        # placing tasks in priority queue
        strategies = []
        for reaction_plugin in self.reaction_plugins:
            if selected is not None:
                if reaction_plugin.name != selected:
                    continue

            self.crr_tasks.put(
                Task(
                    self,
                    f=self._query_reaction,
                    kwargs={"reaction_plugin": reaction_plugin},
                    out=reaction_plugin.name,
                )
            )

            # get supported kmc algorithm
            strategies.append(reaction_plugin.config.kmc)

        # find decision strategy
        if len(kmc := self.config.kmc) > 0:
            # algorithm overwrite
            self.kmc_algorithm = kmc

        else:
            # get algorithm from reaction
            if len(set(strategies)) > 1:
                raise RuntimeError(
                    "Incompatible kmc algorithms chosen in the same reaction.\n"
                    "Split the reactions in separate steps or choose different algorithms\n"
                    "Attempted to combine:\n"
                    f"{ {rp.name:rp.config.kmc for rp in self.reaction_plugins} }"
                )
            self.kmc_algorithm = strategies[0]

        logger.info(f"Queued {len(self.reaction_plugins)} reaction plugin(s)")
        return None

    def _query_reaction(
        self, reaction_plugin: ReactionPlugin, files: TaskFiles
    ) -> TaskFiles:
        logger = files.logger
        logger.info(f"Start query {reaction_plugin.name}")

        new_recipes = reaction_plugin.get_recipe_collection(files).recipes
        self.recipe_collection.recipes.extend(new_recipes)

        logger.info(
            f"Done with Query reactions, {len(new_recipes)} "
            f"recipes recived from {reaction_plugin.name}"
        )

        return files

    def _decide_recipe(
        self,
        files: TaskFiles,
    ) -> TaskFiles:
        logger = files.logger
        logger.info(
            f"Start Decide recipe using {self.kmc_algorithm}, "
            f"{len(self.recipe_collection.recipes)} recipes available."
        )
        kmc = self.kmc_mapping[self.kmc_algorithm.lower()]

        # FIXME Hotfix for #355 aggregate not working for big systems
        if "rfkmc" != self.kmc_algorithm.lower():
            self.recipe_collection.aggregate_reactions()
        if "extrande" in self.kmc_algorithm.lower():
            kmc = partial(kmc, tau_scale=self.config.tau_scale)
        self.kmcresult = kmc(self.recipe_collection, logger=logger)
        recipe = self.kmcresult.recipe

        if self.config.save_recipes:
            self.recipe_collection.to_csv(files.outputdir / "recipes.csv", recipe)

        try:
            if self.config.plot_rates:
                kwargs = {
                    "outfile": files.outputdir / "reaction_rates.svg",
                    "highlight_r": recipe,
                }
                if (self.kmcresult.time_start is not None) and (
                    self.kmcresult.time_start != 0
                ):
                    kwargs["highlight_t"] = self.kmcresult.time_start
                self.recipe_collection.plot(**kwargs)
        except Exception as e:
            logger.warning(f"Error occured during plotting:\n{e}")

        if self.kmcresult.time_delta:
            self.time += self.kmcresult.time_delta
        logger.info("Done with Decide recipe.")
        if len(recipe.rates) == 0:
            logger.info("No reaction selected")
        elif self.kmcresult.time_delta:
            logger.info(
                f"Overall time {self.time*1e-12:.4e} s, reaction occured after {self.kmcresult.time_delta*1e-12:.4e} s"
            )

        # capture state of radicals
        write_json(
            {
                "overall_time": self.time,
                "residence_time": self.kmcresult.time_delta,
                "radicals": list(self.top.radicals.keys()),
            },
            files.outputdir / "radicals.json",
        )

        return files

    def _apply_recipe(self, files: TaskFiles) -> TaskFiles:
        logger = files.logger

        if self.kmcresult is None:
            m = "Attempting to _apply_recipe without having chosen one with _decide_recipe."
            logger.error(m)
            raise RuntimeError(m)

        recipe = self.kmcresult.recipe
        logger.info(f"Start Recipe in KIMMDY iteration {self.iteration}")
        logger.info(f"Recipe: {recipe.get_recipe_name()}")
        logger.debug(f"Performing recipe steps:\n{pformat(recipe.recipe_steps)}")

        # Set time to chosen 'time_start' of KMCResult
        ttime = self.kmcresult.time_start
        if any([isinstance(step, Place) for step in recipe.recipe_steps]):
            # only first time of interval is valid for placement
            ttime = recipe.timespans[0][0]

        truncate_sim_files(files, ttime)

        top_initial = deepcopy(self.top)
        focus_nrs = set()
        for step in recipe.recipe_steps:
            if isinstance(step, Break):
                self.top.break_bond((step.atom_id_1, step.atom_id_2))
                focus_nrs.update([step.atom_id_1, step.atom_id_2])
                if hasattr(self.config, "plumed"):
                    break_bond_plumed(
                        files,
                        (step.atom_id_1, step.atom_id_2),
                        files.outputdir / self.config.plumed.name.replace(".", "_mod."),
                    )
            elif isinstance(step, Bind):
                self.top.bind_bond((step.atom_id_1, step.atom_id_2))
                focus_nrs.update([step.atom_id_1, step.atom_id_2])
            elif isinstance(step, Place):
                task = Task(
                    self,
                    f=place_atom,
                    kwargs={"step": step, "ttime": None},
                    out="place_atom",
                )
                logger.info(f"Starting task: {task.name} with args: {task.kwargs}")
                place_files = task()
                logger.info(f"Finished task: {task.name}")
                if place_files is not None:
                    self._discover_output_files(task.name, place_files)
                focus_nrs.update([step.id_to_place])

            elif isinstance(step, Relax):
                logger.info("Starting relaxation md as part of reaction..")
                if not hasattr(self.config.changer.coordinates, "md"):
                    logger.warning("Relax task requested but no MD specified for it!")
                    continue

                if self.config.changer.coordinates.slow_growth:
                    # Create a slow growth topology for sub-task run_md, afterwards, top will be reset properly
                    self.top.update_parameters(focus_nrs)
                    top_merge = merge_top_slow_growth(top_initial, deepcopy(self.top))
                    top_merge_path = files.outputdir / self.config.top.name.replace(
                        ".", "_mod."
                    )
                    write_top(top_merge.to_dict(), top_merge_path)
                    self.latest_files["top"] = top_merge_path
                instance = self.config.changer.coordinates.md
                task = Task(
                    self, f=self._run_md, kwargs={"instance": instance}, out=instance
                )
                logger.info(f"Starting task: {task.name} with args: {task.kwargs}")
                md_files = task()
                logger.info(f"Finished task: {task.name}")
                if md_files is not None:
                    self._discover_output_files(task.name, md_files)

            elif isinstance(step, CustomTopMod):
                step.f(self.top)

        self.top.update_partial_charges(recipe.recipe_steps)
        self.top.update_parameters(focus_nrs)

        write_top(self.top.to_dict(), files.outputdir / self.config.top.name)
        files.output["top"] = files.outputdir / self.config.top.name

        # Recipe done, reset runmanger state
        self.kmcresult = None

        logger.info("Done with Apply recipe")
        return files

    def __repr__(self):
        return "Runmanager"
