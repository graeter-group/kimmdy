"""
The Runmanager is the main entry point of the program.

It manages the queue of tasks, communicates with the
rest of the program and keeps track of global state.
"""

from __future__ import annotations

import logging
import queue
import re
import shutil
import time
from copy import copy, deepcopy
from datetime import timedelta
from enum import Enum, auto
from functools import partial
from pathlib import Path
from pprint import pformat
from subprocess import CalledProcessError
from typing import Callable, Optional

from kimmdy.config import Config
from kimmdy.constants import (
    MARK_DONE,
    MARK_FAILED,
    MARK_FINISHED,
    MARK_STARTED,
    MARKERS,
)
from kimmdy.coordinates import break_bond_plumed, merge_top_slow_growth, place_atom
from kimmdy.kmc import (
    KMCError,
    KMCReject,
    KMCAccept,
    KMCResult,
    dummy_first_kmc,
    extrande,
    extrande_mod,
    frm,
    rf_kmc,
    multi_rfkmc,
    total_index_to_index_within_plugin,
)
from kimmdy.parsing import read_top, write_json, write_time_marker, write_top
from kimmdy.plugins import (
    BasicParameterizer,
    ReactionPlugin,
    parameterization_plugins,
    reaction_plugins,
)
from kimmdy.recipe import (
    Bind,
    Break,
    CustomTopMod,
    DeferredRecipeSteps,
    Place,
    RecipeCollection,
    Relax,
)
from kimmdy.tasks import Task, TaskFiles, get_plumed_out
from kimmdy.topology.topology import Topology
from kimmdy.topology.utils import get_is_reactive_predicate_from_config_f
from kimmdy.utils import (
    flatten_recipe_collections,
    get_task_directories,
    run_gmx,
    write_coordinate_files_at_reaction_time,
    write_reaction_time_marker,
)

logger = logging.getLogger(__name__)

# file types of which there will be multiple files per type
AMBIGUOUS_SUFFS = ["dat", "xvg", "log", "itp", "mdp"]
# file strings which to ignore
IGNORE_SUBSTR = [
    "_prev.cpt$",
    r"step\d+[bc]\.pdb$",
    r"\.tail$",
    r"_relax\.top$",
    r"_before_with_solvent_bonds\.top$",
    r"_before\.top$",
    r"_after\.top$",
    r"\.\d+#$",
    r"\.log$",
    "rotref",
    r"^\.",  # all hidden files
] + MARKERS
# are there cases where we have multiple trr files?
TASKS_WITHOUT_DIR = ["place_reaction_task"]


class State(Enum):
    """State of the system.
    one of IDLE, MD, REACTION, SETUP, DONE.
    """

    SETUP = auto()
    IDLE = auto()
    MD = auto()
    REACTION = auto()
    DONE = auto()


def get_existing_files(config: Config, section: str = "config") -> dict:
    """Initialize latest_files with every existing file defined in config"""
    files = {}
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

            files[key] = attr
        elif isinstance(attr, Config):
            files.update(get_existing_files(attr, section=attr_name))

    # for starting reactions from existing trajectories,
    # we need to read the name of the plumed output file
    # from the plumed input file and interpret it
    # relative to the parent directory of the trajectory
    # or coordinate files (xtc, trr, tpr, gro).
    if section == "config" and hasattr(config, "plumed"):
        for f in ["xtc", "tpr", "trr", "gro"]:
            if hasattr(config, f):
                plumed_out = get_plumed_out(files["plumed"])
                trajectory_dir = files[f].parent
                files["plumed_out"] = trajectory_dir / plumed_out
                logger.debug(
                    f"Setting up plumed_out from {f} file. Plumed out: {plumed_out} in {trajectory_dir}"
                )
                break

    return files


class RunManager:
    """The Runmanager is the main entry point of the program.

    Manages the queue of tasks, communicates with the
    rest of the program and keeps track of global state.

    Attributes
    ----------
    config
        The configuration object.
    tasks
        Tasks from config.
    priority_tasks
        Additional tasks added during the run by other tasks.
    iteration
        Current iteration.
    state
        Current state of the system.
    recipe_collections
        Dictionary of recipe collections. Keyed by the name of the reaction plugin.
    latest_files
        Dictionary of latest files.
    histfile
        Path to history file.
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
        self.tasks: queue.Queue[Task] = queue.Queue()
        self.priority_tasks: queue.Queue[Task] = queue.Queue()
        self.iteration: int = -1  # start at -1 to have iteration 0 be the initial setup
        self.state: State = State.IDLE
        self.recipe_collections: dict[str, RecipeCollection] = {}
        self.kmcresult: KMCResult | None = None
        self.time: float = 0.0  # [ps]
        self.latest_files: dict[str, Path] = get_existing_files(config)
        logger.debug(f"Initialized latest files:\n{pformat(self.latest_files)}")
        self.histfile: Path = self.config.out / "kimmdy.history"
        self.cptfile: Path = self.config.out / "kimmdy.cpt"
        self.kmc_algorithm: str

        with open(self.histfile, "w") as f:
            f.write("KIMMDY task file history\n")
            f.write(
                "Filepaths in the output directory are shortened to be relative to the output directory.\n\n"
            )

        logger.info(f"Initialized KIMMDY at cwd: {config.cwd}")
        logger.info(f"with output directory {config.out}")
        try:
            if self.config.changer.topology.parameterization == "basic":
                self.parameterizer = BasicParameterizer()
            else:
                self.parameterizer = parameterization_plugins[
                    self.config.changer.topology.parameterization
                ](**self.config.changer.topology.parameterization_kwargs.__dict__)
        except KeyError as e:
            raise KeyError(
                f"The parameterization tool chosen in the configuration file: "
                f"'{self.config.changer.topology.parameterization}' can not be found in "
                f"the parameterization plugins: {list(parameterization_plugins.keys())}"
            ) from e

        nrexcl = getattr(self.config.topology, "nrexcl", None)
        self.top = Topology(
            top=read_top(self.config.top, self.config.ff),
            parametrizer=self.parameterizer,
            is_reactive_predicate_f=get_is_reactive_predicate_from_config_f(
                self.config.topology.reactive
            ),
            radicals=getattr(self.config, "radicals", None),
            residuetypes_path=getattr(self.config, "residuetypes", None),
            reactive_nrexcl=nrexcl,
            gromacs_alias=self.config.gromacs_alias,
        )
        self.filehist: list[dict[str, TaskFiles]] = [
            {"0_setup": TaskFiles(self.get_latest)}
        ]

        # Initialize reaction plugins used in the sequence
        self.reaction_plugins: list[ReactionPlugin] = []
        for name in self.config.reactions.get_attributes():
            logger.debug(f"Initializing reaction: {name}")
            Plugin = reaction_plugins[name]
            reaction_plugin = Plugin(name=name, runmng=self)
            self.reaction_plugins.append(reaction_plugin)

        self.kmc_mapping: dict[str, Callable[..., KMCResult]] = {
            "extrande": extrande,
            "rfkmc": rf_kmc,
            "frm": frm,
            "extrande_mod": extrande_mod,
            "multi_rfkmc": multi_rfkmc,
            "dummy_first": dummy_first_kmc,
        }

        self.task_mapping = {
            "md": {
                "f": self._run_md,
                "kwargs": {},
                "out": None,
            },  # name of out director is inferred from the instance
            "reactions": [
                {
                    "f": self._place_reaction_tasks,
                    "kwargs": {},
                    "out": None,
                },  # has no output directory
                {
                    "f": self._decide_recipe,
                    "kwargs": {},
                    "out": "decide_recipe",
                },
                {"f": self._apply_recipe, "kwargs": {}, "out": "apply_recipe"},
            ],
        }
        """Mapping of task names to functions and their keyword arguments."""

        self._setup_tasks()

    def run(self):
        logger.info("Start run")
        self.start_time = time.time()

        if self.config.restart:
            logger.info(f"Restarting from previous run in: {self.config.out.name}")
            if (self.config.out / MARK_FINISHED).exists():
                m = f"Run in {self.config.out} already finished. Exiting."
                logger.info(m)
                return
            self._setup_restart()

        while (
            self.state is not State.DONE
            and (self.iteration <= self.config.max_tasks or self.config.max_tasks == 0)
            and (
                (time.time() - self.start_time) / 3600 < self.config.max_hours
                or self.config.max_hours == 0
            )
        ):
            next(self)

        write_time_marker(self.config.out / MARK_FINISHED, "finished")
        logger.info(
            f"Finished running last task, state: {self.state} after "
            f"{timedelta(seconds=(time.time() - self.start_time))} "
            f"In output directory {self.config.out}"
        )

    def _setup_restart(self):
        """Set up RunManager to restart from an existing run directory"""

        task_dirs = get_task_directories(self.config.out)
        if task_dirs == []:
            # no tasks found in the output directory. this is a fresh run
            return

        logger.info(
            f"Found task directories in existing output directory ({self.config.out.name}): {[p.name for p in task_dirs]}"
        )
        logger.info(f"Task queue: {self.tasks}")

        found_restart_point = False
        restart_task_name = None
        restart_from_incomplete = False

        # restarting during or after MD except for relaxation/slow_growth MDs are valid restart points
        if hasattr(self.config.changer, "coordinates") and hasattr(
            self.config.changer.coordinates, "md"
        ):
            relax_md_name = self.config.changer.coordinates.md
        else:
            relax_md_name = None
        if hasattr(self.config, "mds"):
            md_task_names = [
                name
                for name in self.config.mds.get_attributes()
                if name != relax_md_name
            ]
        else:
            md_task_names = []

        # keep track of how often we have completed each md instance
        # such that we can later restart correctly at the latest instance
        md_instance_dir_counter = {}
        for name in md_task_names:
            md_instance_dir_counter[name] = 0

        # discover completed or half completed tasks
        logger.info("Checking for restart point in existing task dirs.")
        for task_dir in task_dirs:
            task_n, task_name = task_dir.name.split(sep="_", maxsplit=1)
            task_n = int(task_n)
            logger.info(f"Checking task: {task_n}_{task_name}")
            if (task_dir / MARK_FAILED).exists():
                m = (
                    f"Task in directory `{task_dir.name}` is indicated to "
                    "have failed. Aborting restart. Remove this task "
                    "directory if you want to restart from before the failed task."
                )
                logger.warning(m)
                inp = input(
                    "Do you want to continue and delete this task directory? [y/n]"
                )
                if inp.lower() != "y":
                    exit(1)
            elif (
                (task_dir / MARK_STARTED).exists()
                and not (task_dir / MARK_DONE).exists()
                and task_name in md_task_names
            ):
                # Continue from started but not finished md task is a valid restart point
                # if it got so far that is has written at least one checkpoint file
                checkpoint_files = list(task_dir.glob("*.cpt"))
                if len(checkpoint_files) == 0:
                    m = f"Last started but not done task is {task_dir.name}, but no checkpoint file found. Using an earlier MD or setup task as restart point instead."
                    logger.warning(m)
                    continue
                logger.info(f"Found started but not finished task {task_dir.name}.")
                logger.info(f"Will continue task {task_dir.name}")
                found_restart_point = True
                self.iteration = task_n
                restart_task_name = task_name
                restart_from_incomplete = True
                md_instance_dir_counter[task_name] += 1
                # no need to search further, as this task is unfinished
                break
            elif (
                (task_dir / MARK_STARTED).exists()
                and (task_dir / MARK_DONE).exists()
                and task_name in md_task_names
            ):
                # Continue after last finished MD task is a valid restart point
                # but continue searching for newer tasks after this
                logger.info(f"Found completed task {task_dir.name}")
                logger.info(f"Will continue after task {task_dir.name}")
                found_restart_point = True
                self.iteration = task_n
                restart_task_name = task_name
                restart_from_incomplete = False
                md_instance_dir_counter[task_name] += 1
            elif (
                (task_dir / MARK_STARTED).exists()
                and (task_dir / MARK_DONE).exists()
                and task_name == "setup"
            ):
                # Continuing just after 0_setup is a valid restart point
                found_restart_point = True
                self.iteration = task_n
                restart_task_name = task_name
                restart_from_incomplete = False
            elif (task_dir / MARK_STARTED).exists() and (task_dir / MARK_DONE).exists():
                # Completed task, but not an MD task
                pass
            elif (task_dir / MARK_STARTED).exists() and not (
                task_dir / MARK_DONE
            ).exists():
                # Started but not done task, not an MD task
                m = f"Last started but not done task is {task_dir.name}, which can not be restarted from. Restarting instead from after the last completed MD task."
                logger.info(m)
                break
            else:
                m = f"Encountered task directory {task_dir.name} but kimmdy does not know how to handle this task. Aborting restart."
                logger.error(m)
                raise RuntimeError(m)

        if not found_restart_point or not restart_task_name:
            m = "No valid restart point found in existing task directories."
            logger.error(m)
            raise RuntimeError(m)

        m = f"Restarting from iteration (task number) {self.iteration} with name {restart_task_name}"
        logger.info(m)

        # pop from the task queue until the restart point
        # all md tasks from which we may restart are in the task queue.
        # Only e.g. relax mds would just show up in the runtime prioroty queue,
        # regular md tasks are known before the run starts from the config.sequence.
        task = None
        instance = None
        found_restart_task = False
        md_instance_task_counter = {k: 0 for k in md_instance_dir_counter.keys()}
        while not self.tasks.empty():
            task = self.tasks.get()
            if task.name == restart_task_name and restart_task_name == "setup":
                # this is the case if setup ends up as the restart point
                instance = "setup"
                found_restart_task = True
                break
            if task.name == "run_md":
                instance = task.kwargs["instance"]
                md_instance_task_counter[instance] += 1
                logger.debug(
                    f"{task}, {md_instance_task_counter[instance]}, {md_instance_dir_counter[instance]}"
                )
                if (
                    instance == restart_task_name
                    and md_instance_task_counter[instance]
                    == md_instance_dir_counter[instance]
                ):
                    # restart from the last completed (or half completed) valid restart point
                    found_restart_task = True
                    # put the task back in the queue
                    if restart_from_incomplete:
                        logger.info(
                            "Restarting from incomplete task. Will continue this task."
                        )
                        # append to the front of the queue
                        task.kwargs.update({"continue_md": True})
                        self.tasks.queue.appendleft(task)
                    else:
                        logger.info("Restarting after completed task.")
                    break

        if not isinstance(task, Task) or not found_restart_task or not instance:
            m = f"Could not find task {restart_task_name} in task queue. Either '{task}' is no 'Task', restart task found state is False: '{found_restart_task}' or MD task instance is None: '{instance}'. Aborting restart."
            logger.error(m)
            raise RuntimeError(m)
        if restart_from_incomplete:
            fragment = "within"
        else:
            fragment = "after"
        logger.info(
            f"Restarting from {fragment} task {task.name} with instance {instance}"
        )

        # clean up old task directories that will be overwritten
        for task_dir in task_dirs[self.iteration + 1 :]:
            shutil.rmtree(task_dir)

        # discover after it is clear which tasks are in queue
        for task_dir in task_dirs[: self.iteration + 1]:
            task_name = "_".join(task_dir.name.split(sep="_")[1:])
            task_files = TaskFiles(
                self.get_latest, {}, {}, self.config.out / task_dir.name
            )
            self._discover_output_files(taskname=task_name, files=task_files)

        # plumed fix
        for md_config in self.config.mds.__dict__.values():
            if getattr(md_config, "use_plumed"):
                plumed_out_name = get_plumed_out(self.latest_files["plumed"]).name
                plumed_out = self.get_latest(plumed_out_name)
                if plumed_out is not None:
                    self.latest_files["plumed_out"] = plumed_out
                    self.latest_files.pop(plumed_out_name)
                else:
                    logger.warning(
                        f"Plumed out file {plumed_out_name} not found. Continuing without it."
                    )

        # use latest top file
        top_path = self.get_latest("top")
        if top_path is None:
            m = "No topology file found in output directory."
            logger.error(m)
            raise FileNotFoundError(m)
        self.top = Topology(
            top=read_top(top_path, self.config.ff),
            parametrizer=self.parameterizer,
            is_reactive_predicate_f=get_is_reactive_predicate_from_config_f(
                self.config.topology.reactive
            ),
            radicals=getattr(self.config, "radicals", None),
            residuetypes_path=getattr(self.config, "residuetypes", None),
        )

        # if we restart from within an imcomplete task,
        # decrement the iteration because if will be incremented
        # when the task starts again
        if restart_from_incomplete:
            self.iteration -= 1

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
        for task_name in self.config.sequence:
            if task_name in self.config.mds.get_attributes():
                # entry is a type of MD
                md = self.task_mapping["md"]
                kwargs: dict = copy(md["kwargs"])
                kwargs.update({"instance": task_name})
                if task_name is None:
                    raise ValueError("MD task name is None")
                task = Task(
                    self,
                    f=md["f"],
                    kwargs=kwargs,
                    out=task_name,
                )
                self.tasks.put(task)

            elif task_name in self.config.reactions.get_attributes():
                # entry is a single reaction
                task_list = copy(self.task_mapping["reactions"])
                # 0 is place_reaction_tasks
                task_list[0] = copy(task_list[0])
                task_list[0]["kwargs"] = {"selected": task_name}

                for task_kwargs in task_list:
                    self.tasks.put(Task(self, **task_kwargs))
            elif task_name == "reactions":
                # check all reactions
                for task_kwargs in self.task_mapping["reactions"]:
                    self.tasks.put(Task(self, **task_kwargs))
            else:
                m = f"Unknown task encountered in the sequence: {task_name}"
                logger.error(m)
                raise ValueError(m)
        logger.info(f"Task list build:\n{pformat(list(self.tasks.queue), indent=8)}")

    def get_latest(self, suffix: str) -> Path | None:
        """Returns path to latest file of given type.

        For .dat files (in general ambiguous extensions) use full file name.
        Return None if file is not found.
        """
        logger.debug("Getting latest suffix: " + suffix)
        try:
            path = self.latest_files[suffix]
            logger.debug("Found: " + str(path))
            return path
        except Exception:
            m = f"File {suffix} requested but not found!"
            logger.warning(m)
            return None

    def __iter__(self):
        return self

    def __next__(self):
        if self.tasks.empty() and self.priority_tasks.empty():
            self.state = State.DONE
            return
        if not self.priority_tasks.empty():
            task = self.priority_tasks.get()
        else:
            task = self.tasks.get()
        if self.config.dryrun:
            logger.info(f"Pretend to run: {task.name} with args: {task.kwargs}")
            return
        files = task()
        if files is not None:
            self._discover_output_files(taskname=task.name, files=files)

    def _discover_output_files(
        self, taskname: str, files: TaskFiles
    ) -> Optional[TaskFiles]:
        """Discover further files written by a task.

        and add those files to the `files` as well as
        the file history and latest files.
        and check if double suffs are properly defined declared files by the task
        """

        if not hasattr(files, "outputdir"):
            logger.debug("No output directory found for task: " + taskname)
            return None

        discovered_files = [
            p
            for p in files.outputdir.iterdir()
            if not any(re.search(s, p.name) for s in IGNORE_SUBSTR)
        ]
        suffs = [p.suffix[1:] for p in discovered_files]
        counts = [suffs.count(s) for s in suffs]
        for suff, c, path in zip(suffs, counts, discovered_files):
            if c != 1 and suff not in AMBIGUOUS_SUFFS:
                if files.output.get(suff) is None:
                    e = (
                        "ERROR: Task produced multiple files with same suffix but "
                        "did not define with which to continue!\n"
                        f"Task {taskname}, Suffix {suff} found {c} times"
                    )
                    logger.error(e)
                    raise RuntimeError(e)

        # register discovered output files
        for path in discovered_files:
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
        self.filehist.append({files.outputdir.name: files})

        shortpaths_input = ""
        for k, v in files.input.items():
            if v is not None:
                shortpaths_input += (
                    f'  {k}: {str(v).removeprefix(str(self.config.out) + "/")}\n'
                )

        shortpaths_output = ""
        for k, v in files.output.items():
            if v is not None:
                shortpaths_output += (
                    f'  {k}: {str(v).removeprefix(str(self.config.out) + "/")}\n'
                )

        with open(self.histfile, "a") as f:
            f.write(f"Task: {files.outputdir.name}\n")
            f.write(f"Input:\n{shortpaths_input}")
            f.write(f"Output:\n{shortpaths_output}\n")

        return files

    def _setup(self, files: TaskFiles) -> TaskFiles:
        """A setup task to collect files processed by kimmdy such as the topology"""
        logger = files.logger
        logger.info("Start setup task")
        self.state = State.SETUP
        logger.info("Writing initial topology after parsing")

        if self.config.parameterize_at_setup:
            self.top.parameterization_focus_ids = set(self.top.atoms.keys())
            self.top.needs_parameterization = True
            self.top.update_parameters()

        write_top(self.top.to_dict(), files.outputdir / self.config.top.name)
        files.output["top"] = files.outputdir / self.config.top.name
        logger.info("Done with setup")

        # TODO: find a better solution for this
        #
        # copy input files that are potentially modified
        # to the setup task directory
        # by applying a recipe (e.g. by trunkate_sim_files)
        for f in ["xtc", "tpr", "trr", "plumed", "gro"]:
            if hasattr(self.config, f):
                if path := self.latest_files.get(f):
                    logger.debug(f"Copying {path} to {files.outputdir}")
                    shutil.copy(
                        path, files.outputdir / path.name, follow_symlinks=False
                    )
                    files.output[f] = files.outputdir / path.name

        return files

    def _run_md(
        self, instance: str, files: TaskFiles, continue_md: bool = False
    ) -> TaskFiles:
        """General MD simulation"""
        logger = files.logger
        logger.info(f"Start MD {instance}")
        self.state = State.MD

        md_config = self.config.mds.attr(instance)
        grompp_prefix = self.config.grompp_prefix
        gmx_alias = self.config.gromacs_alias
        mdrun_prefix = self.config.mdrun_prefix
        gmx_mdrun_flags = self.config.gmx_mdrun_flags
        top = files.input["top"]
        gro = files.input["gro"]
        files.input["mdp"] = md_config.mdp
        mdp = files.input["mdp"]
        ndx = files.input["ndx"]
        logger.debug(f"Using the following input files: top: {top}, gro: {gro}")

        outputdir = files.outputdir

        # to continue MD after timeout
        if continue_md:
            cpt = files.input["cpt"]
            logger.info(f"Restart from checkpoint file: {cpt}")
        else:
            cpt = f"{instance}.cpt"

            # running grompp again fails for pulling MD, skip it for restart because it is not necessary
            grompp_cmd = (
                f"{grompp_prefix + ' ' if grompp_prefix else ''}{gmx_alias} grompp -p {top} -c {gro} "
                f"-f {mdp} -n {ndx} -o {instance}.tpr -maxwarn 5"
            )
            # optional files for grompp:
            if self.latest_files.get("trr") is not None:
                trr = files.input["trr"]
                grompp_cmd += f" -t {trr}"
            if self.latest_files.get("edr") is not None:
                edr = files.input["edr"]
                if edr is not None and edr.exists():
                    grompp_cmd += f" -e {edr}"
                else:
                    logger.warning(f"edr file {edr} not found")
            logger.debug(f"grompp cmd: {grompp_cmd}")

        mdrun_cmd = (
            f"{mdrun_prefix + ' ' if mdrun_prefix else ''}{gmx_alias} mdrun -s {instance}.tpr -cpi {cpt} "
            f"-x {instance}.xtc -o {instance}.trr "
            f"-cpo {instance}.cpt "
            f"-c {instance}.gro -g {instance}.log -e {instance}.edr "
            f"-px {instance}_pullx.xvg -pf {instance}_pullf.xvg "
            f"-ro {instance}-rotation.xvg -ra {instance}-rotangles.log "
            f"-rs {instance}-rotslabs.log -rt {instance}-rottorque.log "
            f"{gmx_mdrun_flags}  "
        )

        if getattr(md_config, "use_plumed"):
            plumed_in = files.input["plumed"]
            if plumed_in is None:
                m = "Plumed input file not found in input files."
                logger.error(m)
                raise FileNotFoundError(m)
            mdrun_cmd += f" -plumed {plumed_in}"

            plumed_out = files.outputdir / get_plumed_out(plumed_in)
            files.output["plumed_out"] = plumed_out

        logger.debug(f"mdrun cmd: {mdrun_cmd}")
        try:
            if continue_md is False:
                run_gmx(grompp_cmd, outputdir)
            run_gmx(mdrun_cmd, outputdir)

            # specify trr to prevent rotref trr getting set as standard trr
            path = files.outputdir / f"{instance}.trr"
            if path.exists():
                files.output["trr"] = path

        except CalledProcessError as e:
            write_time_marker(files.outputdir / MARK_FAILED, "failed")
            logger.error(f"Error occured during MD {instance}:\n{e}")

            # TODO: debug
            # attempt to write out the first frame as gro
            # for easier debugging in VMD
            dump_cmd = rf"echo '0\n' | {gmx_alias} trjconv -s {instance}.tpr -f {instance}.xtc -dump 0 -o {instance}_start.gro"
            run_gmx(dump_cmd, outputdir)

            raise e

        logger.info(f"Done with MD {instance}")
        return files

    def _place_reaction_tasks(self, selected: Optional[str] = None) -> None:
        logger.info("Start query reactions")
        self.state = State.REACTION
        # empty collection for every new round of queries
        self.recipe_collections = {}

        # placing tasks in priority queue
        strategies = []
        for reaction_plugin in self.reaction_plugins:
            if selected is not None:
                if reaction_plugin.name != selected:
                    continue

            self.priority_tasks.put(
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
            self.kmc_algorithm = kmc.lower()

        else:
            # get algorithm from reaction
            if len(set(strategies)) > 1:
                raise RuntimeError(
                    "Incompatible kmc algorithms chosen in the same reaction.\n"
                    "Split the reactions in separate tasks or choose different algorithms\n"
                    "Attempted to combine:\n"
                    f"{ {rp.name:rp.config.kmc for rp in self.reaction_plugins} }"
                )
            self.kmc_algorithm = strategies[0].lower()

        logger.info(f"Queued {len(self.reaction_plugins)} reaction plugin(s)")
        return None

    def _query_reaction(
        self, reaction_plugin: ReactionPlugin, files: TaskFiles
    ) -> TaskFiles:
        logger = files.logger
        logger.info(f"Start query {reaction_plugin.name}")

        new_recipes = reaction_plugin.get_recipe_collection(files).recipes
        self.recipe_collections[reaction_plugin.name] = RecipeCollection(
            recipes=new_recipes
        )

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
            f"{sum(len(collection.recipes) for collection in self.recipe_collections.values())} total recipes available. "
            f"from {len(self.recipe_collections.keys())} reaction plugins."
        )
        kmc = self.kmc_mapping.get(self.kmc_algorithm, None)
        if kmc is None:
            m = f"Unknown KMC algorithm: {self.kmc_algorithm}"
            logger.error(m)
            raise ValueError(m)

        # FIXME: Hotfix for #355 aggregate not working for big systems
        if "rfkmc" != self.kmc_algorithm:
            for collection in self.recipe_collections.values():
                # NOTE: this is making the implicit assumption that
                # different reaction plugins produce different types
                # of reactions.
                # If this changes we can flatten before aggregating
                # WARN: This can have a subtle bug for
                # reactions with DeferredRecipeSteps,
                # if they produce non-unique reactions,
                # which would aggregate and mess up the time_start_index.
                # compared to their internal counting
                collection.aggregate_reactions()
        if "extrande" in self.kmc_algorithm:
            kmc = partial(kmc, tau_scale=self.config.tau_scale)
        elif "multi" in self.kmc_algorithm:
            logger.debug(
                f"Setting {self.kmc_algorithm} up to pick "
                f"{self.config.multi_kmc} reactions."
            )
            kmc = partial(kmc, n=self.config.multi_kmc)

        self.kmcresult = kmc(flatten_recipe_collections(self.recipe_collections))
        if not isinstance(self.kmcresult, KMCAccept):
            # rejection or error, nothing to be done
            if isinstance(self.kmcresult, KMCReject):
                logger.info(f"Rejected recipe: {self.kmcresult.reason}")
            elif isinstance(self.kmcresult, Exception):
                logger.error(f"Error during KMC: {self.kmcresult}")
            return files

        # Correct the offset of time_start_index by the concatenation
        # of all rates of all recipes from all reactions back onto the offset
        # within the one chosen reaction plugin
        n_recipes_per_plugin = [
            sum([len(r.rates) for r in v.recipes])
            for v in self.recipe_collections.values()
        ]
        logger.info(f"Plugins: {[k for k in self.recipe_collections.keys()]}")
        logger.info(f"Number of recipes per reaction plugin: {n_recipes_per_plugin}")

        self.kmcresult.time_start_index_within_plugin = (
            total_index_to_index_within_plugin(
                self.kmcresult.time_start_index, n_recipes_per_plugin
            )
        )

        recipe = self.kmcresult.recipe
        if self.config.save_recipes:
            flatten_recipe_collections(self.recipe_collections).to_csv(
                files.outputdir / "recipes.csv", recipe
            )

        self.time += self.kmcresult.time_delta
        if len(recipe.rates) == 0:
            logger.info("No reaction selected")
        else:
            logger.info(
                f"Reaction jumps ahead in time: {self.kmcresult.time_delta*10e-12:.4e} s. Overall new time {self.time*10e-12:.4e} s"
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

        logger.info("Done with Decide recipe.")
        return files

    def _apply_recipe(self, files: TaskFiles) -> TaskFiles:
        logger = files.logger

        if self.kmcresult is None:
            m = "Attempting to _apply_recipe without having chosen one with _decide_recipe."
            logger.error(m)
            write_time_marker(files.outputdir / MARK_FAILED, "failed")
            raise RuntimeError(m)
        elif isinstance(self.kmcresult, KMCReject):
            m = f"No reaction has been accepted (KMCRejection {self.kmcresult})."
            logger.info(m)
            return files
        elif isinstance(self.kmcresult, KMCError):
            m = f"No reaction has been accepted (KMCError {self.kmcresult})."
            logger.warning(m)
            return files

        recipe = self.kmcresult.recipe
        logger.info(f"Start Recipe in KIMMDY iteration {self.iteration}")
        logger.info(f"Recipe: {recipe.get_recipe_name()}")

        # Set time to chosen 'time_start' of KMCResult
        ttime = self.kmcresult.time_start
        plugin_time_index = self.kmcresult.time_start_index_within_plugin

        shadow_files_binding = None

        logger.info(f"Chosen time_start: {ttime} ps")
        logger.info(f"Time index within plugin: {plugin_time_index}")
        if plugin_time_index is None:
            m = f"Time index within plugin is None, this should not happen. Was _apply_recipe called before _decide_recipe?"
            logger.error(m)
            raise RuntimeError(m)

        if isinstance(recipe.recipe_steps, list):
            recipe.recipe_steps = recipe.recipe_steps
        elif isinstance(recipe.recipe_steps, DeferredRecipeSteps):
            logger.info(
                f"Steps of recipe where deferred, calling callback with key {recipe.recipe_steps.key} and time_index {plugin_time_index}"
            )
            recipe.recipe_steps = recipe.recipe_steps.callback(
                recipe.recipe_steps.key, plugin_time_index, ttime
            )
            logger.info(
                f"Got {len(recipe.recipe_steps)} steps in recipe {recipe.get_recipe_name()}"
            )
        else:
            m = f"Recipe steps of {recipe} are neither a list nor a DeferredRecipeSteps object."
            logger.error(m)
            raise ValueError(m)
        if any([isinstance(step, Place) for step in recipe.recipe_steps]):
            # only first time of interval is valid for placement
            ttime = recipe.timespans[0][0]

        # get vmd selection (after deferred steps are resolved)
        vmd_selection = recipe.get_vmd_selection()
        logger.info(f"VMD selection: {vmd_selection}")
        with open(files.outputdir / "vmd_selection.txt", "w") as f:
            f.write(vmd_selection)

        # write time marker for reaction time
        # in the current task dir (<n>_apply_recipe)
        # but also in the output dir of the MD task
        # onto which the reaction is applied
        write_reaction_time_marker(dir=files.outputdir, time=ttime)
        gro = files.input["gro"]
        if gro is None:
            m = "No gro file found from the previous md run."
            logger.error(m)
        else:
            write_reaction_time_marker(dir=gro.parent, time=ttime)

        # because the gro_reaction file is written to files.output
        # it will be discovered by _discover_output_files
        # and set as the latest gro file for the next tasks
        # but this only happens after the apply_recipe task
        # so we need to set it manually here for intermediate tasks
        # like Relax and Place to have the correct coordinates
        logger.info(f"Writing coordinates (gro and trr) and energy (edr) for reaction.")
        if ttime is not None:
            write_coordinate_files_at_reaction_time(files=files, time=ttime)
            self.latest_files["gro"] = files.output["gro"]
            self.latest_files["trr"] = files.output["trr"]
            self.latest_files["edr"] = files.output["edr"]

        top_initial = deepcopy(self.top)
        for step in recipe.recipe_steps:
            if isinstance(step, Break):
                self.top.break_bond((step.atom_id_1, step.atom_id_2))
                if hasattr(self.config, "plumed"):
                    break_bond_plumed(
                        files,
                        (step.atom_id_1, step.atom_id_2),
                        files.outputdir / self.config.plumed.name.replace(".", "_mod."),
                    )
            elif isinstance(step, Bind):
                self.top.bind_bond((step.atom_id_1, step.atom_id_2))
            elif isinstance(step, Place):
                relax_task = Task(
                    self,
                    f=place_atom,
                    kwargs={"step": step, "ttime": None},
                    out="place_atom",
                )
                place_files = relax_task()
                if place_files is not None:
                    self._discover_output_files(
                        taskname=relax_task.name, files=place_files
                    )
                    shadow_files_binding = place_files
                if step.id_to_place is not None:
                    self.top.parameterization_focus_ids.update([step.id_to_place])
            elif isinstance(step, Relax):
                logger.info("Starting relaxation md as part of reaction..")
                if not hasattr(self.config.changer.coordinates, "md"):
                    logger.warning("Relax task requested but no MD specified for it!")
                    continue

                if self.config.changer.coordinates.slow_growth:
                    # Create a temporary slow growth topology for sub-task run_md, afterwards, top will be reset properly.

                    self.top.update_parameters()

                    # write out original topology (topA) and target (topB) for easier debugging
                    write_top(
                        top_initial.to_dict(),
                        files.outputdir
                        / self.config.top.name.replace(".top", "_before.top"),
                    )
                    write_top(
                        self.top.to_dict(),
                        files.outputdir
                        / self.config.top.name.replace(".top", "_after.top"),
                    )

                    # top_initial is still the topology before the reaction
                    # we need to do some (temporary) changes to it to stabilize
                    # the slow_growth
                    # First we find out if there are solvent atoms among the involved atoms
                    solvent_atoms: set[str] = set()
                    logger.debug(f"Checking for reacting solvent residues..")
                    for ai in self.top.parameterization_focus_ids:
                        if top_initial.atoms[ai].residue == "SOL":
                            solvent_atoms.add(ai)
                            logger.debug(
                                f"Reacting solvent atom: {top_initial.atoms[ai]}"
                            )
                    if len(solvent_atoms) > 0:
                        logger.info(
                            f"{len(solvent_atoms)} solvent atoms are involved "
                            "in the reaction, they will get tempoary bonds for "
                            "the start of the slow growth simulation."
                        )
                        ow = None
                        hw1 = None
                        hw2 = None
                        for ai in solvent_atoms:
                            a = top_initial.atoms[ai]
                            if a.atom == "OW":
                                ow = a
                            if a.atom == "HW1":
                                hw1 = a
                            if a.atom == "HW2":
                                hw2 = a
                        if ow is not None and hw1 is not None and hw2 is not None:
                            logger.info(
                                "Found one complete water molecule that takes part in the reaction."
                            )
                            top_initial.bind_bond((ow.nr, hw1.nr))
                            top_initial.bind_bond((ow.nr, hw2.nr))
                            b1 = top_initial.bonds.get((ow.nr, hw1.nr))
                            b2 = top_initial.bonds.get((ow.nr, hw2.nr))
                            logger.info(f"Added bonds: {b1}, {b2}")

                    write_top(
                        top_initial.to_dict(),
                        files.outputdir
                        / self.config.top.name.replace(
                            ".top", "_before_with_solvent_bonds.top"
                        ),
                    )

                    top_merge = merge_top_slow_growth(
                        # top_a was copied before parameters are updated
                        top_a=top_initial,
                        # top_b is parameterized for after the reaction
                        # top_b is modified and returned as the merged top
                        # hence it must be copied here to not modify self.top
                        top_b=deepcopy(self.top),
                        morse_only=self.config.changer.coordinates.slow_growth
                        == "morse_only",
                    )
                    top_merge_path = files.outputdir / self.config.top.name.replace(
                        ".top", "_relax.top"
                    )
                    write_top(top_merge.to_dict(), top_merge_path)
                    # declare this temporary topolgy the latest topology
                    # such that it will be used for the subsequent relax md task
                    self.latest_files["top"] = top_merge_path
                md_instance = self.config.changer.coordinates.md
                relax_task = Task(
                    self,
                    f=self._run_md,
                    kwargs={"instance": md_instance},
                    out=md_instance,
                )
                relax_task_files = relax_task()
                if relax_task_files is not None:
                    self._discover_output_files(relax_task.name, relax_task_files)
                    shadow_files_binding = relax_task_files

            elif isinstance(step, CustomTopMod):
                step.f(self.top)

        logger.info(f"Updating partial charges")
        self.top.update_partial_charges(recipe.recipe_steps)
        self.top.update_parameters()

        # this is the new topology after the reaction
        # not the temporay topology <name>_relax.top for slow_growth
        write_top(self.top.to_dict(), files.outputdir / self.config.top.name)
        files.output["top"] = files.outputdir / self.config.top.name

        # Recipe done, reset runmanger state
        self.kmcresult = None

        if shadow_files_binding is not None:
            # if a relaxation or placement task was run,
            # we overwrite the coordinate output files of
            # the files object with the files from the relaxation or placement task
            # (whichever was later)
            # such that the next task will use these files
            for ext in ["gro", "trr", "xtc", "edr"]:
                out = shadow_files_binding.output.get(ext)
                if out is not None:
                    files.output[ext] = out
            # but not the `top`, because the top for the relaxation
            # is only temporary and should not be used for the next task

        logger.info("Done with Apply recipe")
        return files

    def __repr__(self):
        return "Runmanager"
