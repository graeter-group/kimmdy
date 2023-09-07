"""
The Runmanager is the main entry point of the program.

It manages the queue of tasks, communicates with the
rest of the program and keeps track of global state.
"""
from __future__ import annotations
import logging
from pathlib import Path
from copy import copy, deepcopy
import dill
import queue
from enum import Enum, auto
from typing import Callable, Union
from kimmdy.config import Config
from kimmdy.parsing import read_top, write_json, read_plumed, write_top
from kimmdy.recipe import RecipeCollection, Recipe, Break, Bind, Place, Relax
from kimmdy.plugins import ReactionPlugin, BasicParameterizer
from kimmdy.coordinates import place_atom, break_bond_plumed, merge_top_slow_growth
from kimmdy.tasks import Task, TaskFiles, TaskMapping
from kimmdy.utils import run_gmx
from pprint import pformat
from kimmdy import reaction_plugins, parameterization_plugins
from kimmdy.topology.topology import Topology
import time
from kimmdy.kmc import rf_kmc, KMCResult

logger = logging.getLogger(__name__)

# file types of which there will be multiple files per type
AMBIGUOUS_SUFFS = ["dat", "xvg", "log", "itp", "mdp"]
# are there cases where we have multiple trr files?


class State(Enum):
    """State of the system.
    one of IDLE, MD, REACTION, DONE.
    """

    IDLE = auto()
    MD = auto()
    REACTION = auto()
    DONE = auto()


def get_plumed_out(plumed: Path) -> str:
    plumed_parsed = read_plumed(plumed)
    plumed_out = None
    for part in plumed_parsed["prints"]:
        if file := part.get("FILE"):
            if not plumed_out:
                plumed_out = file
                logger.debug(f"Found plumed print FILE {plumed_out}")
            else:
                raise RuntimeError("Multiple FILE sections found in plumed dat")
    return plumed_out


def get_existing_files(config: Config, section: str = "root") -> dict:
    """Initialize latest_files with every existing file defined in config"""
    file_d = {}
    attr_names = filter(lambda s: s[0] != "_", config.__dir__())
    for attr_name in attr_names:
        attr = getattr(config, attr_name)
        if isinstance(attr, Path):
            if not attr.exists():
                # File exists checks are in the config
                continue
            key = attr.suffix[1:]  # rm the get_existing_files

            # AMBIGUOUS_SUFFS -> key whole name
            if attr.suffix[1:] in AMBIGUOUS_SUFFS:
                key = attr.name

            if attr_name == "plumed":
                key = "plumed"
                file_d["plumed_out"] = attr.parent / get_plumed_out(attr)

            file_d[key] = attr
        elif isinstance(attr, Config):
            file_d.update(get_existing_files(attr, attr_name))
    return file_d


class RunManager:
    """The Runmanager is the main entry point of the program.

    Manages the queue of tasks, communicates with the
    rest of the program and keeps track of global state.

    Attributes
    ----------
    config :
        The configuration object.
    from_checkpoint :
        Whether the runmanager was initialized from a checkpoint.
    tasks :
        Tasks from config.
    crr_tasks :
        Current tasks.
    iteration :
        Current iteration.
    state :
        Current state of the system.
    recipe_collection :
        Collection of recipes.
    latest_files :
        Dictionary of latest files.
    histfile :
        Path to history file.
    cptfile :
        Path to checkpoint file.
    top :
        Topology object.
    filehist :
        List of dictionaries of TaskFiles.
    task_mapping :
        Mapping of task names to runmanager methods.
    """

    def __init__(self, config: Config):
        self.config: Config = config
        self.from_checkpoint: bool = False
        self.tasks: queue.Queue[Task] = queue.Queue()  # tasks from config
        self.crr_tasks: queue.Queue[Task] = queue.Queue()  # current tasks
        self.iteration: int = 0
        self.state: State = State.IDLE
        self.recipe_collection: RecipeCollection = RecipeCollection([])
        self.recipe: Union[Recipe, None] = None
        self.time: float = 0.0  # [ps]
        self.latest_files: dict[str, Path] = get_existing_files(config)
        logger.debug("Initialized latest files:")
        logger.debug(pformat(self.latest_files))
        self.histfile: Path = self.config.out / "kimmdy.history"
        self.cptfile: Path = self.config.out / "kimmdy.cpt"

        self.top = Topology(read_top(self.config.top, self.config.ff))
        try:
            if self.config.changer.topology.parameterization == "basic":
                self.parameterizer = BasicParameterizer()
            else:
                self.parameterizer = parameterization_plugins[
                    self.config.changer.topology.parameterization
                ]()
        except KeyError as e:
            raise KeyError(
                f"The parameterization tool chosen in the configuration file: "
                f"'{self.config.changer.topology.parameterization}' can not be found in "
                f"the parameterization plugins: {list(parameterization_plugins.keys())}"
            ) from e

        self.filehist: list[dict[str, TaskFiles]] = [
            {"setup": TaskFiles(self.get_latest)}
        ]

        self.task_mapping: TaskMapping = {
            "md": self._run_md,
            "reactions": [
                {"f": self._place_reaction_tasks},
                {"f": self._decide_recipe, "kwargs": {"decision_strategy": rf_kmc}},
                {"f": self._apply_recipe, "out": "apply_recipe"},
            ],
        }

        # decide_recipe only needs a outputdir some times:
        if self.config.plot_rates or self.config.save_recipes:
            self.task_mapping["reactions"][1]["out"] = "decide_recipe"

        # Instantiate reactions
        self.reaction_plugins: list[ReactionPlugin] = []
        react_names = self.config.reactions.get_attributes()
        logger.debug(f"Instantiating Reactions: {react_names}")
        for rp_name in react_names:
            r = reaction_plugins[rp_name]
            reaction_plugin: ReactionPlugin = r(rp_name, self)
            self.reaction_plugins.append(reaction_plugin)

    def run(self):
        logger.info("Start run")
        self.start_time = time.time()
        self.current_time = time.time()

        if self.from_checkpoint:
            print(logger)
            logger.info(f"KIMMDY is starting from a checkpoint.")
        else:
            self._setup_tasks()

        while (
            self.state is not State.DONE
            and (self.iteration <= self.config.max_tasks or self.config.max_tasks == 0)
            and (
                (self.current_time - self.start_time) / 3600 < self.config.max_hours
                or self.config.max_hours == 0
            )
        ):
            logger.info("Writing checkpoint before next task")
            with open(self.cptfile, "wb") as f:
                dill.dump(self, f)
            next(self)
            self.current_time = time.time()
            logger.info("Done with:")
            logger.info(f"task: {self.iteration}, max: {self.config.max_tasks}")
            logger.info(
                f"hours: {(self.current_time - self.start_time) / 3600}, max: {self.config.max_hours}"
            )

        logger.info(f"Finished running tasks, state: {self.state}")
        logger.info(f"KIMMDY is done.")

    def _setup_tasks(self):
        """Populates the tasks queue.
        Allows for mapping one config entry to multiple tasks
        """
        logger.info("Build task list")
        for entry in self.config.sequence:
            # handle md steps
            if hasattr(self.config, "mds"):
                if entry in self.config.mds.get_attributes():
                    task = Task(
                        self,
                        f=self.task_mapping["md"],
                        kwargs={"instance": entry},
                        out=entry,
                    )
                    self.tasks.put(task)
                    continue

            # handle single reaction steps
            if entry in self.config.reactions.get_attributes():
                task_list = copy(self.task_mapping["reactions"])
                task_list[0] = {
                    "f": self._place_reaction_tasks,
                    "kwargs": {"selected": entry},
                }

                for task_kargs in task_list:
                    self.tasks.put(Task(self, **task_kargs))
                continue

            # handle combined reaction steps
            if entry == "reactions":
                for task_kargs in self.task_mapping[entry]:
                    self.tasks.put(Task(self, **task_kargs))
                continue
            m = f"Unknown task encounterd: {entry}"
            logger.error(m)
            raise ValueError(m)

    def write_one_checkoint(self):
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
        logger.info("Starting task: " + pformat(task))
        files = task()
        self._discover_output_files(task.name, files)

    def _discover_output_files(self, taskname, files: TaskFiles):
        """Discover further files written by a task.

        and add those files to the `files` as well as
        the file history and latest files.
        """
        # discover other files written by the task
        if hasattr(files, "outputdir"):
            # check whether double suffs are properly defined in files by the task
            suffs = [p.suffix[1:] for p in files.outputdir.iterdir()]
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

            logger.debug("Update latest files with: ")
            logger.debug(pformat(files.output))
            self.latest_files.update(files.output)
            logger.debug("Append to file history")
            self.filehist.append({taskname: files})

            logger.info(f"Current task files:")
            m = f"""
            Task: {taskname} with output directory: {files.outputdir}
            Task: {taskname}, input:\n{pformat(files.input)}
            Task: {taskname}, output:\n{pformat(files.output)}
            """
            logger.info(m)
            with open(self.histfile, "a") as f:
                f.write(m)

    def _run_md(self, instance, files) -> TaskFiles:
        """General MD simulation"""
        logger = files.logger
        logger.info(f"Start MD {instance}")
        self.state = State.MD

        md_config = self.config.mds.attr(instance)
        gmx_alias = self.config.gromacs_alias
        gmx_mdrun_flags = self.config.gmx_mdrun_flags

        logger.debug(self.latest_files)
        top = files.input["top"]
        gro = files.input["gro"]
        mdp = md_config.mdp
        files.input["mdp"] = mdp
        ndx = files.input["ndx"]

        outputdir = files.outputdir

        grompp_cmd = (
            f"{gmx_alias} grompp -p {top} -c {gro} "
            f"-f {mdp} -n {ndx} -o {instance}.tpr -maxwarn 5"
        )
        # only appends these lines if there are trr and edr files
        try:
            trr = files.input["trr"]
            edr = files.input["edr"]
            grompp_cmd += f" -t {trr} -e {edr}"
        except FileNotFoundError:
            pass

        mdrun_cmd = (
            f"{gmx_alias} mdrun -s {instance}.tpr -cpi {instance}.cpt "
            f"-x {instance}.xtc -o {instance}.trr -cpo {instance}.cpt "
            f"-c {instance}.gro -g {instance}.log -e {instance}.edr "
            f"-px {instance}_pullx.xvg -pf {instance}_pullf.xvg "
            f"-ro {instance}-rotation.xvg -ra {instance}-rotangles.log "
            f"-rs {instance}-rotslabs.log -rt {instance}-rottorque.log "
            f"{gmx_mdrun_flags}  "
        )
        # -ntomp {ntomp} removed for now
        # like this, the previous checkpoint file would not be used,
        # -t and -e options from grompp
        # replace the checkpoint file if gen_vel = no in the mdp file

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

    def _place_reaction_tasks(self, selected: Union[str, None] = None):
        logger.info("Query reactions")
        self.state = State.REACTION
        # empty list for every new round of queries
        self.recipe_collection: RecipeCollection = RecipeCollection([])

        # placing tasks in priority queue
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

        logger.info(f"Queued {len(self.reaction_plugins)} reaction plugin(s)")
        return

    def _query_reaction(self, reaction_plugin, files):
        logger = files.logger
        logger.info(f"Start query {reaction_plugin.name}")

        self.recipe_collection.recipes.extend(
            reaction_plugin.get_recipe_collection(files).recipes
        )
        self.recipe_collection.aggregate_reactions()

        logger.info(f"Recipes recived from {reaction_plugin.name}")
        return files

    def _decide_recipe(
        self, decision_strategy: Callable[[RecipeCollection], KMCResult], files=None
    ):
        if files is None:
            logger = logging.getLogger(__name__)
        else:
            logger = files.logger

        logger.info(
            f"Decide on a recipe from {len(self.recipe_collection.recipes)} available"
        )
        decision_d = decision_strategy(self.recipe_collection)
        self.recipe = decision_d.recipe

        if self.config.save_recipes:
            self.recipe_collection.to_csv(files.outputdir / "recipes.csv", self.recipe)

        try:
            if self.config.plot_rates:
                kwargs = {
                    "outfile": files.outputdir / "reaction_rates.svg",
                    "highlight_r": self.recipe,
                }
                if (decision_d.time_start is not None) and (decision_d.time_start != 0):
                    kwargs["highlight_t"] = decision_d.time_start
                self.recipe_collection.plot(**kwargs)
        except Exception as e:
            logger.warning(f"Error occured during plotting:\n{e}")

        if decision_d.time_delta:
            self.time += decision_d.time_delta
        logger.info(
            f"Chosen recipe is: {self.recipe.get_recipe_name()} at time {self.time}"
        )
        return

    def _apply_recipe(self, files: TaskFiles) -> TaskFiles:
        logger = files.logger
        logger.info(f"Start Recipe in KIMMDY iteration {self.iteration}")
        logger.info(f"Recipe: {self.recipe.get_recipe_name()}")

        logger.debug(f"Chose recipe steps: {self.recipe.recipe_steps}")

        top_initial = deepcopy(self.top)
        for step in self.recipe.recipe_steps:
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
                task = Task(
                    self,
                    f=place_atom,
                    kwargs={"step": step, "timespan": self.recipe.timespans},
                    out="place_atom",
                )
                place_files = task()
                self._discover_output_files(task.name, place_files)

            elif isinstance(step, Relax):
                logger.info("Starting relaxation md as part of reaction..")
                if not hasattr(self.config.changer.coordinates, "md"):
                    logger.warning("Relax task requested but no MD specified for it!")
                    continue
                self.parameterizer.parameterize_topology(
                    self.top
                )  # not optimal that this is here

                if self.config.changer.coordinates.slow_growth:
                    # Create a slow growth topology for sub-task run_md, afterwards, top will be set back properly
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
                md_files = task()
                self._discover_output_files(task.name, md_files)

        # parameterize topology
        self.parameterizer.parameterize_topology(self.top)
        write_top(self.top.to_dict(), files.outputdir / self.config.top.name)
        files.output["top"] = files.outputdir / self.config.top.name

        # capture state of radicals
        write_json(
            {"time": self.time, "radicals": list(self.top.radicals.keys())},
            files.outputdir / "radicals.json",
        )

        # Recipe done, reset runmanger state
        self.recipe = None

        logger.info("Reaction done")
        return files
