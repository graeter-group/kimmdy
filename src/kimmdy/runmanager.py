from __future__ import annotations
import logging
import queue
from pathlib import Path
from enum import Enum, auto
from typing import Callable
from kimmdy.config import Config
from kimmdy.reactions.homolysis import Homolysis
from kimmdy.reaction import ReactionResult, ConversionRecipe, ConversionType
import kimmdy.mdmanager as md
import kimmdy.changemanager as changer
from pprint import pformat
import random


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


class Task:
    """A task to be performed as as a step in the RunManager.
    consists of a function and it's keyword arguments and is
    itself callable.
    """

    def __init__(self, f, kwargs={}):
        self.f = f
        self.kwargs = kwargs
        self.name = self.f.__name__

    def __call__(self):
        return self.f(**self.kwargs)

    def __repr__(self) -> str:
        return str(self.f) + " args: " + str(self.kwargs)


class RunManager:
    """The RunManager, a central piece.
    Manages the queue of tasks, communicates with the
    rest of the program and keeps track of global state.
    """

    reaction_results: list[ReactionResult]

    def __init__(self, input_file: Path):
        self.config = Config(Path(input_file))
        self.tasks = queue.Queue()
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
        self.filehist: list[dict[str, dict[str, Path]]] = []
        self.ambiguos_suffs = ["dat", "xvg", "log", "trr"]
        # index initilial files
        self.filehist.append({"in": dict(), "out": dict()})
        self.filehist[-1]["in"]["top"] = self.config.top
        self.filehist[-1]["in"]["gro"] = self.config.gro
        self.filehist[-1]["in"]["idx"] = self.config.idx
        self.filehist[-1]["in"]["plumed_dat"] = self.config.plumed.dat
        self.filehist[-1]["in"]["distances_dat"] = self.config.plumed.distances

    def get_latest(self, f_type: str) -> Path:
        """Returns path to latest file of given type, raises an error if the file is not found.
        For .dat files (in general ambiuos extensions) use full file name.
        """
        f_type = f_type.replace(".", "_")
        for step in self.filehist[::-1]:
            for io in (step["out"], step["in"]):
                if f_type in io.keys():
                    # if f_type in self.ambiguos_suffs: # suffix was ambiguos
                    #     logging.warn(f"{f_type} ambiguos! Please specify full file name!")
                    return io[f_type]
                # elif f_type[-3:] in io.keys():
                #     logging.debug(f"Full file name {f_type} given, but only suffix found.")
                #     return io[f_type[-3:]]
        else:
            m = f"File {f_type} requested but not found!"
            logging.error(m)
            raise FileNotFoundError(m)

    def run(self):
        logging.info("Start run")
        # TODO: Make task order dynamic: iterate over seq in config
        # golbal step counter for dirs
        # LATER: counter for groups as defined initially

        self.tasks.put(Task(self._run_md_equil))

        for _ in range(self.iterations):
            self.tasks.put(Task(self._run_md_prod))
            self.tasks.put(Task(self._query_reactions))
            self.tasks.put(Task(self._decide_reaction))
            self.tasks.put(Task(self._run_recipe))

            # self.tasks.put(Task(self._dummy, {}))
            # self.tasks.put(Task(self._run_md_minim, {}))

        while (self.state is not State.DONE) or (self.iteration < self.iterations):
            next(self)

        logging.info(
            f"Stop running tasks, state: {self.state}, iteration:{self.iteration}, max:{self.iterations}"
        )
        logging.info(f"History:\n{pformat(self.filehist)}")

    def __iter__(self):
        return self

    def __next__(self):
        if self.tasks.empty():
            self.state = State.DONE
            return
        task = self.tasks.get()
        if self.config.dryrun:
            logging.info(f"Pretending to run: {task.name} with args: {task.kwargs}")
            return
        in_d, out_dir = task()

        out_d = {}
        if out_dir is not None:
            for p in out_dir.iterdir():
                p_k = p.suffix[1:]  # strip dot
                if p_k in self.ambiguos_suffs:
                    p_k = p.name.replace(".", "_")
                out_d[p_k] = p

        if in_d is None:
            in_d = {}
        self.filehist.append({"in": in_d, "out": out_d})
        self.iteration += 1

    def _dummy(self):
        logging.info("Start dummy task")
        out_dir = self.config.out / f"dummy_{self.iteration}"
        out_dir.mkdir()
        in_d = {
            "top": self.get_latest("top"),
            "gro": self.get_latest("gro"),
        }
        # perform step
        md.dummy_step(out_dir, **in_d)

        return in_d, out_dir

    def _run_md_equil(self):
        logging.info("Start equilibration MD")
        self.state = State.MD
        out_dir = self.config.out / f"equil_{self.iteration}"
        out_dir.mkdir()
        (out_dir / self.config.ff.name).symlink_to(self.config.ff)
        in_d = {
            "top": self.get_latest("top"),
            "gro": self.get_latest("gro"),
            "mdp": self.config.equilibrium.mdp,
            "idx": self.config.idx,
        }
        # perform step
        md.equilibrium(out_dir, **in_d)
        logging.info("Done equilibrating")
        return in_d, out_dir

    def _run_md_minim(self):
        logging.info("Start minimization md")
        self.state = State.MD
        out_dir = self.config.out / f"min_{self.iteration}"
        out_dir.mkdir()
        (out_dir / self.config.ff.name).symlink_to(self.config.ff)
        in_d = {
            "top": self.get_latest("top"),
            "gro": self.get_latest("gro"),
            "mdp": self.config.minimization.mdp,
        }
        # perform step
        md.minimzation(out_dir, **in_d)
        logging.info("Done minimizing")
        return in_d, out_dir

    def _run_md_eq(self):
        # TODO: combine w/ other equilibration?
        logging.info("Start _run_md_eq MD")
        self.state = State.MD
        out_dir = self.config.out / f"equil_{self.iteration}"
        out_dir.mkdir()
        (out_dir / self.config.ff.name).symlink_to(self.config.ff)
        in_d = {
            "top": self.get_latest("top"),
            "mdp": self.config.equilibrium.mdp,
            "gro": self.get_latest("gro"),
        }
        # perform step
        md.equilibration(out_dir, **in_d)
        logging.info("Done equilibrating")
        return in_d, out_dir

    def _run_md_prod(self):
        logging.info("Start production MD")
        self.state = State.MD
        out_dir = self.config.out / f"prod_{self.iteration}"
        out_dir.mkdir()
        (out_dir / self.config.ff.name).symlink_to(self.config.ff)
        in_d = {
            "top": self.get_latest("top"),
            "gro": self.get_latest("gro"),
            "mdp": self.config.prod.mdp,
            "idx": self.config.idx,
            "cpt": self.get_latest("cpt"),
            "plumed_dat": self.get_latest("plumed_dat"),
        }
        # perform step
        md.production(out_dir, **in_d)
        logging.info("Done minimizing")
        return in_d, out_dir

    def _run_md_relax(self):
        logging.info("Start relaxation MD")
        self.state = State.MD
        out_dir = self.config.out / f"prod_{self.iteration}"
        out_dir.mkdir()
        (out_dir / self.config.ff.name).symlink_to(self.config.ff)
        in_d = {
            "top": self.get_latest("top"),
            "gro": self.get_latest("gro"),
            "mdp": self.config.changer.coordinates.md.mdp,
            "idx": self.config.idx,
            "cpt": self.get_latest("cpt"),
        }
        # perform step
        md.relaxation(out_dir, **in_d)
        logging.info("Done simulating")
        return in_d, out_dir

    def _query_reactions(self):
        logging.info("Query reactions")
        self.state = State.REACTION
        in_d = None
        out_dir = None

        # TODO: this could be more elegant
        if hasattr(self.config.reactions, "homolysis"):
            self.reaction = Homolysis()

            in_d = {
                "plumed_dat": self.get_latest("plumed.dat"),
                "distances_dat": self.get_latest("distances.dat"),
                "top": self.get_latest("top"),
                "ffbonded_itp": self.config.reactions.homolysis.bonds,
                "edissoc_dat": self.config.reactions.homolysis.edis,
            }
            self.reaction_results.append(self.reaction.get_reaction_result(**in_d))

        logging.info("Rates and Recipes:")
        # TODO: this would fail with a out of range error,
        # we need a way of writing out results
        # logging.info(self.rates[0:10])
        # logging.info(self.recipe.type)
        # logging.info(self.recipe.atom_idx[0:10])
        # logging.info("Reaction done")
        return in_d, out_dir

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
        return None, None

    def _run_recipe(self):
        logging.info(f"Start Recipe in step {self.iteration}")
        logging.info(f"Breakpair: {self.chosen_recipe.atom_idx}")
        out_dir = self.config.out / f"recipe_{self.iteration}"
        out_dir.mkdir()
        (out_dir / self.config.ff.name).symlink_to(self.config.ff)

        in_d = {
            "top": self.get_latest("top"),
            "plumed_dat": self.get_latest("plumed.dat"),
        }

        newtop = out_dir / "topol_mod.top"
        newplumeddat = out_dir / "plumed.dat"
        newplumeddist = out_dir / "distances.dat"
        changer.modify_top(self.chosen_recipe, in_d["top"], newtop)
        logging.info(f"Wrote new topology to {newtop.parts[-3:]}")
        changer.modify_plumed(
            self.chosen_recipe, in_d["plumed_dat"], newplumeddat, newplumeddist
        )
        logging.info(f"Wrote new plumedfile to {newplumeddat.parts[-3:]}")
        logging.info(f"Looking for md in {self.config.changer.coordinates.__dict__}")
        # TODO: clean this up, maybe make function for this in config
        if hasattr(self.config, "changer"):
            if hasattr(self.config.changer, "coordinates"):
                if hasattr(self.config.changer.coordinates, "md"):
                    self.tasks.put(Task(self._run_md_relax))
        return in_d, out_dir
