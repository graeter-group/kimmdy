from __future__ import annotations
import logging
import queue
from pathlib import Path
from enum import Enum, auto
from dataclasses import dataclass
from typing import Callable
from kimmdy.utils import run_shell_cmd, identify_atomtypes
from kimmdy.config import Config
from kimmdy.reaction import ConversionType, ConversionRecipe
from kimmdy.reactions.homolysis import Homolysis
from kimmdy.mdmanager import MDManager
from kimmdy.changemanager import ChangeManager
from pprint import pformat

import os
import numpy as np
import random


def default_decision_strategy(list_of_rates, list_of_recipes):
    # compare e.g. https://en.wikipedia.org/wiki/Kinetic_Monte_Carlo#Rejection-free_KMC

    total_rate = sum(list_of_rates)  # sum all rates
    random.seed()
    t = random.random()  # t in [0.0,1.0)
    rate_running_sum = 0
    for i in range(len(list_of_rates)):
        rate_running_sum += list_of_rates[i]
        if (t * total_rate) <= rate_running_sum:
            idx = i
            # breakpair = list_of_recipes[i][1]
            # atomtypes = list_of_recipes[i][2]
            # rate = list_of_rates[i]
            break
    u = random.random()
    delta_t = np.log(1 / u) / total_rate

    return list_of_recipes[idx]

    # print ('Rate = ' + str(rate) + ', breakpair = ' + str (breakpair) + ', atomtypes = ' + str(atomtypes) + 'jump [ps] = ' + str(delta_t))
    # logging.info(('Rate = ' + str(rate) + ', breakpair = ' + str (breakpair) + ', atomtypes = ' + str(atomtypes) + 'jump [ps] = ' + str(delta_t)))


class State(Enum):
    IDLE = auto()
    MD = auto()
    REACTION = auto()
    DONE = auto()


class Task:
    def __init__(self, f, kwargs):
        self.f = f
        self.kwargs = kwargs
        self.name = self.f.__name__

    def __call__(self):
        return self.f(**self.kwargs)

    def __repr__(self) -> str:
        return str(self.f) + " args: " + str(self.kwargs)


class RunManager:
    def __init__(self, input_file: Path):
        self.config = Config(Path(input_file))
        self.tasks = queue.Queue()
        self.iteration = 0
        self.iterations = self.config.iterations
        self.state = State.IDLE
        self.rates = []
        self.md = MDManager()
        # self.measurements = Path("measurements")
        # self.structure = self.config.gro
        # self.top = self.config.top
        # self.trj = self.config.cwd / ("prod_" + str(self.iteration) + ".trj")
        # self.plumeddat = self.config.plumed.dat
        # self.plumeddist = self.config.plumed.distances

        self.filehist: list[dict[str, dict[str, Path]]] = []

        # index initilial files
        self.filehist.append({"in": dict(), "out": dict()})
        self.filehist[-1]["in"]["top"] = self.config.top
        self.filehist[-1]["in"]["gro"] = self.config.gro
        self.filehist[-1]["in"]["idx"] = self.config.idx

    def get_latest(self, f_type: str):
        """Returns path to latest file of given type, or None if not found."""
        for step in self.filehist[::-1]:
            for io in (step["out"], step["in"]):
                if f_type in io.keys():
                    if f_type + "#" in io.keys():
                        logging.warn(f"Multiple files {f_type} found!")
                    return io[f_type]
        else:
            logging.error(f"File {f_type} requested but not found!")

    def run(self):
        logging.info("Start run")
        # TODO: Make task order dynamic

        self.tasks.put(Task(self._run_md_equil, {}))

        for i in range(self.iterations):
            # self.tasks.put(Task(self._run_md_prod, {"it": i}))
            # self.tasks.put(Task(self._query_reactions, {}))
            # self.tasks.put(Task(self._decide_reaction, {}))
            # self.tasks.put(Task(self._run_recipe, {"it": i}))

            self.tasks.put(Task(self._dummy, {}))
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
        for p in out_dir.iterdir():
            p_k = p.suffix[1:]  # strip dot
            while p_k in out_d.keys():
                p_k = p_k + "#"
            out_d[p_k] = p

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
        self.md.dummy_step(out_dir, **in_d)

        return in_d, out_dir
    
    
    def _run_md_equil(self):
        logging.info("Start equilibration MD")
        self.state = State.MD
        out_dir = self.config.out / f"equil_{self.iteration}"
        out_dir.mkdir()
        in_d = {
            "top": self.get_latest("top"),
            "gro": self.get_latest("gro"),
            "mdp": self.config.equilibrium.mdp,
            "idx": self.config.idx,
        }
        # perform step
        self.md.equilibrium(out_dir, **in_d)
        logging.info("Done equilibrating")
        return in_d, out_dir

    def _run_md_minim(self):
        logging.info("Start minimization md")
        self.state = State.MD
        out_dir = self.config.out / f"min_{self.iteration}"
        out_dir.mkdir()
        in_d = {
            "top": self.get_latest("top"),
            "gro": self.get_latest("gro"),
            "mdp": self.config.minimization.mdp,
        }
        # perform step
        self.md.minimzation(out_dir, **in_d)
        logging.info("Done minimizing")
        return in_d, out_dir

    def _run_md_eq(self):
        #TODO: combine w/ other equilibration?
        logging.info("Start _run_md_eq MD")
        self.state = State.MD
        out_dir = self.config.out / f"equil_{self.iteration}"
        out_dir.mkdir()
        in_d = {
            "top": self.get_latest("top"),
            "mdp": self.config.equilibrium.mdp,
            "gro": self.get_latest("gro"),
        }
        # perform step
        self.md.equilibration(out_dir, **in_d)
        logging.info("Done equilibrating")
        return in_d, out_dir


    def _run_md_prod(self):
        logging.info("Start production MD")
        self.state = State.MD
        out_dir = self.config.out / f"prod_{self.iteration}"
        out_dir.mkdir()
        in_d = {
            "top": self.get_latest("top"),
            "gro": self.get_latest("gro"),
            "mdp": self.config.prod.mdp,
            "idx": self.config.idx,
            "cpt": self.get_latest("cpt"),
            "dat": self.config.plumed.dat,
        }
        # perform step
        self.md.production(out_dir, **in_d)
        logging.info("Done minimizing")
        return in_d, out_dir

    def _run_md_relax(self):
        logging.info("Start relaxation MD")
        self.state = State.MD
        out_dir = self.config.out / f"prod_{self.iteration}"
        out_dir.mkdir()
        in_d = {
            "top": self.get_latest("top"),
            "gro": self.get_latest("gro"),
            "mdp": self.config.changer.coordinates.md.mdp,
            "idx": self.config.idx,
            "cpt": self.get_latest("cpt"),
        }
        # perform step
        self.md.production(out_dir, **in_d)
        logging.info("Done simulating")
        return in_d, out_dir
        

    def _query_reactions(self):
        logging.info("Query reactions")
        self.state = State.REACTION
        if hasattr(self.config.reactions, "homolysis"):
            self.reaction = Homolysis()
            self.rates, self.recipes = self.reaction.rupturerates(
                self.plumeddat,
                self.plumeddist,
                self.top,
                self.config.reactions.homolysis.bonds,
                self.config.reactions.homolysis.edis,
            )
        logging.info("Rates and Recipes:")
        logging.info(self.rates[0:10])
        logging.info(self.recipes[0:10])
        logging.info("Done")

    def _decide_reaction(self, decision_strategy=default_decision_strategy):
        logging.info("Decide on a reaction")
        self.chosen_recipe = decision_strategy(self.rates, self.recipes)
        logging.info("Chosen recipe is:")
        logging.info(self.chosen_recipe)

    def _run_recipe(self, it):
        logging.info(f"Start Recipe {it}")
        logging.info(f"Breakpair: {self.chosen_recipe.atom_idx}")
        self.changer = ChangeManager(it, self.chosen_recipe.atom_idx)
        newtop = (
            os.path.splitext(self.top)[0] + str(it) + "_.top"
        )  # there must be a better way to do this :/
        self.top = self.changer.modify_top(self.top, newtop)
        logging.info(f"Wrote new topology to {self.top}")
        newplumeddat = os.path.splitext(self.plumeddat)[0] + str(it) + "_.dat"
        newplumeddist = os.path.splitext(self.plumeddist)[0] + str(it) + "_.dat"
        self.plumeddat, self.plumeddist = self.changer.modify_plumed(
            self.plumeddat, newplumeddat, newplumeddist
        )
        logging.info(f"Wrote new plumedfile to {self.plumeddat}")
        logging.info(f"Looking for md in {self.config.changer.coordinates.__dict__}")
        if hasattr(self.config.changer.coordinates, "md"):
            self.tasks.put(Task(self._run_md_relax, {"it": it}))
        self.state = State.IDLE
