from __future__ import annotations
import logging
import queue
from pathlib import Path
from enum import Enum, auto
from dataclasses import dataclass
from typing import Callable
from kimmdy.utils import run_shell_cmd,identify_atomtypes
from kimmdy.config import Config
from kimmdy.reaction import ConversionType, ConversionRecipe
from kimmdy.reactions.homolysis import Homolysis
from kimmdy.mdmanager import MDManager
from kimmdy.changemanager import ChangeManager

import os
import numpy as np
import random

def default_decision_strategy(list_of_rates,list_of_recipes):
    #compare e.g. https://en.wikipedia.org/wiki/Kinetic_Monte_Carlo#Rejection-free_KMC

    total_rate = sum(list_of_rates)  #sum all rates
    random.seed()
    t = random.random() #t in [0.0,1.0)
    rate_running_sum = 0
    for i in range(len(list_of_rates)): 
        rate_running_sum += list_of_rates[i]
        if (t*total_rate) <= rate_running_sum:
            idx = i
            #breakpair = list_of_recipes[i][1]
            #atomtypes = list_of_recipes[i][2]
            #rate = list_of_rates[i]                 
            break
    u = random.random()
    delta_t = np.log(1/u)/total_rate

    return list_of_recipes[idx]

    #print ('Rate = ' + str(rate) + ', breakpair = ' + str (breakpair) + ', atomtypes = ' + str(atomtypes) + 'jump [ps] = ' + str(delta_t))
    #logging.info(('Rate = ' + str(rate) + ', breakpair = ' + str (breakpair) + ', atomtypes = ' + str(atomtypes) + 'jump [ps] = ' + str(delta_t)))

    #return breakpair, atomtypes, delta_t

    ### moved to changer
    #rupture_time = delta_t / 1000.0 #convert ps to ns                          
    #newtop = '../broken_topol.top'
    # modify_top (topfile, newtop, breakpair)  
    #topfile = newtop
    
class State(Enum):
    IDLE = auto()
    MD = auto()
    REACTION = auto()
    DONE = auto()


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
    plumeddat: Path
    plumeddist: Path
    measurements: Path
    rates: list
    md: MDManager

    def __init__(self, input_file: Path):
        self.config = Config(input_file)
        self.tasks = queue.Queue()
        self.iteration = 0
        self.iterations=self.config.iterations
        self.state = State.IDLE
        #self.trj = self.config.cwd / Path("prod_" + str(self.iteration) + ".trj")
        self.top = self.config.top
        self.structure=self.config.gro
        self.measurements = Path("measurements")
        self.rates = []
        self.config.cwd.mkdir(parents=True,exist_ok=True)
        os.chdir(self.config.cwd)

        self.plumeddat = self.config.plumed['dat']
        self.plumeddist = self.config.plumed['distances']
        # self.md = MDManager()

    def run(self):
        logging.info("Start run")

        self.tasks.put(Task(self._run_md_equil,{'it':0}))

        for i in range(self.iterations):
            self.tasks.put(Task(self._run_md_prod,{'it':i}))
            self.tasks.put(Task(self._query_reactions,{}))
            self.tasks.put(Task(self._decide_reaction,{}))
            self.tasks.put(Task(self._run_recipe,{'it':i}))


        #self.tasks.put(Task(self._run_md_eq, {"ensemble": "npt"}))

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
        self.md = MDManager(self.top,self.config.gro,self.iteration,self.config.dryrun,self.config.minimization["mdp"])
        self.md.minimzation()
        self.state = State.IDLE
        logging.info("Done minimizing")

    def _run_md_eq(self, ensemble):
        logging.info("Start equilibration md")
        self.state = State.MD
        self.md = MDManager(self.top,self.config.gro,self.iteration,self.config.dryrun,self.config.equilibration[ensemble]["mdp"])
        self.md.equilibration(ensemble)
        self.state = State.IDLE
        logging.info("Done equilibrating")

    def _run_md_equil(self,it):
        logging.info("Start equilibration MD")
        self.state = State.MD
        self.md = MDManager(self.top,self.structure,self.config.idx,self.config.equilibrium['mdp'],it,'equilibrium',self.config.dryrun)
        self.structure = self.md.equilibrium()
        self.state = State.IDLE
        logging.info("Done equilibrating")

    def _run_md_prod(self,it):
        logging.info("Start production MD")
        self.state = State.MD
        self.md = MDManager(self.top,self.structure,self.config.idx,self.config.prod['mdp'],it,'pull',self.config.dryrun)
        self.structure = self.md.production('state.cpt', self.plumeddat)
        self.state = State.IDLE
        logging.info("Done simulating")

    def _run_md_relax(self,it):
        logging.info("Start relaxation MD")
        self.state = State.MD
        self.md = MDManager(self.top,self.structure,self.config.idx,self.config.changer['coordinates']['md']['mdp'],it,'relax',self.config.dryrun)
        self.structure = self.md.relaxation('state.cpt')
        self.state = State.IDLE
        logging.info("Done simulating")

    def _query_reactions(self):
        logging.info("Query reactions")
        self.state = State.REACTION
        if 'homolysis' in self.config.reactions:
            self.reaction = Homolysis()
            self.rates,self.recipes = self.reaction.rupturerates(self.plumeddat,self.plumeddist,self.top,self.config.reactions['homolysis']['bonds'],self.config.reactions['homolysis']['edis'])
        logging.info("Rates and Recipes:")
        logging.info(self.rates[0:10])
        logging.info(self.recipes[0:10])
        logging.info("Done")       

    def _decide_reaction(self, decision_strategy=default_decision_strategy):
        logging.info("Decide on a reaction")
        self.chosen_recipe = decision_strategy(self.rates,self.recipes)
        logging.info("Chosen recipe is:")
        logging.info(self.chosen_recipe)
 
    def _run_recipe(self,it):
        logging.info(f"Start Recipe {it}")
        logging.info(f"Breakpair: {self.chosen_recipe.atom_idx}")
        self.changer = ChangeManager(it,self.chosen_recipe.atom_idx)  
        self.top = self.changer.modify_top(self.top,f"topol_{it}.top")  #self.top.replace(".top","_0.top")
        logging.info(f"Wrote new topology to {self.top}")
        self.plumeddat, self.plumeddist = self.changer.modify_plumed(self.plumeddat,self.plumeddat.replace(".dat",f"_{it}.dat"),self.plumeddist.replace(".dat",f"_{it}.dat"))
        logging.info(f"Wrote new plumedfile to {self.plumeddat}")
        logging.info(f"Looking for md in {self.config.changer['coordinates'].keys()}")
        if 'md' in self.config.changer['coordinates'].keys():
            self.tasks.put(Task(self._run_md_relax,{'it':it}))
        self.state = State.IDLE
        





