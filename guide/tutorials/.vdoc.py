# type: ignore
# flake8: noqa
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
<<<<<<< HEAD
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#| code-fold: true
#| code-summary: "KIMMDY config file (kimmdy.yml)"

# yaml-language-server: $schema=../../src/kimmdy/kimmdy-yaml-schema.json

dryrun: false
name: 'hexalanine_homolysis_000'
max_tasks: 100
gromacs_alias: 'gmx'
gmx_mdrun_flags: -maxh 24 -dlb yes -nt 8
ff: 'amber99sb-star-ildnp.ff' # optional, dir endinng with .ff by default 
top: 'hexala_out.top'
gro: 'npt.gro'
ndx: 'index.ndx'
plumed: 'plumed.dat'
mds:
  equilibrium:
    mdp: 'md.mdp'
  pull:
    mdp: 'md.mdp'
    use_plumed: true
  relax:
    mdp: 'md_slow_growth.mdp'
changer:
  coordinates:
    md: 'relax'
    slow_growth: True
  topology:
    parameterization: 'basic' 
reactions:
  homolysis:
    edis: 'edissoc.dat'
    itp: 'ffbonded.itp'
sequence:
  - equilibrium
  - pull  
  - reactions
  - equilibrium
  - pull

=======
#| code-fold: true
#| code-summary: "Plugin main code (reaction.py)"

from kimmdy.recipe import (
    Bind,
    Recipe,
    RecipeCollection,
)
from kimmdy.plugins import ReactionPlugin
from kimmdy.tasks import TaskFiles
from kimmdy.utils import (
    get_atomnrs_from_plumedid,
)
from kimmdy.parsing import (
    read_plumed,
    read_distances_dat,
)


class BindReaction(ReactionPlugin):
    """Reaction to bind two particles if they are in proximity
    """

    def get_recipe_collection(self, files: TaskFiles):
        logger = files.logger
        logger.debug("Getting recipe for reaction: Bind")

        # get cutoff distance from config, unit is [nm]
        cutoff = self.config.distance

        # Initialization of objects from files
        distances = read_distances_dat(files.input["plumed_out"])
        plumed = read_plumed(files.input["plumed"])

        recipes = []
        # check distances file for values lower than the cutoff
        for plumedid, dists in distances.items():
            if plumedid == "time":
                continue
            # checks only last frame
            if dists[-1] < cutoff:
                # get atomnrs from plumedid 
                atomnrs = get_atomnrs_from_plumedid(plumedid, plumed)
                recipes.append(
                    Recipe(
                        recipe_steps=[
                            Bind(atom_id_1=atomnrs[0], atom_id_2=atomnrs[1]),
                        ],
                        rates=[1],
                        timespans=[(distances["time"][0], distances["time"][-1])],
                    )
                )

        return RecipeCollection(recipes)
>>>>>>> main
#
#
#
#
#
#
#
#
<<<<<<< HEAD
=======
#| code-fold: true
#| code-summary: "Plugin schema (kimmdy-yaml-schema.json)"

{
  "$schema": "https://json-schema.org/draft/2020-12/schema",
  "type": "object",
  "$id": "bind-config",
  "description": "Settings for bind reactions",
  "additionalProperties": false,
  "properties": {
    "distance": {
      "description": "Cutoff distance for two particles to bind [nm]",
      "type": "float",
      "pytype": "float"
    },
    "kmc": {
      "description": "KMC algorithm for this reaction.",
      "type": "string",
      "pytype": "str",
      "enum": ["rfkmc", "frm", "extrande"],
      "default": "rfkmc"
    }
  },
  "required": ["distance"]
}
>>>>>>> main
#
#
#
#
#
#
#
#| code-fold: true
<<<<<<< HEAD
#| code-summary: "KIMMDY config file (kimmdy.yml)"

# yaml-language-server: $schema=../../src/kimmdy/kimmdy-yaml-schema.json

dryrun: false
name: 'fibril_000'
max_tasks: 100
gromacs_alias: 'gmx'
gmx_mdrun_flags: -maxh 24 -dlb yes -ntomp 4 -ntmpi 10  -pme gpu -bonded gpu -npme 1
ff: 'amber99sb-star-ildnp.ff' # optional, dir endinng with .ff by default 
top: 'hexala_out.top'
gro: 'npt.gro'
ndx: 'index.ndx'
plumed: 'plumed.dat'
mds:
  equilibrium:
    mdp: 'md.mdp'
  pull:
    mdp: 'md.mdp'
    use_plumed: true
  relax:
    mdp: 'md_slow_growth.mdp'
changer:
  coordinates:
    md: 'relax'
    slow_growth: True
  topology:
    parameterization: 'basic' 
reactions:
  homolysis:
    edis: 'edissoc.dat'
    itp: 'ffbonded.itp'
sequence:
  - equilibrium
  - pull  
  - reactions
  - equilibrium
  - pull

#
#
#
#
=======
#| code-summary: "setup.py"

from setuptools import setup

setup()
#
#
#
#| code-fold: true
#| code-summary: "setup.cfg"

[metadata]
name = kimmdy-reactions
version = 0.1
license = GPL-3.0 
description = Reaction template for KIMMDY
long_description = file: README.md
author = Eric Hartmann
author_email = eric.Hartmann@h-its.org
classifiers=
        Programming Language :: Python :: 3
        License :: OSI Approved :: MIT License
        Operating System :: OS Independent

[options]
packages = find:
package_dir =
    =src
include_package_data = True
install_requires =
    MDAnalysis

python_requires = >= 3.9

[options.packages.find]
where=src

[options.entry_points]
kimmdy.reaction_plugins =
    bind = bind.reaction:BindReaction
>>>>>>> main
#
#
#
