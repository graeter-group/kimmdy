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
#
#
#
#
#
#
#
#
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
#
#
#
#
#
#
#
#| code-fold: true
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
#
#
#
