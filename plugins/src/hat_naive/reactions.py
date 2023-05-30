from kimmdy.topology.atomic import Atom
from kimmdy.constants import ATOMTYPE_BONDORDER_FLAT
from kimmdy.reaction import (
    Break,
    Bind,
    Recipe,
    RecipeCollection,
    ReactionPlugin,
)
import MDAnalysis as mda
import random as rng
import logging


def find_radical(atoms: list[Atom]):
    for atom in atoms:
        if atom.is_radical:
            return atom
        bo = ATOMTYPE_BONDORDER_FLAT.get(atom.type)
        if bo and bo > len(atom.bound_to_nrs):
            return atom

    return None


class HAT_naive(ReactionPlugin):
    """Naive HAT reaction, selects hydrogens at random"""

    def get_recipe_collection(self, files) -> RecipeCollection:
        logging.info("Starting naive HAT reaction")
        top = self.runmng.top

        tpr = files.input["tpr"]
        trr = files.input["trr"]
        u = mda.Universe(str(tpr), str(trr), topology_format="tpr", format="trr")

        if not top.radicals:
            radical = find_radical(list(top.atoms.values()))
            if radical:
                top.radicals[radical.nr] = radical

        if top.radicals:
            rad = rng.sample(list(top.radicals.values()), 1)[0]
            hs = []
            froms = []
            for nr in rad.bound_to_nrs:
                atom = top.atoms[nr]
                for nr2 in atom.bound_to_nrs:
                    atom2 = top.atoms[nr2]
                    if atom2.type.startswith("H"):
                        froms.append(atom.nr)
                        hs.append(atom2.nr)
            logging.info(f"hs: {hs}")
            logging.info(f"froms: {froms}")
            i = rng.randint(0, len(hs) - 1)
            r = rad.nr
            h = hs[i]
            f = froms[i]
            logging.info(f"i: {i}")
            logging.info(f"radical: {rad}")
            logging.info(f"h: {top.atoms[h]}")
            logging.info(f"from: {top.atoms[f]}")

            recipe = Recipe(
                recipe_steps=[Break(f, h), Bind(h, r)],
                rates=[1],
                times=[u.trajectory[-1].time],
            )
            return RecipeCollection([recipe])

        return RecipeCollection([])
