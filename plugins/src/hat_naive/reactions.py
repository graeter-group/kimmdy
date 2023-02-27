from kimmdy.topology.atomic import Atom
from kimmdy.constants import ATOMTYPE_BONDORDER_FLAT
from kimmdy.reaction import (
    Conversion,
    ConversionRecipe,
    Reaction,
    ReactionOutcome,
    ReactionResult,
    ConversionType,
)
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


class HAT_naive(Reaction):
    """HAT reaction"""

    def get_reaction_result(self, files) -> ReactionResult:
        logging.info("Starting naive HAT reaction")
        top = self.runmng.top

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
            outcome = ReactionOutcome(
                recipe=[
                    Conversion(ConversionType.BREAK, (f, h)),
                    Conversion(ConversionType.BIND, (h, r)),
                ],
                rate=1,
            )
            return [outcome]

        return [ReactionOutcome([], 0)]
