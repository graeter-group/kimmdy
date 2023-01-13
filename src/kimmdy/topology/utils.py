from __future__ import annotations  # for 3.7 <= Python version < 3.10
from itertools import takewhile, permutations
from typing import Optional, Any
from pathlib import Path
from xml.etree.ElementTree import Element
import re
from typing import TYPE_CHECKING
from kimmdy.topology.atomic import *

if TYPE_CHECKING:
    from kimmdy.topology.topology import Topology


def field_or_none(l: list[str], i) -> Optional[str]:
    try:
        return l[i]
    except IndexError as _:
        return None


def attributes_to_list(obj) -> list[str]:
    return list(takewhile(lambda x: x is not None, obj.__dict__.values()))


def is_not_none(x) -> bool:
    return x is None


def match_attr(patches: list[Element], attr: str, m: str) -> Optional[Element]:
    matches = []
    for p in patches:
        if value := p.get(attr):
            if value == m:
                return p
            pattern = value.replace("*", r".*").replace("+", r".+")
            if re.match(pattern, m):
                matches.append(p)
    if matches:
        matches.sort(key=lambda x: x.get(attr))
        return matches[0]
    else:
        return None


def match_multi_attr(
    patches: list[Element], attrs: list[str], m: list[str]
) -> Optional[Element]:
    raise NotImplementedError("WIP")


def get_by_permutations(d: dict, key) -> Optional[Any]:
    for k in permutations(key):
        value = d.get(k, None)
        if value is not None:
            return value
    return None


def get_element_id(e: Element) -> Optional[str]:
    id = None
    if e.tag == "Atom":
        id = ""
    elif e.tag == "Bond":
        id = ""
    return id


def generate_topology_from_bound_to(
    atoms: list[Atom], ffdir: Path, ffpatch: Path
) -> Topology:
    top = Topology({}, ffdir, ffpatch)
    for atom in atoms:
        top.atoms[atom.nr] = atom

    # bonds
    keys = []
    for atom in top.atoms.values():
        keys = top._get_atom_bonds(atom.nr)
        for key in keys:
            top.bonds[key] = Bond(key[0], key[1], "1")

    # angles
    for atom in top.atoms.values():
        keys = top._get_atom_angles(atom.nr)
        for key in keys:
            top.angles[key] = Angle(key[0], key[1], key[2], "1")

    # dihedrals and pass
    for atom in top.atoms.values():
        keys = top._get_atom_proper_dihedrals(atom.nr)
        for key in keys:
            top.proper_dihedrals[key] = Dihedral(key[0], key[1], key[2], key[3], "9")
            pairkey = tuple(str(x) for x in sorted([key[0], key[3]], key=int))
            if top.pairs.get(pairkey) is None:
                top.pairs[pairkey] = Pair(pairkey[0], pairkey[1], "1")

    for atom in top.atoms.values():
        impropers = top._get_atom_improper_dihedrals(atom.nr)
        for key, improper in impropers:
            top.improper_dihedrals[key] = Dihedral(
                improper.atom1,
                improper.atom2,
                improper.atom3,
                improper.atom4,
                "4",
                improper.cq,
            )

    return top


def match_id_to_patch_id(id: list[str], keys: list[str]) -> Optional[str]:
    result = None
    longest_match = 0
    for key in keys:
        if key == id:
            # early return exact match
            return key
        # escape special regex characters
        # that can appear in a forcecield
        # use X as the wildcard
        id = [s.replace("*", "STAR").replace("+", "PLUS") for s in id]
        key_re = key.replace("*", "STAR").replace("+", "PLUS").replace("X", ".*")
        # construct id permutations
        for perm in permutations(id):
            id_perm = "---".join(perm)
            match = re.match(key_re, id_perm)
            if match is not None:
                # favor longer (=more specific) and later matches
                if len(key) >= longest_match:
                    result = key

    return result

