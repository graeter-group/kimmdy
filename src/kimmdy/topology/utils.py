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
    from kimmdy.topology.ff import Patch, Patches


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


def match_id_to_patch(id: list[str], patches: Patches) -> Optional[Patch]:
    id_str = "---".join([s.replace("*", "STAR").replace("+", "PLUS") for s in id])
    result = None
    longest_match = 0
    for _, patch in patches.items():
        # escape special regex characters
        # that can appear in a forcecield
        # use X as the wildcard
        for key in [patch.id, patch.id_sym]:
            key = key.replace("*", "STAR").replace("+", "PLUS")
            # early return exact match
            if key == id_str:
                return patch
            key_re = key.replace("X", ".*")
            match = re.match(key_re, id_str)
            if match is not None:
                # favor longer (=more specific) and later matches
                if len(key) >= longest_match:
                    longest_match = len(key)
                    result = patch

    return result

