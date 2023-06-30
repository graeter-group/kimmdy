from __future__ import annotations  # for 3.7 <= Python version < 3.10
from itertools import takewhile, permutations
from typing import Optional, Any
from xml.etree.ElementTree import Element
import re
from typing import TYPE_CHECKING

from dill.logger import logging

if TYPE_CHECKING:
    from kimmdy.topology.ff import Patch, Patches
    from kimmdy.topology.atomic import AtomicType, AtomicTypes


def get_top_section(
    top: dict, name: str, moleculetype: Optional[int] = None
) -> Optional[list[list]]:
    """Get content of a section from a topology dict.
    By resolving any `#ifdef` statements by check in the top['define'] dict
    and choosing the 'content' or 'else_content' depending on the result.
    """
    if moleculetype is not None:
        parent_name = f"moleculetype_{moleculetype}"
        parent_section = top.get(parent_name)
        if parent_section is None:
            logging.warning(f"topology does not contain moleculetype {moleculetype}")
            return None
        section = parent_section["subsections"].get(name)
    else:
        section = top.get(name)

    if section is None:
        logging.warning(f"topology does not contain section {name}")
        return None
    condition = section.get("condition")
    if condition is not None:
        condition_type = condition.get("type")
        condition_value = condition.get("value")
        if condition_type == "ifdef":
            if condition_value in top["define"].keys():
                return section.get("content")
            else:
                return section.get("else_content")
        elif condition_type == "ifndef":
            if condition_value not in top["define"].keys():
                return section.get("content")
            else:
                return section.get("else_content")
        else:
            raise NotImplementedError(
                f"condition type {condition_type} is not supported"
            )
    return section.get("content")


def get_protein_section(top: dict, name: str) -> Optional[list[list]]:
    """
    Get content of a section in the first moleculetype (protein) from a topology dict.
    """
    return get_top_section(top, name, moleculetype=0)


def set_top_section(
    top: dict, name: str, value: list, moleculetype: Optional[int] = None
) -> Optional[list[list]]:
    """Set content of a section from a topology dict.
    By resolving any `#ifdef` statements by check in the top['define'] dict
    and choosing the 'content' or 'else_content' depending on the result.
    """
    if moleculetype is not None:
        parent_name = f"moleculetype_{moleculetype}"
        parent_section = top.get(parent_name)
        if parent_section is None:
            raise ValueError(f"topology does not contain moleculetype {moleculetype}")
        section = parent_section["subsections"].get(name)
    else:
        section = top.get(name)

    if section is None:
        raise ValueError(f"topology does not contain section {name}")
    condition = section.get("condition")
    if condition is not None:
        condition_type = condition.get("type")
        condition_value = condition.get("value")
        if condition_type == "ifdef":
            if condition_value in top["define"].keys():
                section["content"] = value
            else:
                section["else_content"] = value
        elif condition_type == "ifndef":
            if condition_value not in top["define"].keys():
                section["content"] = value
            else:
                section["else_content"] = value
        else:
            raise NotImplementedError(
                f"condition type {condition_type} is not supported"
            )
    section["content"] = value


def set_protein_section(top: dict, name: str, value: list) -> Optional[list[list]]:
    """
    Set content of a section in the first moleculetype (protein) from a topology dict.
    """
    set_top_section(top, name, value, moleculetype=0)


def field_or_none(l: list[str], i) -> Optional[str]:
    try:
        return l[i]
    except IndexError as _:
        return None


def attributes_to_list(obj) -> list[str]:
    attrs = []
    for k, v in obj.__dict__.items():
        if k in ["bound_to_nrs", "is_radical", "id", "id_sym"]:
            continue
        if v in [None, ""]:
            continue
        attrs.append(v)
    return attrs


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
    id = [s.replace("*", "STAR").replace("+", "PLUS") for s in id]
    id_str = "---".join(id)
    id_sym = reversed(id)
    id_str = "---".join(id)
    id_sym_str = "---".join(id_sym)
    result = None
    longest_match = 0
    for _, patch in patches.items():
        # escape special regex characters
        # that can appear in a forcecield
        # use X as the wildcard
        for key in [patch.id, patch.id_sym]:
            key = key.replace("*", "STAR").replace("+", "PLUS")
            # early return exact match
            if key == id_str or key == id_sym_str:
                return patch
            key_re = key.replace("X", ".*")
            match = re.match(key_re, id_str)
            if match is not None:
                # favor longer (=more specific) and later matches
                if len(key) >= longest_match:
                    longest_match = len(key)
                    result = patch

    return result


def match_atomic_item_to_atomic_type(
    id: list[str], types: AtomicTypes
) -> Optional[AtomicType]:
    def escape_re_atomtypes(s: str) -> str:
        """
        escape special regex characters
        that can appear in a forcecield
        use X as the wildcard
        """
        return s.replace("*", "STAR").replace("+", "PLUS")

    id = [escape_re_atomtypes(s) for s in id]
    id_sym = reversed(id)
    id_str = "---".join(id)
    id_sym_str = "---".join(id_sym)
    result = None
    longest_match = 0
    for _, atomic_type in types.items():
        for key in [atomic_type.id, atomic_type.id_sym]:
            key = escape_re_atomtypes(key).replace("X", ".*")
            keys = key.split("---")
            # early return exact match
            if key == id_str or key == id_sym_str:
                return atomic_type
            matches = [re.match(pattern, s) for pattern, s in zip(keys, id)]
            if all(matches):
                # favor longer (=more specific) and later matches
                if len(key) >= longest_match:
                    longest_match = len(key)
                    result = atomic_type

    return result
