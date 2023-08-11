"""
Use quartodoc renderer from <https://github.com/machow/quartodoc/issues/135>
until the feature is upstreamed into quartodoc
to automatically link type annotations.
"""
from __future__ import annotations

import base64
import html
import re
from importlib.resources import files
from pathlib import Path
from typing import Literal, TypedDict, Union

import quartodoc.ast as qast
from griffe import dataclasses as dc
from griffe.expressions import Expression
from griffe.docstrings import dataclasses as ds
from plum import dispatch
from quartodoc import MdRenderer
from quartodoc.renderers.base import convert_rst_link_to_md, sanitize

KIMMDY_PATH = Path(files("kimmdy").joinpath())


DOCSTRING_TEMPLATE = """\
{rendered}

{header} Examples

{examples}
"""


# This is the same as the FileContentJson type in TypeScript.
class FileContentJson(TypedDict):
    name: str
    content: str
    type: Literal["text", "binary"]


class Renderer(MdRenderer):
    style = "shiny"

    @dispatch
    def render(self, el: qast.DocstringSectionSeeAlso):
        # The See Also section in the Shiny docs has bare function references, ones that
        # lack a leading :func: and backticks. This function fixes them. In the future,
        # we can fix the docstrings in Shiny, once we decide on a standard. Then we can
        # remove this function.
        return prefix_bare_functions_with_func(el.value)

    @dispatch
    def render(self, el: Union[dc.Object, dc.Alias]):
        rendered = super().render(el)

        converted = convert_rst_link_to_md(rendered)

        p_example_dir = KIMMDY_PATH / "examples" / el.name
        if (p_example_dir / "app.py").exists():
            example = ""

            files = list(p_example_dir.glob("**/*"))

            # Sort, and then move app.py to first position.
            files.sort()
            app_py_idx = files.index(p_example_dir / "app.py")
            files = [files[app_py_idx]] + files[:app_py_idx] + files[app_py_idx + 1 :]

            for f in files:
                file_info = read_file(f, p_example_dir)
                if file_info["type"] == "text":
                    example += f"\n## file: {file_info['name']}\n{file_info['content']}"
                else:
                    example += f"\n## file: {file_info['name']}\n## type: binary\n{file_info['content']}"

            return DOCSTRING_TEMPLATE.format(
                rendered=converted,
                examples=example,
                header="#" * (self.crnt_header_level + 1),
            )

        return converted

    @dispatch
    def render(self, el: ds.DocstringSectionText):
        # functions like shiny.ui.tags.b have html in their docstrings, so
        # we escape them. Note that we are only escaping text sections, but
        # since these cover the top text of the docstring, it should solve
        # the immediate problem.
        rendered = super().render(el)
        return html_escape_except_backticks(rendered)

    @dispatch
    def render_annotation(self, el: str):
        return sanitize(el)

    @dispatch
    def render_annotation(self, el: None):
        return ""

    @dispatch
    def render_annotation(self, el: Expression):
        # an expression is essentially a list[dc.Name | str]
        # e.g. Optional[TagList]
        #   -> [Name(source="Optional", ...), "[", Name(...), "]"]

        return "".join(map(self.render_annotation, el))

    @dispatch
    def render_annotation(self, el: dc.Name):
        # e.g. Name(source="Optional", full="typing.Optional")
        return f"[{el.source}](`{el.full}`)"

    @dispatch
    def summarize(self, el: dc.Object | dc.Alias):
        result = super().summarize(el)
        return html.escape(result)


def html_escape_except_backticks(s: str) -> str:
    """
    HTML-escape a string, except for content inside of backticks.

    Examples
    --------
        s = "This is a <b>test</b> string with `backticks <i>unescaped</i>`."
        print(html_escape_except_backticks(s))
        #> This is a &lt;b&gt;test&lt;/b&gt; string with `backticks <i>unescaped</i>`.
    """
    # Split the string using backticks as delimiters
    parts = re.split(r"(`[^`]*`)", s)

    # Iterate over the parts, escaping the non-backtick parts, and preserving backticks in the backtick parts
    escaped_parts = [
        html.escape(part) if i % 2 == 0 else part for i, part in enumerate(parts)
    ]

    # Join the escaped parts back together
    escaped_string = "".join(escaped_parts)
    return escaped_string


def prefix_bare_functions_with_func(s: str) -> str:
    """
    The See Also section in the Shiny docs has bare function references, ones that lack
    a leading :func: and backticks. This function fixes them.

    If there are bare function references, like "~shiny.ui.panel_sidebar", this will
    prepend with :func: and wrap in backticks.

    For example, if the input is this:
        "~shiny.ui.panel_sidebar  :func:`~shiny.ui.panel_sidebar`"
    This function will return:
        ":func:`~shiny.ui.panel_sidebar`  :func:`~shiny.ui.panel_sidebar`"
    """

    def replacement(match: re.Match[str]) -> str:
        return f":func:`{match.group(0)}`"

    pattern = r"(?<!:func:`)~\w+(\.\w+)*"
    return re.sub(pattern, replacement, s)


def read_file(file: str | Path, root_dir: str | Path | None = None) -> FileContentJson:
    file = Path(file)
    if root_dir is None:
        root_dir = Path("/")
    root_dir = Path(root_dir)

    type: Literal["text", "binary"] = "text"

    try:
        with open(file, "r") as f:
            file_content = f.read()
            type = "text"
    except UnicodeDecodeError:
        # If text failed, try binary.
        with open(file, "rb") as f:
            file_content_bin = f.read()
            file_content = base64.b64encode(file_content_bin).decode("utf-8")
            type = "binary"

    return {
        "name": str(file.relative_to(root_dir)),
        "content": file_content,
        "type": type,
    }
