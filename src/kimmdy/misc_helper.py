"""
Miscelaneous utilitiies that didn't fit anywhere else for now.
"""
from pathlib import Path
import shutil
import argparse
import subprocess
from typing import Union
from kimmdy.topology.topology import Topology


def topology_to_edgelist(top: Topology):
    ls = []
    for b in top.bonds:
        ai = b[0]
        aj = b[1]
        x = top.atoms[ai]
        y = top.atoms[aj]
        left = f'"{ai} {x.type}" '
        right = f'"{aj} {y.type}" '
        ls.append(f"{left} -- {right};")
        if x.is_radical:
            ls.append(f'{left} [color="red"]')
        if y.is_radical:
            ls.append(f'{right} [color="red"]')
    return ls


def edgelist_to_dot_graph(ls: list[str], overlap: str = "true"):
    header = f"""
        graph G {{
          layout=neato
          overlap={overlap}
          node [shape="circle"]
    """
    tail = """
        }
    """
    body = "\n".join(ls)
    return header + body + tail


def top_to_graph(top: Topology, overlap: str = "true"):
    return edgelist_to_dot_graph(topology_to_edgelist(top), overlap)


def _build_examples(args: argparse.Namespace):
    overwrite = True if args.restore else False

    rename = {
        "hat_naive": "alanine",
        "homolysis": "hexalanine",
        "whole_run": "charged_peptide",
        "pull": "tripelhelix",
    }
    basedir = Path(__file__).parents[2]
    testpath = basedir / "tests" / "test_files" / "test_integration"
    examplepath = basedir / "example"
    relassetpath = Path(".") / ".." / ".." / "tests" / "test_files" / "assets"

    for k, v in rename.items():
        try:
            print("Building example", v, "...", end="")
            srcpath = testpath / k
            dstpath = examplepath / v
            if args.restore == "hard":
                shutil.rmtree(dstpath, ignore_errors=True)
            shutil.copytree(srcpath, dstpath, dirs_exist_ok=overwrite)
            (dstpath / "amber99sb-star-ildnp.ff").unlink(missing_ok=True)
            (dstpath / "amber99sb-star-ildnp.ff").symlink_to(
                relassetpath / "amber99sb-star-ildnp.ff",
                target_is_directory=True,
            )
            # fix scheme path
            yml = dstpath / "kimmdy.yml"
            if yml.exists():
                scheme_p = (
                    Path("..") / ".." / "src" / "kimmdy" / "kimmdy-yaml-schema.json"
                )
                first_line = f"# yaml-language-server: $schema={scheme_p}\n"
                with open(yml, "r+") as f:
                    lines = f.readlines()
                    if lines[0][:22] == "# yaml-language-server":
                        lines[0] = first_line
                    else:
                        lines = [first_line] + lines
                    f.seek(0)
                    f.writelines(lines)
                    f.truncate()

            print("done")
        except FileExistsError as e:
            raise FileExistsError(
                f"Could not build example directory {v} because it already exists. Try the --restore option to still build it."
            ) from e
