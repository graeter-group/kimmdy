"""
Standalone tools that are complementary to KIMMDY.
"""

from pathlib import Path
import shutil

from kimmdy.topology.topology import Topology
from kimmdy.parsing import read_top, write_top
from kimmdy.utils import run_gmx, run_shell_cmd


def build_examples(restore: str):
    print(f"build_examples arguments: restore {restore}")
    overwrite = True if restore else False

    rename = {
        "hat_naive": "alanine",
        "homolysis": "hexalanine",
        "whole_run": "charged_peptide",
        "pull": "tripelhelix",
        "single_reaction": "single_reaction",
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
            if restore == "hard":
                shutil.rmtree(dstpath, ignore_errors=True)
            print(srcpath, dstpath)
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


def remove_hydrogen(
    gro: str,
    top: str,
    nr: str,
    parameterize: bool,
    equilibrate: bool,
    gmx_mdrun_flags: str,
):
    """remove hydrogen from a gro and top file"""
    gro_path = Path(gro)
    top_path = Path(top)

    topology = Topology(read_top(top_path))

    ## check for input validity
    assert (atom_type := topology.atoms[nr].type).startswith(
        "H"
    ), f"Wrong atom type {atom_type} with nr {nr} for remove hydrogen, should start with 'H'"

    ## deal with top file, order is important here
    [heavy_nr] = topology.atoms[nr].bound_to_nrs

    del topology.atoms[nr]

    topology.atoms[heavy_nr].is_radical = True
    topology.radicals[heavy_nr] = topology.atoms[heavy_nr]

    update_map = {
        atom_nr: str(i + 1) for i, atom_nr in enumerate(topology.atoms.keys())
    }

    topology.reindex_atomnrs()

    ## parameterize with grappa
    if parameterize:
        from kimmdy import parameterization_plugins

        if "grappa" in parameterization_plugins.keys():
            grappa = parameterization_plugins["grappa"]()
            grappa.parameterize_topology(topology)

        else:
            raise KeyError(
                "No grappa in parameterization plugins. Can't continue to parameterize molecule"
            )

    ## write top file
    top_stem = top_path.stem
    top_outpath = top_path.with_stem(top_stem + f"_d{nr}")
    write_top(topology.to_dict(), top_outpath)

    ## deal with gro file
    with open(gro_path, "r") as f:
        gro_raw = f.readlines()

    if len(gro_raw[3].split()) > 6:
        print("gro file likely contains velocities. These will be discarded.")

    gro_stem = gro_path.stem
    gro_outpath = gro_path.with_stem(gro_stem + f"_d{nr}")
    with open(gro_outpath, "w") as f:
        f.write(gro_raw[0])
        f.write(f"   {str(int(gro_raw[1])-1)}\n")
        for line in gro_raw[2:-1]:
            linesplit = line.split()
            if val := update_map.get(linesplit[2]):
                linesplit[2] = val
            else:
                continue
            f.write("{:>8s}   {:>4s} {:>4s} {:>7s} {:>7s} {:>7s}\n".format(*linesplit))
            # format not exactly as defined by GROMACS but should be good enough
        f.write(gro_raw[-1])

    ## minimize and equilibrate system using GROMACS
    if equilibrate:
        cwd = top_outpath.parent
        run_shell_cmd(
            f"cp {str(Path(__file__).parents[2] / 'tests'/'test_files'/'assets'/'md'/'*')} .",
            cwd=cwd,
        )
        run_gmx(
            f"gmx editconf -f {gro_outpath.name} -o {gro_outpath.stem}_box.gro -c -d 1.0 -bt dodecahedron",
            cwd=cwd,
        )
        run_gmx(
            f"gmx solvate -cp {gro_outpath.stem}_box.gro -p {top_outpath.name} -o {gro_outpath.stem}_solv.gro",
            cwd=cwd,
        )
        run_gmx(
            f"gmx grompp -f ions.mdp -c {gro_outpath.stem}_solv.gro -p {top_outpath.name} -o {gro_outpath.stem}_genion.tpr",
            cwd=cwd,
        )
        run_gmx(
            f"echo 'SOL' | gmx genion -s {gro_outpath.stem}_genion.tpr -p {top_outpath.name} -o {gro_outpath.stem}_ion.gro -conc 0.15 -neutral",
            cwd=cwd,
        )
        run_gmx(
            f"gmx grompp -f minim.mdp -c {gro_outpath.stem}_ion.gro -p {top_outpath.name} -o {gro_outpath.stem}_min.tpr -maxwarn 2",
            cwd=cwd,
        )
        run_gmx(f"gmx mdrun -deffnm {gro_outpath.stem}_min" + gmx_mdrun_flags, cwd=cwd)
        run_gmx(
            f"gmx grompp -f nvt.mdp -c {gro_outpath.stem}_min.gro -p {top_outpath.name} -o {gro_outpath.stem}_nvt.tpr -maxwarn 2",
            cwd=cwd,
        )
        run_gmx(f"gmx mdrun -deffnm {gro_outpath.stem}_nvt" + gmx_mdrun_flags, cwd=cwd)
        run_gmx(
            f"gmx grompp -f npt.mdp -c {gro_outpath.stem}_nvt.gro -p {top_outpath.name} -o {gro_outpath.stem}_npt.tpr -maxwarn 2",
            cwd=cwd,
        )
        run_gmx(
            f"gmx mdrun -v -deffnm {gro_outpath.stem}_npt" + gmx_mdrun_flags, cwd=cwd
        )


## dot graphs


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
