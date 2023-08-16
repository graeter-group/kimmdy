from pathlib import Path

from kimmdy.topology.topology import Topology
from kimmdy.parsing import read_top, write_top
from kimmdy.utils import run_gmx, run_shell_cmd
from kimmdy.parameterize import Parameterizer


def remove_hydrogen(args):
    """remove hydrogen from a gro and top file"""
    gro_path = Path(args.gro)
    top_path = Path(args.top)
    nr = args.nr
    parameterize = args.parameterize
    equilibrate = args.equilibrate

    top = Topology(read_top(top_path))

    ## check for input validity
    assert (atom_type := top.atoms[nr].type).startswith(
        "H"
    ), f"Wrong atom type {atom_type} with nr {nr} for remove hydrogen, should start with 'H'"

    ## deal with top file, order is important here
    [heavy_nr] = top.atoms[nr].bound_to_nrs

    del top.atoms[nr]

    top.atoms[heavy_nr].is_radical = True
    top.radicals[heavy_nr] = top.atoms[heavy_nr]

    update_map = {atom_nr: str(i + 1) for i, atom_nr in enumerate(top.atoms.keys())}

    top.reindex_atomnrs()

    ## parameterize with grappa
    if parameterize:
        import sys

        if sys.version_info > (3, 10):
            from importlib_metadata import entry_points

            discovered_parameterization_plugins = entry_points(
                group="kimmdy.parameterization_plugins"
            )
        else:
            from importlib.metadata import entry_points

            discovered_parameterization_plugins = entry_points()[
                "kimmdy.parameterization_plugins"
            ]

        parameterization_plugins: dict[str, Parameterizer | Exception] = {}
        for _ep in discovered_parameterization_plugins:
            try:
                parameterization_plugins[_ep.name] = _ep.load()
            except Exception as _e:
                parameterization_plugins[_ep.name] = _e

        if "grappa" in parameterization_plugins.keys():
            grappa = parameterization_plugins["grappa"]()
            grappa.parameterize_topology(top)

        else:
            raise KeyError(
                "No grappa in parameterization plugins. Can't continue to parameterize molecule"
            )

    ## write top file
    top_stem = top_path.stem
    top_outpath = top_path.with_stem(top_stem + f"_d{nr}")
    write_top(top.to_dict(), top_outpath)

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
        run_gmx(f"gmx mdrun -deffnm {gro_outpath.stem}_min", cwd=cwd)
        run_gmx(
            f"gmx grompp -f nvt.mdp -c {gro_outpath.stem}_min.gro -p {top_outpath.name} -o {gro_outpath.stem}_nvt.tpr -maxwarn 2",
            cwd=cwd,
        )
        run_gmx(f"gmx mdrun -deffnm {gro_outpath.stem}_nvt", cwd=cwd)
        run_gmx(
            f"gmx grompp -f npt.mdp -c {gro_outpath.stem}_nvt.gro -p {top_outpath.name} -o {gro_outpath.stem}_npt.tpr -maxwarn 2",
            cwd=cwd,
        )
        run_gmx(f"gmx mdrun -v -deffnm {gro_outpath.stem}_npt", cwd=cwd)
