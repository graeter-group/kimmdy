"""Reaction plugin building blocks"""

import logging
from pathlib import Path
from typing import Optional, TypeAlias, TypedDict

import MDAnalysis as mda
import numpy as np
from MDAnalysis.lib.distances import calc_bonds

from kimmdy.constants import R
from kimmdy.parsing import Plumed_dict, read_distances_dat, read_edissoc, read_plumed
from kimmdy.topology.atomic import BondId
from kimmdy.topology.topology import Topology

logger = logging.getLogger(__name__)


class Stats(TypedDict):
    plumed_id: str
    mean_d: float
    mean_f: float
    delta_d: float
    b0: float
    harmonic_f: Optional[float]  # Optional, only if calculated


BondStats: TypeAlias = dict[BondId, Stats]
PlumedId: TypeAlias = str
Time: TypeAlias = float
Distance: TypeAlias = float
BONDSTATS_COLUMNS = "ai,aj,plumed_id,mean_d,mean_f,delta_d,b0,harmonic_f"
BondToPlumedID: TypeAlias = dict[BondId, PlumedId]
PlumedDistances: TypeAlias = dict[Time, dict[PlumedId, Distance]]


def bondstats_to_csv(stats: BondStats, path: str | Path):
    ls = []
    ls.append(BONDSTATS_COLUMNS)
    for k, s in stats.items():
        l = f"{k[0]},{k[1]},{s['plumed_id']},{s['mean_d']:.6f},{s['mean_f']:.6f},{s['delta_d']:.6f},{s['b0']:.6f}"
        if s.get("harmonic_f") is not None:
            l += f",{s['harmonic_f']:.6f}"
        ls.append(l)

    with open(path, "w") as f:
        f.write("\n".join(ls))


def bondstats_from_csv(path: str | Path) -> BondStats:
    stats: BondStats = {}
    with open(path, "r") as f:
        next(f)
        for line in f:
            line = line.strip().split(",")
            harmonic_f = None
            if len(line) == 7:
                (ai, aj, plumed_id, mean_d, mean_f, delta_d, b0) = line
            elif len(line) == 8:
                (ai, aj, plumed_id, mean_d, mean_f, delta_d, b0, harmonic_f) = line
            else:
                m = f"Line in bondstats csv file does not have 7 or 8 columns: {line}"
                logger.error(m)
                raise ValueError(m)

            stats[(ai, aj)] = {
                "plumed_id": plumed_id,
                "mean_d": float(mean_d),
                "mean_f": float(mean_f),
                "delta_d": float(delta_d),
                "b0": float(b0),
                "harmonic_f": (float(harmonic_f) if harmonic_f is not None else None),
            }

    return stats


def read_plumed_input(path: str | Path) -> BondToPlumedID:
    if not isinstance(path, Path):
        path = Path(path)
    plumed = read_plumed(path)
    d = {}
    for k, v in plumed["labeled_action"].items():
        if v["keyword"] != "DISTANCE":
            continue
        atoms = v["atoms"]
        bondkey = tuple(sorted(atoms, key=int))
        d[bondkey] = k
    return d


def select_bonds(
    top: Topology,
    atom_nrs: Optional[set[str]] = None,
    exclude: Optional[list[str]] = None,
) -> list[BondId]:
    """Select bonds of a topology for bond distance sampling.

    Parameters
    ----------
    top
        Topology to take the bonds from.
    atom_nrs
        Only keep bonds with both atoms in this set of atom numbers.
        All bonds are considered if None.
    exclude
        Atom names or bonds to leave out. An atom name (e.g. "O") excludes
        bonds containing that atom, a trailing `*` (e.g. "H*") excludes
        all atom names starting with the prefix and two atom names joined
        by a hyphen (e.g. "C-N") exclude bonds between those atoms.

    Returns
    -------
    Selected bonds as tuples of atom numbers.
    """
    atoms_to_remove = []
    atoms_start_to_remove = []
    bonds_to_remove = []
    for entry in exclude or []:
        entrysplit = entry.split(sep="-")
        if len(entrysplit) == 1:
            if entrysplit[0].endswith("*"):
                atoms_start_to_remove.append(entrysplit[0][:-1])
            else:
                atoms_to_remove.append(entrysplit[0])
        elif len(entrysplit) == 2:
            bonds_to_remove.append(set(entrysplit))
        else:
            logger.warning(f"Don't know how to deal with entry {entry}, skipping it!")

    bonds = []
    for bond in top.bonds.keys():
        if atom_nrs is not None and not (bond[0] in atom_nrs and bond[1] in atom_nrs):
            continue
        ai_name = top.atoms[bond[0]].atom
        aj_name = top.atoms[bond[1]].atom
        if set([ai_name, aj_name]) in bonds_to_remove:
            continue
        if any([b in atoms_to_remove for b in [ai_name, aj_name]]):
            continue
        if any(
            [
                b.startswith(prefix)
                for b in [ai_name, aj_name]
                for prefix in atoms_start_to_remove
            ]
        ):
            continue
        bonds.append(bond)
    return bonds


def calculate_bondstats(
    top: Topology,
    plumed_in: Path,
    plumed_out: Path,
    dt: float = 0.0,
    edissoc_dat: Optional[Path] = None,
) -> BondStats:
    distances = read_distances_dat(path=plumed_out, dt=dt)
    bond_to_plumed_id = read_plumed_input(plumed_in)
    return get_bondstats(
        top=top,
        distances=distances,
        bond_to_plumed_id=bond_to_plumed_id,
        edissoc_dat=edissoc_dat,
    )


def get_bondstats(
    top: Topology,
    distances: PlumedDistances,
    bond_to_plumed_id: BondToPlumedID,
    edissoc_dat: Path | None = None,
) -> BondStats:
    if edissoc_dat is not None:
        edissoc = read_edissoc(edissoc_dat)
    else:
        edissoc = None

    stats: BondStats = {}
    for bondkey, plumed_id in bond_to_plumed_id.items():
        ds = np.asarray([values[plumed_id] for values in distances.values()])
        stats[bondkey] = _stats_for_bond(
            top=top, bondkey=bondkey, bond_id=plumed_id, ds=ds, edissoc=edissoc
        )
    return stats


def _stats_for_bond(
    top: Topology,
    bondkey: BondId,
    bond_id: str,
    ds: np.ndarray,
    edissoc: Optional[dict],
) -> Stats:
    """Mean distance and forces of one bond from its distances in nm."""
    ai = top.atoms[bondkey[0]]
    aj = top.atoms[bondkey[1]]
    typekey = (ai.type, aj.type)
    bondtype = top.ff.bondtypes.get(typekey)
    if bondtype is None:
        # attempt the reverse key
        typekey = (aj.type, ai.type)
        bondtype = top.ff.bondtypes.get(typekey)
    if not bondtype or bondtype.c0 is None or bondtype.c1 is None:
        m = f"Could not find bondtype of atoms with ids {ai.nr} and {aj.nr}"
        logger.error(m)
        raise ValueError(m)
    b0 = float(bondtype.c0)
    kb = float(bondtype.c1)

    if edissoc is not None:
        edis = get_edissoc_from_atomnames(
            atomnames=[ai.atom, aj.atom], edissoc=edissoc, residue=ai.residue
        )
    else:
        m = f"Could not find dissociation energy for atoms with ids {ai.nr} and {aj.nr}. Using default."
        logger.debug(m)
        edis = 500

    mean_d = float(np.mean(ds))

    beta = calculate_beta(kb=kb, edis=edis)
    forces = calculate_forces(ds=ds, b0=b0, edis=edis, beta=beta)
    mean_f = float(np.mean(forces))
    harmonic_f = float(np.mean(calculate_harmonic_forces(ds=ds, b0=b0, k=kb)))

    return {
        "plumed_id": bond_id,
        "mean_d": mean_d,
        "mean_f": mean_f,
        "delta_d": mean_d - b0,
        "b0": b0,
        "harmonic_f": harmonic_f,
    }


# Residue names that are not part of a "Protein" compressed-x-grps output.
# NOTE: duplicated from kimmdy-hat (kimmdy_hat/reaction.py, which additionally
# excludes dummy atoms). Both could be merged into a shared util in the future.
SOL_RESNAMES = ["SOL", "WAT", "TIP3", "TIP4", "CL", "NA", "K", "CLA", "POT"]


def load_trajectory_universe(
    structure: Path,
    trajectory: Path,
    compressed_group: Optional[str] = None,
) -> mda.Universe:
    """Load a trajectory onto a structure, accounting for compressed-x-grps.

    xtc files may only contain the atoms of the `compressed-x-grps` group of the
    mdp file. Following kimmdy-hat, the group "System" (or no group) uses all
    atoms of the structure, while "Protein" and unknown groups use all atoms
    that are not solvent or ions.
    NOTE: duplicated from kimmdy-hat (kimmdy_hat/reaction.py), could be merged
    into a shared util in the future.

    Parameters
    ----------
    structure
        Structure file with all atoms of the system, e.g. a gro file.
    trajectory
        xtc or trr trajectory file.
    compressed_group
        Value of `compressed-x-grps` in the mdp file of the MD run.
        Only used for xtc files.

    Returns
    -------
    MDAnalysis Universe with the trajectory loaded. `u.atoms.ids` are the
    1-based atom numbers of the full system.
    """
    u = mda.Universe(str(structure))
    system_indices = u.atoms.indices
    if (
        trajectory.suffix == ".xtc"
        and compressed_group
        and str(compressed_group).lower() != "system"
    ):
        if str(compressed_group).lower() != "protein":
            logger.warning(
                f"Unknown compressed-x-grps {compressed_group}, "
                "assuming it contains all non-solvent atoms."
            )
        selection = u.select_atoms(f"not resname {' '.join(SOL_RESNAMES)}")
        system_indices = selection.indices
        u = mda.Merge(selection)

    try:
        u.load_new(str(trajectory))
    except ValueError as e:
        m = (
            f"Number of atoms in {trajectory.name} does not match the "
            f"{len(u.atoms)} atoms taken from {structure.name}. Check that "
            "compressed-x-grps is set to System or Protein in the mdp file "
            f"or use the trr trajectory instead. ({e})"
        )
        logger.error(m)
        raise ValueError(m) from e

    u.atoms.ids = system_indices + 1
    return u


def atom_nrs_from_selection(u: mda.Universe, selection: str) -> set[str]:
    """Atom numbers of an MDAnalysis selection in a Universe from
    `load_trajectory_universe`."""
    return {str(i) for i in u.select_atoms(selection).ids}


def calculate_bond_distances(
    u: mda.Universe, bonds: list[BondId], dt: float = 0.0
) -> tuple[np.ndarray, np.ndarray]:
    """Calculate bond distances over a trajectory.

    Distances use the minimum image convention.

    Parameters
    ----------
    u
        Universe from `load_trajectory_universe`.
    bonds
        Bonds as tuples of atom numbers.
    dt
        Minimum time in ps between used frames. 0 uses every frame.

    Returns
    -------
    Times in ps with shape (n_frames,) and distances in nm with shape
    (n_frames, n_bonds).
    """
    nr_to_index = {str(nr): i for i, nr in enumerate(u.atoms.ids)}
    missing = {nr for bond in bonds for nr in bond if nr not in nr_to_index}
    if missing:
        m = f"Atoms {sorted(missing, key=int)} of selected bonds are not in the trajectory."
        logger.error(m)
        raise ValueError(m)
    idx_i = np.array([nr_to_index[b[0]] for b in bonds], dtype=int)
    idx_j = np.array([nr_to_index[b[1]] for b in bonds], dtype=int)

    times = []
    distances = []
    t_last = None
    for ts in u.trajectory:
        # truncate to 1 fs to avoid floating point errors
        time = round(float(ts.time), 3)
        if dt > 0 and t_last is not None and time - t_last < dt - 1e-6:
            continue
        t_last = time
        times.append(time)
        # MDAnalysis works in Angstrom, GROMACS in nm
        ds = calc_bonds(ts.positions[idx_i], ts.positions[idx_j], box=ts.dimensions)
        distances.append(ds / 10.0)

    return np.asarray(times), np.asarray(distances).reshape(len(times), len(bonds))


def calculate_bondstats_from_trajectory(
    top: Topology,
    u: mda.Universe,
    bonds: list[BondId],
    dt: float = 0.0,
    edissoc_dat: Optional[Path] = None,
) -> BondStats:
    """Mean bond distances and forces from a trajectory.

    Parameters
    ----------
    top
        Topology with the bond parameters.
    u
        Universe from `load_trajectory_universe`.
    bonds
        Bonds as tuples of atom numbers, e.g. from `select_bonds`.
    dt
        Minimum time in ps between used frames. 0 uses every frame.
    edissoc_dat
        File with dissociation energies.
    """
    times, distances = calculate_bond_distances(u=u, bonds=bonds, dt=dt)
    if len(times) == 0:
        m = "No trajectory frames found for calculating bond statistics."
        logger.error(m)
        raise ValueError(m)
    edissoc = read_edissoc(edissoc_dat) if edissoc_dat is not None else None

    stats: BondStats = {}
    for i, bondkey in enumerate(bonds):
        stats[bondkey] = _stats_for_bond(
            top=top,
            bondkey=bondkey,
            bond_id=f"{bondkey[0]}-{bondkey[1]}",
            ds=distances[:, i],
            edissoc=edissoc,
        )
    return stats


def calculate_beta(kb: float, edis: float) -> float:
    """Calculate the beta parameter for the Morse potential.

    parameters are in gromacs units with kb as kJ/mol/nm^2 and edis in kJ/mol.
    """
    return np.sqrt(kb / (2 * edis))


def calculate_x_dagger(kb: float, edis: float) -> float:
    return np.sqrt(2 * edis / kb)


def calculate_harmonic_forces(ds: np.ndarray, b0: float, k: float) -> np.ndarray:
    """Calculate harmonic forces in a bond from a distances using the harmonic potential.

    Forces are returned in gromacs units kJ/mol/nm.
    The force is returned without a sign, i.e. it is positive for a stretched bond and negative for a compressed bond.
    """
    dds = ds - b0
    return k * dds


def harmonic_transition_rate(
    r_curr: list[float],
    r_0: float,
    dissociation_energy: float,
    k_f: float,
    frequency_factor: float = 0.288,
    temperature: float = 300,
) -> tuple[list[float], list[float]]:
    rs = np.asarray(r_curr)
    fs = calculate_harmonic_forces(ds=rs, b0=r_0, k=k_f)
    ks = harmonic_rates_from_forces(
        fs=fs,
        edis=dissociation_energy,
        kb=k_f,
        frequency_factor=frequency_factor,
        temperature=temperature,
    )
    return list(ks), list(fs)


def harmonic_rates_from_forces(
    fs: np.ndarray,
    edis: float,
    kb: float,
    frequency_factor: float = 0.288,
    temperature: float = 300,
) -> np.ndarray:
    """Calculate reaction rate constant from forces using the harmonic potential."""
    x_dagger = calculate_x_dagger(kb=kb, edis=edis)
    delta_v = edis - fs * x_dagger
    return frequency_factor * np.exp(-delta_v / (R * temperature))  # [1/ps]


def calculate_forces(
    ds: np.ndarray, b0: float, edis: float, beta: float, max_f: bool = True
) -> np.ndarray:
    """Calculate force in a bond from a distances using the Morse potential.

    Forces are returned in gromacs units kJ/mol/nm.
    """
    d_inflection = (beta * b0 + np.log(2)) / beta
    dds = ds - b0
    if max_f:
        # if the bond is stretched beyond the inflection point,
        # take the inflection point force because this force must have acted on the bond at some point
        ds_mask = ds > d_inflection
        dds[ds_mask] = d_inflection - b0
    return 2 * beta * edis * np.exp(-beta * dds) * (1 - np.exp(-beta * dds))


def morse_transition_rate(
    r_curr: list[float],
    r_0: float,
    dissociation_energy: float,
    k_f: float,
    frequency_factor: float = 0.288,
    temperature: float = 300,
) -> tuple[list[float], list[float]]:
    """Calculates reaction rate constant for a bond breaking event.

    Uses the Morse potential model for this calculation. For an array of bond distances of the same bond,
    first calculates the forces on the bond, then the minima and maxima of the shifted Morse potential
    to get an energy barrier and finally a reaction rate constant using the Arrhenius equation.
    For intramolecular reactions, the reaction rate constant is equal to the reaction rate.

    The calculation should be according to the derivation in the original KIMMDY paper: DOI: 10.1021/acs.jctc.9b00786

    Parameters
    ----------
    r_curr:
        Bond distances for a single bond, typically from a time series.
    r_0:
        Equilibrium bond length of the bond.
    dissociation energy:
        Dissociation energy of the bond.
    k_f:
        Spring constant of the bond.
    frequency_factor:
        Prefactor of the Arrhenius equation in [1/ps]. Default value from fitting averaged C_a - N data to gromacs data, see original KIMMDY paper
        Alternatively 1/2pi sqrt(k/m).
    temperature:
        Temperature for the Arrhenius equation in GROMACS units.
    """
    rs = np.asarray(r_curr)
    beta = calculate_beta(kb=k_f, edis=dissociation_energy)
    fs = calculate_forces(ds=rs, b0=r_0, edis=dissociation_energy, beta=beta)
    ks = morse_rates_from_forces(
        fs=fs,
        b0=r_0,
        edis=dissociation_energy,
        beta=beta,
        frequency_factor=frequency_factor,
        temperature=temperature,
    )

    return list(ks), list(fs)


def morse_rates_from_forces(
    fs: np.ndarray,
    b0: float,
    edis: float,
    beta: float,
    frequency_factor: float = 0.288,
    temperature: float = 300,
) -> np.ndarray:
    # calculate extrema of shifted potential i.e. get barrier height of V_eff = V_morse - F*X
    s = np.sqrt(
        (beta**2 * edis**2 - 2 * edis * beta * fs)
        + 1e-7  # prevent rounding issue close to zero
    )
    r_min = b0 - 1 / beta * np.log((beta * edis + s) / (2 * beta * edis))
    r_max = b0 - 1 / beta * np.log((beta * edis - s) / (2 * beta * edis))
    r_max = np.where(
        ~np.isfinite(r_max), 10 * b0, r_max
    )  # set rmax to r0 * 10 where no rmax can be found

    v_max = edis * (1 - np.exp(-beta * (r_max - b0))) ** 2 - fs * (r_max - b0)
    v_min = edis * (1 - np.exp(-beta * (r_min - b0))) ** 2 - fs * (r_min - b0)
    # Note: F*r should lead to same result as F*(r-r_0) since the shifts in Vmax-Vmin adds up to zero
    delta_v = v_max - v_min

    # calculate reaction rate constant from barrier height
    return frequency_factor * np.exp(-delta_v / (R * temperature))  # [1/ps]


def get_atomnrs_from_plumedid(
    plumedid: str,
    plumed: Plumed_dict,
) -> list[str]:
    """
    Convert from plumedid to atomnr, information from the plumed file is used.

    Parameters
    ----------
    plumedid:
        Identifier from a plumed input file (e.g d0).
    plumed:
        Parsed plumed input file
    """
    # lookup_atomnr_plumedid = {k: frozenset(v["atoms"])
    plumed_action = plumed["labeled_action"][plumedid]
    if a := plumed_action.get("atoms"):
        atomnrs = sorted(a, key=int)
        return atomnrs
    else:
        raise NotImplementedError(
            f"Can't get atomnrs for {plumedid}, is this for an unexpected plumed action?"
        )


def get_atominfo_from_atomnrs(
    atomnrs: (
        list[str]
        | tuple[str]
        | tuple[str, str]
        | tuple[str, str, str]
        | tuple[str, str, str, str]
    ),
    top: Topology,
) -> tuple[list[str], list[str]]:
    """Use topology atoms section to convert from atomnr to atomtype"""
    atomtypes = []
    atomnames = []
    for atomnr in atomnrs:
        atomtypes.append(top.atoms[atomnr].type)
        atomnames.append(top.atoms[atomnr].atom)
    return atomtypes, atomnames


def get_bondprm_from_atomtypes(
    atomtypes: list[str],
    ffbonded: dict,
) -> tuple[float, float]:
    """Returns bond parameters (b0, kb) for a set of two atomtypes.

    Parameters
    ----------
    atomtypes:
        Two atomtypes as defined in the respective force field
    ffbonded:
        Force field ffbonded.itp file parsed through the rtp parser
    """
    # search for b0 and kb for the given atomtypes in ffbonded bondtypes
    for bondtype in ffbonded["bondtypes"]["content"]:
        if set(atomtypes) == set(bondtype[:2]):
            b0, kb = [float(x) for x in bondtype[3:5]]
            break
    else:
        raise KeyError(
            f"Did not find bond parameters for atomtypes {atomtypes} in ffbonded file"
        )

    return b0, kb


def get_edissoc_from_atomnames(
    atomnames: list[str], edissoc: dict, residue: str = "_"
) -> float:
    """Returns dissociation energy E_dissoc for a set of two atomnames.

    Parameters
    ----------
    atomnames:
        Two atomnames as defined in the respective force field
    edissoc:
        Parsed file with dissociation energies per bond between two atomtypes or elements
    residue:
        Residue for which the atomnames are defined
    """
    if residue not in edissoc.keys():
        if "general" in edissoc.keys():
            logger.debug(
                f"residue {residue} not in edissoc keys: {edissoc.keys()}, using 'general' as residue."
            )
            residue = "general"
        else:
            raise KeyError(f"Did not find residue {residue} in edissoc file")

    try:
        interaction_key = tuple(sorted(atomnames))
        E_dis = edissoc[residue][interaction_key]
    except KeyError:
        # continue with guessed edissoc
        logger.warning(
            f"Did not find dissociation energy for atomtypes {atomnames}, residue {residue} in edissoc file, using standard value of 400.0"
        )
        E_dis = 400.0

    return E_dis
