"""Reaction plugin building blocks"""

import logging
from pathlib import Path
from typing import Optional, TypeAlias, TypedDict

import numpy as np

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
            if len(line) == 8:
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

        ds = np.asarray([values[plumed_id] for values in distances.values()])
        mean_d = float(np.mean(ds))

        beta = calculate_beta(kb=kb, edis=edis)
        forces = calculate_forces(ds=ds, b0=b0, edis=edis, beta=beta)
        mean_f = float(np.mean(forces))
        harmonic_f = float(np.mean(calculate_harmonic_forces(ds=ds, b0=b0, k=kb)))

        stats[bondkey] = {
            "plumed_id": plumed_id,
            "mean_d": mean_d,
            "mean_f": mean_f,
            "delta_d": mean_d - b0,
            "b0": b0,
            "harmonic_f": harmonic_f,
        }
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
