# %%
from itertools import combinations
from pathlib import Path
from re import A
from typing import Union

import MDAnalysis as MDA
import numpy as np
import random
from scipy.spatial.transform import Rotation
from tqdm.autonotebook import tqdm

# %% Functions

def check_cylinderempty(a, b, t, r_min=0.8, verbose=False):
    """Checks for atoms in a cylinder around a possible HAT reaction path.
    Cylinder begins/ends 10% after/before the end points.
    Ref.: https://geomalgorithms.com/a02-_lines.html

    Parameters
    ----------
    a : Union[np.ndarray, list]
        Center of cylinder base
    b : Union[np.ndarray, list]
        Center of cylinder top
    t : Union[np.ndarray, list]
        Testpoint, or list of testpoints
    r_min : int, optional
        Radius of cylinder, by default 0.8

    Returns
    -------
    bool
        True if no points are inside the cylinder, False otherwise.
    """

    def _check_point(a, b, t, r_min, verbose):
        v = b - a  # apex to base
        w = t - a  # apex to testpoint
        c1 = np.dot(v, w)
        c2 = np.dot(v, v)
        x = c1 / c2  # percentage along axis

        tx = a + (x * v)  # projection of t on x
        r = np.linalg.norm(t - tx)

        if x < 0.1 or x > 0.9 or r > r_min:
            # point is outside cylinder
            return True
        if verbose:
            print(f"a = {a}, b = {b}, t = {t}")
            print(f"x = {x}, r = {r}")
        return False

    if isinstance(t, np.ndarray) and t.shape == (3,):
        return _check_point(a, b, t, r_min, verbose)

    elif (
        (isinstance(t, np.ndarray) and t.shape[1] == 3)
        or isinstance(t, list)
        and len(t[0]) == 3
    ):
        return all(
            [
                _check_point(np.array(a), np.array(b), np.array(p), r_min, verbose)
                for p in t
            ]
        )

    else:
        raise ValueError(f"Type/Shape of t wrong: {t}")


def find_radical_pos(
    center: MDA.core.groups.Atom, bonded: MDA.core.groups.AtomGroup, tetrahedral=False
):
    """Calculates possible radical positions of a given radical atom

    Parameters
    ----------
    center : MDA.core.groups.Atom
        Radical atom
    bonded : MDA.core.groups.AtomGroup
        Atom group of bonded atoms. From its length the geometry is predicted.
    tetrahedral : bool
        Whether to assume a tetrahedral conformation around C and N

    Returns
    -------
    list
        List of radical positions, three dimensional arrays
    """
    scale_C = 1.10
    scale_N = 1.04
    scale_O = 0.97
    scale_S = 1.41

    if center.element == "C" and len(bonded) == 2:
        c_alphas = center.bonded_atoms
        c_o = c_alphas[np.nonzero(c_alphas.elements == "O")]
        assert len(c_o) in [0, 1], "Carboxyl radical?"
        if len(c_o) > 0 and len(c_o[0].bonded_atoms) > 1:
            tetrahedral = True

    if tetrahedral:
        # MAKE SP3 from C with 2 bonds
        # invert bonded positions
        inv_1 = (center.position - bonded[0].position) + center.position
        inv_2 = (center.position - bonded[1].position) + center.position
        # construct rotation axis
        midpoint = (inv_2 - inv_1) * 0.5 + inv_1
        rot_ax = midpoint - center.position
        # 90 degree rotation
        r = Rotation.from_rotvec((np.pi / 2) * (rot_ax / np.linalg.norm(rot_ax)))
        # rotated bonds relative to center
        rad_1 = r.apply(inv_1 - center.position)
        rad_2 = r.apply(inv_2 - center.position)

        # scale to correct bond length, make absolute
        if center.element == "N":
            scale = scale_N
        elif center.element == "C":
            scale = scale_C
        else:
            raise NotImplementedError("H Bondlength to central atom missing")

        rad_1 = (rad_1 / np.linalg.norm(rad_1)) * scale + center.position
        rad_2 = (rad_2 / np.linalg.norm(rad_2)) * scale + center.position
        return [rad_1, rad_2]

    if len(bonded) in [2, 3]:
        assert center.element in ["C", "N"], "Element element does not match bond number"

        # prediction: inverse midpoint between bonded
        # scale to correct bond length, make absolute
        if center.element == "N":
            scale = scale_N
        elif center.element == "C":
            scale = scale_C

        b_normed = []
        for b in bonded:
            b_vec = b.position - center.position
            b_vec_norm = b_vec / np.linalg.norm(b_vec)
            b_normed.append(b_vec_norm)

        midpoint = sum(b_normed)

        v = midpoint / np.linalg.norm(midpoint)
        rad = center.position + (-1 * v * scale)
        return [rad]

    # Radicals w/ only one bond:
    elif len(bonded) == 1:
        assert center.element in ["O", "S"], "Element element does not match bond number"
        if center.element == "O":
            scale = scale_O
        elif center.element == "S":
            scale = scale_S

        b = bonded[0]
        b_vec = b.position - center.position
        b_vec = b_vec / np.linalg.norm(b_vec)
        rnd_vec = [1, 1, 1]  # to find a vector perpendicular to b_vec

        rnd_rot_ax = np.cross(b_vec, rnd_vec)
        rnd_rot_ax = rnd_rot_ax / np.linalg.norm(rnd_rot_ax)

        r1 = Rotation.from_rotvec(1.911 * rnd_rot_ax)  # 109.5 degree
        r2 = Rotation.from_rotvec(0.785 * b_vec)  # 45 degree

        ends = [r1.apply(b_vec)]  # up to 109.5

        for i in range(8):
            ends.append(r2.apply(ends[-1]))  # turn in 45d steps

        # norm and vec --> position
        ends = [(e / np.linalg.norm(e)) * scale + center.position for e in ends]

        return ends

    else:
        raise ValueError(f"Weired count of bonds: {list(bonded)}\nCorrect radicals?")

def get_residue(atm):
    """Builds iteratively an atomgroup containing the residue of the given atom.
    Relies on the resid attribute.

    Parameters
    ----------
    atm : MDA.core.groups.Atom
        Starting atom

    Returns
    -------
    MDA.AtomGroup
        Residue
    """
    resid = atm.resid
    atm_env = atm.universe.select_atoms(f"point { str(atm.position).strip('[ ]') } 15")
    to_check = atm_env.select_atoms(f"bonded index {atm.ix}")
    checked = atm_env.select_atoms(f"index {atm.ix}")
    resid_group = atm_env.select_atoms(f"index {atm.ix}")

    while len(to_check) > 0:
        for c_atm in to_check:
            if c_atm.resid == resid:
                resid_group = resid_group + c_atm
                to_check = to_check + atm_env.select_atoms(f"bonded index {c_atm.ix}")
            checked = checked + c_atm
            to_check = to_check - checked
    assert (
        len(np.nonzero(resid_group.names == "CA")) == 1
    ), "ERROR: Multiple CA found in one residue!"
    return resid_group

def get_res_union(atms):
    """Builds atomgroupe containing union of all residues of given atoms.
    Avoids unnecessary calls to get_residue

    Parameters
    ----------
    atms : MDA.core.groups.Atom
        Atoms which residues will be unionized.
    """
    res = atms[0].universe.atoms[[]]
    for atm in atms:
        if atm not in res:
            res = res | get_residue(atm)
    return res

def cap_aa(atms):
    """Caps aminoacid(s) with amin or amid group.

    Parameters
    ----------
    atms : MDA.AtomGroup
        All atoms of the aminoacid(s) to cap.

    Returns
    -------
    MDA.AtomGroup
        Capping atoms. Should be used as union with the
        aminoacid residue as capping occures outside this residue.
    List[int]
        Indices of atoms used as cappes in the old universe.
    """

    possible_ends = atms.select_atoms(f"name C or name N")

    # Handle special cases aka. crosslinks
    for res in filter(lambda r: r.resname == "L5Y", atms.residues):
        possible_ends += res.atoms.select_atoms(f"name NZ")
    for res in filter(lambda r: r.resname == "L4Y", atms.residues):
        possible_ends += res.atoms.select_atoms(f"name CE")

    assert len(possible_ends) > 0, "ERROR: No possible ends to cap found!"
    env = atms.universe.select_atoms(
        f"point { str(possible_ends[0].position).strip('[ ]') } 20"
    )
    cap_atms = []
    cap_ix = []
    cap = MDA.Universe.empty(0).atoms
    return cap, cap_ix              # dont need caps at the moment

    h_bondlength = 1.01

    # print(f"!! possible ends: {possible_ends}")
    for pe in possible_ends:
        if pe.element == "C":
            # build cap from atoms of next AA
            cap_d = {
                "N": env.select_atoms(
                    f"(bonded index {pe.ix}) and not resid {pe.resid}"
                )
            }

            # check cap selection
            assert (
                len(cap_d["N"]) <= 1
            ), f"ERROR: Wrong cap atom selection at index {pe.ix}, cap: {list(cap_d['N'])}"
            if len(cap_d["N"]) == 0:
                continue  # chain end reached (on radical?)
            if cap_d["N"].residues in atms.residues:
                continue  # next aminoacid included, no need to cap
            if cap_d["N"].resnames[0] == "NME":
                continue  # already capped

            # Special treatment if single AA is in between:
            # Build Gly linker
            linker = False
            if (cap_d["N"][0].resindex - pe.resindex) + cap_d["N"][
                0
            ].resindex in atms.resindices:
                if cap_d["N"][0].resname not in ["L4Y", "L5Y"]:
                    linker = True

            N_alphas = env.select_atoms(
                f"(bonded index {cap_d['N'][0].ix}) and (resid {cap_d['N'][0].resid})"
            )

            # C_a --> CH3,
            # everything w/ more or less than 1 H attached --> H

            if "H" in N_alphas.elements:
                cap_d["N_C"] = N_alphas[np.nonzero(N_alphas.elements == "C")[0]][0]
                cap_d["N_H"] = N_alphas[np.nonzero(N_alphas.elements == "H")[0]][0]
            else:  # two C atoms bond to N, only in PRO, HYP
                cap_d["N_C"] = N_alphas[np.nonzero(N_alphas.names == "CA")[0]][0]
                cap_d["N_H"] = N_alphas[np.nonzero(N_alphas.names != "CA")[0]][0]

            assert all(
                [k in cap_d.keys() for k in ["N_C", "N_H"]]
            ), f"ERROR while building capping group on C-term!\nAtom:{cap_d['N'][0]}"

            if not linker:
                cap_d["NC_3H"] = cap_d["N_C"].bonded_atoms - cap_d["N"]
                assert len(cap_d["NC_3H"]) == 3, f"CAPPING ERROR: Atom {cap_d['N'][0]}"

                cap = sum([cap_d[k] for k in ["N", "N_H", "N_C", "NC_3H"]])
                h_idxs = (1, 3, 4, 5)

                # make new universe with fixed order
                u_cap = MDA.core.universe.Merge(cap)
                u_cap.residues[0].resname = "NME"

            else:  # LINKER:
                cap_d["bb"] = cap_d["N"][0].residue.atoms.select_atoms("backbone")
                assert (
                    len(cap_d["bb"]) == 4
                ), f"CAPPING ERROR in linker: backbone {cap_d['N'][0]}"
                cap_d["NC_2H"] = cap_d["N_C"].bonded_atoms - cap_d["bb"]

                cap = sum([cap_d[k] for k in ["N_H", "NC_2H", "bb"]])
                h_idxs = (0, 1, 2)

                # make new universe with fixed order
                u_cap = MDA.core.universe.Merge(cap)
                u_cap.residues[0].resname = "GLY"

            # scale positions and mutate elements to H
            for h in [u_cap.atoms[i] for i in h_idxs]:
                h_alpha = h.bonded_atoms[0]
                bond = h.position - h_alpha.position
                h.position = h_alpha.position + (
                    (bond / np.linalg.norm(bond)) * h_bondlength
                )

                h.element = "H"
                h.name = "H"
                h.element = "H"

            [cap_ix.append(i) for i in cap.ix]
            cap_atms.append(u_cap.atoms)

        if pe.element == "N":
            # build cap from atoms of next AA
            cap_d = {
                "C": env.select_atoms(
                    f"(bonded index {pe.ix}) and not resid {pe.resid}"
                )
            }

            # check cap selection
            assert (
                len(cap_d["C"]) <= 1
            ), f"ERROR: Wrong cap atom selection at index {pe.ix}, cap: {list(cap_d['C'])}"
            if len(cap_d["C"]) == 0:
                continue  # chain end reached
            if cap_d["C"].residues in atms.residues:
                continue  # next aminoacid included, no need to cap
            if cap_d["C"].resnames[0] == "ACE":
                continue  # already capped

            # skip if linker AA, always build from C term
            if (cap_d["C"][0].resindex - pe.resindex) + cap_d["C"][
                0
            ].resindex in atms.resindices:
                if cap_d["C"][0].resname not in ["L4Y", "L5Y"]:
                    continue

            C_alphas = env.select_atoms(
                f"(bonded index {cap_d['C'][0].ix}) and (resid {cap_d['C'][0].resid})"
            )
            assert (
                len(C_alphas) == 2
            ), f"ERROR while building capping group on N-term!\nAtom:{cap_d['C'][0]}"

            cap_d["O"] = filter(lambda a: a.element == "O", C_alphas).__next__()
            cap_d["CC"] = (C_alphas - cap_d["O"])[0]

            cap_d["CC_H3"] = cap_d["CC"].bonded_atoms - cap_d["C"]

            cap = sum([cap_d[k] for k in ["C", "O", "CC", "CC_H3"]])
            [cap_ix.append(i) for i in cap.ix]

            # make new universe with fixed order
            u_cap = MDA.core.universe.Merge(cap)
            u_cap.residues[0].resname = "ACE"

            # scale positions and mutate elements
            for h in [u_cap.atoms[i] for i in (3, 4, 5)]:
                h_alpha = h.bonded_atoms[0]
                bond = h.position - h_alpha.position
                h.position = h_alpha.position + (
                    (bond / np.linalg.norm(bond)) * h_bondlength
                )

                h.element = "H"
                h.name = "H"
                h.element = "H"

            cap_atms.append(u_cap.atoms)

    if len(cap_atms) == 0:
        cap = MDA.Universe.empty(0).atoms
    else:
        cap = MDA.core.universe.Merge(*cap_atms).atoms

    #print(f"cap: {cap},{len(cap_atms)}")
    return cap, cap_ix

def _get_charge(atm):
    aa_charge_dict = {
        "ala": 0,
        "arg": 1,
        "asn": 0,
        "asp": -1,
        "cys": 0,
        "dop": 0,  # ?
        "gln": 0,
        "glu": -1,
        "gly": 0,
        "his": 0,  # ?
        "hyp": 0,
        "ile": 0,
        "leu": 0,
        "lys": 1,
        "met": 0,
        "phe": 0,
        "pro": 0,
        "ser": 0,
        "thr": 0,
        "trp": 0,
        "tyr": 0,
        "val": 0,
        "l4y": 0,
        "l5y": 0,
        "ace": 0,
        "nme": 0,
    }

    return aa_charge_dict[atm.resname.lower()]


def cap_single_rad(u, ts, rad, bonded_rad, h_cutoff=3, env_cutoff=7):
    """Builds capped systems around a single radical in a single frame.
    Aminoacids are capped at peptide bonds resulting in amines and amides.
    Subsystems contain the reactive hydrogen at index 0 followed by the
    radical atom.

    Parameters
    ----------
    u : MDA.Universe
        Main universe
    ts : MDAnalysis.coordinates.base.Timestep
        On which timestep to operate
    rad : MDA.AtomGroup
        Radical atom in its own group
    bonded_rad : MDA.AtomGroup
        AtomGroup containing all atoms bonded to the radical
    h_cutoff : float, optional
        Cutoff radius for hydrogen search around radical, by default 3
    env_cutoff : float, optional
        Cutoff radius for local env used for better performance, by default 7

    Returns
    -------
    List
        List of capped subsystems
    """
    # selecting in a smaller env is faster than in whole universe
    #print('Starting cap_single_rad')
    env = u.atoms.select_atoms(
        f"point { str(rad.positions).strip('[ ]') } {env_cutoff}"
    )

    end_poss = find_radical_pos(rad[0], bonded_rad)
    #print(f" end pos {end_poss}")

    hs = []
    for end_pos in end_poss:
        hs.append(
            env.select_atoms(
                f"point { str(end_pos).strip('[ ]') } {h_cutoff} and element H and not resname SOL"
            )
        )
    hs = sum(hs)
    #print(f"!!!! hs: {hs}")

    # hs = (
    #     env.select_atoms(
    #         f"point { str(rad.positions).strip('[ ]') } {h_cutoff} and element H"
    #     )
    #     - bonded_rad
    # )

    clashes = np.empty((len(hs), len(end_poss)), dtype=bool)
    for h_idx, h in enumerate(hs):
        for end_idx, end_pos in enumerate(end_poss):
            clashes[h_idx, end_idx] = check_cylinderempty(
                end_pos, h.position, env.positions, r_min=0.8
            )
    #print(f"clashes: {clashes}")
    # get whole residues near radical
    rad_alphas = env.select_atoms(f"bonded index {rad[0].ix}")
    # rad_betas = sum([env.select_atoms(f'bonded index {alpha.ix}') for alpha in rad_alphas]) - rad

    rad_aa = get_res_union(rad_alphas)

    capped_systems = np.zeros((len(hs),), dtype=object)
    min_translations = np.ones((len(hs),)) * 99

    # iterate over defined HAT reactions
    for h_idx, end_idx in zip(*np.nonzero(clashes)):
        end_pos = end_poss[end_idx]
        h = env.select_atoms(f"index {hs[h_idx].ix}")

        translation = np.linalg.norm(end_pos - h.positions)
        # only keep reaction w/ smallest translation
        if translation > min_translations[h_idx]:
            continue
        min_translations[h_idx] = translation

        # get whole residues near reacting H
        h_alpha = env.select_atoms(f"bonded index {h[0].ix}")[0]
        h_betas = sum(env.select_atoms(f"bonded index {h_alpha.ix}")) - h
        # h_gammas = sum(env.select_atoms(f'bonded index {" ".join([str(i) for i in h_betas.ix])}')) - h_alpha

        h_aa = get_res_union(h_betas)

        core = h_aa | rad_aa
        caps, caps_ix = cap_aa(core)

        # core can have more residues than just h-res and rad-res, important for charge!
        core = core - h - rad

        # N terminal end capped w/ RNH3+, fix charge:
        charge_correction = 0
        if "OC1" in core.names:  # OC1 and OC2 form COO- end
            charge_correction += -1
        if "H1" in core.names or "H2" in core.names:  # H1 H2 H3 form NH3+ end
            charge_correction += 1

        #print(h,rad,caps,core)
        ags = []
        for ag in [h,rad,caps,core]:
            if len(ag) > 0:
                ags.append(ag)
        #print(f"!!! {ags}")
        capped_systems[h_idx] = {
            "start_u": MDA.core.universe.Merge(*ags),
            "end_u": MDA.core.universe.Merge(*ags),
            "meta": {
                "translation": translation,
                "u1_name": rad[0].resname.lower() + "-sim",
                "u2_name": h[0].resname.lower() + "-sim",
                "charge_u1": _get_charge(rad[0]),
                "charge_u2": _get_charge(h[0]),
                "trajectory": u._trajectory.filename,
                "frame": ts.frame,
                "indices": (*h.ix, *rad.ix, *caps_ix, *h_aa.ix, *rad_aa.ix),
                "intramol": rad[0].residue == h[0].residue,
                "charge": sum([_get_charge(res.atoms[0]) for res in core.residues])
                + charge_correction,
            },
        }

        # change H position in end universe
        capped_systems[h_idx]["end_u"].atoms[0].position = end_pos

        # hashes based on systems rather than subgroups, subgroubs would collide
        capped_systems[h_idx]["meta"]["hash_u1"] = abs(
            hash(capped_systems[h_idx]["start_u"])
        )
        capped_systems[h_idx]["meta"]["hash_u2"] = abs(
            hash(capped_systems[h_idx]["end_u"])
        )

    return capped_systems[np.nonzero(capped_systems)[0]]

def find_radical(u):
    """
    assumes there is a radical
    """
    nbonds = {'H':1,'HC':1,'H1':1,'O':1,'N':3,'C':3,'CT':4}
    for atom in u.atoms:
        if len(atom.bonded_atoms) < nbonds[atom.type]:
            return MDA.AtomGroup([atom])  