import logging

from kimmdy import plugin_utils
from kimmdy.parsing import read_top
from kimmdy.topology.topology import Topology
from pathlib import Path
import MDAnalysis as mda
from MDAnalysis.lib.distances import calc_bonds
import numpy as np
import pytest

# backbone bonds of hexalanine_homolysis as in its former plumed.dat
# (backbone minus "C-N H* O*", the kimmdy-create-plumed default)
HEXALA_BACKBONE_BONDS = {
    ("7", "9"),
    ("9", "15"),
    ("17", "19"),
    ("19", "25"),
    ("27", "29"),
    ("29", "35"),
    ("37", "39"),
    ("39", "45"),
    ("47", "49"),
    ("49", "55"),
    ("57", "59"),
    ("59", "65"),
}


def test_bondstate_csv_io(arranged_tmp_path):
    bondstat = plugin_utils.bondstats_from_csv(".kimmdy.bondstats")
    assert isinstance(bondstat, dict)
    assert all(
        [
            k in d.keys()
            for k in ["b0", "delta_d", "harmonic_f", "mean_d", "mean_f", "plumed_id"]
            for d in bondstat.values()
        ]
    )

    out = arranged_tmp_path / "bondstats_out"
    plugin_utils.bondstats_to_csv(bondstat, out)

    with open(out, "r") as f:
        written_lines = f.readlines()
    with open(".kimmdy.bondstats", "r") as f:
        original_lines = f.readlines()
    assert all([w == r for w, r in zip(written_lines, original_lines)])


def test_bondstats_csv_without_harmonic_force(tmp_path):
    path = tmp_path / "bondstats_7col"
    path.write_text(
        "ai,aj,plumed_id,mean_d,mean_f,delta_d,b0,harmonic_f\n"
        "7,9,d0,0.146791,422.989854,0.001891,0.144900\n"
    )
    bondstats = plugin_utils.bondstats_from_csv(path)
    assert bondstats[("7", "9")]["harmonic_f"] is None
    assert bondstats[("7", "9")]["mean_d"] == pytest.approx(0.146791)

    out = tmp_path / "bondstats_out"
    plugin_utils.bondstats_to_csv(bondstats, out)
    assert out.read_text() == path.read_text().rstrip("\n")


@pytest.mark.parametrize(
    "arranged_tmp_path", ["test_integration/hexalanine_homolysis"], indirect=True
)
def test_select_bonds_matches_create_plumed_defaults(arranged_tmp_path):
    top = Topology(read_top(Path("hexala_out.top")))
    backbone = {
        nr for nr, atom in top.atoms.items() if atom.atom in {"N", "CA", "C", "O"}
    }
    bonds = plugin_utils.select_bonds(
        top=top, atom_nrs=backbone, exclude=["C-N", "H*", "O*"]
    )
    assert set(bonds) == HEXALA_BACKBONE_BONDS


@pytest.mark.parametrize(
    "arranged_tmp_path", ["test_integration/hexalanine_homolysis"], indirect=True
)
def test_select_bonds_exclusions(arranged_tmp_path):
    top = Topology(read_top(Path("hexala_out.top")))
    assert len(plugin_utils.select_bonds(top=top)) == len(top.bonds)

    bonds = plugin_utils.select_bonds(top=top, exclude=["H*", "O"])
    names = [{top.atoms[a].atom, top.atoms[b].atom} for a, b in bonds]
    assert names
    assert not any(n.startswith("H") or n == "O" for pair in names for n in pair)
    assert {"C", "N"} in names

    bonds = plugin_utils.select_bonds(top=top, exclude=["C-N"])
    names = [{top.atoms[a].atom, top.atoms[b].atom} for a, b in bonds]
    assert {"C", "N"} not in names


def _write_trajectory(tmp_path, frames, resnames, xtc_atoms=None, dt=2.0):
    """Write a gro of all atoms and an xtc (optionally of a subset of atoms)."""
    n = len(resnames)
    box = [30.0, 30.0, 30.0, 90.0, 90.0, 90.0]
    u = mda.Universe.empty(
        n, n_residues=n, atom_resindex=list(range(n)), trajectory=True
    )
    u.add_TopologyAttr("names", [f"A{i}" for i in range(n)])
    u.add_TopologyAttr("resnames", list(resnames))
    u.add_TopologyAttr("resids", list(range(1, n + 1)))
    u.dimensions = box
    u.atoms.positions = frames[0]
    gro = tmp_path / "structure.gro"
    u.atoms.write(str(gro))

    u.load_new(np.asarray(frames, dtype=np.float32), order="fac", dimensions=box, dt=dt)
    atoms = u.atoms if xtc_atoms is None else u.atoms[xtc_atoms]
    xtc = tmp_path / "traj.xtc"
    with mda.Writer(str(xtc), atoms.n_atoms) as w:
        for _ in u.trajectory:
            w.write(atoms)
    return gro, xtc


def _pbc_frames():
    # atoms 1 and 2 are bonded across the periodic boundary in x
    # (1.0, 1.2 and 1.4 Angstrom apart), atoms 3 and 4 are solvent
    frames = []
    for d in [1.0, 1.2, 1.4]:
        frames.append(
            [
                [0.5, 5.0, 5.0],
                [30.5 - d, 5.0, 5.0],
                [10.0, 10.0, 10.0],
                [11.0, 10.0, 10.0],
            ]
        )
    return frames


def test_eyring_prefactor():
    assert plugin_utils.eyring_prefactor(300) == pytest.approx(6.251, rel=1e-3)
    assert plugin_utils.eyring_prefactor(300, 0.5) == pytest.approx(3.1255, rel=1e-3)


def test_bond_distances_from_trajectory_use_minimum_image(tmp_path):
    gro, xtc = _write_trajectory(tmp_path, _pbc_frames(), ["ALA", "ALA", "SOL", "SOL"])
    u = plugin_utils.load_trajectory_universe(gro, xtc)
    times, distances = plugin_utils.calculate_bond_distances(u, [("1", "2")])
    assert times == pytest.approx([0.0, 2.0, 4.0])
    assert distances.shape == (3, 1)
    assert distances[:, 0] == pytest.approx([0.10, 0.12, 0.14], abs=1e-3)


def test_bond_distances_from_trajectory_stride(tmp_path):
    gro, xtc = _write_trajectory(tmp_path, _pbc_frames(), ["ALA", "ALA", "SOL", "SOL"])
    u = plugin_utils.load_trajectory_universe(gro, xtc)
    times, distances = plugin_utils.calculate_bond_distances(u, [("1", "2")], stride=2)
    assert times == pytest.approx([0.0, 4.0])
    assert distances[:, 0] == pytest.approx([0.10, 0.14], abs=1e-3)


def test_bond_distances_from_trajectory_logs_frames(tmp_path, caplog, monkeypatch):
    monkeypatch.setattr(logging.getLogger("kimmdy"), "propagate", True)
    caplog.set_level(logging.INFO, logger="kimmdy.plugin_utils")
    gro, xtc = _write_trajectory(tmp_path, _pbc_frames(), ["ALA", "ALA", "SOL", "SOL"])
    u = plugin_utils.load_trajectory_universe(gro, xtc)
    plugin_utils.calculate_bond_distances(u, [("1", "2")], stride=2)
    assert "Analysed 2 of 3 trajectory frames (0 to 4 ps, stride 2)" in caplog.text

    with pytest.raises(ValueError, match="stride must be at least 1"):
        plugin_utils.calculate_bond_distances(u, [("1", "2")], stride=0)


def test_trajectory_with_protein_compressed_group(tmp_path):
    gro, xtc = _write_trajectory(
        tmp_path, _pbc_frames(), ["ALA", "ALA", "SOL", "SOL"], xtc_atoms=[0, 1]
    )
    u = plugin_utils.load_trajectory_universe(gro, xtc, compressed_group="Protein")
    assert list(u.atoms.ids) == [1, 2]
    _, distances = plugin_utils.calculate_bond_distances(u, [("1", "2")])
    assert distances[:, 0] == pytest.approx([0.10, 0.12, 0.14], abs=1e-3)

    with pytest.raises(ValueError, match="not in the trajectory"):
        plugin_utils.calculate_bond_distances(u, [("2", "3")])

    for group in [None, "System"]:
        with pytest.raises(ValueError, match="compressed-x-grps"):
            plugin_utils.load_trajectory_universe(gro, xtc, compressed_group=group)


@pytest.mark.parametrize(
    "arranged_tmp_path", ["test_integration/hexalanine_homolysis"], indirect=True
)
def test_bondstats_from_trajectory_hexala_backbone(arranged_tmp_path):
    top = Topology(read_top(Path("hexala_out.top")))
    gro = Path("npt.gro")
    structure = mda.Universe(str(gro))
    positions = structure.atoms.positions
    xtc = arranged_tmp_path / "traj.xtc"
    with mda.Writer(str(xtc), structure.atoms.n_atoms) as w:
        for _ in range(2):
            w.write(structure.atoms)

    u = plugin_utils.load_trajectory_universe(gro, xtc, compressed_group="System")
    atom_nrs = plugin_utils.atom_nrs_from_selection(u, "backbone")
    bonds = plugin_utils.select_bonds(
        top=top, atom_nrs=atom_nrs, exclude=["C-N", "H*", "O*"]
    )
    assert set(bonds) == HEXALA_BACKBONE_BONDS

    stats = plugin_utils.calculate_bondstats_from_trajectory(top, u, bonds)
    assert set(stats.keys()) == HEXALA_BACKBONE_BONDS
    for (ai, aj), s in stats.items():
        expected = calc_bonds(
            positions[int(ai) - 1], positions[int(aj) - 1], box=structure.dimensions
        )
        assert s["mean_d"] == pytest.approx(expected / 10, abs=1e-3)
        assert s["plumed_id"] == f"{ai}-{aj}"
        assert s["delta_d"] == pytest.approx(s["mean_d"] - s["b0"])
