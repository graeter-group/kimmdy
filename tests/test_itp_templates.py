from pathlib import Path

import pytest

from kimmdy.parsing import read_top
from kimmdy.topology.topology import Topology


@pytest.fixture
def dopa_topology() -> Topology:
    filedir = Path(__file__).parent / "test_files" / "test_itp_templates"
    raw_top = read_top(filedir / "topol.top", ffdir=filedir)
    return Topology(
        raw_top,
        radicals="92",
        product_template_sources=("itp",),
        product_template_itp_glob="*.itp",
    )


def _atoms_in_residue(topology: Topology, resnr: str):
    return [atom for atom in topology.atoms.values() if atom.resnr == resnr]


def _hat_atoms(topology: Topology):
    """Return donor H37, its O2 donor, and radical O3 from the fixture."""
    donor_h = next(
        atom
        for atom in topology.atoms.values()
        if atom.atom == "H37"
        and any(
            candidate.atom == "H38"
            for candidate in _atoms_in_residue(topology, atom.resnr)
        )
    )
    donor = topology.atoms[donor_h.bound_to_nrs[0]]
    acceptor = next(
        atom
        for atom in topology.atoms.values()
        if atom.atom == "O3"
        and not any(
            candidate.atom == "H38"
            for candidate in _atoms_in_residue(topology, atom.resnr)
        )
    )
    return donor_h, donor, acceptor


def test_loads_closed_and_radical_dopa_templates(dopa_topology):
    templates = dopa_topology.ff.product_templates.templates["DOPA"]
    templates_by_moleculetype = {
        template.moleculetype: template for template in templates
    }

    assert set(templates_by_moleculetype) == {"DOPA", "DOPA_rad"}

    closed = templates_by_moleculetype["DOPA"]
    radical = templates_by_moleculetype["DOPA_rad"]

    assert set(closed.atoms) - set(radical.atoms) == {"H38"}
    assert set(radical.atoms) < set(closed.atoms)
    assert set(closed.bonds) - set(radical.bonds) == {("H38", "O3")}


def test_matches_hydrogenated_dopa_radical(dopa_topology):
    donor_h, donor, acceptor = _hat_atoms(dopa_topology)

    dopa_topology.break_bond((donor.nr, donor_h.nr))
    donor_h.residue = acceptor.residue
    donor_h.resnr = acceptor.resnr

    match = dopa_topology.ff.product_templates.match_hydrogenated_residue(
        topology=dopa_topology,
        acceptor=acceptor,
        transferred_h=donor_h,
    )

    assert match is not None
    assert match.template.moleculetype == "DOPA"
    assert match.transferred_h_id == donor_h.nr
    assert match.acceptor_id == acceptor.nr
    assert match.transferred_h_template_name == "H38"
    assert match.atom_mapping[donor_h.nr] == match.template.atoms["H38"]
    assert match.atom_mapping[acceptor.nr] == match.template.atoms["O3"]


def test_rejects_template_when_existing_atom_type_differs(dopa_topology):
    donor_h, donor, acceptor = _hat_atoms(dopa_topology)
    dopa_topology.break_bond((donor.nr, donor_h.nr))
    donor_h.residue = acceptor.residue
    donor_h.resnr = acceptor.resnr
    acceptor.type = "CG321"

    match = dopa_topology.ff.product_templates.match_hydrogenated_residue(
        topology=dopa_topology,
        acceptor=acceptor,
        transferred_h=donor_h,
    )

    assert match is None


def test_rejects_template_with_unknown_atomtype(dopa_topology):
    known_atomtypes = dict(dopa_topology.ff.atomtypes)
    known_atomtypes.pop("HGP1")

    with pytest.raises(ValueError, match="atom type 'HGP1'.*missing"):
        dopa_topology.ff.product_templates.validate_atomtypes(known_atomtypes)


def test_rejects_unknown_product_template_source():
    filedir = Path(__file__).parent / "test_files" / "test_itp_templates"
    raw_top = read_top(filedir / "topol.top", ffdir=filedir)

    with pytest.raises(ValueError, match="Unknown product template sources.*banana"):
        Topology(raw_top, product_template_sources=("itp", "banana"))
