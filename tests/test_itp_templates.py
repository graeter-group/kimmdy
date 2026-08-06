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


@pytest.fixture
def ohr_topology() -> Topology:
    filedir = Path(__file__).parent / "test_files" / "test_itp_templates"
    raw_top = read_top(filedir / "ohr_topol.top", ffdir=filedir)
    return Topology(
        raw_top,
        product_template_sources=("itp",),
        product_template_itp_glob="*.itp",
        product_template_transitions={"OHR": "OHRH"},
    )


def _ohr_transition_atoms(topology: Topology):
    acceptor = next(atom for atom in topology.atoms.values() if atom.type == "ORA")
    donor = next(atom for atom in topology.atoms.values() if atom.type == "OT")
    donor_h = next(
        topology.atoms[atom_nr]
        for atom_nr in donor.bound_to_nrs
        if topology.atoms[atom_nr].atom == "H2"
    )
    return donor_h, donor, acceptor


def _prepare_ohr_transition(topology: Topology):
    donor_h, donor, acceptor = _ohr_transition_atoms(topology)
    topology.break_bond((donor.nr, donor_h.nr))
    donor_h.residue = acceptor.residue
    donor_h.resnr = acceptor.resnr
    return donor_h, donor, acceptor


def test_configured_transition_allows_exact_ohr_type_conversion(ohr_topology):
    donor_h, _donor, acceptor = _prepare_ohr_transition(ohr_topology)

    match = ohr_topology.ff.product_templates.match_hydrogenated_residue(
        topology=ohr_topology,
        acceptor=acceptor,
        transferred_h=donor_h,
    )

    assert match is not None
    assert match.template.moleculetype == "OHRH"
    assert match.transferred_h_template_name == "H2"
    mapped_by_name = {
        template_atom.name: template_atom.atom_type
        for template_atom in match.atom_mapping.values()
    }
    assert mapped_by_name == {"O1": "OT", "H1": "HT", "H2": "HT"}


def test_ohr_type_conversion_remains_disabled_without_transition():
    filedir = Path(__file__).parent / "test_files" / "test_itp_templates"
    raw_top = read_top(filedir / "ohr_topol.top", ffdir=filedir)
    topology = Topology(
        raw_top,
        product_template_sources=("itp",),
        product_template_itp_glob="*.itp",
    )
    donor_h, _donor, acceptor = _prepare_ohr_transition(topology)

    match = topology.ff.product_templates.match_hydrogenated_residue(
        topology=topology,
        acceptor=acceptor,
        transferred_h=donor_h,
    )

    assert match is None


def test_configured_transition_requires_exact_precursor(ohr_topology):
    donor_h, _donor, acceptor = _prepare_ohr_transition(ohr_topology)
    acceptor.type = "OT"

    match = ohr_topology.ff.product_templates.match_hydrogenated_residue(
        topology=ohr_topology,
        acceptor=acceptor,
        transferred_h=donor_h,
    )

    assert match is None


def test_rejects_unknown_product_template_transition():
    filedir = Path(__file__).parent / "test_files" / "test_itp_templates"
    raw_top = read_top(filedir / "ohr_topol.top", ffdir=filedir)

    with pytest.raises(ValueError, match="Unknown moleculetype.*MISSING"):
        Topology(
            raw_top,
            product_template_sources=("itp",),
            product_template_itp_glob="*.itp",
            product_template_transitions={"OHR": "MISSING"},
        )
