from copy import deepcopy
from pathlib import Path

import pytest

from kimmdy.parsing import read_top
from kimmdy.topology.topology import Topology
from kimmdy.topology.utils import get_residue_by_bonding


@pytest.fixture
def plpc_topology() -> Topology:
    filedir = Path(__file__).parent / "test_files" / "test_itp_templates"
    raw_top = read_top(filedir / "topol.top", ffdir=filedir)
    return Topology(
        raw_top,
        radicals="222",  # Radical C33: 132 closed-shell atoms + local atom 90.
        product_template_sources=("itp",),
        product_template_itp_glob="*.itp",
    )


def _atoms_in_residue(topology: Topology, resnr: str):
    return [atom for atom in topology.atoms.values() if atom.resnr == resnr]


def _hat_atoms(topology: Topology):
    """Return donor H3Y, its C33 donor, and radical C33 from the fixture."""
    donor_h = next(
        atom
        for atom in topology.atoms.values()
        if atom.atom == "H3Y"
        and any(
            candidate.atom == "H3X"
            for candidate in _atoms_in_residue(topology, atom.resnr)
        )
    )
    donor = topology.atoms[donor_h.bound_to_nrs[0]]
    acceptor = next(
        atom
        for atom in topology.atoms.values()
        if atom.atom == "C33"
        and not any(
            candidate.atom == "H3X"
            for candidate in _atoms_in_residue(topology, atom.resnr)
        )
    )
    return donor_h, donor, acceptor


def _transfer(topology: Topology, donor_h, donor, acceptor) -> Topology:
    """Move donor_h from donor to acceptor; returns the topology before (A)."""
    top_a = deepcopy(topology)
    topology.break_bond((donor.nr, donor_h.nr))
    topology.bind_bond((acceptor.nr, donor_h.nr))
    return top_a


def _match_acceptor(topology: Topology, top_a: Topology, acceptor, donor_h):
    """The ITP template of the acceptor's residue, which donor_h joined."""
    residue = get_residue_by_bonding(acceptor, topology.atoms)
    return topology.ff.templates.match_itp(
        residue, {donor_h.nr}, topology.atoms, top_a.atoms
    )


def test_loads_closed_and_radical_plpc_templates(plpc_topology):
    templates = plpc_topology.ff.templates.itp["PLPC"]
    templates_by_moleculetype = {template.name: template for template in templates}

    assert set(templates_by_moleculetype) == {"PLPC", "PLPC_rad"}

    closed = templates_by_moleculetype["PLPC"]
    radical = templates_by_moleculetype["PLPC_rad"]

    assert set(closed.atoms) - set(radical.atoms) == {"H3X"}
    assert set(radical.atoms) < set(closed.atoms)
    assert set(closed.bonds) - set(radical.bonds) == {("C33", "H3X")}


def test_matches_hydrogenated_plpc_radical(plpc_topology):
    donor_h, donor, acceptor = _hat_atoms(plpc_topology)
    top_a = _transfer(plpc_topology, donor_h, donor, acceptor)

    match = _match_acceptor(plpc_topology, top_a, acceptor, donor_h)

    assert match is not None
    assert match.template.name == "PLPC"
    assert match.names[donor_h.nr] == "H3X"
    assert match.names[acceptor.nr] == "C33"


def test_rejects_template_when_existing_atom_type_differs(plpc_topology):
    donor_h, donor, acceptor = _hat_atoms(plpc_topology)
    top_a = _transfer(plpc_topology, donor_h, donor, acceptor)
    acceptor.type = "CTL3"

    match = _match_acceptor(plpc_topology, top_a, acceptor, donor_h)

    assert match is None


def test_rejects_template_with_unknown_atomtype(plpc_topology):
    known_atomtypes = dict(plpc_topology.ff.atomtypes)
    known_atomtypes.pop("HAL2")

    with pytest.raises(ValueError, match="atom type 'HAL2'.*missing"):
        plpc_topology.ff.templates.validate_atomtypes(known_atomtypes)


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
    """Move H2 from the water-like OHRH to the OHR radical; returns A and the atoms."""
    donor_h, donor, acceptor = _ohr_transition_atoms(topology)
    top_a = _transfer(topology, donor_h, donor, acceptor)
    return top_a, donor_h, donor, acceptor


def test_configured_transition_allows_exact_ohr_type_conversion(ohr_topology):
    top_a, donor_h, _donor, acceptor = _prepare_ohr_transition(ohr_topology)

    match = _match_acceptor(ohr_topology, top_a, acceptor, donor_h)

    assert match is not None
    assert match.template.name == "OHRH"
    assert match.names[donor_h.nr] == "H2"
    mapped_by_name = {
        name: match.template.atoms[name].type for name in match.names.values()
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
    top_a, donor_h, _donor, acceptor = _prepare_ohr_transition(topology)

    match = _match_acceptor(topology, top_a, acceptor, donor_h)

    assert match is None


def test_configured_transition_requires_exact_precursor(ohr_topology):
    top_a, donor_h, _donor, acceptor = _prepare_ohr_transition(ohr_topology)
    # the precursor is the residue before the transfer (A)
    top_a.atoms[acceptor.nr].type = "OT"

    match = _match_acceptor(ohr_topology, top_a, acceptor, donor_h)

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
