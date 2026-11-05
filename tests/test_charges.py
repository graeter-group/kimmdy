"""Tests for the partial charges of reaction recipes (kimmdy.topology.charges)."""

import logging
from copy import deepcopy
from pathlib import Path
from types import SimpleNamespace
from typing import Optional

import numpy as np
import pytest

from kimmdy import schema
from kimmdy.parsing import read_top
from kimmdy.recipe import Bind, Break
from kimmdy.topology import charges
from kimmdy.topology.atomic import Atom
from kimmdy.topology.topology import Topology

FILEDIR = Path(__file__).parent / "test_files" / "test_topology"


@pytest.fixture(scope="module")
def imree_raw():
    return read_top(FILEDIR / "IMREE.top")


@pytest.fixture(scope="module")
def amber_ff(imree_raw):
    """Force field (with the amber99sb-star-ildnp residue templates) of IMREE.top."""
    return Topology(deepcopy(imree_raw)).ff


@pytest.fixture
def imree_hat(imree_raw):
    """IMREE with a Cα radical on GLU 5, as topologies A and B of a HAT from MET 4 Hα.

    del_atom moves the removed Hα's charge to its Cα (the radical convention).
    Returns (top_a, top_b, recipe_steps), with the steps applied to top_b only.
    """
    top = Topology(deepcopy(imree_raw))
    top.del_atom("70")  # GLU 5 HA, atom numbers after it shift down by one
    top_a = deepcopy(top)
    steps = [Break(atom_id_1="52", atom_id_2="53"), Bind(atom_id_1="69", atom_id_2="53")]
    top.break_bond(("52", "53"))
    top.bind_bond(("69", "53"))
    return top_a, top, steps


def _total_charge(top) -> float:
    return sum(float(atom.charge) for atom in top.atoms.values())


def _residue_charge(top, resnr: str) -> float:
    return sum(float(atom.charge) for atom in top.atoms.values() if atom.resnr == resnr)


# hydroxyl H charge of Hyp HD1, carried by a hydrogen that moves to a Cα
HYDROXYL_H_CHARGE = "0.3851"


class _RecipeTopology:
    """The parts of a Topology that assign_recipe_charges uses."""

    def __init__(self, atoms: dict[str, Atom], ff):
        self.atoms = atoms
        self.ff = ff
        self.radicals: dict[str, Atom] = {}
        self.pending_product_template_matches: dict = {}


def _template_residue(
    ff, template: str, resnr: str, first_nr: int, residue: Optional[str] = None
) -> dict[str, Atom]:
    """Atoms of an rtp entry, with its charges and its bonds inside the residue.

    residue is the residue name in the topology, the entry's name by default.
    """
    residuetype = ff.residuetypes[template]
    nrs = {name: str(first_nr + i) for i, name in enumerate(residuetype.atoms)}
    atoms = {
        nrs[name]: Atom(
            nr=nrs[name],
            type=spec.type,
            resnr=resnr,
            residue=residue or template,
            atom=name,
            cgnr=nrs[name],
            charge=spec.charge,
        )
        for name, spec in residuetype.atoms.items()
    }
    for name_1, name_2 in residuetype.bonds:
        if name_1 in nrs and name_2 in nrs:
            atoms[nrs[name_1]].bound_to_nrs.append(nrs[name_2])
            atoms[nrs[name_2]].bound_to_nrs.append(nrs[name_1])
    return atoms


def _bond(atoms: dict[str, Atom], nr_1: str, nr_2: str, bound: bool) -> None:
    for nr, other in ((nr_1, nr_2), (nr_2, nr_1)):
        if bound:
            atoms[nr].bound_to_nrs.append(other)
        else:
            atoms[nr].bound_to_nrs.remove(other)


def _hat_onto_gly(
    ff, template: str, residue: str = "GLY", h_name: str = "HX", h_type: str = "HO"
):
    """A hydroxyl hydrogen moves from a donor residue to the Cα radical of a Gly.

    Residue 1 has the atoms and charges of the rtp entry template, without
    HA2; its Cα carries HA2's charge (radical convention). Residue 2 (DON,
    no template) is an O-H with the hydroxyl charges. In B, the hydrogen is
    bound to the Cα and bind_bond named it h_name with type h_type.
    Returns (top_a, top_b, recipe_steps, Cα number).
    """
    atoms = _template_residue(ff, template, "1", 1, residue)
    removed = next(nr for nr, atom in atoms.items() if atom.atom == "HA2")
    ca = next(nr for nr, atom in atoms.items() if atom.atom == "CA")
    atoms[ca].bound_to_nrs.remove(removed)
    atoms[ca].charge = f"{float(atoms[ca].charge) + float(atoms.pop(removed).charge):7.4f}"
    atoms["31"] = Atom(
        nr="31", type="OH", resnr="2", residue="DON", atom="O", cgnr="31",
        charge=f"-{HYDROXYL_H_CHARGE}",
    )
    atoms["32"] = Atom(
        nr="32", type="HO", resnr="2", residue="DON", atom="HO", cgnr="32",
        charge=HYDROXYL_H_CHARGE,
    )
    _bond(atoms, "31", "32", True)
    top_a = _RecipeTopology(atoms, ff)

    top_b = _RecipeTopology(deepcopy(atoms), ff)
    _bond(top_b.atoms, "31", "32", False)
    _bond(top_b.atoms, ca, "32", True)
    hydrogen = top_b.atoms["32"]
    hydrogen.resnr, hydrogen.residue, hydrogen.atom, hydrogen.type = (
        "1", residue, h_name, h_type,
    )
    top_b.radicals = {"31": top_b.atoms["31"]}
    steps = [Break(atom_id_1="31", atom_id_2="32"), Bind(atom_id_1=ca, atom_id_2="32")]
    return top_a, top_b, steps, ca


def _assert_template_charges(top, resnr: str, template: str) -> None:
    residuetype = top.ff.residuetypes[template]
    for atom in top.atoms.values():
        if atom.resnr == resnr:
            assert np.isclose(float(atom.charge), float(residuetype.atoms[atom.atom].charge))


@pytest.mark.parametrize(
    "template, expected_type", [("GLY", "H1"), ("NGLY", "HP"), ("CGLY", "H1")]
)
def test_rtp_selects_terminal_variant(amber_ff, caplog, template, expected_type):
    """A Gly named GLY gets the charges of the variant its atoms come from; the hydrogen is renamed."""
    top_a, top_b, steps, _ = _hat_onto_gly(amber_ff, template)

    with caplog.at_level(logging.INFO, logger="kimmdy.topology.charges"):
        charges.assign_recipe_charges(top_a, top_b, steps)

    assert f"Charges of residue 1 from rtp:{template} (renamed hydrogen)" in caplog.text
    assert "Charges of residue 2 from local" in caplog.text
    hydrogen = top_b.atoms["32"]
    assert (hydrogen.atom, hydrogen.type) == ("HA2", expected_type)
    _assert_template_charges(top_b, "1", template)
    assert np.isclose(float(top_b.atoms["31"].charge), 0.0)
    assert np.isclose(_total_charge(top_b), _total_charge(top_a))


def test_rtp_exact_match_keeps_the_hydrogen_name(amber_ff, caplog):
    top_a, top_b, steps, _ = _hat_onto_gly(amber_ff, "GLY", h_name="HA2", h_type="H1")

    with caplog.at_level(logging.INFO, logger="kimmdy.topology.charges"):
        charges.assign_recipe_charges(top_a, top_b, steps)

    assert "Charges of residue 1 from rtp:GLY (exact)" in caplog.text
    _assert_template_charges(top_b, "1", "GLY")


def test_local_rules_when_templates_would_change_the_residue_charge(amber_ff, caplog):
    """A residue whose charges differ from its template gets the local rules: the hydrogen keeps its charge."""
    top_a, top_b, steps, ca = _hat_onto_gly(amber_ff, "CGLY")
    c_atom = next(nr for nr, atom in top_a.atoms.items() if atom.atom == "C")
    for top in (top_a, top_b):
        top.atoms[c_atom].charge = "0.7731"  # CGLY C is 0.7231

    with caplog.at_level(logging.DEBUG, logger="kimmdy.topology.charges"):
        charges.assign_recipe_charges(top_a, top_b, steps)

    assert "rtp:CGLY (renamed hydrogen) changes the charge of residue 1 by -0.0500 e" in caplog.text
    assert "Charges of residue 1 from local" in caplog.text
    assert np.isclose(float(top_b.atoms["32"].charge), float(HYDROXYL_H_CHARGE))
    assert np.isclose(
        float(top_b.atoms[ca].charge),
        float(top_a.atoms[ca].charge) - float(HYDROXYL_H_CHARGE),
    )
    assert top_b.atoms[c_atom].charge == "0.7731"
    for resnr in ("1", "2"):
        assert np.isclose(_residue_charge(top_b, resnr), _residue_charge(top_a, resnr))


def test_local_rules_without_template(amber_ff, caplog):
    top_a, top_b, steps, ca = _hat_onto_gly(amber_ff, "GLY", residue="XYZ")

    with caplog.at_level(logging.INFO, logger="kimmdy.topology.charges"):
        charges.assign_recipe_charges(top_a, top_b, steps)

    assert "Charges of residue 1 from local" in caplog.text
    assert top_b.atoms["32"].atom == "HX"
    assert np.isclose(float(top_b.atoms["32"].charge), float(HYDROXYL_H_CHARGE))
    assert np.isclose(_total_charge(top_b), _total_charge(top_a))


def test_hydrogen_without_new_partner_leaves_its_charge_on_the_old_one(amber_ff, caplog):
    top_a, _, _, _ = _hat_onto_gly(amber_ff, "GLY")
    top_b = _RecipeTopology(deepcopy(top_a.atoms), amber_ff)
    _bond(top_b.atoms, "31", "32", False)
    steps = [Break(atom_id_1="31", atom_id_2="32")]

    with caplog.at_level(logging.INFO, logger="kimmdy.topology.charges"):
        charges.assign_recipe_charges(top_a, top_b, steps)

    assert "Charges of residue 2 from local" in caplog.text
    assert np.isclose(float(top_b.atoms["31"].charge), 0.0)
    assert np.isclose(float(top_b.atoms["32"].charge), 0.0)


def _itp_match(top_b, template: str, moleculetype: str = "FAKE"):
    """A stand-in for the exact ITP match of residue 1, with the charges of an rtp entry."""
    residuetype = top_b.ff.residuetypes[template]
    names = {nr: atom.atom for nr, atom in top_b.atoms.items() if atom.resnr == "1"}
    names["32"] = "HA2"
    mapping = {
        nr: SimpleNamespace(
            name=name,
            atom_type=residuetype.atoms[name].type,
            charge=residuetype.atoms[name].charge,
            mass=None,
        )
        for nr, name in names.items()
    }
    return SimpleNamespace(
        atom_mapping=mapping, template=SimpleNamespace(moleculetype=moleculetype)
    )


def _add_itp_match(top_b, ca: str, match) -> None:
    top_b.pending_product_template_matches[tuple(sorted((ca, "32"), key=int))] = match


def test_itp_match_precedes_rtp(amber_ff, caplog):
    top_a, top_b, steps, ca = _hat_onto_gly(amber_ff, "CGLY")
    _add_itp_match(top_b, ca, _itp_match(top_b, "CGLY"))

    with caplog.at_level(logging.INFO, logger="kimmdy.topology.charges"):
        charges.assign_recipe_charges(top_a, top_b, steps)

    assert "Charges of residue 1 from itp:FAKE" in caplog.text
    assert top_b.atoms["32"].atom == "HA2"
    _assert_template_charges(top_b, "1", "CGLY")
    assert top_b.pending_product_template_matches == {}


def test_itp_match_changing_the_residue_charge_hands_over_to_rtp(amber_ff, caplog):
    top_a, top_b, steps, ca = _hat_onto_gly(amber_ff, "CGLY")
    match = _itp_match(top_b, "CGLY")
    match.atom_mapping[ca].charge = "-0.1493"  # CGLY CA is -0.2493
    _add_itp_match(top_b, ca, match)

    with caplog.at_level(logging.DEBUG, logger="kimmdy.topology.charges"):
        charges.assign_recipe_charges(top_a, top_b, steps)

    assert "itp:FAKE changes the charge of residue 1 by 0.1000 e" in caplog.text
    assert "Charges of residue 1 from rtp:CGLY (renamed hydrogen)" in caplog.text
    _assert_template_charges(top_b, "1", "CGLY")


def test_failing_scheme_is_skipped(amber_ff, caplog, monkeypatch):
    top_a, top_b, steps, ca = _hat_onto_gly(amber_ff, "CGLY")
    _add_itp_match(top_b, ca, _itp_match(top_b, "CGLY"))

    def broken_itp_scheme(residue, match):
        raise RuntimeError("broken template")

    monkeypatch.setattr(charges, "_itp_scheme", broken_itp_scheme)
    with caplog.at_level(logging.INFO, logger="kimmdy.topology.charges"):
        charges.assign_recipe_charges(top_a, top_b, steps)

    assert "Charge scheme itp failed for residue 1; trying the next one." in caplog.text
    assert "Charges of residue 1 from rtp:CGLY (renamed hydrogen)" in caplog.text
    _assert_template_charges(top_b, "1", "CGLY")


def test_assign_recipe_charges_hat_between_residues(imree_hat, caplog):
    """A HAT between residues conserves both residue charges and uses the template H charge."""
    top_a, top_b, steps = imree_hat
    donor_charge_a = float(top_a.atoms["52"].charge)
    h_charge_a = float(top_a.atoms["53"].charge)

    with caplog.at_level(logging.INFO, logger="kimmdy.topology.charges"):
        charges.assign_recipe_charges(top_a, top_b, steps)

    assert "Charges of residue 5 from rtp:GLU (exact)" in caplog.text
    glu = top_b.ff.residuetypes["GLU"]
    assert np.isclose(float(top_b.atoms["53"].charge), float(glu.atoms["HA"].charge))
    assert np.isclose(float(top_b.atoms["69"].charge), float(glu.atoms["CA"].charge))
    assert np.isclose(float(top_b.atoms["52"].charge), donor_charge_a + h_charge_a)
    for resnr in ("4", "5"):
        assert np.isclose(_residue_charge(top_b, resnr), _residue_charge(top_a, resnr))
    assert np.isclose(_total_charge(top_b), _total_charge(top_a))


@pytest.mark.parametrize("flexible_water", [False, True])
def test_residue_templates_set_hydrolysis_termini(amber_ff, caplog, flexible_water):
    """Hydrolysis of Gly-Gly: the termini get CGLY and NGLY charges without the plugin setting them."""
    # A: GLY 1 (atoms 1-7) - GLY 2 (atoms 11-17) and a water (21-23)
    atoms = _template_residue(amber_ff, "GLY", "1", 1)
    atoms |= _template_residue(amber_ff, "GLY", "2", 11)
    for nr, name, atom_type, charge in (
        ("21", "OW", "OW", "-0.834"),
        ("22", "HW1", "HW", "0.417"),
        ("23", "HW2", "HW", "0.417"),
    ):
        atoms[nr] = Atom(
            nr=nr, type=atom_type, resnr="3", residue="SOL", atom=name, cgnr=nr, charge=charge
        )
    _bond(atoms, "6", "11", True)  # peptide bond C-N
    if flexible_water:
        _bond(atoms, "21", "22", True)
        _bond(atoms, "21", "23", True)
    top_a = _RecipeTopology(atoms, amber_ff)

    # B: what the hydrolysis recipe does to the bonds, names and types, without charges
    top_b = _RecipeTopology(deepcopy(atoms), amber_ff)
    steps = [Break(atom_id_1="6", atom_id_2="11")]
    if flexible_water:
        steps += [Break(atom_id_1="21", atom_id_2="22"), Break(atom_id_1="21", atom_id_2="23")]
    steps += [
        Bind(atom_id_1="21", atom_id_2="6"),
        Bind(atom_id_1="22", atom_id_2="11"),
        Bind(atom_id_1="23", atom_id_2="11"),
    ]
    for step in steps:
        _bond(top_b.atoms, step.atom_id_1, step.atom_id_2, isinstance(step, Bind))
    renames = {
        "7": ("1", "CGLY", "OC1", "O2"),
        "21": ("1", "CGLY", "OC2", "O2"),
        "12": ("2", "NGLY", "H1", "H"),
        "22": ("2", "NGLY", "H2", "H"),
        "23": ("2", "NGLY", "H3", "H"),
        "11": ("2", "NGLY", "N", "N3"),
    }
    for atom in top_b.atoms.values():
        if atom.resnr in ("1", "2"):
            atom.residue = {"1": "CGLY", "2": "NGLY"}[atom.resnr]
    for nr, (resnr, residue, name, atom_type) in renames.items():
        atom = top_b.atoms[nr]
        atom.resnr, atom.residue, atom.atom, atom.type = resnr, residue, name, atom_type

    with caplog.at_level(logging.INFO, logger="kimmdy.topology.charges"):
        charges.assign_recipe_charges(top_a, top_b, steps)

    assert "Charges of residue 1 from rtp:CGLY (exact)" in caplog.text
    assert "Charges of residue 2 from rtp:NGLY (exact)" in caplog.text
    for atom in top_b.atoms.values():
        template = amber_ff.residuetypes[atom.residue]
        assert np.isclose(float(atom.charge), float(template.atoms[atom.atom].charge))
    assert np.isclose(_total_charge(top_b), _total_charge(top_a))


def test_residue_templates_set_hat_donor_radical_and_acceptor(imree_hat, caplog):
    """After a HAT the donor residue gets its template as a radical, the acceptor as closed shell."""
    top_a, top_b, steps = imree_hat

    with caplog.at_level(logging.INFO, logger="kimmdy.topology.charges"):
        charges.assign_recipe_charges(top_a, top_b, steps)

    assert "Charges of residue 4 from rtp:MET (radical)" in caplog.text
    assert "Charges of residue 5 from rtp:GLU (exact)" in caplog.text
    met = top_b.ff.residuetypes["MET"]
    glu = top_b.ff.residuetypes["GLU"]
    assert np.isclose(
        float(top_b.atoms["52"].charge),
        float(met.atoms["CA"].charge) + float(met.atoms["HA"].charge),
    )
    assert top_b.atoms["53"].atom == "HA"
    assert np.isclose(float(top_b.atoms["53"].charge), float(glu.atoms["HA"].charge))
    assert np.isclose(_total_charge(top_b), _total_charge(top_a))


def test_residue_templates_are_rejected_when_charges_would_change(imree_hat, caplog):
    """A residue whose charges differ from its template gets the local rules."""
    top_a, top_b, steps = imree_hat
    # GLU 5 CB (atom 70 after the deleted HA) no longer carries its template charge
    for top in (top_a, top_b):
        top.atoms["70"].charge = "0.1060"

    with caplog.at_level(logging.INFO, logger="kimmdy.topology.charges"):
        charges.assign_recipe_charges(top_a, top_b, steps)

    assert "Charges of residue 5 from local" in caplog.text
    assert top_b.atoms["70"].charge == "0.1060"
    for resnr in ("4", "5"):
        assert np.isclose(_residue_charge(top_b, resnr), _residue_charge(top_a, resnr))


def test_check_charge_conservation_actions(imree_raw, caplog):
    top_a = Topology(deepcopy(imree_raw))
    top_b = deepcopy(top_a)
    top_b.atoms["1"].charge = f"{float(top_b.atoms['1'].charge) + 0.1:7.4f}"

    with caplog.at_level(logging.WARNING):
        delta = charges.check_charge_conservation(top_a, top_b, action="warn")
    assert np.isclose(delta, 0.1)
    assert "changed the net charge" in caplog.text

    with pytest.raises(ValueError, match="changed the net charge"):
        charges.check_charge_conservation(top_a, top_b, action="error")

    assert charges.check_charge_conservation(top_a, top_b, action="off") == 0.0
    assert np.isclose(
        charges.check_charge_conservation(top_a, deepcopy(top_a), action="error"), 0.0
    )


def test_prepare_product_parameters_parameterizes_once(imree_raw):
    """The product is parameterized once, on the reaction atoms, then only radicals stay in focus."""

    class RecordingParameterizer:
        def __init__(self):
            self.calls = []

        def parameterize_topology(self, current_topology, focus_nrs=None):
            self.calls.append(set(focus_nrs or set()))
            return current_topology

    top = Topology(deepcopy(imree_raw))
    recorder = RecordingParameterizer()
    top.parametrizer = recorder
    top.break_bond(("52", "53"))

    report = top.prepare_product_parameters()
    top.prepare_product_parameters()

    assert recorder.calls == [{"52", "53"}]
    assert report.event_focus_ids == {"52", "53"}
    assert report.active_radical_ids == set(top.radicals)
    assert top.parameterization_focus_ids == set(top.radicals)


def test_update_charges_bonded_sets_charges_before_parameterizing(imree_hat):
    """The parameterizer runs once and sees the final charges of B."""
    top_a, top_b, steps = imree_hat
    seen_h_charges = []

    class RecordingParameterizer:
        def parameterize_topology(self, current_topology, focus_nrs=None):
            seen_h_charges.append(current_topology.atoms["53"].charge)
            return current_topology

    top_b.parametrizer = RecordingParameterizer()
    report = top_b.update_charges_bonded(top_a, steps, conservation_check="error")

    glu = top_b.ff.residuetypes["GLU"]
    assert len(seen_h_charges) == 1
    assert np.isclose(float(seen_h_charges[0]), float(glu.atoms["HA"].charge))
    assert {"52", "53", "69"} <= report.event_focus_ids
    assert top_b.parameterization_focus_ids == set(top_b.radicals)


def test_update_charges_bonded_without_topology_change_only_prepares_parameters(
    imree_raw, monkeypatch
):
    """Without A (no Break, Bind or CustomTopMod), no charges are assigned."""

    def fail(*args, **kwargs):
        raise AssertionError("charges assigned without a topology change")

    monkeypatch.setattr("kimmdy.topology.topology.assign_recipe_charges", fail)
    monkeypatch.setattr("kimmdy.topology.topology.check_charge_conservation", fail)
    top = Topology(deepcopy(imree_raw))
    top.parameterization_focus_ids = {"1"}  # e.g. from a Place

    report = top.update_charges_bonded(None, [])

    assert report.event_focus_ids == {"1"}
    assert top.parameterization_focus_ids == set(top.radicals)


def test_residue_variants_use_amber_terminal_names():
    assert charges.residue_variants("GLY") == ["GLY", "NGLY", "CGLY"]


def test_apply_recipe_schema_options():
    changer = schema.load_kimmdy_schema()["properties"]["changer"]["properties"]
    charge_options = changer["apply_recipe"]["properties"]["charges"]["properties"]
    assert "mode" not in charge_options
    assert charge_options["conservation_check"]["enum"] == list(
        charges.CONSERVATION_CHECK_ACTIONS
    )
    conserve = changer["topology"]["properties"]["conserve_transferred_h_charge"]
    assert conserve["deprecated"] is True
    assert "default" not in conserve
