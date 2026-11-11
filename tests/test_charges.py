"""Tests for the partial charges of reaction recipes (kimmdy.topology.charges)."""

import logging
from copy import deepcopy
from dataclasses import replace
from pathlib import Path
from typing import Optional

import numpy as np
import pytest

from kimmdy import schema
from kimmdy.parsing import read_top
from kimmdy.recipe import Bind, Break
from kimmdy.topology import charges
from kimmdy.topology.atomic import Atom
from kimmdy.topology.templates import ResidueTemplate, assign_recipe_names
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
    assign_recipe_names(top_a, top_b, steps)

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


def _itp_template(
    ff, template: str, moleculetype: str = "FAKE", **atom_charges: str
) -> ResidueTemplate:
    """A stand-in for an ITP template of the GLY residue, with the atoms of an rtp entry.

    atom_charges replace the entry's charges of some atoms, by atom name.
    """
    rtp = ResidueTemplate.from_rtp(ff.residuetypes[template])
    atoms = {
        name: replace(atom, charge=atom_charges.get(name, atom.charge))
        for name, atom in rtp.atoms.items()
    }
    return replace(
        rtp,
        source="itp",
        name=moleculetype,
        residue="GLY",
        atoms=atoms,
        boundary_bonds=frozenset(),
    )


def _add_itp_template(monkeypatch, ff, template: ResidueTemplate) -> None:
    monkeypatch.setitem(ff.templates.itp, template.residue, [template])


def test_itp_match_precedes_rtp(amber_ff, caplog, monkeypatch):
    top_a, top_b, steps, ca = _hat_onto_gly(amber_ff, "CGLY")
    _add_itp_template(monkeypatch, amber_ff, _itp_template(amber_ff, "CGLY"))

    with caplog.at_level(logging.INFO, logger="kimmdy.topology.charges"):
        charges.assign_recipe_charges(top_a, top_b, steps)
    assign_recipe_names(top_a, top_b, steps)

    assert "Charges of residue 1 from itp:FAKE" in caplog.text
    assert top_b.atoms["32"].atom == "HA2"
    _assert_template_charges(top_b, "1", "CGLY")


def test_itp_match_changing_the_residue_charge_hands_over_to_rtp(amber_ff, caplog, monkeypatch):
    top_a, top_b, steps, ca = _hat_onto_gly(amber_ff, "CGLY")
    # CGLY CA is -0.2493
    _add_itp_template(monkeypatch, amber_ff, _itp_template(amber_ff, "CGLY", CA="-0.1493"))

    with caplog.at_level(logging.DEBUG, logger="kimmdy.topology.charges"):
        charges.assign_recipe_charges(top_a, top_b, steps)
    assign_recipe_names(top_a, top_b, steps)

    assert "itp:FAKE changes the charge of residue 1 by 0.1000 e" in caplog.text
    assert "Charges of residue 1 from rtp:CGLY (renamed hydrogen)" in caplog.text
    _assert_template_charges(top_b, "1", "CGLY")


def test_failing_scheme_is_skipped(amber_ff, caplog, monkeypatch):
    top_a, top_b, steps, ca = _hat_onto_gly(amber_ff, "CGLY")
    _add_itp_template(monkeypatch, amber_ff, _itp_template(amber_ff, "CGLY"))

    def broken_itp_scheme(*args):
        raise RuntimeError("broken template")

    monkeypatch.setattr(charges, "_itp_scheme", broken_itp_scheme)
    with caplog.at_level(logging.INFO, logger="kimmdy.topology.charges"):
        charges.assign_recipe_charges(top_a, top_b, steps)
    assign_recipe_names(top_a, top_b, steps)

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

    report = top.update_all(None, [])
    top.update_all(None, [])

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
    report = top_b.update_all(top_a, steps, conservation_check="error")

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
    monkeypatch.setattr("kimmdy.topology.topology.assign_recipe_names", fail)
    top = Topology(deepcopy(imree_raw))
    top.parameterization_focus_ids = {"1"}  # e.g. from a Place

    report = top.update_all(None, [])

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


# --- residues with the same residue number (merged topologies) ---------------------------------


def _nr_of(atoms: dict[str, Atom], first_nr: int, last_nr: int, name: str) -> str:
    return next(
        nr for nr, atom in atoms.items() if first_nr <= int(nr) <= last_nr and atom.atom == name
    )


def _tyr_oh_radical(ff, resnr: str, first_nr: int) -> dict[str, Atom]:
    """A Tyr with the rtp charges whose HH is abstracted: OH carries HH's charge (radical convention)."""
    atoms = _template_residue(ff, "TYR", resnr, first_nr)
    last_nr = first_nr + len(atoms) - 1
    oh = _nr_of(atoms, first_nr, last_nr, "OH")
    hh = _nr_of(atoms, first_nr, last_nr, "HH")
    atoms[oh].bound_to_nrs.remove(hh)
    atoms[oh].charge = f"{float(atoms[oh].charge) + float(atoms.pop(hh).charge):7.4f}"
    return atoms


def _hat_onto_tyr_oh(ff, donor: dict[str, Atom], donor_nr: str, h_nr: str, first_nr: int):
    """Hydrogen h_nr moves from donor_nr to the OH radical of a Tyr with the same residue number.

    The Tyr starts at atom first_nr. In B, bind_bond gave the hydrogen the
    Tyr's residue (here the same number) and named it HX.
    Returns (top_a, top_b, recipe_steps, Tyr OH number).
    """
    tyr = _tyr_oh_radical(ff, donor[donor_nr].resnr, first_nr)
    oh = next(nr for nr, atom in tyr.items() if atom.atom == "OH")
    top_a = _RecipeTopology({**donor, **tyr}, ff)
    top_b = _RecipeTopology(deepcopy(top_a.atoms), ff)
    _bond(top_b.atoms, donor_nr, h_nr, False)
    _bond(top_b.atoms, oh, h_nr, True)
    hydrogen = top_b.atoms[h_nr]
    hydrogen.residue, hydrogen.atom = "TYR", "HX"
    top_b.radicals = {donor_nr: top_b.atoms[donor_nr]}
    steps = [Break(atom_id_1=donor_nr, atom_id_2=h_nr), Bind(atom_id_1=oh, atom_id_2=h_nr)]
    return top_a, top_b, steps, oh


def test_hat_between_residues_with_the_same_number_uses_both_templates(amber_ff, caplog):
    """Two Tyr 5 on different chains stay two residues: donor and acceptor get their templates.

    Grouped by residue number, they were one residue with two fragments, and
    only the local rules ran: the hydrogen kept its HB2 charge on the OH.
    """
    donor = _template_residue(amber_ff, "TYR", "5", 1)
    cb = _nr_of(donor, 1, 21, "CB")
    hb2 = _nr_of(donor, 1, 21, "HB2")
    top_a, top_b, steps, oh = _hat_onto_tyr_oh(amber_ff, donor, cb, hb2, 101)
    template = amber_ff.residuetypes["TYR"].atoms

    with caplog.at_level(logging.INFO, logger="kimmdy.topology.charges"):
        charges.assign_recipe_charges(top_a, top_b, steps)
    assign_recipe_names(top_a, top_b, steps)

    assert "Charges of residue 5 from rtp:TYR (radical)" in caplog.text
    assert "Charges of residue 5 from rtp:TYR (renamed hydrogen)" in caplog.text
    assert "from local" not in caplog.text
    hydrogen = top_b.atoms[hb2]
    assert (hydrogen.atom, hydrogen.type) == ("HH", template["HH"].type)
    assert np.isclose(float(hydrogen.charge), float(template["HH"].charge))
    assert np.isclose(float(top_b.atoms[oh].charge), float(template["OH"].charge))
    assert np.isclose(
        float(top_b.atoms[cb].charge),
        float(template["CB"].charge) + float(template["HB2"].charge),
    )
    assert np.isclose(_total_charge(top_b), _total_charge(top_a))


def test_hat_from_non_integer_residue_with_the_same_number_keeps_its_charge(amber_ff, caplog):
    """A donor that is not integer on its own (like a crosslink half) keeps its charge.

    Grouped by residue number with the acceptor, the local rules made the
    donor fragment integer, and the scheme was then rejected.
    """
    donor = {
        "201": Atom(nr="201", type="OH", resnr="5", residue="DON", atom="O", cgnr="201",
                    charge="-0.6000"),
        "202": Atom(nr="202", type="HO", resnr="5", residue="DON", atom="HO", cgnr="202",
                    charge="0.4000"),
    }
    _bond(donor, "201", "202", True)
    top_a, top_b, steps, oh = _hat_onto_tyr_oh(amber_ff, donor, "201", "202", 1)

    with caplog.at_level(logging.INFO, logger="kimmdy.topology.charges"):
        charges.assign_recipe_charges(top_a, top_b, steps)

    assert "Charges of residue 5 from local" in caplog.text
    assert "Charges of residue 5 from rtp:TYR (renamed hydrogen)" in caplog.text
    assert "No charge scheme" not in caplog.text
    assert np.isclose(float(top_b.atoms["201"].charge), -0.2)
    assert np.isclose(_total_charge(top_b), _total_charge(top_a))


def test_backbone_homolysis_keeps_the_split_residue_together(amber_ff, caplog):
    """CA-C homolysis inside a Gly: one residue with two fragments, each made integer by the local rules."""
    atoms = _template_residue(amber_ff, "GLY", "1", 1)
    ca = _nr_of(atoms, 1, 99, "CA")
    c = _nr_of(atoms, 1, 99, "C")
    o = _nr_of(atoms, 1, 99, "O")
    top_a = _RecipeTopology(atoms, amber_ff)
    top_b = _RecipeTopology(deepcopy(atoms), amber_ff)
    _bond(top_b.atoms, ca, c, False)
    top_b.radicals = {ca: top_b.atoms[ca], c: top_b.atoms[c]}
    steps = [Break(atom_id_1=ca, atom_id_2=c)]

    with caplog.at_level(logging.INFO, logger="kimmdy.topology.charges"):
        charges.assign_recipe_charges(top_a, top_b, steps)

    assert caplog.text.count("Charges of residue 1 from local") == 1
    assert np.isclose(float(top_b.atoms[c].charge) + float(top_b.atoms[o].charge), 0.0)
    assert np.isclose(_residue_charge(top_b, "1"), _residue_charge(top_a, "1"))


def test_recombination_of_a_split_residue_is_one_residue(amber_ff, caplog):
    """CA and C of a split Gly bind again: two residues in A, one in B, which fits the template."""
    atoms = _template_residue(amber_ff, "GLY", "1", 1)
    ca = _nr_of(atoms, 1, 99, "CA")
    c = _nr_of(atoms, 1, 99, "C")
    _bond(atoms, ca, c, False)
    top_a = _RecipeTopology(atoms, amber_ff)
    top_b = _RecipeTopology(deepcopy(atoms), amber_ff)
    _bond(top_b.atoms, ca, c, True)
    steps = [Bind(atom_id_1=ca, atom_id_2=c)]

    with caplog.at_level(logging.INFO, logger="kimmdy.topology.charges"):
        charges.assign_recipe_charges(top_a, top_b, steps)

    assert caplog.text.count("Charges of residue 1 from") == 1
    assert "Charges of residue 1 from rtp:GLY (exact)" in caplog.text
    _assert_template_charges(top_b, "1", "GLY")


# --- protonation states (HIS = HID/HIE/HIP, ...) and radicals the rtp entry has to explain ----


def test_residue_variants_include_protonation_states():
    assert charges.residue_variants("HIS") == [
        "HIS", "NHIS", "CHIS",
        "HID", "NHID", "CHID",
        "HIE", "NHIE", "CHIE",
        "HIP", "NHIP", "CHIP",
    ]
    assert charges.residue_variants("LYN") == ["LYN", "NLYN", "CLYN", "LYS", "NLYS", "CLYS"]
    assert charges.residue_variants("CYS")[3:6] == ["CYX", "NCYX", "CCYX"]


def _template_charge(ff, template: str) -> float:
    return sum(float(atom.charge) for atom in ff.residuetypes[template].atoms.values())


def _hat_onto_cb(ff, template: str, residue: str):
    """A hydroxyl hydrogen of DON 2 moves to the CB radical (HB2 abstracted) of residue 1.

    Residue 1 is named residue and has the atoms and charges of the rtp entry
    template; its CB carries HB2's charge (radical convention).
    Returns (top_a, top_b, recipe_steps, hydrogen number).
    """
    atoms = _template_residue(ff, template, "1", 1, residue=residue)
    cb = _nr_of(atoms, 1, 99, "CB")
    hb2 = _nr_of(atoms, 1, 99, "HB2")
    atoms[cb].bound_to_nrs.remove(hb2)
    atoms[cb].charge = f"{float(atoms[cb].charge) + float(atoms.pop(hb2).charge):7.4f}"
    atoms["101"] = Atom(nr="101", type="OH", resnr="2", residue="DON", atom="O", cgnr="101",
                        charge=f"-{HYDROXYL_H_CHARGE}")
    atoms["102"] = Atom(nr="102", type="HO", resnr="2", residue="DON", atom="HO", cgnr="102",
                        charge=HYDROXYL_H_CHARGE)
    _bond(atoms, "101", "102", True)
    top_a = _RecipeTopology(atoms, ff)
    top_b = _RecipeTopology(deepcopy(atoms), ff)
    _bond(top_b.atoms, "101", "102", False)
    _bond(top_b.atoms, cb, "102", True)
    hydrogen = top_b.atoms["102"]
    hydrogen.resnr, hydrogen.residue, hydrogen.atom = "1", residue, "HX"
    top_b.radicals = {"101": top_b.atoms["101"]}
    steps = [Break(atom_id_1="101", atom_id_2="102"), Bind(atom_id_1=cb, atom_id_2="102")]
    return top_a, top_b, steps, "102"


def test_hat_onto_his_uses_the_hie_entry(amber_ff, caplog):
    """A HIS with the atoms of HIE gets the HIE charges; the hydrogen takes the free HB2 position."""
    top_a, top_b, steps, h = _hat_onto_cb(amber_ff, "HIE", "HIS")

    with caplog.at_level(logging.INFO, logger="kimmdy.topology.charges"):
        charges.assign_recipe_charges(top_a, top_b, steps)
    assign_recipe_names(top_a, top_b, steps)

    assert "Charges of residue 1 from rtp:HIE (renamed hydrogen)" in caplog.text
    hie = amber_ff.residuetypes["HIE"].atoms
    assert (top_b.atoms[h].atom, top_b.atoms[h].type) == ("HB2", hie["HB2"].type)
    _assert_template_charges(top_b, "1", "HIE")
    assert np.isclose(_total_charge(top_b), _total_charge(top_a))


@pytest.mark.parametrize("template", ["CYX", "CYM"])
def test_cys_states_with_the_same_atoms_keep_their_charge(amber_ff, caplog, template):
    """CYX (0) and CYM (-1) have the same atom names; each residue named CYS keeps its own state.

    Both entries fit at the same level, the one that keeps the residue
    charge is taken (a change of +-1 alone would also pass).
    """
    top_a, top_b, steps, _ = _hat_onto_cb(amber_ff, template, "CYS")

    with caplog.at_level(logging.INFO, logger="kimmdy.topology.charges"):
        charges.assign_recipe_charges(top_a, top_b, steps)

    assert f"Charges of residue 1 from rtp:{template} (renamed hydrogen)" in caplog.text
    assert np.isclose(_residue_charge(top_b, "1"), _template_charge(amber_ff, template))
    assert np.isclose(_total_charge(top_b), _total_charge(top_a))


@pytest.mark.parametrize(
    "template, residue, heavy, hydrogen, other_state",
    [
        ("HIP", "HIS", "NE2", "HE2", "HID"),
        ("LYS", "LYS", "NZ", "HZ3", "LYN"),
        ("ASH", "ASP", "OD2", "HD2", "ASP"),
        ("GLH", "GLU", "OE2", "HE2", "GLU"),
        ("CYS", "CYS", "SG", "HG", "CYX"),
    ],
)
def test_abstraction_keeps_the_protonation_state(
    amber_ff, caplog, template, residue, heavy, hydrogen, other_state
):
    """A residue that lost a hydrogen has the atom names of another protonation state.

    E.g. HIP minus HE2 has the names of HID, but its NE2 radical is not a
    radical in HID. The residue becomes a radical of its own state and keeps
    its charge (the hydrogen's charge on the radical).
    """
    atoms = _template_residue(amber_ff, template, "1", 1, residue=residue)
    heavy_nr = _nr_of(atoms, 1, 99, heavy)
    h_nr = _nr_of(atoms, 1, 99, hydrogen)
    atoms["101"] = Atom(nr="101", type="OH", resnr="2", residue="DON", atom="O", cgnr="101",
                        charge="0.0000")
    top_a = _RecipeTopology(atoms, amber_ff)
    top_b = _RecipeTopology(deepcopy(atoms), amber_ff)
    _bond(top_b.atoms, heavy_nr, h_nr, False)
    _bond(top_b.atoms, "101", h_nr, True)
    moved = top_b.atoms[h_nr]
    moved.resnr, moved.residue, moved.atom = "2", "DON", "HX"
    top_b.radicals = {heavy_nr: top_b.atoms[heavy_nr]}
    steps = [Break(atom_id_1=heavy_nr, atom_id_2=h_nr), Bind(atom_id_1="101", atom_id_2=h_nr)]

    with caplog.at_level(logging.DEBUG, logger="kimmdy.topology.charges"):
        charges.assign_recipe_charges(top_a, top_b, steps)

    assert f"rtp:{other_state} (exact) does not fit residue 1: radical {heavy}" in caplog.text
    assert f"Charges of residue 1 from rtp:{template} (radical)" in caplog.text
    entry = amber_ff.residuetypes[template].atoms
    assert np.isclose(
        float(top_b.atoms[heavy_nr].charge),
        float(entry[heavy].charge) + float(entry[hydrogen].charge),
    )
    assert np.isclose(_residue_charge(top_b, "1"), _template_charge(amber_ff, template))
    assert np.isclose(_total_charge(top_b), _total_charge(top_a))


def test_radicals_of_a_broken_bond_between_residues_keep_the_exact_entry(amber_ff, caplog):
    """C and N radicals after a C-N break between two Gly are explained by -C N in the entry."""
    atoms = _template_residue(amber_ff, "GLY", "1", 1)
    atoms.update(_template_residue(amber_ff, "GLY", "2", 11))
    c = _nr_of(atoms, 1, 9, "C")
    n = _nr_of(atoms, 11, 19, "N")
    _bond(atoms, c, n, True)
    top_a = _RecipeTopology(atoms, amber_ff)
    top_b = _RecipeTopology(deepcopy(atoms), amber_ff)
    _bond(top_b.atoms, c, n, False)
    top_b.radicals = {c: top_b.atoms[c], n: top_b.atoms[n]}
    steps = [Break(atom_id_1=c, atom_id_2=n)]

    with caplog.at_level(logging.INFO, logger="kimmdy.topology.charges"):
        charges.assign_recipe_charges(top_a, top_b, steps)

    assert "Charges of residue 1 from rtp:GLY (exact)" in caplog.text
    assert "Charges of residue 2 from rtp:GLY (exact)" in caplog.text
    _assert_template_charges(top_b, "1", "GLY")
    _assert_template_charges(top_b, "2", "GLY")
