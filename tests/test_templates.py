"""Tests for the residue templates (kimmdy.topology.templates) and their use in the topology."""

import logging
from copy import deepcopy
from pathlib import Path
from typing import Optional

import pytest

from kimmdy.parsing import read_top
from kimmdy.recipe import Bind, Break
from kimmdy.topology.atomic import Atom
from kimmdy.topology.templates import (
    ResidueTemplate,
    TemplateAtom,
    TemplateLibrary,
    assign_recipe_names,
    match_residue,
)
from kimmdy.topology.topology import Topology
from kimmdy.topology.utils import get_residue_by_bonding

FILEDIR = Path(__file__).parent / "test_files" / "test_topology"


@pytest.fixture(scope="module")
def imree_raw():
    return read_top(FILEDIR / "IMREE.top")


@pytest.fixture(scope="module")
def amber_ff(imree_raw):
    """Force field (with the amber99sb-star-ildnp residue templates) of IMREE.top."""
    return Topology(deepcopy(imree_raw)).ff


@pytest.fixture
def hexala() -> Topology:
    return Topology(read_top(FILEDIR / "hexala.top"))


@pytest.fixture(scope="module")
def r10() -> Topology:
    """One SSGP-[GRGDSPYS]10-GY chain from pdb2gmx with a99SB-disp (rtp NSER ... CTYR)."""
    return Topology(
        read_top(FILEDIR / "R10_a99SBdisp.top", ffdir=FILEDIR / "a99SBdisp.ff")
    )


def _r10_residue(top: Topology, resnr: str) -> dict[str, Atom]:
    atom = next(a for a in top.atoms.values() if a.resnr == resnr)
    return get_residue_by_bonding(atom, top.atoms)


def _residue(ff, template: str, residue: Optional[str] = None) -> dict[str, Atom]:
    """Atoms 1, 2, ... of an rtp entry, with its types, charges and bonds; residue is the name in the topology."""
    residuetype = ff.residuetypes[template]
    nrs = {name: str(i + 1) for i, name in enumerate(residuetype.atoms)}
    atoms = {
        nrs[name]: Atom(
            nr=nrs[name],
            type=spec.type,
            resnr="1",
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


def _nr(atoms: dict[str, Atom], name: str) -> str:
    return next(nr for nr, atom in atoms.items() if atom.atom == name)


def _abstract(atoms: dict[str, Atom], name: str) -> str:
    """Remove a hydrogen; its partner keeps its charge (radical convention). Returns the partner."""
    h = _nr(atoms, name)
    partner = atoms[h].bound_to_nrs[0]
    atoms[partner].bound_to_nrs.remove(h)
    atoms[partner].charge = f"{float(atoms[partner].charge) + float(atoms.pop(h).charge):7.4f}"
    return partner


def _add_hydrogen(atoms: dict[str, Atom], partner: str, name: str, nr: str = "99") -> str:
    atoms[nr] = Atom(
        nr=nr, type="HO", resnr="1", residue=atoms[partner].residue, atom=name, cgnr=nr,
        charge="0.4000",
    )
    atoms[nr].bound_to_nrs.append(partner)
    atoms[partner].bound_to_nrs.append(nr)
    return nr


def test_match_levels(amber_ff):
    """A Gly is an exact match; without HA2 only a radical match, with CA as the radical."""
    template = amber_ff.templates.rtp["GLY"]
    atoms = _residue(amber_ff, "GLY")
    match = match_residue(atoms, template)
    assert match is not None
    assert match.names == {nr: atom.atom for nr, atom in atoms.items()}
    assert match.missing == {}

    ca = _abstract(atoms, "HA2")
    assert match_residue(atoms, template, "exact") is None
    assert match_residue(atoms, template, "renamed hydrogen") is None
    match = match_residue(atoms, template, "radical")
    assert match.missing == {"HA2": ca}
    assert match.radicals == {ca}


def test_match_needs_the_bonds_of_the_template(amber_ff):
    """The atom names alone do not match: HA2 bonded to N instead of CA does not fit."""
    template = amber_ff.templates.rtp["GLY"]
    atoms = _residue(amber_ff, "GLY")
    ha2 = _nr(atoms, "HA2")
    _abstract(atoms, "HA2")
    _add_hydrogen(atoms, _nr(atoms, "N"), "HA2", nr=ha2)
    assert match_residue(atoms, template, "radical") is None


def test_renamed_hydrogens_take_a_free_position_on_their_partner(amber_ff):
    template = amber_ff.templates.rtp["GLY"]
    atoms = _residue(amber_ff, "GLY")
    ca = _abstract(atoms, "HA2")
    h = _add_hydrogen(atoms, ca, "HX")
    assert match_residue(atoms, template, "exact", renamed={h}) is None
    match = match_residue(atoms, template, "renamed hydrogen", renamed={h})
    assert match.names[h] == "HA2"

    # no free position on C
    atoms = _residue(amber_ff, "GLY")
    _abstract(atoms, "HA2")
    h = _add_hydrogen(atoms, _nr(atoms, "C"), "HX")
    assert match_residue(atoms, template, "radical", renamed={h}) is None


def test_match_compares_types_when_asked(amber_ff):
    template = amber_ff.templates.rtp["GLY"]
    atoms = _residue(amber_ff, "GLY")
    atoms[_nr(atoms, "CA")].type = "CX"
    assert match_residue(atoms, template) is not None
    assert match_residue(atoms, template, same_types=True) is None


def test_rtp_template_by_residue_name_and_atom_names(amber_ff):
    """A residue named HIS has no entry of its own; its atom names choose HIE or HIP."""
    templates = amber_ff.templates
    assert templates.rtp_template(_residue(amber_ff, "ALA")).name == "ALA"
    assert templates.rtp_template(_residue(amber_ff, "HIE", "HIS")).name == "HIE"
    assert templates.rtp_template(_residue(amber_ff, "HIP", "HIS")).name == "HIP"
    # a hydrogen without a template position does not decide
    atoms = _residue(amber_ff, "HIE", "HIS")
    _add_hydrogen(atoms, _nr(atoms, "CB"), "HX")
    assert templates.rtp_template(atoms).name == "HIE"
    assert templates.rtp_template(_residue(amber_ff, "GLY", "XYZ")) is None


def test_hydrogen_identity_from_the_partner_position(amber_ff):
    """A hydrogen binding to CA takes the free HA position, or is HX with the HA type."""
    templates = amber_ff.templates
    hydrogen = Atom(nr="99", type="HO", resnr="2", residue="DON", atom="HO", cgnr="99", charge="0.4")

    atoms = _residue(amber_ff, "ALA")
    ca = _nr(atoms, "CA")
    assert templates.hydrogen_identity(atoms, atoms[ca], hydrogen) == ("HX", "H1")

    _abstract(atoms, "HA")
    assert templates.hydrogen_identity(atoms, atoms[ca], hydrogen) == ("HA", "H1")
    # no hydrogen position on C: HX, the type from C (bonded to O: H1)
    c = _nr(atoms, "C")
    assert templates.hydrogen_identity(atoms, atoms[c], hydrogen) == ("HX", "H1")
    # no rtp entry
    atoms = _residue(amber_ff, "ALA", "XYZ")
    assert templates.hydrogen_identity(atoms, atoms[ca], hydrogen) is None


@pytest.mark.parametrize(
    "template, residue, abstracted, atom, radical",
    [
        ("LYS", "LYS", "HZ3", "NZ", True),
        ("LYN", "LYS", None, "NZ", False),
        ("CYS", "CYS", "HG", "SG", True),
        ("CYM", "CYS", None, "SG", False),
    ],
)
def test_radical_by_template_uses_the_residue_charge(
    amber_ff, template, residue, abstracted, atom, radical
):
    """LYS without HZ3 has the atoms of LYN: its charge (+1, radical convention) makes NZ a radical."""
    atoms = _residue(amber_ff, template, residue)
    if abstracted is not None:
        _abstract(atoms, abstracted)
    nr = _nr(atoms, atom)
    assert amber_ff.templates.radical_by_template(atoms, atoms[nr]) is radical


def test_radical_by_template_without_a_fitting_entry(amber_ff):
    atoms = _residue(amber_ff, "LYS")
    nz = _abstract(atoms, "HZ3")
    # a residue charge of 1.5: neither LYS (+1) nor LYN (0)
    atoms[nz].charge = f"{float(atoms[nz].charge) + 0.5:7.4f}"
    assert amber_ff.templates.radical_by_template(atoms, atoms[nz]) is None


def _itp_residue(name: str, atoms: list[tuple[str, str, str]]) -> ResidueTemplate:
    """An ITP template of residue OHR: (name, type, charge) of O1 and its hydrogens."""
    return ResidueTemplate(
        source="itp",
        name=name,
        residue="OHR",
        atoms={n: TemplateAtom(name=n, type=t, charge=q) for n, t, q in atoms},
        bonds=frozenset(("O1", n) if "O1" < n else (n, "O1") for n, _, _ in atoms[1:]),
        linked=frozenset(),
        boundary_bonds=frozenset(),
        local_resnr="1",
    )


WATER = [("O1", "OT", "-0.834"), ("H1", "HT", "0.417"), ("H2", "HT", "0.417")]


def test_transition_to_a_copy_of_a_template_uses_the_kept_one(monkeypatch, caplog):
    """A moleculetype dropped as a copy of another can still be named in a transition, with a warning."""
    # earlier tests may configure the kimmdy logger with propagate=False
    monkeypatch.setattr(logging.getLogger("kimmdy"), "propagate", True)
    library = TemplateLibrary()
    library.add_itp(_itp_residue("OHR", [("O1", "ORA", "-0.436"), ("H1", "HRA", "0.436")]))
    library.add_itp(_itp_residue("OHRH", WATER))
    library.add_itp(_itp_residue("OHRH_COPY", WATER))
    assert [template.name for template in library.itp["OHR"]] == ["OHR", "OHRH"]

    with caplog.at_level(logging.WARNING, logger="kimmdy.topology.templates"):
        library.set_transitions({"OHR": "OHRH_COPY"})

    assert library.transitions == {"OHR": "OHRH"}
    assert "OHRH_COPY has the same residues as OHRH" in caplog.text
    with pytest.raises(ValueError, match="Unknown moleculetype.*MISSING"):
        library.set_transitions({"OHR": "MISSING"})


def _water_residue() -> dict[str, Atom]:
    """Atoms 1-3 of an OHR residue with the atoms, types and bonds of WATER."""
    atoms = {
        str(i): Atom(nr=str(i), type=t, resnr="1", residue="OHR", atom=n, cgnr=str(i), charge=q)
        for i, (n, t, q) in enumerate(WATER, start=1)
    }
    for h in ("2", "3"):
        atoms["1"].bound_to_nrs.append(h)
        atoms[h].bound_to_nrs.append("1")
    return atoms


def test_itp_copies_with_differently_written_numbers_are_copies():
    library = TemplateLibrary()
    library.add_itp(_itp_residue("OHRH", WATER))
    library.add_itp(
        _itp_residue("OHRH_COPY", [("O1", "OT", "-0.834000"), ("H1", "HT", "0.4170"), ("H2", "HT", ".417")])
    )
    assert [template.name for template in library.itp["OHR"]] == ["OHRH"]


def test_several_fitting_itp_templates_use_the_first(monkeypatch, caplog):
    """Two charge sets of the same atoms both fit; the first is used, with a warning."""
    # earlier tests may configure the kimmdy logger with propagate=False
    monkeypatch.setattr(logging.getLogger("kimmdy"), "propagate", True)
    library = TemplateLibrary()
    library.add_itp(_itp_residue("OHRH", WATER))
    library.add_itp(
        _itp_residue("OHRH_ALT", [("O1", "OT", "-0.800"), ("H1", "HT", "0.400"), ("H2", "HT", "0.400")])
    )
    residue = _water_residue()

    with caplog.at_level(logging.WARNING, logger="kimmdy.topology.templates"):
        match = library.match_itp(residue, set(), residue)

    assert match.template.name == "OHRH"
    assert "Several ITP templates fit residue 1 OHR" in caplog.text
    assert "the first, itp:OHRH, is used" in caplog.text


def test_find_radicals_decides_ambiguous_atoms_by_template(imree_raw):
    """Arg NE (type N2, 2 or 3 bonds) without HE has 2 bonds: a radical by the ARG entry, not by its type."""
    top = Topology(deepcopy(imree_raw))
    he = next(nr for nr, a in top.atoms.items() if a.residue == "ARG" and a.atom == "HE")
    top.del_atom(he)
    ne = next(nr for nr, a in top.atoms.items() if a.residue == "ARG" and a.atom == "NE")
    mol = top.reactive_molecule
    assert len(mol.atoms[ne].bound_to_nrs) == 2

    mol.find_radicals()
    assert ne not in mol.radicals

    mol.find_radicals(top.ff)
    assert ne in mol.radicals
    assert mol.atoms[ne].is_radical


def test_update_residue_impropers_restores_the_template_impropers(hexala):
    mol = hexala.reactive_molecule
    before = deepcopy(mol.improper_dihedrals)
    with_h = [key for key in mol.improper_dihedrals if "8" in key]
    assert with_h  # -C CA N H of residue 2
    for key in with_h:
        mol.improper_dihedrals.pop(key)

    residue = get_residue_by_bonding(mol.atoms["8"], mol.atoms)
    mol.update_residue_impropers([residue], hexala.ff)

    assert mol.improper_dihedrals == before


def test_update_all_names_swapped_hydrogens_and_sets_their_impropers(hexala):
    """Amide H and HB1 of Ala 2 swap places: the templates name them by their new positions.

    bind_bond leaves the names and the impropers; update_all renames both
    and sets the -C CA N H improper on the new amide H.
    """
    top = hexala
    top_a = deepcopy(top)
    steps = [
        Break(atom_id_1="7", atom_id_2="8"),
        Break(atom_id_1="11", atom_id_2="12"),
        Bind(atom_id_1="7", atom_id_2="12"),
        Bind(atom_id_1="11", atom_id_2="8"),
    ]
    for step in steps:
        if isinstance(step, Break):
            top.break_bond((step.atom_id_1, step.atom_id_2))
        else:
            top.bind_bond((step.atom_id_1, step.atom_id_2))
    assert (top.atoms["12"].atom, top.atoms["8"].atom) == ("HB1", "H")

    top.update_all(top_a, steps, conservation_check="error")

    assert (top.atoms["12"].atom, top.atoms["12"].type) == ("H", "H")
    assert (top.atoms["8"].atom, top.atoms["8"].type) == ("HB1", "HC")
    assert ("5", "9", "7", "12") in top.improper_dihedrals
    assert ("5", "9", "7", "8") not in top.improper_dihedrals


class _Top:
    """The parts of a Topology that assign_recipe_names uses."""

    def __init__(self, atoms: dict[str, Atom], ff, radicals: tuple[str, ...] = ()):
        self.atoms = atoms
        self.ff = ff
        self.radicals = {nr: atoms[nr] for nr in radicals}


def test_names_come_from_the_template_whatever_the_charges(amber_ff):
    """A hydroxyl H that moves onto the CA radical of a Gly becomes HA2, also if the Gly's charges are not the template's."""
    atoms = _residue(amber_ff, "GLY")
    ca = _abstract(atoms, "HA2")
    atoms[_nr(atoms, "C")].charge = "0.7731"  # GLY C is 0.5973
    atoms["31"] = Atom(nr="31", type="OH", resnr="2", residue="DON", atom="O", cgnr="31",
                       charge="-0.4000")
    h = _add_hydrogen(atoms, "31", "HO", nr="32")
    atoms[h].resnr, atoms[h].residue = "2", "DON"
    top_a = _Top(atoms, amber_ff, radicals=(ca,))

    top_b = _Top(deepcopy(atoms), amber_ff)
    b = top_b.atoms
    b["31"].bound_to_nrs.remove(h)
    b[h].bound_to_nrs = [ca]
    b[ca].bound_to_nrs.append(h)
    b[h].resnr, b[h].residue = "1", "GLY"  # as bind_bond leaves it
    top_b.radicals = {"31": b["31"]}
    steps = [Break(atom_id_1="31", atom_id_2=h), Bind(atom_id_1=ca, atom_id_2=h)]

    residues = assign_recipe_names(top_a, top_b, steps)

    assert sorted(residue.resnr for residue in residues) == ["1", "2"]
    assert top_a.atoms[h].atom != b[h].atom  # h got another name
    assert (b[h].atom, b[h].type) == ("HA2", "H1")
    assert b[h].mass == amber_ff.atomtypes["H1"].mass


def _hydroxyl_hydrogen() -> Atom:
    return Atom(nr="99", type="HO", resnr="2", residue="DON", atom="HO", cgnr="99", charge="0.4")


@pytest.mark.parametrize(
    "residue, partner, h_type",
    [
        ("ALA", "C", "H1"),  # acyl C after a backbone break: bonded to O
        ("PHE", "CG", "HA"),  # aromatic CA: the rtp entries put HA on it
        ("PRO", "N", "H"),  # N: the rtp entries put H on it
    ],
)
def test_hydrogen_without_a_template_position_gets_the_type_of_its_partner(
    amber_ff, residue, partner, h_type
):
    """A hydrogen where the rtp entry has none is HX, with the type for its partner, not its old one."""
    atoms = _residue(amber_ff, residue)
    nr = _nr(atoms, partner)
    identity = amber_ff.templates.hydrogen_identity(atoms, atoms[nr], _hydroxyl_hydrogen())
    assert identity == ("HX", h_type)


def test_hydrogen_type_on_an_unknown_partner_type_is_none(amber_ff):
    """A partner type the force field does not have gives no type; the hydrogen keeps its own."""
    templates = amber_ff.templates
    atoms = _residue(amber_ff, "ALA")
    c = _nr(atoms, "C")
    atoms[c].type = "XX"
    hydrogen = _hydroxyl_hydrogen()
    assert templates.hydrogen_type_on(atoms[c], hydrogen, atoms) is None
    assert templates.hydrogen_identity(atoms, atoms[c], hydrogen) == ("HX", "HO")


def test_hydrogen_moving_into_a_residue_without_rtp_entry_gets_the_type_of_its_partner(amber_ff):
    """A hydroxyl H that moves onto the CB radical of an Ala named XYZ keeps its name and becomes HC."""
    atoms = _residue(amber_ff, "ALA", "XYZ")
    cb = _abstract(atoms, "HB3")
    atoms["31"] = Atom(nr="31", type="OH", resnr="2", residue="DON", atom="O", cgnr="31",
                       charge="-0.4000")
    h = _add_hydrogen(atoms, "31", "HO", nr="32")
    atoms[h].resnr, atoms[h].residue = "2", "DON"
    top_a = _Top(atoms, amber_ff, radicals=(cb,))

    top_b = _Top(deepcopy(atoms), amber_ff)
    b = top_b.atoms
    b["31"].bound_to_nrs.remove(h)
    b[h].bound_to_nrs = [cb]
    b[cb].bound_to_nrs.append(h)
    b[h].resnr, b[h].residue = "1", "XYZ"  # as bind_bond leaves it
    top_b.radicals = {"31": b["31"]}
    steps = [Break(atom_id_1="31", atom_id_2=h), Bind(atom_id_1=cb, atom_id_2=h)]

    residues = assign_recipe_names(top_a, top_b, steps)

    assert sorted(residue.resnr for residue in residues) == ["1", "2"]
    assert top_a.atoms[h].atom == b[h].atom  # h keeps its name
    assert (b[h].atom, b[h].type) == ("HO", "HC")
    assert b[h].mass == amber_ff.atomtypes["HC"].mass


def test_rtp_template_picks_terminal_entries_by_atoms(r10):
    """pdb2gmx names the terminal residues SER and TYR; their atoms and bonds pick NSER and CTYR."""
    templates = r10.ff.templates
    assert templates.rtp_template(_r10_residue(r10, "1")).name == "NSER"
    assert templates.rtp_template(_r10_residue(r10, "2")).name == "SER"
    assert templates.rtp_template(_r10_residue(r10, "86")).name == "CTYR"


def test_rtp_template_matches_a_radical_without_its_hydrogen(r10):
    """A SER without HA (abstracted) matches SER at the radical level."""
    residue = _r10_residue(r10, "2")
    radical = {nr: a for nr, a in residue.items() if a.atom != "HA"}
    assert r10.ff.templates.rtp_template(radical).name == "SER"


def test_rtp_template_is_none_if_no_entry_fits(r10):
    """A SER without OG (a heavy atom) fits no entry."""
    residue = _r10_residue(r10, "2")
    broken = {nr: a for nr, a in residue.items() if a.atom != "OG"}
    assert r10.ff.templates.rtp_template(broken) is None


def test_update_names_LJ_sets_the_impropers_of_a_bind_without_renaming(hexala):
    """Break and Bind of the amide H 8 of Ala 2 rename nothing; update_names_LJ still sets its impropers."""
    top = hexala
    top_a = deepcopy(top)
    steps = [Break(atom_id_1="7", atom_id_2="8"), Bind(atom_id_1="7", atom_id_2="8")]
    top.break_bond(("7", "8"))
    top.bind_bond(("7", "8"))
    assert ("5", "9", "7", "8") not in top.improper_dihedrals  # left to update_names_LJ

    top.update_names_LJ(top_a, steps)

    assert top.improper_dihedrals == top_a.improper_dihedrals
