"""Installed-package contracts exercised without accessing plugin internals.

Only the external engine launch/version probe is replaced. Input writing,
format selection, parsing and the ASE calculator handoff remain real.
"""

from pathlib import Path
import runpy
import shutil

import numpy as np
import pytest
from ase import Atoms

from abacuslite import Abacus, AbacusProfile
from abacuslite.io import latestio, legacyio
from abacuslite.io.generalio import read_input, read_kpt, read_stru
from abacuslite.io.generalio import write_input, write_kpt, write_stru


FIXTURES = Path(__file__).resolve().parents[1] / "fixtures"


class ReplayProfile(AbacusProfile):
    """A controlled launcher which returns recorded engine output."""

    def __init__(self, version, output):
        super().__init__("unused-engine")
        self.engine_version = version
        self.output = output

    def version(self):
        return self.engine_version

    def run(self, directory, inputfile, outputfile, errorfile=None, append=False):
        assert inputfile is None
        assert read_input(Path(directory) / "INPUT")
        destination = Path(directory) / "OUT.ABACUS"
        destination.mkdir(exist_ok=True)
        for path in self.output.iterdir():
            if path.name in {"running_scf.log", "running_md.log", "eig_occ.txt", "MD_dump"}:
                shutil.copy2(path, destination / path.name)


@pytest.mark.parametrize("order", [("legacy", "latest"), ("latest", "legacy")])
def test_independent_calculators_select_formats_through_public_calls(tmp_path, order):
    versions = {"legacy": "v3.10.1", "latest": "v3.11.0-beta10"}
    expected = {"legacy": -1940.5509086572, "latest": -1940.6311063727}
    atoms = {}
    for name in order:
        atoms[name] = Atoms("AsGa", positions=[[1.4, 1.4, 1.4], [0, 0, 0]],
                            cell=[6, 6, 6], pbc=True)
        atoms[name].calc = Abacus(
            profile=ReplayProfile(versions[name], FIXTURES / f"multiframe_md_{name}"),
            directory=tmp_path / name,
            pseudopotentials={"As": "As.upf", "Ga": "Ga.upf"},
            inp={"calculation": "md", "basis_type": "pw", "gamma_only": True},
        )
    for name in reversed(order):
        assert atoms[name].get_potential_energy() == pytest.approx(expected[name])
        assert atoms[name].get_potential_energy(force_consistent=True) == pytest.approx(expected[name])
        assert atoms[name].get_forces().shape == (2, 3)
    assert atoms["legacy"].calc.results["energy_without_entropy"] == pytest.approx(-1940.3851061129)


def test_scf_matches_input_structure_and_restores_interleaved_atom_order(tmp_path):
    atoms = Atoms("HAuH", positions=[[4, 5, 5], [5, 5, 1], [6, 5, 5]],
                  cell=[10, 10, 10], pbc=True)
    atoms.calc = Abacus(
        profile=ReplayProfile("v3.10.1", FIXTURES / "multiframe_scf_trial_last"),
        directory=tmp_path,
        pseudopotentials={"H": "H.upf", "Au": "Au.upf"},
        inp={"basis_type": "pw", "gamma_only": True, "cal_force": True},
    )
    np.testing.assert_allclose(atoms.get_forces(),
                               [[-0.5, -0.01, 0], [0, 0, 0], [0.5, 0.01, 0]])
    assert atoms.calc.results["efermi"] == pytest.approx(-2.04)
    assert atoms.calc.results["fermi_level"] == pytest.approx(-2.04)
    assert "energy_without_entropy" not in atoms.calc.results


def test_scf_refuses_output_from_a_different_structure(tmp_path):
    atoms = Atoms("HHAu", positions=[[1, 1, 1], [2, 2, 2], [3, 3, 3]],
                  cell=[10, 10, 10], pbc=True)
    atoms.calc = Abacus(
        profile=ReplayProfile("v3.10.1", FIXTURES / "multiframe_scf_trial_last"),
        directory=tmp_path, pseudopotentials={"H": "H.upf", "Au": "Au.upf"},
        inp={"basis_type": "pw", "gamma_only": True},
    )
    with pytest.raises(RuntimeError):
        atoms.get_potential_energy()
    assert not atoms.calc.results


def test_core_script_entry_can_collect_scf_results(tmp_path):
    """Keep the existing `python core.py` test entry usable."""
    core_file = Path(__file__).resolve().parents[2] / "abacuslite" / "core.py"
    namespace = runpy.run_path(str(core_file), run_name="script_entry_contract")
    shutil.copy(FIXTURES / "multiframe_scf_trial_last" / "STRU", tmp_path)
    output = tmp_path / "OUT.ABACUS"
    output.mkdir()
    shutil.copy(
        FIXTURES / "multiframe_scf_trial_last" / "running_scf.log", output
    )
    template = namespace["AbacusTemplate"]()
    template.suffix = "ABACUS"
    template.calculation = "scf"
    template.atomorder = [0, 1, 2]
    np.testing.assert_allclose(
        template.read_results(tmp_path)["forces"],
        [[-0.5, -0.01, 0], [0.5, 0.01, 0], [0, 0, 0]],
    )


@pytest.mark.parametrize("name,reader,energy", [
    ("latest", latestio.read_abacus_out, -1940.6311063727),
    ("legacy", legacyio.read_abacus_out, -1940.5509086572),
])
def test_output_readers_are_repeatable_indexed_and_read_only(tmp_path, name, reader, energy):
    output = tmp_path / name
    shutil.copytree(FIXTURES / f"multiframe_md_{name}", output)
    before = {path.name: path.read_bytes() for path in output.iterdir()}
    log = output / "running_md.log"
    frames = reader(log)
    selected = reader(log, index=-1)
    assert isinstance(frames, list) and len(frames) == (2 if name == "latest" else 1)
    assert isinstance(selected, Atoms)
    assert selected.get_potential_energy() == pytest.approx(energy)
    assert selected.get_potential_energy(force_consistent=True) == pytest.approx(energy)
    np.testing.assert_allclose(selected.get_forces(), frames[-1].get_forces())
    reordered = reader(log, index=-1, sort_atoms_with=[1, 0])
    assert reordered.get_chemical_symbols() == selected.get_chemical_symbols()[::-1]
    np.testing.assert_allclose(reordered.get_forces(), selected.get_forces()[::-1])
    assert {path.name: path.read_bytes() for path in output.iterdir()} == before


@pytest.mark.parametrize("name,version,reader", [
    ("latest", "v3.11.0-beta10", latestio.read_abacus_out),
    ("legacy", "v3.10.1", legacyio.read_abacus_out),
])
def test_noncollinear_results_need_no_global_ase_patch(tmp_path, name, version, reader):
    source = FIXTURES / "noncollinear_scf" / name
    # Compare the public reader's values with the file-backed handoff.
    # Independent raw-table assertions belong to the real-engine harness.
    frame = reader(source / "running_scf.log", index=-1)
    moments = frame.get_magnetic_moments()
    assert moments.shape == (2, 3)
    assert np.all(np.linalg.norm(moments, axis=1) > 1e-3)
    energy = frame.get_potential_energy()
    atoms = Atoms("Fe2", positions=[[0, 0, 0], [2.6, 2.6, 2.6]],
                  cell=[5.2, 5.2, 5.2], pbc=True, magmoms=[[1, 0, 0], [0, 1, 0]])
    atoms.calc = Abacus(
        profile=ReplayProfile(version, source), directory=tmp_path,
        pseudopotentials={"Fe": "Fe.upf"},
        inp={"basis_type": "pw", "gamma_only": True, "nspin": 4},
    )
    assert atoms.get_potential_energy() == pytest.approx(energy)
    np.testing.assert_allclose(atoms.calc.results["magmoms"], moments)
    assert atoms.calc.results["eigenvalues"].shape[-1] > 0
    exported = atoms.calc.export_properties()
    assert exported["energy"] == pytest.approx(energy)
    assert "fermi_level" in exported
    assert "magmoms" not in exported
    assert "energy_without_entropy" not in exported
    assert atoms.get_properties(["energy"])["energy"] == pytest.approx(energy)
    # A noncollinear calculation must not make later collinear ASE result
    # conversion reject its scalar per-atom moments.
    collinear = legacyio.read_abacus_out(
        FIXTURES / "multiframe_scf_trial_last" / "running_scf.log", index=-1)
    assert np.asarray(collinear.calc.properties()["magmoms"]).shape == (3,)


def test_scf_handoff_accepts_negative_noncollinear_component(tmp_path):
    source = FIXTURES / "noncollinear_scf" / "latest"

    class NegativeMomentProfile(ReplayProfile):
        def run(self, directory, inputfile, outputfile, errorfile=None, append=False):
            super().run(directory, inputfile, outputfile, errorfile, append)
            path = Path(directory) / "STRU"
            path.write_text(path.read_text().replace(
                "mag 0.0 1.0 0.0", "mag 0.0 -1.0 0.0"
            ))

    atoms = Atoms("Fe2", positions=[[0, 0, 0], [2.6, 2.6, 2.6]],
                  cell=[5.2, 5.2, 5.2], pbc=True)
    atoms.calc = Abacus(
        profile=NegativeMomentProfile("v3.11.0-beta10", source),
        directory=tmp_path,
        pseudopotentials={"Fe": "Fe.upf"},
        inp={"basis_type": "pw", "gamma_only": True, "nspin": 4},
    )
    atoms.get_potential_energy()
    assert atoms.calc.results["magmoms"].shape == (2, 3)


def test_input_io_without_an_engine_preserves_declared_values(tmp_path):
    write_input({"suffix": "public", "scf_thr": 1e-8}, tmp_path / "INPUT")
    assert read_input(tmp_path / "INPUT") == {"suffix": "public", "scf_thr": "1e-08"}
    sampling = {"mode": "mp-sampling", "gamma-centered": True,
                "nk": (2, 3, 4), "kshift": (0, 0, 0)}
    write_kpt(sampling, tmp_path / "KPT")
    assert read_kpt(tmp_path / "KPT") == sampling
    atoms = Atoms("HHeH", positions=[[0, 0, 0], [1, 1, 1], [2, 2, 2]],
                  cell=[5, 5, 5], pbc=True)
    atoms.set_masses([2.014, 4.003, 2.014])
    write_stru(atoms, tmp_path, {"H": "H.upf", "He": "He.upf"})
    data = read_stru(tmp_path / "STRU")
    assert [item["symbol"] for item in data["species"]] == ["H", "He"]
    assert [item["mass"] for item in data["species"]] == pytest.approx([2.014, 4.003])


@pytest.mark.parametrize("module,header,marker", [
    (latestio, "--> #ION MOVE# 1 #ELEC ITER#", "#"),
    (legacyio, "LCAO ALGORITHM ---- ION=1 ELEC=", ""),
])
@pytest.mark.parametrize("source_kind", ["path", "string", "lines"])
def test_selected_log_fields_keep_units_shapes_and_missing_values(
        tmp_path, module, header, marker, source_kind):
    # Literal independent tables: two electronic iterations, one ionic step.
    content = f"""{header} 1
ENERGY Rydberg eV
E_KohnSham -1.0 -13.6
end
{header} 2
ENERGY Rydberg eV
E_KohnSham -2.0 -27.2
E_Fermi 0.1 1.36
end
{marker}TOTAL-FORCE (eV/Angstrom){marker}
Si1 1.0 -2.0 3.0
Si2 -1.0 2.0 -3.0
end
{marker}TOTAL-STRESS (kbar){marker}
1.0 0.0 0.0
0.0 2.0 0.0
0.0 0.0 -3.0
end
"""
    log = tmp_path / "running_scf.log"
    log.write_text(content)
    src = {"path": log, "string": str(log), "lines": content.splitlines()}[source_kind]
    rydberg, ev = module.read_energies_from_running_log(src)
    assert rydberg[-1] == {"E_KohnSham": -2.0, "E_Fermi": 0.1}
    assert ev == [{"E_KohnSham": -13.6}, {"E_KohnSham": -27.2, "E_Fermi": 1.36}]
    headers = module.read_iter_header_from_running_log(src)
    assert headers == [(1, 1), (1, 2)]
    assert module.find_final_info_with_iter_header(ev, headers) == [ev[-1]]
    np.testing.assert_allclose(module.read_forces_from_running_log(src),
                               [[[1, -2, 3], [-1, 2, -3]]])
    np.testing.assert_allclose(module.read_stress_from_running_log(src),
                               [np.diag([-0.0006241509, -0.0012483018, 0.0018724527])],
                               rtol=1e-6, atol=0)
    assert module.read_energies_from_running_log(["no energy table"]) == ([], [])
    assert module.read_forces_from_running_log(["no force table"]) == []
    assert module.read_stress_from_running_log(["no stress table"]) == []
