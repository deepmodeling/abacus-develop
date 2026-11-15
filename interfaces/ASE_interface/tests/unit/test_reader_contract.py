"""Regression tests for abacuslite output reader atom reordering."""

from __future__ import annotations

import numpy as np
import pytest

from abacuslite.io import latestio, legacyio


def _mock_reader_data():
    frame = {
        "elem": ["Na", "Na", "Cl"],
        "coords": np.array(
            [
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [2.0, 0.0, 0.0],
            ]
        ),
        "cell": np.eye(3),
    }
    elecstate = [
        {
            "k": np.zeros((1, 1, 3)),
            "e": np.zeros((1, 1, 1)),
            "occ": np.ones((1, 1, 1)),
        }
    ]
    energies = [{"E_KohnSham": -1.0, "E_Fermi": 0.0}]
    kpoints = ((np.zeros((1, 3)), np.ones(1), None), None, None, None, None)
    magmoms = [np.array([10.0, 20.0, 30.0])]
    forces = [np.array([[1.0, 0.0, 0.0], [2.0, 0.0, 0.0], [3.0, 0.0, 0.0]])]
    return frame, elecstate, energies, kpoints, magmoms, forces


def _patch_reader(module, monkeypatch, frame, elecstate, energies, kpoints, magmoms, forces):
    monkeypatch.setattr(module, "read_esolver_type_from_running_log", lambda lines: "ksdft")
    monkeypatch.setattr(module, "read_traj_from_running_log", lambda lines: [frame])
    monkeypatch.setattr(module, "read_forces_from_running_log", lambda lines: forces)
    monkeypatch.setattr(module, "read_stress_from_running_log", lambda lines: [])
    monkeypatch.setattr(module, "read_kpoints_from_running_log", lambda lines: kpoints)
    monkeypatch.setattr(module, "read_energies_from_running_log", lambda lines: ([], energies))
    monkeypatch.setattr(module, "read_iter_header_from_running_log", lambda lines: [])
    monkeypatch.setattr(module, "find_final_info_with_iter_header", lambda rows, headers: energies)
    monkeypatch.setattr(module, "read_magmom_from_running_log", lambda lines: magmoms)
    if module is latestio:
        monkeypatch.setattr(module, "read_band_from_eig_occ", lambda path: elecstate)
    else:
        monkeypatch.setattr(module, "read_band_from_running_log", lambda lines: elecstate)


@pytest.mark.parametrize("module", [latestio, legacyio])
@pytest.mark.parametrize("missing_fermi", [False, True])
def test_read_abacus_out_reorders_calculator_magmoms_with_atoms(module, missing_fermi, tmp_path, monkeypatch):
    """Reader result magmoms should follow reordered Atoms indices."""
    running_log = tmp_path / "running_scf.log"
    running_log.write_text("", encoding="utf-8")
    (tmp_path / "eig_occ.txt").write_text("", encoding="utf-8")
    frame, elecstate, energies, kpoints, magmoms, forces = _mock_reader_data()
    if missing_fermi:
        energies[0].pop("E_Fermi")
    _patch_reader(module, monkeypatch, frame, elecstate, energies, kpoints, magmoms, forces)

    atoms = module.read_abacus_out(running_log, sort_atoms_with=[0, 2, 1])[0]

    assert atoms.get_chemical_symbols() == ["Na", "Cl", "Na"]
    np.testing.assert_allclose(atoms.get_initial_magnetic_moments(), [10.0, 30.0, 20.0])
    np.testing.assert_allclose(atoms.calc.results["magmoms"], [10.0, 30.0, 20.0])
    # 力与原子重排一致（帧选择依赖 results 与 atoms 同序，spec P3 ①）
    np.testing.assert_allclose(
        atoms.calc.results["forces"],
        [[1.0, 0.0, 0.0], [3.0, 0.0, 0.0], [2.0, 0.0, 0.0]],
    )

    assert atoms.calc.get_fermi_level() == (None if missing_fermi else 0.0)
