#!/usr/bin/env python3
"""
Example: LCAO Basis - Gamma-only

Key differences:
- gamma_only = 1
- mp_grid = [1, 1, 1]
- Much faster for large systems
"""

import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from abacuswannier import ABACUSWannier90

def main():
    print("=== ABACUS Wannier90 Example: LCAO Gamma-only ===")
    job = ABACUSWannier90(work_dir="./Bi2Se3_gamma", scf_dir="./Bi2Se3_scf_gamma")

    # Structure
    lattice = [[-2.069, -3.583614, 0.0], [2.069, -3.583614, 0.0], [0.0, 2.389075, 9.546667]]
    atoms = [
        {"name": "Bi", "pos": [0.399, 0.399, 0.697]}, {"name": "Bi", "pos": [0.601, 0.601, 0.303]},
        {"name": "Se", "pos": [0.000, 0.000, 0.500]}, {"name": "Se", "pos": [0.206, 0.206, 0.118]},
        {"name": "Se", "pos": [0.794, 0.794, 0.882]},
    ]
    job.set_structure(lattice, atoms)

    job.pp_orbitals = {"Bi": "Bi.upf", "Se": "Se.upf"}
    job.orbital_files = ["Bi.orb", "Se.orb"]

    # Wannier90 Settings for Gamma
    job.set_wannier_parameters(
        num_wann=30, num_bands=100,
        projections=["Bi : pz; px; py", "Se : pz; px; py"],
        dis_win_min=3.0, dis_win_max=18.0,
        dis_froz_min=3.0, dis_froz_max=14.8,
        mp_grid=[1, 1, 1]  # Gamma only
    )

    # ABACUS Settings for Gamma
    job.set_abacus_parameters(
        ecutwfc=100, nbands=100,
        basis_type="lcao", ks_solver="genelpa",
        gamma_only=1,        # Enable Gamma-only
        nspin=4, lspinorb=1
    )

    job.run()

if __name__ == "__main__":
    main()
