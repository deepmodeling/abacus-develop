#!/usr/bin/env python3
"""
Basic Example: LCAO Basis with Multi-k (Standard Workflow)

This represents the most common use case:
- LCAO basis (requires orbital files)
- Multi-k point mesh (mp_grid > 1x1x1)
- Automated execution via `run()`
"""

import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from abacusw90 import ABACUSWannier90

def main():
    print("=== ABACUS Wannier90 Example: Basic (LCAO Multi-k) ===")

    # 1. Initialization
    job = ABACUSWannier90(
        work_dir="./Bi2Se3_basic",
        scf_dir="./Bi2Se3_scf"
    )

    # 2. Structure
    lattice = [
        [-2.069, -3.583614, 0.0], [2.069, -3.583614, 0.0], [0.0, 2.389075, 9.546667]
    ]
    atoms = [
        {"name": "Bi", "pos": [0.399, 0.399, 0.697]}, {"name": "Bi", "pos": [0.601, 0.601, 0.303]},
        {"name": "Se", "pos": [0.000, 0.000, 0.500]}, {"name": "Se", "pos": [0.206, 0.206, 0.118]},
        {"name": "Se", "pos": [0.794, 0.794, 0.882]},
    ]
    job.set_structure(lattice, atoms)

    # 3. Dependencies
    job.pp_orbitals = {"Bi": "Bi.upf", "Se": "Se.upf"}
    job.orbital_files = ["Bi.orb", "Se.orb"]

    # 4. Wannier90 Parameters
    job.set_wannier_parameters(
        num_wann=30, num_bands=100,
        projections=["Bi : pz; px; py", "Se : pz; px; py"],
        dis_win_min=3.0, dis_win_max=18.0,
        dis_froz_min=3.0, dis_froz_max=14.8,
        mp_grid=[4, 4, 4]
    )

    # 5. ABACUS Parameters (LCAO Multi-k)
    job.set_abacus_parameters(
        ecutwfc=100, nbands=100,
        basis_type="lcao",
        ks_solver="genelpa",
        gamma_only=0,        # Explicitly False for Multi-k
        nspin=4, lspinorb=1
    )

    # 6. Run
    job.run()

if __name__ == "__main__":
    main()
