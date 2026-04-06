#!/usr/bin/env python3
"""
Example: PW (Plane Wave) Basis

Key differences:
- basis_type = 'pw'
- ks_solver = 'cg'
- No orbital files needed
"""

import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from abacusw90 import ABACUSWannier90

def main():
    print("=== ABACUS Wannier90 Example: PW Basis ===")
    job = ABACUSWannier90(work_dir="./Bi2Se3_pw", scf_dir="./Bi2Se3_scf_pw")

    # Structure
    lattice = [[-2.069, -3.583614, 0.0], [2.069, -3.583614, 0.0], [0.0, 2.389075, 9.546667]]
    atoms = [
        {"name": "Bi", "pos": [0.399, 0.399, 0.697]}, {"name": "Bi", "pos": [0.601, 0.601, 0.303]},
        {"name": "Se", "pos": [0.000, 0.000, 0.500]}, {"name": "Se", "pos": [0.206, 0.206, 0.118]},
        {"name": "Se", "pos": [0.794, 0.794, 0.882]},
    ]
    job.set_structure(lattice, atoms)

    # Only PP needed
    job.pp_orbitals = {"Bi": "Bi.upf", "Se": "Se.upf"}

    job.set_wannier_parameters(
        num_wann=30, num_bands=100, mp_grid=[4, 4, 4],
        projections=["Bi : pz; px; py", "Se : pz; px; py"],
        dis_win_min=3.0, dis_win_max=18.0, dis_froz_min=3.0, dis_froz_max=14.8
    )

    # PW Settings
    job.set_abacus_parameters(
        ecutwfc=100, nbands=100,
        basis_type="pw", ks_solver="cg",
        nspin=4, lspinorb=1
    )

    job.run()

if __name__ == "__main__":
    main()
