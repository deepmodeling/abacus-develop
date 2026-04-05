#!/usr/bin/env python3
"""
Basic Example: Bi2Se3 Wannierization

This script demonstrates the standard workflow to generate Wannier functions
using the ABACUS-Wannier90 interface. It assumes the user has already 
completed the SCF calculation and obtained the charge density files.
"""

import os
import sys

# Add parent directory to path to allow importing the package
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from abacusw90 import ABACUSWannier90

def main():
    print("="*60)
    print(" ABACUS Wannier90 Basic Example: Bi2Se3")
    print("="*60)

    # ------------------------------------------------------------------
    # 1. Initialization
    # ------------------------------------------------------------------
    # Define working directory and SCF directory
    # scf_dir should contain SPIN1~4_CHG.cube and ABACUS output files.
    job = ABACUSWannier90(
        work_dir="./Bi2Se3_wannier_basic",
        scf_dir="./Bi2Se3_scf", 
        abacus_exe="abacus",            # Or "mpirun -np 4 abacus"
        wannier90_exe="wannier90.x"
    )

    # ------------------------------------------------------------------
    # 2. Define Structure (Bi2Se3)
    # ------------------------------------------------------------------
    lattice = [
        [-2.069, -3.583614, 0.000000],
        [ 2.069, -3.583614, 0.000000],
        [ 0.000,  2.389075, 9.546667]
    ]
    
    atoms = [
        {"name": "Bi", "pos": [0.399, 0.399, 0.697]},
        {"name": "Bi", "pos": [0.601, 0.601, 0.303]},
        {"name": "Se", "pos": [0.000, 0.000, 0.500]},
        {"name": "Se", "pos": [0.206, 0.206, 0.118]},
        {"name": "Se", "pos": [0.794, 0.794, 0.882]},
    ]
    job.set_structure(lattice, atoms)

    # Setup pseudopotentials and orbitals (modify paths to your actual files)
    # These files must exist in the current directory or specified paths
    job.pp_orbitals = {
        "Bi": "./Bi.oncvpsp.upf",
        "Se": "./Se.oncvpsp.upf"
    }
    job.orbital_files = [
        "./Bi_gga_8.0au_100Ry_4s2p2d.orb",
        "./Se_gga_8.0au_100Ry_4s2p2d.orb"
    ]

    # ------------------------------------------------------------------
    # 3. Configure Wannier90 Parameters
    # ------------------------------------------------------------------
    # These parameters are determined by analyzing the band structure (Step 2 in tutorial)
    job.set_wannier_parameters(
        num_wann=30,
        num_bands=100,
        projections=[
            "Bi : pz; px; py", 
            "Se : pz; px; py"
        ],
        # Energy windows (from tutorial)
        dis_win_min=3.0,
        dis_win_max=18.0,
        dis_froz_min=3.0,
        dis_froz_max=14.8,
        # K-grid
        mp_grid=[4, 4, 4],
        # Band path for plotting
        kpath=[
            {"start_label": "G", "start_pos": [0,0,0], "end_label": "Z", "end_pos": [0,0,0.5]},
            {"start_label": "Z", "start_pos": [0,0,0.5], "end_label": "F", "end_pos": [0.5,0.5,0.0]},
            {"start_label": "F", "start_pos": [0.5,0.5,0.0], "end_label": "G", "end_pos": [0,0,0]},
            {"start_label": "G", "start_pos": [0,0,0], "end_label": "L", "end_pos": [0.5,0.0,0.0]}
        ],
        spinors=True
    )

    # ------------------------------------------------------------------
    # 4. Configure ABACUS Parameters
    # ------------------------------------------------------------------
    # Setup parameters for the NSCF calculation (Interface run)
    job.set_abacus_parameters(
        ecutwfc=100,
        nbands=100,
        nspin=4,        # SOC calculation
        lspinorb=1,     # Turn on SOC
        scf_thr=1e-8,
        scf_nmax=200
    )

    # ------------------------------------------------------------------
    # 5. Run Workflow
    # ------------------------------------------------------------------
    try:
        job.run()
        print("\nSUCCESS: Workflow completed.")
        print(f"Results are located in: {job.work_dir}")
        print(f"Tight-binding model: {job.work_dir}/wannier90_hr.dat")
    except Exception as e:
        print(f"\nERROR: Workflow failed with exception:\n{e}")

if __name__ == "__main__":
    main()
