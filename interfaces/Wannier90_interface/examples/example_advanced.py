#!/usr/bin/env python3
"""
Advanced Example: Customized Workflow and Step-by-Step Execution

This script demonstrates:
1. Executing steps manually (Step 1 -> Step 2 -> Step 3 -> Step 4).
2. Customizing advanced input parameters.
3. Handling specific file dependencies.
"""

import os
import sys
import time

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from abacusw90 import ABACUSWannier90
from pathlib import Path

def check_input_files():
    """Helper function to ensure necessary input files exist."""
    required_files = ["Bi.oncvpsp.upf", "Se.oncvpsp.upf"] # Add orbitals if needed
    for f in required_files:
        if not Path(f).exists():
            print(f"Warning: Input file '{f}' not found in current directory.")
            print("Please ensure pseudo-potential and orbital files are available.")
            return False
    return True

def main():
    print("="*60)
    print(" ABACUS Wannier90 Advanced Example")
    print("="*60)

    if not check_input_files():
        return

    # ------------------------------------------------------------------
    # 1. Initialization
    # ------------------------------------------------------------------
    job = ABACUSWannier90(
        work_dir="./Bi2Se3_wannier_advanced",
        scf_dir="./Bi2Se3_scf"
    )

    # Structure Setup
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
    job.pp_orbitals = {"Bi": "Bi.oncvpsp.upf", "Se": "Se.oncvpsp.upf"}

    # ------------------------------------------------------------------
    # 2. Advanced Parameter Configuration
    # ------------------------------------------------------------------
    print("\n[Setup] Configuring advanced parameters...")
    
    # Customizing Wannier90 parameters
    job.set_wannier_parameters(
        num_wann=30,
        num_bands=100,
        projections=["Bi : pz; px; py", "Se : pz; px; py"],
        dis_win_min=3.0,
        dis_win_max=18.0,
        dis_froz_min=3.0,
        dis_froz_max=14.8,
        mp_grid=[4, 4, 4],
        # Advanced kwargs passed directly to .win
        dis_num_iter=500,         # Increase disentanglement iterations
        dis_mix_ratio=0.5,        # Mixing ratio
        num_iter=200,             # Wannierisation iterations
        guiding_centres=True,     # Use guiding centres
        write_xyz=True            # Output WF centres in xyz format
    )

    # Customizing ABACUS parameters
    job.set_abacus_parameters(
        ecutwfc=100,
        nbands=100,
        nspin=4,
        lspinorb=1,
        # Advanced kwargs passed directly to INPUT
        scf_thr=1e-9,
        mixing_beta=0.4,
        mixing_type='pulay'
    )

    # ------------------------------------------------------------------
    # 3. Manual Step-by-Step Execution
    # ------------------------------------------------------------------
    
    # Step 1: Generate .win and .nnkp
    try:
        print("\n>>> Step 1: Generating Wannier90 inputs...")
        job.step1_generate_wannier_win()
        print("    Success: wannier90.nnkp created.")
    except Exception as e:
        print(f"    Failed: {e}")
        return

    # Step 2: Prepare ABACUS inputs
    # Here we can inject logic, e.g., checking if nnkp parsing is correct
    try:
        print("\n>>> Step 2: Preparing ABACUS inputs...")
        job.step2_prepare_abacus_input()
        print("    Success: INPUT, KPT, STRU created.")
    except Exception as e:
        print(f"    Failed: {e}")
        return

    # Step 3: Run ABACUS
    # This is the most time-consuming step
    try:
        print("\n>>> Step 3: Running ABACUS NSCF calculation...")
        start_time = time.time()
        job.step3_run_abacus()
        end_time = time.time()
        print(f"    Success: ABACUS finished in {end_time - start_time:.2f} seconds.")
    except Exception as e:
        print(f"    Failed: {e}")
        # Check log file for details
        log_file = job.work_dir / "abacus_nscf.log"
        if log_file.exists():
            print(f"    Check log file: {log_file}")
        return

    # Step 4: Run Wannier90
    try:
        print("\n>>> Step 4: Running Wannier90 minimization...")
        job.step4_run_wannier90()
        print("    Success: Wannier90 finished.")
    except Exception as e:
        print(f"    Failed: {e}")
        return

    # ------------------------------------------------------------------
    # 4. Post-Processing Check
    # ------------------------------------------------------------------
    print("\n[Post-Processing] Checking output files...")
    hr_file = job.work_dir / "wannier90_hr.dat"
    wout_file = job.work_dir / "wannier90.wout"

    if hr_file.exists():
        print(f"  [+] Tight-binding model found: {hr_file}")
        print(f"  [+] File size: {hr_file.stat().st_size / 1024:.2f} KB")
    else:
        print("  [-] Error: wannier90_hr.dat not found.")

    if wout_file.exists():
        # Basic check for convergence
        with open(wout_file, 'r') as f:
            content = f.read()
            if "All done: wannier90 exiting" in content:
                print("  [+] Wannier90 exited normally.")
            else:
                print("  [!] Warning: Wannier90 may not have exited normally.")

if __name__ == "__main__":
    main()
