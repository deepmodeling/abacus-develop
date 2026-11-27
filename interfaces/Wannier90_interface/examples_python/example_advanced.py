#!/usr/bin/env python3
"""
Advanced Example: Step-by-Step Control & Error Handling
=========================================================
Demonstrates:
  - Manual execution of individual steps (step0 → step1 → ... → step4)
  - Customizing advanced parameters (iter nums, mixings, guiding centres)
  - Post-processing checks
  - Dry run mode for input validation

Directory layout:
  Bi2Se3_advanced/
  ├── scf/        ← Step 0: ABACUS SCF output
  └── wannier/    ← Step 1~4: Wannier90 workflow
"""

import os
import sys
import time

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from abacusw90.interface import ABACUSWannier90

# ============================================================
# Config
# ============================================================
DRY_RUN = False  # True: stop after generating inputs; False: run everything
BASE_DIR = "./Bi2Se3_advanced"
SCF_DIR = f"{BASE_DIR}/scf"
WORK_DIR = f"{BASE_DIR}/wannier"


def main():
    print("=== ABACUS Wannier90 Example: Advanced Usage ===")
    print(f"    SCF dir  : {SCF_DIR}")
    print(f"    Work dir : {WORK_DIR}")
    print(f"    Dry run  : {DRY_RUN}")
    print()

    # ----------------------------------------------------------
    # 1. Initialization
    # ----------------------------------------------------------
    job = ABACUSWannier90(
        work_dir=WORK_DIR,
        scf_dir=SCF_DIR,
    )

    # ----------------------------------------------------------
    # 2. Structure: Bi2Se3 (Rhombohedral, R-3m)
    # ----------------------------------------------------------
    lattice = [
        [-2.069, -3.583614, 0.0],
        [2.069, -3.583614, 0.0],
        [0.000, 2.389075, 9.546667],
    ]
    atoms = [
        {"name": "Bi", "pos": [0.399, 0.399, 0.697]},
        {"name": "Bi", "pos": [0.601, 0.601, 0.303]},
        {"name": "Se", "pos": [0.000, 0.000, 0.500]},
        {"name": "Se", "pos": [0.206, 0.206, 0.118]},
        {"name": "Se", "pos": [0.794, 0.794, 0.882]},
    ]
    job.set_structure(lattice, atoms)

    # ----------------------------------------------------------
    # 3. Dependency files
    # ----------------------------------------------------------
    job.pp_orbitals = {"Bi": "../../../tests/PP_ORB/Bi_pbe_fr.upf", "Se": "../../../tests/PP_ORB/Se_pbe_fr.upf"}
    job.orbital_files = ["../../../tests/PP_ORB/Bi_gga_10au_100Ry_2s2p2d.orb", "../../../tests/PP_ORB/Se_gga_10au_100Ry_2s2p2d.orb"]

    # ----------------------------------------------------------
    # 4. Advanced Parameter Configuration
    # ----------------------------------------------------------
    print("[Setup] Configuring advanced parameters...")
    job.set_wannier_parameters(
        num_wann=30,
        num_bands=100,
        projections=["Bi : pz; px; py", "Se : pz; px; py"],
        dis_win_min=3.0,
        dis_win_max=18.0,
        dis_froz_min=3.0,
        dis_froz_max=14.8,
        mp_grid=[4, 4, 4],
        # --- Advanced Wannier90 params ---
        dis_num_iter=500,  # more disentanglement iterations
        num_iter=200,  # more minimization iterations
        guiding_centres=True,  # use guiding centres for initial guess
    )
    job.set_abacus_parameters(
        ecutwfc=100,
        nbands=100,
        basis_type="lcao",
        ks_solver="genelpa",
        nspin=4,
        lspinorb=1,
        # --- Advanced ABACUS params ---
        scf_thr=1e-9,  # tighter convergence
        mixing_beta=0.4,  # lower mixing for stability
        mixing_type="pulay",  # Pulay mixing
    )

    # ----------------------------------------------------------
    # 5. Step-by-step Execution
    # ----------------------------------------------------------
    try:
        job._validate_inputs()

        # === Step 0: SCF ===
        print("\n>>> Step 0: Running ABACUS SCF...")
        t0 = time.time()
        job.step0_run_scf(
            scf_mp_grid=[6, 6, 6],  # denser k-mesh for SCF
            scf_nmax=200,
            scf_thr=1e-9,
        )
        print(f"    Step 0 elapsed: {time.time() - t0:.1f}s")

        # === Step 1: wannier90 -pp ===
        print("\n>>> Step 1: Generating wannier90.win & preprocessing...")
        job.step1_generate_wannier_win()

        # === Dry run gate: stop here in dry mode ===
        if DRY_RUN:
            # Still generate NSCF inputs so user can inspect them
            print("\n>>> Step 2: Preparing ABACUS NSCF inputs...")
            job.step2_prepare_abacus_input()

            print()
            print("=" * 60)
            print("DRY RUN COMPLETE")
            print(f"  SCF inputs  : {job.scf_dir}")
            print(f"  NSCF inputs : {job.work_dir}")
            print()
            print("Generated files to inspect:")
            print(f"  {job.scf_dir}/INPUT       (SCF parameters)")
            print(f"  {job.scf_dir}/KPT         (SCF k-points)")
            print(f"  {job.scf_dir}/STRU        (crystal structure)")
            print(f"  {job.work_dir}/wannier90.win  (Wannier90 input)")
            print(f"  {job.work_dir}/wannier90.nnkp (k-point mapping)")
            print(f"  {job.work_dir}/INPUT       (NSCF parameters)")
            print(f"  {job.work_dir}/KPT         (NSCF k-points)")
            print(f"  {job.work_dir}/STRU        (crystal structure)")
            print()
            print("To continue: set DRY_RUN = False, or run manually:")
            print("  job.step3_run_abacus()")
            print("  job.step4_run_wannier90()")
            print("=" * 60)
            return

        # === Step 2: NSCF input preparation ===
        print("\n>>> Step 2: Preparing ABACUS NSCF inputs...")
        job.step2_prepare_abacus_input()

        # === Step 3: ABACUS NSCF ===
        print("\n>>> Step 3: Running ABACUS NSCF (overlap matrices)...")
        t3 = time.time()
        job.step3_run_abacus()
        print(f"    Step 3 elapsed: {time.time() - t3:.1f}s")

        # === Step 4: Wannier90 minimization ===
        print("\n>>> Step 4: Running Wannier90 minimization...")
        t4 = time.time()
        job.step4_run_wannier90()
        print(f"    Step 4 elapsed: {time.time() - t4:.1f}s")

        # === Post-processing check ===
        print("\n>>> Post-processing check:")
        hr = job.work_dir / "wannier90_hr.dat"
        wout = job.work_dir / "wannier90.wout"
        if hr.exists():
            size_kb = hr.stat().st_size / 1024
            print(f"  [OK] {hr}  ({size_kb:.1f} KB)")
        else:
            print(f"  [MISSING] {hr}")
        if wout.exists():
            # Check spread
            with open(wout) as f:
                for line in f:
                    if "Final" in line and "Spread" in line:
                        print(f"  [OK] {line.strip()}")
                        break
        print()

    except FileNotFoundError as e:
        print(f"\n[FILE ERROR] {e}")
    except RuntimeError as e:
        print(f"\n[RUNTIME ERROR] {e}")
    except ValueError as e:
        print(f"\n[INPUT ERROR] {e}")


if __name__ == "__main__":
    main()
