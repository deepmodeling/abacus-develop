#!/usr/bin/env python3
"""
Example: PW (Plane Wave) Basis
================================
Key differences from basic.py:
  - basis_type = 'pw'
  - ks_solver  = 'cg'  (or 'dav', 'cg', 'bpcg')
  - No orbital files needed

Directory layout:
  Bi2Se3_pw/
  ├── scf/        ← Step 0: ABACUS SCF output
  └── wannier/    ← Step 1~4: Wannier90 workflow
"""

import os
import sys

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from abacusw90.interface import ABACUSWannier90

# ============================================================
# Config
# ============================================================
DRY_RUN = False
BASE_DIR = "./Bi2Se3_pw"
SCF_DIR = f"{BASE_DIR}/scf"
WORK_DIR = f"{BASE_DIR}/wannier"


def main():
    print("=== ABACUS Wannier90 Example: PW Basis ===")
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
    # 3. Dependency files (PW only needs pseudopotentials)
    # ----------------------------------------------------------
    job.pp_orbitals = {"Bi": "../../../tests/PP_ORB/Bi_pbe_fr.upf", "Se": "../../../tests/PP_ORB/Se_pbe_fr.upf"}
    # ← PW基组不需要轨道文件，不设置 orbital_files

    # ----------------------------------------------------------
    # 4. Wannier90 Parameters
    # ----------------------------------------------------------
    job.set_wannier_parameters(
        num_wann=30,
        num_bands=100,
        projections=["Bi : pz; px; py", "Se : pz; px; py"],
        dis_win_min=3.0,
        dis_win_max=18.0,
        dis_froz_min=3.0,
        dis_froz_max=13.5,
        mp_grid=[4, 4, 4],
    )

    # ----------------------------------------------------------
    # 5. ABACUS Parameters (PW basis)
    # ----------------------------------------------------------
    job.set_abacus_parameters(
        ecutwfc=100,
        nbands=100,
        basis_type="pw",  # ← Plane wave
        ks_solver="cg",  # ← Conjugate gradient
        nspin=4,
        lspinorb=1,
    )

    # ----------------------------------------------------------
    # 6. Run
    # ----------------------------------------------------------
    try:
        if DRY_RUN:
            job._validate_inputs()

            # Step 0: SCF (PW基组)
            job.step0_run_scf(
                scf_mp_grid=[4, 4, 4],
            )

            # Step 1: wannier90 -pp
            job.step1_generate_wannier_win()

            # Step 2: NSCF input files
            job.step2_prepare_abacus_input()

            print()
            print("=" * 60)
            print("DRY RUN COMPLETE — Input files generated:")
            print(f"  SCF   : {job.scf_dir}")
            print(f"  Work  : {job.work_dir}")
            print()
            print("To run the full workflow, set DRY_RUN = False")
            print("=" * 60)
        else:
            job.run(run_scf=True)

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
                with open(wout) as f:
                    for line in f:
                        if "Final" in line and "Spread" in line:
                            print(f"  [OK] {line.strip()}")
                            break

    except FileNotFoundError as e:
        print(f"\n[FILE ERROR] {e}")
        print("  [1] wannier90.x in PATH?  →  which wannier90.x")
        print("  [2] abacus in PATH?       →  which abacus")
        print("  [3] PP files?             →  ls ../../../tests/PP_ORB/Bi*.upf")
    except RuntimeError as e:
        print(f"\n[RUNTIME ERROR] {e}")
    except ValueError as e:
        print(f"\n[INPUT ERROR] {e}")


if __name__ == "__main__":
    main()

