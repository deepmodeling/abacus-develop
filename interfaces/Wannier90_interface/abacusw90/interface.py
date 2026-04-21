"""Core Interface Class for ABACUS Wannier90 workflow."""

import os
import shutil
from pathlib import Path
from typing import List, Dict, Optional
from dataclasses import dataclass, field

from . import io_utils
from . import runner


@dataclass
class ABACUSWannier90:
    """Main interface class to handle ABACUS-Wannier90 workflow."""

    # Directory settings
    work_dir: str = "./wannier_work"
    scf_dir: str = "./scf_out"

    # Executables
    abacus_exe: str = "abacus"
    wannier90_exe: str = "wannier90.x"

    # Structure and files
    structure: Dict = field(default_factory=dict)
    pp_orbitals: Dict = field(default_factory=dict)
    orbital_files: List[str] = field(default_factory=list)

    # Internal objects
    _wannier_input: Optional[io_utils.Wannier90Input] = field(default=None, init=False)
    _abacus_input: Optional[io_utils.AbacusInput] = field(default=None, init=False)
    _scf_input: Optional[io_utils.AbacusInput] = field(
        default=None, init=False
    )

    def __post_init__(self):
        self.work_dir = Path(self.work_dir)
        self.scf_dir = Path(self.scf_dir)

    # ----------------------------------------------------------------
    # Public API: Configuration
    # ----------------------------------------------------------------
    def set_structure(self, lattice: List[List[float]], atoms: List[Dict]):
        self.structure = {"lattice": lattice, "atoms": atoms}

    def set_wannier_parameters(
        self,
        num_wann,
        num_bands,
        projections,
        dis_win_min,
        dis_win_max,
        dis_froz_min,
        dis_froz_max,
        mp_grid,
        kpath=None,
        spinors=True,
        write_hr=True,
        **kwargs,
    ):
        self._wannier_input = io_utils.Wannier90Input(
            num_wann=num_wann,
            num_bands=num_bands,
            projections=projections,
            dis_win_min=dis_win_min,
            dis_win_max=dis_win_max,
            dis_froz_min=dis_froz_min,
            dis_froz_max=dis_froz_max,
            mp_grid=mp_grid,
            spinors=spinors,
            write_hr=write_hr,
            kpath=kpath,
            **kwargs,
        )

    def set_abacus_parameters(
        self,
        ecutwfc,
        nbands,
        nspin=4,
        lspinorb=1,
        noncolin=0,
        scf_thr=1e-8,
        scf_nmax=200,
        **kwargs,
    ):
        if self._wannier_input and nbands != self._wannier_input.params["num_bands"]:
            print(
                f"Warning: ABACUS nbands ({nbands}) overriden by Wannier90 "
                f"num_bands ({self._wannier_input.params['num_bands']})"
            )
            nbands = self._wannier_input.params["num_bands"]
        defaults = {
            "calculation": "nscf",
            "towannier90": 1,
            "wannier_method": 2,
            "nnkpfile": "wannier90.nnkp",
            "symmetry": -1,
            "init_chg": "file",
            "scf_nmax": scf_nmax,
            "scf_thr": scf_thr,
            "ecutwfc": ecutwfc,
            "nbands": nbands,
            "nspin": nspin,
            "lspinorb": lspinorb,
            "noncolin": noncolin,
        }
        defaults.update(kwargs)
        self._abacus_input = io_utils.AbacusInput(**defaults)

    # ----------------------------------------------------------------
    # ★ Step 0: SCF
    # ----------------------------------------------------------------
    def step0_run_scf(
        self,
        scf_mp_grid: List[int] = None,
        scf_nmax: int = 100,
        scf_thr: float = 1e-8,
        smearing_method: str = "gaussian",
        smearing_sigma: float = 0.02,
        mixing_type: str = "broyden",
        mixing_beta: float = 0.7,
    ):
        """
        Step 0: Run ABACUS SCF calculation.

        Creates scf_dir, generates INPUT/KPT/STRU, copies dependency files,
        then executes ABACUS SCF. The resulting charge densities are used
        by subsequent Wannier90 interface steps.

        Args:
            scf_mp_grid:  k-point mesh for SCF. Defaults to Wannier90 mp_grid.
            scf_nmax:     max SCF iterations.
            scf_thr:      energy convergence threshold (Ry).
            smearing_method: smearing type for metals (gaussian/mp/gauss).
            smearing_sigma:  smearing width (Ry).
            mixing_type:  charge mixing method.
            mixing_beta:  mixing parameter (0~1).
        """
        if not self._abacus_input:
            raise ValueError(
                "ABACUS parameters not set. Call set_abacus_parameters() first."
            )
        if not self.structure:
            raise ValueError("Structure not set. Call set_structure() first.")

        print(">>> Step 0: Running ABACUS SCF calculation...")
        self.scf_dir.mkdir(parents=True, exist_ok=True)

        # --- 1. Build SCF INPUT params (reuse shared params, remove wannier90-specific) ---
        skip_keys = {
            "calculation",
            "towannier90",
            "wannier_method",
            "nnkpfile",
            "init_chg",
            "symmetry",
            "out_wannier_mmn",
            "out_wannier_amn",
            "out_wannier_eig",
            "out_wannier_unk",
        }
        scf_params = {
            k: v for k, v in self._abacus_input.params.items() if k not in skip_keys
        }
        scf_params.update(
            {
                "calculation": "scf",
                "scf_nmax": scf_nmax,
                "scf_thr": scf_thr,
                "smearing_method": smearing_method,
                "smearing_sigma": smearing_sigma,
                "mixing_type": mixing_type,
                "mixing_beta": mixing_beta,
                "out_chg": 1,
            }
        )
        self._scf_input = io_utils.AbacusInput(**scf_params)

        # --- 2. Write INPUT ---
        input_file = self.scf_dir / "INPUT"
        self._scf_input.write(input_file)
        print(f"    Written: {input_file}")

        # --- 3. Write KPT (Gamma-centered MP grid) ---
        mp = scf_mp_grid or self._wannier_input.params.get("mp_grid", [4, 4, 4])
        kpt_file = self.scf_dir / "KPT"
        self._write_scf_kpt(kpt_file, mp)
        print(f"    Written: {kpt_file} (Gamma-centered {mp[0]}x{mp[1]}x{mp[2]})")

        # --- 4. Write STRU ---
        stru_file = self.scf_dir / "STRU"
        self._write_stru(stru_file)
        print(f"    Written: {stru_file}")

        # --- 5. Copy dependency files ---
        self._copy_dep_files(self.scf_dir)

        # --- 6. Run ABACUS SCF ---
        cmd = self.abacus_exe
        log_file = self.scf_dir / "abacus_scf.log"
        runner.run_command(
            cmd,
            cwd=str(self.scf_dir),
            log_file=str(log_file),
            label="abacus_scf",
        )

        # --- 7. Verify output ---
        converged = False
        if log_file.exists():
            with open(log_file, "r") as f:
                for line in f:
                    if (
                        "convergence has been achieved" in line.lower()
                        or "converged!" in line.lower()
                        or "scf is converged" in line.lower()
                    ):
                        converged = True
                        break

        if converged:
            print("    SCF converged successfully.")
        else:
            print(f"    Warning: SCF may not have converged. Check {log_file}")

        # Check charge density files
        chg_found = False
        for d in [self.scf_dir, self.scf_dir / "OUT.ABACUS"]:
            if not d.exists():
                continue
            for pat in ["SPIN*_CHG.cube", "SPIN*_CHG", "*CHARGE-DENSITY.restart"]:
                if list(d.glob(pat)):
                    chg_found = True
                    break
            if chg_found:
                break

        if chg_found:
            print("    Charge density files generated.")
        else:
            print(f"    Warning: No charge density files found in {self.scf_dir}")
            print("    → ABACUS may need more iterations or different parameters.")

    def _write_scf_kpt(self, filepath, mp_grid: List[int]):
        """Write ABACUS KPT file for SCF (Gamma-centered Monkhorst-Pack mesh)."""
        with open(filepath, "w") as f:
            f.write("K_POINTS\n")
            f.write("0\n")
            f.write("Gamma\n")
            f.write(f"{mp_grid[0]} {mp_grid[1]} {mp_grid[2]}\n")
            f.write("0 0 0\n")

    def _copy_dep_files(self, target_dir: Path):
        """Copy pseudopotential and orbital files to target directory."""
        for orb in self.orbital_files:
            src = Path(orb)
            if src.exists():
                shutil.copy2(src, target_dir / src.name)
                print(f"    Copied orbital: {src.name}")
            else:
                print(f"    Warning: Orbital file not found: {orb}")
        for pp_name, pp_file in self.pp_orbitals.items():
            src = Path(pp_file)
            if src.exists():
                shutil.copy2(src, target_dir / src.name)
                print(f"    Copied PP: {src.name}")
            else:
                print(f"    Warning: PP file not found: {pp_file}")

    # ----------------------------------------------------------------
    # Step 1~4 (existing, with minor fix in _prepare_scf_files)
    # ----------------------------------------------------------------
    def _prepare_scf_files(self):
        """Copy charge density from scf_dir to work_dir/OUT.ABACUS/."""
        if not self.scf_dir.exists():
            print(f"    Warning: SCF directory not found: {self.scf_dir}")
            return

        target_dir = self.work_dir / "OUT.ABACUS"
        target_dir.mkdir(parents=True, exist_ok=True)

        found = set()
        patterns = [
            "SPIN*_CHG.cube",
            "SPIN*_CHG",
            "chg*.cube",
            "CHG*.cube",
            "CHG*",
            "*CHARGE-DENSITY.restart",
        ]
        for pat in patterns:
            for p in self.scf_dir.glob(pat):
                if p.is_file():
                    found.add(p)
        out_dir = self.scf_dir / "OUT.ABACUS"
        if out_dir.exists():
            for pat in patterns:
                for p in out_dir.glob(pat):
                    if p.is_file():
                        found.add(p)

        if not found:
            print(f"    Warning: No charge density files found in {self.scf_dir}")
            return

        copied = []
        for src in sorted(found):
            dst = target_dir / src.name
            if not dst.exists():
                shutil.copy2(src, dst)
                copied.append(src.name)

        if copied:
            print(f"    Copied {len(copied)} file(s) to {target_dir}/:")
            for name in copied:
                print(f"      {name}")
        else:
            print(f"    Charge density files already present in {target_dir}/")

    def step1_generate_wannier_win(self):
        """Generate wannier90.win and run -pp to generate nnkp."""
        print(">>> Step 1: Generating wannier90.win and preprocessing...")
        self.work_dir.mkdir(parents=True, exist_ok=True)
        win_file = self.work_dir / "wannier90.win"
        self._wannier_input.write(win_file, self.structure)
        print(f"    Written: {win_file}")

        cmd = f"{self.wannier90_exe} -pp wannier90"
        log_file = self.work_dir / "wannier90_pp.log"
        runner.run_command(
            cmd,
            cwd=str(self.work_dir),
            log_file=str(log_file),
            label="wannier90_pp",
            check=True,
        )

        nnkp_file = self.work_dir / "wannier90.nnkp"
        runner.check_file_exists(
            nnkp_file,
            hint="wannier90 -pp did not produce .nnkp. Check wannier90_pp.log.",
        )
        print("    wannier90.nnkp generated successfully.")
        return nnkp_file

    def step2_prepare_abacus_input(self):
        """Prepare ABACUS input files for the interface run."""
        print(">>> Step 2: Preparing ABACUS NSCF input for Wannier90...")
        self.work_dir.mkdir(parents=True, exist_ok=True)
        self._prepare_scf_files()

        nnkp_path = self.work_dir / "wannier90.nnkp"
        runner.check_file_exists(nnkp_path, hint="Run step1 first.")

        kpoints = io_utils.parse_nnkp(nnkp_path)

        kpt_file = self.work_dir / "KPT"
        with open(kpt_file, "w") as f:
            f.write("K_POINTS\n")
            f.write(f"{len(kpoints)}\n")
            f.write("Direct\n")
            for k in kpoints:
                f.write(f"{k[0]:.8f} {k[1]:.8f} {k[2]:.8f} 1.0\n")
        print(f"    Written: {kpt_file} ({len(kpoints)} k-points)")

        input_file = self.work_dir / "INPUT"
        self._abacus_input.write(input_file)
        print(f"    Written: {input_file}")

        stru_file = self.work_dir / "STRU"
        self._write_stru(stru_file)
        print(f"    Written: {stru_file}")

        self._copy_dep_files(self.work_dir)

    def step3_run_abacus(self):
        """Run ABACUS to generate mmn, amn, eig."""
        print(">>> Step 3: Running ABACUS to generate overlap matrices...")
        cmd = self.abacus_exe
        log_file = self.work_dir / "abacus_nscf.log"
        runner.run_command(
            cmd, cwd=str(self.work_dir), log_file=str(log_file), label="abacus_nscf"
        )

        out_dir = self.work_dir / "OUT.ABACUS"
        search_dirs = [self.work_dir, out_dir]

        for fname in ["wannier90.mmn", "wannier90.amn", "wannier90.eig"]:
            found = False
            for d in search_dirs:
                src = d / fname
                if src.exists():
                    if d != self.work_dir:
                        dst = self.work_dir / fname
                        shutil.copy2(src, dst)
                        print(f"    Copied {fname} from OUT.ABACUS/ to work_dir/")
                    found = True
                    break
            if not found:
                runner.check_file_exists(
                    self.work_dir / fname,
                    hint=f"ABACUS did not produce {fname}. Check abacus_nscf.log.",
                )
        print("    ABACUS calculation finished. Matrix elements generated.")

    def step4_run_wannier90(self):
        """Run Wannier90 minimization."""
        print(">>> Step 4: Running Wannier90 minimization...")
        cmd = f"{self.wannier90_exe} wannier90"
        log_file = self.work_dir / "wannier90_min.log"
        runner.run_command(
            cmd, cwd=str(self.work_dir), log_file=str(log_file), label="wannier90_min"
        )

        hr_file = self.work_dir / "wannier90_hr.dat"
        runner.check_file_exists(
            hr_file,
            hint="Wannier90 did not produce wannier90_hr.dat. Check wannier90.wout.",
        )
        print("    Wannier90 finished. Output: wannier90_hr.dat")

    # ----------------------------------------------------------------
    # run() — now includes step0
    # ----------------------------------------------------------------
    def _validate_inputs(self):
        errors = []
        if not self.structure:
            errors.append("Structure not set. Call set_structure() first.")
        if self._wannier_input is None:
            errors.append(
                "Wannier90 parameters not set. Call set_wannier_parameters() first."
            )
        if self._abacus_input is None:
            errors.append(
                "ABACUS parameters not set. Call set_abacus_parameters() first."
            )
        if not self.pp_orbitals:
            errors.append("Pseudopotentials not set. Assign job.pp_orbitals.")
        basis = (
            self._abacus_input.params.get("basis_type", "")
            if self._abacus_input
            else ""
        )
        if basis == "lcao" and not self.orbital_files:
            errors.append(
                "LCAO basis requires orbital files. Assign job.orbital_files."
            )
        if errors:
            raise ValueError("Input validation failed:\n  - " + "\n  - ".join(errors))
        print("    Input validation passed.")

    def run(self, run_scf: bool = True):
        """
        Execute the full workflow.

        Args:
            run_scf: If True, run Step 0 (SCF) when charge densities are missing.
                     If scf_dir already has charge density files, Step 0 is skipped.
                     Set to False to skip Step 0 entirely.
        """
        self._validate_inputs()

        # --- Step 0: SCF (auto-skip if charge densities already exist) ---
        if run_scf:
            has_chg = False
            for d in [self.scf_dir, self.scf_dir / "OUT.ABACUS"]:
                if not d.exists():
                    continue
                for pat in ["SPIN*_CHG.cube", "SPIN*_CHG", "*CHARGE-DENSITY.restart"]:
                    if list(d.glob(pat)):
                        has_chg = True
                        break
                if has_chg:
                    break
            if has_chg:
                print(">>> Step 0: SKIPPED (charge densities already exist in scf_dir)")
            else:
                self.step0_run_scf()

        # --- Step 1~4: Wannier90 interface ---
        self.work_dir.mkdir(parents=True, exist_ok=True)
        self.step1_generate_wannier_win()
        self.step2_prepare_abacus_input()
        self.step3_run_abacus()
        self.step4_run_wannier90()

        print("=" * 60)
        print("All steps completed successfully.")
        print(f"  SCF output : {self.scf_dir}")
        print(f"  Wannier90  : {self.work_dir}")
        print("=" * 60)

    def _write_stru(self, filename):
        """Write ABACUS STRU file.

        Strictly follows official block order to avoid parser state machine errors:
          ATOMIC_SPECIES -> NUMERICAL_ORBITAL -> LATTICE_CONSTANT -> LATTICE_VECTORS -> ATOMIC_POSITIONS
        """
        nspin = self._abacus_input.params.get("nspin", 1)

        with open(filename, "w") as f:

            # --- 1. ATOMIC_SPECIES ---
            f.write("ATOMIC_SPECIES\n")
            seen = set()
            for atom in self.structure["atoms"]:
                name = atom["name"]
                if name not in seen:
                    seen.add(name)
                    pp = self.pp_orbitals.get(name, f"{name}.upf")
                    mass = atom.get("mass", 1.0)
                    f.write(f"{name} {mass:.3f} {os.path.basename(pp)}\n")
            f.write("\n")

            # --- 2. NUMERICAL_ORBITAL ---
            if self.orbital_files:
                f.write("NUMERICAL_ORBITAL\n")
                for orb in self.orbital_files:
                    f.write(f"{os.path.basename(orb)}\n")
                f.write("\n")

            # --- 3. LATTICE_CONSTANT ---
            f.write("LATTICE_CONSTANT\n")
            f.write("1.8897162\n")
            f.write("\n")

            # --- 4. LATTICE_VECTORS ---
            f.write("LATTICE_VECTORS\n")
            for vec in self.structure["lattice"]:
                f.write(f"{vec[0]:12.6f} {vec[1]:12.6f} {vec[2]:12.6f}\n")
            f.write("\n")

            # --- 5. ATOMIC_POSITIONS ---
            f.write("ATOMIC_POSITIONS\n")
            f.write("Direct\n")

            atom_types = {}
            for atom in self.structure["atoms"]:
                name = atom["name"]
                if name not in atom_types:
                    atom_types[name] = []
                atom_types[name].append(atom)

            for i, (name, atom_list) in enumerate(atom_types.items()):
                if i > 0:
                    f.write("\n")
                f.write(f"{name}\n")
                f.write(f"{atom_list[0].get('mag', 0.0)}\n")
                f.write(f"{len(atom_list)}\n")

                for atom in atom_list:
                    pos = atom["pos"]
                    if nspin == 4:
                        m = atom.get("mag", 0.0)
                        mf = "0" if m == 0 else f"{m}"
                        f.write(
                            f"{pos[0]:12.6f} {pos[1]:12.6f} {pos[2]:12.6f} "
                            f"{mf} {mf} {mf}\n"
                        )
                    elif nspin == 2:
                        m = atom.get("mag", 0.0)
                        mf = "0" if m == 0 else f"{m}"
                        f.write(
                            f"{pos[0]:12.6f} {pos[1]:12.6f} {pos[2]:12.6f} " f"{mf}\n"
                        )
                    else:
                        f.write(f"{pos[0]:12.6f} {pos[1]:12.6f} {pos[2]:12.6f}\n")
