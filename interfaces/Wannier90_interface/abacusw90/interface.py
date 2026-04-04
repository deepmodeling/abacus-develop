"""
Core Interface Class for ABACUS Wannier90 workflow.
"""

import os
import shutil
from pathlib import Path
from typing import List, Dict, Optional
from dataclasses import dataclass, field

from . import io_utils
from . import runner

@dataclass
class ABACUSWannier90:
    """
    Main interface class to handle ABACUS-Wannier90 workflow.
    """
    # Directory settings
    work_dir: str = "./wannier_work"
    scf_dir: str = "./scf_out"  # Directory containing ABACUS SCF results
    
    # Executables
    abacus_exe: str = "abacus"
    wannier90_exe: str = "wannier90.x"
    
    # Structure and files
    structure: Dict = field(default_factory=dict)
    pp_orbitals: Dict = field(default_factory=dict) # {'Bi': 'Bi.upf', ...}
    orbital_files: List[str] = field(default_factory=list) # ['Bi.orb', ...]
    
    # Internal objects
    _wannier_input: Optional[io_utils.Wannier90Input] = field(default=None, init=False)
    _abacus_input: Optional[io_utils.AbacusInput] = field(default=None, init=False)

    def __post_init__(self):
        self.work_dir = Path(self.work_dir)
        self.scf_dir = Path(self.scf_dir)
        self.work_dir.mkdir(parents=True, exist_ok=True)

    def set_structure(self, lattice: List[List[float]], atoms: List[Dict]):
        """
        Set the crystal structure.
        atoms format: [{'name': 'Bi', 'pos': [0.0, 0.0, 0.0]}, ...]
        """
        self.structure = {
            'lattice': lattice,
            'atoms': atoms
        }

    def set_wannier_parameters(
        self, 
        num_wann: int, 
        num_bands: int,
        projections: List[str],
        dis_win_min: float, 
        dis_win_max: float,
        dis_froz_min: float, 
        dis_froz_max: float,
        mp_grid: List[int],
        kpath: Optional[List[Dict]] = None,
        spinors: bool = True,
        write_hr: bool = True,
        **kwargs
    ):
        """
        Set parameters for Wannier90 calculation.
        """
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
            **kwargs
        )

    def set_abacus_parameters(
        self, 
        ecutwfc: float, 
        nbands: int,
        nspin: int = 4, 
        lspinorb: int = 1, 
        noncolin: int = 0,
        scf_thr: float = 1e-8,
        scf_nmax: int = 200,
        **kwargs
    ):
        """
        Set ABACUS input parameters for NSCF run.
        """
        # Sync nbands with Wannier90 input
        if self._wannier_input and nbands != self._wannier_input.params['num_bands']:
            print(f"Warning: ABACUS nbands ({nbands}) overriden by Wannier90 num_bands ({self._wannier_input.params['num_bands']})")
            nbands = self._wannier_input.params['num_bands']

        # Default interface parameters
        defaults = {
            "calculation": "nscf",
            "towannier90": 1,
            "wannier_method": 2,
            "nnkpfile": "wannier90.nnkp",
            "symmetry": -1,
            "init_chg": "file",
            "scf_nmax": scf_nmax,
            "scf_thr": scf_thr,
            # User parameters
            "ecutwfc": ecutwfc,
            "nbands": nbands,
            "nspin": nspin,
            "lspinorb": lspinorb,
            "noncolin": noncolin
        }
        
        # Update with user kwargs (allows override of defaults)
        defaults.update(kwargs)
        
        self._abacus_input = io_utils.AbacusInput(**defaults)

    def _prepare_scf_files(self):
        """Prepare necessary files from SCF calculation."""
        # Copy charge density files (SPIN1~4_CHG.cube)
        for spin in range(4):
            src = self.scf_dir / f"SPIN{spin+1}_CHG.cube"
            if src.exists():
                shutil.copy(src, self.work_dir)
        
    def step1_generate_wannier_win(self):
        """Generate wannier90.win and run -pp to generate nnkp."""
        print(">>> Step 1: Generating wannier90.win and preprocessing...")
        
        # Write .win file
        win_file = self.work_dir / "wannier90.win"
        self._wannier_input.write(win_file, self.structure)
        
        # Run wannier90 -pp
        cmd = f"{self.wannier90_exe} -pp wannier90"
        runner.run_command(cmd, cwd=self.work_dir)
        
        nnkp_file = self.work_dir / "wannier90.nnkp"
        runner.check_file_exists(nnkp_file)
        
        print(">>> wannier90.nnkp generated successfully.")
        return nnkp_file

    def step2_prepare_abacus_input(self):
        """Prepare ABACUS input files for the interface run."""
        print(">>> Step 2: Preparing ABACUS NSCF input for Wannier90...")
        
        # Copy SCF results
        self._prepare_scf_files()
        
        # Parse nnkp to get KPOINTS for ABACUS
        nnkp_path = self.work_dir / "wannier90.nnkp"
        kpoints = io_utils.parse_nnkp(nnkp_path)
        
        # Write KPT file
        kpt_file = self.work_dir / "KPT"
        with open(kpt_file, 'w') as f:
            f.write("K_POINTS\n")
            f.write(f"{len(kpoints)}\n")
            f.write("Direct\n")
            for k in kpoints:
                f.write(f"{k[0]:.8f} {k[1]:.8f} {k[2]:.8f} 1.0\n")
                
        # Write INPUT file
        input_file = self.work_dir / "INPUT"
        self._abacus_input.write(input_file)
        
        # Write STRU file
        stru_file = self.work_dir / "STRU"
        self._write_stru(stru_file)
        
        # Link/Copy Orbitals and PP
        for orb in self.orbital_files:
            if os.path.exists(orb):
                shutil.copy(orb, self.work_dir)
        for pp_name, pp_file in self.pp_orbitals.items():
            if os.path.exists(pp_file):
                shutil.copy(pp_file, self.work_dir)

    def step3_run_abacus(self):
        """Run ABACUS to generate mmn, amn, eig."""
        print(">>> Step 3: Running ABACUS to generate overlap matrices...")
        
        cmd = self.abacus_exe
        log_file = self.work_dir / "abacus_nscf.log"
        runner.run_command(cmd, cwd=self.work_dir, log_file=log_file)
        
        # Check output files
        required_files = ["wannier90.mmn", "wannier90.amn", "wannier90.eig"]
        for f in required_files:
            runner.check_file_exists(self.work_dir / f)
            
        print(">>> ABACUS calculation finished. Matrix elements generated.")

    def step4_run_wannier90(self):
        """Run Wannier90 minimization."""
        print(">>> Step 4: Running Wannier90 minimization...")
        
        cmd = f"{self.wannier90_exe} wannier90"
        log_file = self.work_dir / "wannier90.log"
        runner.run_command(cmd, cwd=self.work_dir, log_file=log_file)
        
        print(">>> Wannier90 finished. Check wannier90.wout and wannier90_hr.dat")

    def run(self):
        """Execute the full workflow."""
        self.step1_generate_wannier_win()
        self.step2_prepare_abacus_input()
        self.step3_run_abacus()
        self.step4_run_wannier90()
        print("All steps completed successfully.")

    def _write_stru(self, filename):
        """Write ABACUS STRU file."""
        with open(filename, 'w') as f:
            f.write("ATOMIC_SPECIES\n")
            for atom_name, pp_file in self.pp_orbitals.items():
                f.write(f"{atom_name} 1.0 {os.path.basename(pp_file)}\n")
            
            f.write("\nLATTICE_VECTORS\n")
            for vec in self.structure['lattice']:
                f.write(f"{vec[0]:.10f} {vec[1]:.10f} {vec[2]:.10f}\n")
            
            f.write("\nATOMIC_POSITIONS\n")
            f.write("Direct\n")
            
            # Group atoms by type
            atom_types = {}
            for atom in self.structure['atoms']:
                name = atom['name']
                if name not in atom_types:
                    atom_types[name] = []
                atom_types[name].append(atom['pos'])
            
            for name, positions in atom_types.items():
                f.write(f"{name}\n")
                f.write("0.0\n") # Magnetic moment dummy
                for pos in positions:
                    f.write(f"{pos[0]:.10f} {pos[1]:.10f} {pos[2]:.10f} 1 1 1\n")
