"""
High-level runner interface for ABACUS calculations.

This module provides the `abacus()` function - the main entry point
for running ABACUS DFT calculations from Python.

Two implementations are available:
1. C++ bindings (if available) - direct library calls
2. Subprocess fallback - calls the abacus executable
"""

from dataclasses import dataclass, field
from typing import Optional, Dict, Any, Union, List
from pathlib import Path
import numpy as np
import os
import subprocess
import re


# Unit conversion constants
RY_TO_EV = 13.605693122994  # 1 Ry = 13.6057 eV
BOHR_TO_ANG = 0.529177249   # 1 Bohr = 0.529 Angstrom


@dataclass
class CalculationResult:
    """
    Container for ABACUS calculation results.

    All energies are stored in eV units.

    Attributes
    ----------
    converged : bool
        Whether SCF converged
    niter : int
        Number of SCF iterations
    etot : float
        Total energy in eV
    forces : np.ndarray, optional
        Forces on atoms (nat, 3) in eV/Angstrom
    stress : np.ndarray, optional
        Stress tensor (3, 3) in kbar
    energies : dict
        Dictionary of energy components (all in eV)
    fermi_energy : float
        Fermi energy in eV
    bandgap : float
        Band gap in eV
    nat : int
        Number of atoms
    ntype : int
        Number of atom types
    nbands : int
        Number of bands
    nks : int
        Number of k-points
    """
    # Convergence info
    converged: bool = False
    niter: int = 0
    drho: float = 0.0

    # Energies (all in eV)
    etot: float = 0.0
    eband: float = 0.0
    hartree_energy: float = 0.0
    etxc: float = 0.0
    ewald_energy: float = 0.0
    demet: float = 0.0
    exx: float = 0.0
    evdw: float = 0.0

    # Forces (in eV/Angstrom) and stress (in kbar)
    forces: Optional[np.ndarray] = None
    stress: Optional[np.ndarray] = None

    # Electronic structure info
    fermi_energy: float = 0.0  # in eV
    bandgap: float = 0.0       # in eV

    # System info
    nat: int = 0
    ntype: int = 0
    nbands: int = 0
    nks: int = 0

    # Output file tracking
    output_dir: str = ""  # Path to OUT.$suffix folder
    log_file: str = ""    # Path to the main log file (running_*.log)
    output_files: Dict[str, str] = field(default_factory=dict)  # filename -> full path

    @property
    def etot_ev(self) -> float:
        """Total energy in eV (same as etot, for compatibility)."""
        return self.etot

    @property
    def energies(self) -> Dict[str, float]:
        """Dictionary of all energy components (all in eV)."""
        return {
            'etot': self.etot,
            'eband': self.eband,
            'hartree_energy': self.hartree_energy,
            'etxc': self.etxc,
            'ewald_energy': self.ewald_energy,
            'demet': self.demet,
            'exx': self.exx,
            'evdw': self.evdw,
        }

    @property
    def forces_ev_ang(self) -> Optional[np.ndarray]:
        """Forces in eV/Angstrom (same as forces, for compatibility)."""
        return self.forces

    @property
    def has_forces(self) -> bool:
        """Whether forces are available."""
        return self.forces is not None

    @property
    def has_stress(self) -> bool:
        """Whether stress is available."""
        return self.stress is not None

    @property
    def has_output_dir(self) -> bool:
        """Whether output directory exists and is set."""
        return bool(self.output_dir) and os.path.isdir(self.output_dir)

    def get_output_file(self, filename: str) -> Optional[str]:
        """
        Get full path to a specific output file.

        Parameters
        ----------
        filename : str
            Name of the output file (e.g., 'running_scf.log', 'BANDS_1.dat')

        Returns
        -------
        str or None
            Full path to the file if it exists, None otherwise
        """
        return self.output_files.get(filename)

    def list_output_files(self) -> List[str]:
        """
        List all output file names.

        Returns
        -------
        list of str
            List of output file names
        """
        return list(self.output_files.keys())

    def summary(self) -> str:
        """Return a summary string of the calculation result."""
        lines = [
            "=== ABACUS Calculation Result ===",
            f"Converged: {'Yes' if self.converged else 'No'}",
            f"SCF iterations: {self.niter}",
            f"Final drho: {self.drho:.2e}",
            "",
            "Energies (eV):",
            f"  Total energy: {self.etot:.8f}",
            f"  Band energy:  {self.eband:.8f}",
            f"  Hartree:      {self.hartree_energy:.8f}",
            f"  XC energy:    {self.etxc:.8f}",
            f"  Ewald:        {self.ewald_energy:.8f}",
            f"  Entropy(-TS): {self.demet:.8f}",
            f"  EXX:          {self.exx:.8f}",
            f"  VdW:          {self.evdw:.8f}",
        ]

        lines.extend([
            "",
            "System info:",
            f"  Atoms: {self.nat}, Types: {self.ntype}",
            f"  Bands: {self.nbands}, K-points: {self.nks}",
            f"  Fermi energy: {self.fermi_energy:.6f} eV",
            f"  Band gap: {self.bandgap:.6f} eV",
        ])

        lines.append("")
        lines.append("Forces (eV/Angstrom):")
        if self.has_forces and self.forces is not None:
            max_force = np.max(np.abs(self.forces))
            lines.append(f"  Calculated ({self.nat} atoms), Max force: {max_force:.6f}")
            for i, f in enumerate(self.forces):
                lines.append(f"    Atom {i+1}: [{f[0]:12.8f}, {f[1]:12.8f}, {f[2]:12.8f}]")
        else:
            lines.append("  Not calculated")

        lines.append("")
        lines.append("Stress (kbar):")
        if self.has_stress and self.stress is not None:
            lines.append("  Calculated:")
            for i, row in enumerate(self.stress):
                lines.append(f"    [{row[0]:12.6f}, {row[1]:12.6f}, {row[2]:12.6f}]")
        else:
            lines.append("  Not calculated")

        # Output file tracking
        lines.extend([
            "",
            "Output:",
            f"  Directory: {self.output_dir if self.output_dir else 'N/A'}",
            f"  Log file: {os.path.basename(self.log_file) if self.log_file else 'N/A'}",
            f"  Files: {len(self.output_files)} output files",
        ])

        return "\n".join(lines)

    def __repr__(self) -> str:
        return (
            f"<CalculationResult converged={self.converged} "
            f"etot={self.etot:.6f} eV>"
        )


def _find_abacus_executable() -> Optional[str]:
    """Find the abacus executable in PATH or common locations."""
    import shutil

    # Check PATH first
    abacus_path = shutil.which("abacus")
    if abacus_path:
        return abacus_path

    # Check common locations
    common_paths = [
        "/usr/local/bin/abacus",
        "/usr/bin/abacus",
        os.path.expanduser("~/abacus/build/abacus"),
        os.path.expanduser("~/.local/bin/abacus"),
    ]

    for path in common_paths:
        if os.path.isfile(path) and os.access(path, os.X_OK):
            return path

    return None


def _parse_running_log(log_path: str) -> CalculationResult:
    """Parse the running log file to extract calculation results (all energies in eV)."""
    result = CalculationResult()

    if not os.path.exists(log_path):
        return result

    with open(log_path, 'r') as f:
        content = f.read()

    # Parse convergence - check multiple patterns
    convergence_patterns = [
        r"#SCF IS CONVERGED#",
        r"charge density convergence is achieved",
        r"convergence is achieved",
        r"SCF CONVERGED",
    ]
    for pattern in convergence_patterns:
        if re.search(pattern, content, re.IGNORECASE):
            result.converged = True
            break

    # Parse total energy (look for final energy first)
    # Pattern 1: !FINAL_ETOT_IS -57.02190809937956 eV
    final_match = re.search(r"!FINAL_ETOT_IS\s+([-\d.]+)\s+eV", content)
    if final_match:
        result.etot = float(final_match.group(1))  # Already in eV
    else:
        # Pattern 2: E_KohnSham lines - get the last one (Ry, eV)
        ks_matches = re.findall(r"E_KohnSham\s+([-\d.]+)\s+([-\d.]+)", content)
        if ks_matches:
            result.etot = float(ks_matches[-1][1])  # Second value is in eV

    # Parse number of SCF iterations - count ELEC ITER lines
    iter_matches = re.findall(r"#ELEC ITER#\s+(\d+)", content)
    if iter_matches:
        result.niter = int(iter_matches[-1])

    # Parse drho - get the last value (Electron density deviation)
    drho_matches = re.findall(r"Electron density deviation\s+([\d.eE+-]+)", content)
    if drho_matches:
        result.drho = float(drho_matches[-1])

    # Parse number of atoms
    nat_match = re.search(r"TOTAL ATOM NUMBER\s*=\s*(\d+)", content)
    if nat_match:
        result.nat = int(nat_match.group(1))

    # Parse number of types - count "READING ATOM TYPE" lines
    ntype_matches = re.findall(r"READING ATOM TYPE\s+(\d+)", content)
    if ntype_matches:
        result.ntype = int(ntype_matches[-1])

    # Parse number of bands
    nbands_match = re.search(r"Number of electronic states \(NBANDS\)\s*=\s*(\d+)", content)
    if nbands_match:
        result.nbands = int(nbands_match.group(1))

    # Parse number of k-points (nkstot now = X after reduction)
    nkstot_match = re.search(r"nkstot now\s*=\s*(\d+)", content)
    if nkstot_match:
        result.nks = int(nkstot_match.group(1))
    else:
        # Fallback to original nkstot
        nkstot_match = re.search(r"nkstot\s*=\s*(\d+)", content)
        if nkstot_match:
            result.nks = int(nkstot_match.group(1))

    # Parse Fermi energy - format: E_Fermi  0.4657978215  6.3375044881 (Ry, eV)
    fermi_match = re.search(r"E_Fermi\s+([-\d.]+)\s+([-\d.]+)", content)
    if fermi_match:
        result.fermi_energy = float(fermi_match.group(2))  # Second value is in eV

    # Parse band gap - format: E_gap(k)  0.1070873708  1.4569984261 (Ry, eV)
    gap_match = re.search(r"E_gap\(k\)\s+([-\d.]+)\s+([-\d.]+)", content)
    if gap_match:
        result.bandgap = float(gap_match.group(2))  # Second value is in eV

    # Parse energy components from the final SCF iteration
    # Format: E_xxx  value_Ry  value_eV - we take the eV value (second column)

    # E_band (band energy)
    eband_match = re.search(r"E_band\s+([-\d.]+)\s+([-\d.]+)", content)
    if eband_match:
        result.eband = float(eband_match.group(2))  # eV

    # E_Hartree
    hartree_match = re.search(r"E_Hartree\s+([-\d.]+)\s+([-\d.]+)", content)
    if hartree_match:
        result.hartree_energy = float(hartree_match.group(2))  # eV

    # E_xc (exchange-correlation)
    etxc_match = re.search(r"E_xc\s+([-\d.]+)\s+([-\d.]+)", content)
    if etxc_match:
        result.etxc = float(etxc_match.group(2))  # eV

    # E_Ewald
    ewald_match = re.search(r"E_Ewald\s+([-\d.]+)\s+([-\d.]+)", content)
    if ewald_match:
        result.ewald_energy = float(ewald_match.group(2))  # eV

    # E_entropy(-TS) for metals
    demet_match = re.search(r"E_entropy\(-TS\)\s+([-\d.]+)\s+([-\d.]+)", content)
    if demet_match:
        result.demet = float(demet_match.group(2))  # eV

    # E_exx (exact exchange)
    exx_match = re.search(r"E_exx\s+([-\d.]+)\s+([-\d.]+)", content)
    if exx_match:
        result.exx = float(exx_match.group(2))  # eV

    return result


def _get_suffix_from_input(input_dir: str) -> str:
    """Parse the suffix from INPUT file."""
    input_file = os.path.join(input_dir, "INPUT")
    suffix = "ABACUS"  # default suffix

    if os.path.exists(input_file):
        with open(input_file, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('#') or not line:
                    continue
                # Parse suffix parameter
                if 'suffix' in line.lower():
                    parts = line.split()
                    if len(parts) >= 2:
                        suffix = parts[1]
                        break
    return suffix


def _collect_output_files(output_dir: str) -> Dict[str, str]:
    """
    Collect all output files from the output directory.

    Parameters
    ----------
    output_dir : str
        Path to the output directory (OUT.$suffix)

    Returns
    -------
    dict
        Dictionary mapping filename to full path
    """
    output_files = {}
    if output_dir and os.path.isdir(output_dir):
        try:
            for entry in os.listdir(output_dir):
                full_path = os.path.join(output_dir, entry)
                if os.path.isfile(full_path):
                    output_files[entry] = full_path
        except OSError:
            pass  # Ignore errors during directory iteration
    return output_files


def _parse_forces_from_log(log_path: str, nat: int) -> Optional[np.ndarray]:
    """Parse forces from the running log file (returns forces in eV/Angstrom)."""
    if not os.path.exists(log_path) or nat <= 0:
        return None

    with open(log_path, 'r') as f:
        content = f.read()

    # Try multiple force block formats
    # Format 1: #TOTAL-FORCE (eV/Angstrom)#
    #           -------------------------------------------------------------------------
    #               Atoms              Force_x              Force_y              Force_z
    #           -------------------------------------------------------------------------
    #                 Al1         0.0000000000         0.0000000000         0.0000000000
    force_pattern1 = r"#TOTAL-FORCE \(eV/Angstrom\)#.*?-{10,}\s*\n\s*Atoms\s+Force_x\s+Force_y\s+Force_z\s*\n\s*-{10,}\s*\n((?:\s*\S+\s+[-\d.eE+]+\s+[-\d.eE+]+\s+[-\d.eE+]+\s*\n)+)"
    match = re.search(force_pattern1, content, re.DOTALL)

    if not match:
        # Format 2: TOTAL-FORCE (eV/Angstrom)
        #           ----------------------------
        #            atom    x       y       z
        #             Si1  0.001   0.002   0.003
        force_pattern2 = r"TOTAL-FORCE \(eV/Angstrom\).*?-{10,}\s*\n\s*atom.*?\n((?:\s*\S+\s+[-\d.eE+]+\s+[-\d.eE+]+\s+[-\d.eE+]+\s*\n)+)"
        match = re.search(force_pattern2, content, re.DOTALL)

    if match:
        force_lines = match.group(1).strip().split('\n')
        forces = []
        for line in force_lines:
            parts = line.split()
            if len(parts) >= 4:
                # parts[0] is atom label, parts[1:4] are fx, fy, fz in eV/Angstrom
                try:
                    fx, fy, fz = float(parts[1]), float(parts[2]), float(parts[3])
                    forces.append([fx, fy, fz])  # Already in eV/Angstrom
                except (ValueError, IndexError):
                    continue

        if len(forces) == nat:
            return np.array(forces)

    return None


def _parse_stress_from_log(log_path: str) -> Optional[np.ndarray]:
    """Parse stress tensor from the running log file (returns stress in kbar)."""
    if not os.path.exists(log_path):
        return None

    with open(log_path, 'r') as f:
        content = f.read()

    # Try multiple stress block formats
    # Format 1: #TOTAL-STRESS (kbar)#
    #           ----------------------------------------------------------------
    #                    Stress_x             Stress_y             Stress_z
    #           ----------------------------------------------------------------
    #               15.7976835472         0.0000000000         0.0000000000
    stress_pattern1 = r"#TOTAL-STRESS \(kbar\)#.*?-{10,}\s*\n\s*Stress_x\s+Stress_y\s+Stress_z\s*\n\s*-{10,}\s*\n((?:\s*[-\d.eE+]+\s+[-\d.eE+]+\s+[-\d.eE+]+\s*\n){3})"
    match = re.search(stress_pattern1, content, re.DOTALL | re.IGNORECASE)

    if not match:
        # Format 2: TOTAL-STRESS (KBAR)
        #           ----------------------------
        #             1.234   0.000   0.000
        stress_pattern2 = r"TOTAL-STRESS \(KBAR\).*?-{10,}\s*\n((?:\s*[-\d.eE+]+\s+[-\d.eE+]+\s+[-\d.eE+]+\s*\n){3})"
        match = re.search(stress_pattern2, content, re.DOTALL | re.IGNORECASE)

    if match:
        stress_lines = match.group(1).strip().split('\n')
        stress = []
        for line in stress_lines:
            parts = line.split()
            if len(parts) >= 3:
                try:
                    stress.append([float(parts[0]), float(parts[1]), float(parts[2])])
                except (ValueError, IndexError):
                    continue

        if len(stress) == 3:
            return np.array(stress)

    return None


def _modify_input_file(input_dir: str, calculate_force: bool, calculate_stress: bool) -> Optional[str]:
    """
    Modify INPUT file to add cal_force and cal_stress parameters.

    Returns the path to the backup file if modifications were made, None otherwise.
    """
    input_file = os.path.join(input_dir, "INPUT")
    if not os.path.exists(input_file):
        return None

    # Read original content
    with open(input_file, 'r') as f:
        lines = f.readlines()

    # Check if cal_force/cal_stress already exist
    has_cal_force = False
    has_cal_stress = False
    for line in lines:
        line_lower = line.lower().strip()
        if line_lower.startswith('cal_force'):
            has_cal_force = True
        if line_lower.startswith('cal_stress'):
            has_cal_stress = True

    # If both already exist, no need to modify
    if has_cal_force and has_cal_stress:
        return None

    # Create backup
    backup_file = input_file + ".pyabacus_backup"
    with open(backup_file, 'w') as f:
        f.writelines(lines)

    # Add missing parameters
    new_lines = lines.copy()
    if not has_cal_force:
        new_lines.append(f"cal_force {1 if calculate_force else 0}\n")
    if not has_cal_stress:
        new_lines.append(f"cal_stress {1 if calculate_stress else 0}\n")

    # Write modified file
    with open(input_file, 'w') as f:
        f.writelines(new_lines)

    return backup_file


def _restore_input_file(input_dir: str, backup_file: Optional[str]):
    """Restore INPUT file from backup."""
    if backup_file is None:
        return

    input_file = os.path.join(input_dir, "INPUT")
    if os.path.exists(backup_file):
        # Restore original
        with open(backup_file, 'r') as f:
            content = f.read()
        with open(input_file, 'w') as f:
            f.write(content)
        # Remove backup
        os.remove(backup_file)


def _run_abacus_subprocess(
    input_dir: str,
    output_dir: str,
    verbosity: int,
    calculate_force: bool = True,
    calculate_stress: bool = False,
    nprocs: int = 1,
    nthreads: int = 1,
) -> CalculationResult:
    """Run ABACUS using subprocess and parse results."""
    import shutil

    # Find abacus executable
    abacus_exe = _find_abacus_executable()
    if abacus_exe is None:
        raise RuntimeError(
            "ABACUS executable not found. Please ensure 'abacus' is in your PATH "
            "or install ABACUS from https://github.com/deepmodeling/abacus-develop"
        )

    # Convert to absolute path
    input_dir = os.path.abspath(input_dir)

    # Modify INPUT file to add cal_force and cal_stress
    backup_file = _modify_input_file(input_dir, calculate_force, calculate_stress)

    try:
        # Get the suffix from INPUT file to know where output will be
        suffix = _get_suffix_from_input(input_dir)
        expected_out_dir = os.path.join(input_dir, f"OUT.{suffix}")

        # Set up environment
        env = os.environ.copy()
        env["OMP_NUM_THREADS"] = str(nthreads)

        # Build command
        if nprocs > 1:
            # Find mpirun or mpiexec
            mpirun = shutil.which("mpirun") or shutil.which("mpiexec")
            if mpirun is None:
                raise RuntimeError(
                    f"MPI requested (nprocs={nprocs}) but mpirun/mpiexec not found. "
                    "Please install MPI or set nprocs=1."
                )
            cmd = [mpirun, "-np", str(nprocs), abacus_exe]
        else:
            cmd = [abacus_exe]

        # Set up stdout/stderr based on verbosity
        if verbosity >= 2:
            stdout = None
            stderr = None
        elif verbosity == 1:
            stdout = subprocess.PIPE
            stderr = subprocess.PIPE
        else:
            stdout = subprocess.DEVNULL
            stderr = subprocess.DEVNULL

        try:
            proc = subprocess.run(
                cmd,
                cwd=input_dir,
                env=env,
                stdout=stdout,
                stderr=stderr,
                timeout=None,  # No timeout
            )
        except subprocess.TimeoutExpired:
            raise RuntimeError("ABACUS calculation timed out")
        except Exception as e:
            raise RuntimeError(f"Failed to run ABACUS: {e}")

        # Find and parse the output
        # First try the expected output directory based on suffix
        log_path = None
        if os.path.exists(expected_out_dir):
            for log_name in ["running_scf.log", "running_relax.log", "running_cell-relax.log", "running_nscf.log"]:
                candidate = os.path.join(expected_out_dir, log_name)
                if os.path.exists(candidate):
                    log_path = candidate
                    break

        # Fallback: find the most recently modified OUT.* directory
        if log_path is None:
            out_dirs = [d for d in os.listdir(input_dir) if d.startswith("OUT.") and os.path.isdir(os.path.join(input_dir, d))]
            if out_dirs:
                latest_out = max(out_dirs, key=lambda d: os.path.getmtime(os.path.join(input_dir, d)))
                out_dir_path = os.path.join(input_dir, latest_out)
                for log_name in ["running_scf.log", "running_relax.log", "running_cell-relax.log", "running_nscf.log"]:
                    candidate = os.path.join(out_dir_path, log_name)
                    if os.path.exists(candidate):
                        log_path = candidate
                        break

        if log_path and os.path.exists(log_path):
            result = _parse_running_log(log_path)
            # Set output tracking fields
            result.log_file = os.path.abspath(log_path)
            result.output_dir = os.path.abspath(os.path.dirname(log_path))
            result.output_files = _collect_output_files(result.output_dir)

            # Parse forces if requested
            if calculate_force and result.nat > 0:
                forces = _parse_forces_from_log(log_path, result.nat)
                if forces is not None:
                    result.forces = forces

            # Parse stress if requested
            if calculate_stress:
                stress = _parse_stress_from_log(log_path)
                if stress is not None:
                    result.stress = stress
        else:
            result = CalculationResult()
            # Try to find output directory even if log file wasn't found
            if os.path.exists(expected_out_dir):
                result.output_dir = os.path.abspath(expected_out_dir)
                result.output_files = _collect_output_files(result.output_dir)

    finally:
        # Restore original INPUT file
        _restore_input_file(input_dir, backup_file)

    return result


def abacus(
    input_dir: Optional[str] = None,
    *,
    input_file: Optional[str] = None,
    stru_file: Optional[str] = None,
    kpt_file: Optional[str] = None,
    pseudo_dir: Optional[str] = None,
    orbital_dir: Optional[str] = None,
    output_dir: Optional[str] = None,
    calculate_force: bool = True,
    calculate_stress: bool = False,
    verbosity: int = 1,
    nprocs: int = 1,
    nthreads: int = 1,
) -> CalculationResult:
    """
    Run an ABACUS DFT calculation.

    This is the main entry point for running ABACUS calculations from Python.
    It provides the same functionality as the ABACUS command-line program.

    Parameters
    ----------
    input_dir : str, optional
        Directory containing INPUT, STRU, KPT files.
        If not specified, uses current directory.
    input_file : str, optional
        Explicit path to INPUT file. Overrides input_dir/INPUT.
    stru_file : str, optional
        Explicit path to STRU file. Overrides value in INPUT.
    kpt_file : str, optional
        Explicit path to KPT file. Overrides value in INPUT.
    pseudo_dir : str, optional
        Directory containing pseudopotential files.
        Overrides pseudo_dir in INPUT.
    orbital_dir : str, optional
        Directory containing orbital files (for LCAO).
        Overrides orbital_dir in INPUT.
    output_dir : str, optional
        Directory for output files. Default: "OUT.PYABACUS"
    calculate_force : bool, optional
        Whether to calculate forces. Default: True
    calculate_stress : bool, optional
        Whether to calculate stress tensor. Default: False
    verbosity : int, optional
        Output verbosity level:
        - 0: Silent (no output)
        - 1: Normal (default)
        - 2: Verbose (detailed output)
    nprocs : int, optional
        Number of MPI processes. Default: 1
        Equivalent to: mpirun -np nprocs abacus
    nthreads : int, optional
        Number of OpenMP threads. Default: 1
        Equivalent to: OMP_NUM_THREADS=nthreads

    Returns
    -------
    CalculationResult
        Object containing all calculation results including:
        - converged: Whether SCF converged
        - etot: Total energy (Ry)
        - etot_ev: Total energy (eV)
        - forces: Forces on atoms (if calculate_force=True)
        - stress: Stress tensor (if calculate_stress=True)
        - energies: Dictionary of energy components

    Raises
    ------
    FileNotFoundError
        If INPUT file is not found
    RuntimeError
        If calculation fails or ABACUS is not installed

    Examples
    --------
    Basic SCF calculation:

    >>> result = pyabacus.abacus("./Si_scf/")
    >>> print(f"Energy: {result.etot_ev:.6f} eV")
    >>> print(f"Converged: {result.converged}")

    Calculate forces and stress:

    >>> result = pyabacus.abacus(
    ...     "./Si_relax/",
    ...     calculate_force=True,
    ...     calculate_stress=True,
    ... )
    >>> print(f"Max force: {np.max(np.abs(result.forces_ev_ang)):.4f} eV/Ang")

    Parallel calculation with MPI and OpenMP:

    >>> result = pyabacus.abacus(
    ...     "./Si_scf/",
    ...     nprocs=4,      # 4 MPI processes
    ...     nthreads=2,    # 2 OpenMP threads per process
    ... )

    Silent mode:

    >>> result = pyabacus.abacus("./Si_scf/", verbosity=0)
    """
    # Try to use C++ driver first
    try:
        from ._driver_pack import PyDriver
        _HAS_CPP_DRIVER = True
    except ImportError:
        _HAS_CPP_DRIVER = False

    # Determine input directory
    if input_dir is None:
        if input_file is not None:
            input_dir = str(Path(input_file).parent)
        else:
            input_dir = "."

    # Validate input directory exists
    input_path = Path(input_dir)
    if not input_path.exists():
        raise FileNotFoundError(f"Input directory not found: {input_dir}")

    # Check for INPUT file
    if input_file is None:
        input_file_path = input_path / "INPUT"
        if not input_file_path.exists():
            raise FileNotFoundError(
                f"INPUT file not found in {input_dir}. "
                "Please provide input_file parameter or ensure INPUT exists."
            )

    # Set default output directory
    if output_dir is None:
        output_dir = "OUT.PYABACUS"

    if _HAS_CPP_DRIVER:
        # Use C++ driver
        driver = PyDriver()

        cpp_result = driver.run(
            input_dir=str(input_dir),
            input_file=input_file or "",
            stru_file=stru_file or "",
            kpt_file=kpt_file or "",
            pseudo_dir=pseudo_dir or "",
            orbital_dir=orbital_dir or "",
            output_dir=output_dir or "",
            calculate_force=calculate_force,
            calculate_stress=calculate_stress,
            verbosity=verbosity,
        )

        # Convert C++ result to Python dataclass
        result = CalculationResult(
            converged=cpp_result.converged,
            niter=cpp_result.niter,
            drho=cpp_result.drho,
            etot=cpp_result.etot,
            eband=cpp_result.eband,
            hartree_energy=cpp_result.hartree_energy,
            etxc=cpp_result.etxc,
            ewald_energy=cpp_result.ewald_energy,
            demet=cpp_result.demet,
            exx=cpp_result.exx,
            evdw=cpp_result.evdw,
            fermi_energy=cpp_result.fermi_energy,
            bandgap=cpp_result.bandgap,
            nat=cpp_result.nat,
            ntype=cpp_result.ntype,
            nbands=cpp_result.nbands,
            nks=cpp_result.nks,
        )

        # Copy forces if available
        if cpp_result.has_forces:
            result.forces = np.array(cpp_result.forces)

        # Copy stress if available
        if cpp_result.has_stress:
            result.stress = np.array(cpp_result.stress)

        # Copy output tracking fields
        result.output_dir = cpp_result.output_dir
        result.log_file = cpp_result.log_file
        result.output_files = dict(cpp_result.output_files)
    else:
        # Use subprocess fallback
        result = _run_abacus_subprocess(
            input_dir=str(input_dir),
            output_dir=output_dir,
            verbosity=verbosity,
            calculate_force=calculate_force,
            calculate_stress=calculate_stress,
            nprocs=nprocs,
            nthreads=nthreads,
        )

    return result


def run_scf(
    input_dir: str,
    **kwargs
) -> CalculationResult:
    """
    Convenience function for running SCF calculation.

    This is an alias for `abacus()` with default parameters
    suitable for single-point SCF calculations.

    Parameters
    ----------
    input_dir : str
        Directory containing input files
    **kwargs
        Additional arguments passed to `abacus()`

    Returns
    -------
    CalculationResult
        Calculation results
    """
    return abacus(input_dir, **kwargs)


def run_relax(
    input_dir: str,
    **kwargs
) -> CalculationResult:
    """
    Convenience function for running geometry optimization.

    This is an alias for `abacus()` with force calculation enabled.

    Parameters
    ----------
    input_dir : str
        Directory containing input files
    **kwargs
        Additional arguments passed to `abacus()`

    Returns
    -------
    CalculationResult
        Calculation results
    """
    kwargs.setdefault('calculate_force', True)
    return abacus(input_dir, **kwargs)
