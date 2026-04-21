"""Input/Output utilities for ABACUS and Wannier90 files."""

from typing import List, Dict


# ==================================================================
# Monkhorst-Pack k-point grid generator
# ==================================================================
def generate_mp_grid(mp_grid: List[int]) -> List[List[float]]:
    """
    Generate Monkhorst-Pack k-point grid.

    Args:
        mp_grid: [nk1, nk2, nk3] grid dimensions

    Returns:
        List of [kx, ky, kz, weight] for each k-point

    Example:
        mp_grid=[4, 4, 4] → 64 k-points, each weight = 1/64 = 0.015625
    """
    nk1, nk2, nk3 = mp_grid
    ntotal = nk1 * nk2 * nk3
    weight = 1.0 / ntotal

    kpoints = []
    for i in range(nk1):
        for j in range(nk2):
            for k in range(nk3):
                kx = i / nk1
                ky = j / nk2
                kz = k / nk3
                kpoints.append([kx, ky, kz, weight])

    return kpoints


# ==================================================================
# Wannier90 .win Generator
# ==================================================================
class Wannier90Input:
    """Helper class to construct wannier90.win file."""

    def __init__(self, **kwargs):
        self.params = kwargs

    def write(self, filename, structure: Dict):
        with open(filename, "w") as f:

            # ===== 1. Core Parameters =====
            f.write(f"num_wann = {self.params['num_wann']}\n")
            f.write(f"num_bands = {self.params['num_bands']}\n")
            f.write("\n")

            # ===== 2. Disentanglement Parameters =====
            f.write(f"dis_num_iter = {self.params.get('dis_num_iter', 200)}\n")
            f.write("\n")
            f.write("! outer window\n")
            f.write(f"dis_win_min = {self.params['dis_win_min']}\n")
            f.write(f"dis_win_max = {self.params['dis_win_max']}\n")
            f.write("\n")
            f.write("! inner window\n")
            f.write(f"dis_froz_min = {self.params['dis_froz_min']}\n")
            f.write(f"dis_froz_max = {self.params['dis_froz_max']}\n")
            f.write("\n")

            # ===== 3. Global Control =====
            f.write("write_hr = .true.\n")
            f.write(
                f"spinors = {'.true.' if self.params.get('spinors', True) else '.false.'}\n"
            )
            f.write("\n")

            # ===== 4. Projection =====
            f.write("begin projections\n")
            for proj in self.params["projections"]:
                f.write(f"{proj}\n")
            f.write("end projections\n")
            f.write("\n")

            # ===== 5. Unit Cell =====
            f.write("begin unit_cell_cart\n")
            for vec in structure["lattice"]:
                f.write(f"{vec[0]:16.10f} {vec[1]:16.10f} {vec[2]:16.10f}\n")
            f.write("end unit_cell_cart\n")
            f.write("\n")

            # ===== 6. Atomic Coordinate =====
            f.write("begin atoms_frac\n")
            for atom in structure["atoms"]:
                pos = atom["pos"]
                f.write(
                    f"{atom['name']}  {pos[0]:12.6f} {pos[1]:12.6f} {pos[2]:12.6f}\n"
                )
            f.write("end atoms_frac\n")
            f.write("\n")

            # ===== 7. K-point path =====
            if self.params.get("kpath"):
                f.write("bands_plot = true\n")
                f.write(
                    f"bands_num_points {self.params.get('bands_num_points', 101)}\n"
                )
                f.write("begin kpoint_path\n")
                for kp in self.params["kpath"]:
                    sp = kp["start_pos"]
                    ep = kp["end_pos"]
                    f.write(
                        f"{kp['start_label']} " f"{sp[0]:.8f} {sp[1]:.8f} {sp[2]:.8f} "
                    )
                    f.write(
                        f"{kp['end_label']} " f"{ep[0]:.8f} {ep[1]:.8f} {ep[2]:.8f}\n"
                    )
                f.write("end kpoint_path\n")
                f.write("\n")

            # ===== 8. K-point mesh =====
            mp = self.params["mp_grid"]
            kpoints = generate_mp_grid(mp)

            f.write(f"mp_grid : {mp[0]} {mp[1]} {mp[2]}\n")
            f.write("begin kpoints\n")
            for kp in kpoints:
                f.write(f"{kp[0]:.8f} {kp[1]:.8f} {kp[2]:.8f} {kp[3]:.8f}\n")
            f.write("end kpoints\n")


# ==================================================================
# ABACUS INPUT Generator
# ==================================================================
class AbacusInput:
    """Helper class to construct ABACUS INPUT file."""

    def __init__(self, **kwargs):
        self.params = kwargs

    def write(self, filename):
        with open(filename, "w") as f:
            f.write("INPUT_PARAMETERS\n")
            for key, value in self.params.items():
                f.write(f"{key}    {value}\n")


# ==================================================================
# .nnkp Parser
# ==================================================================
def parse_nnkp(filename) -> List[List[float]]:
    """Parse wannier90.nnkp to extract k-points."""
    kpoints = []
    with open(filename, "r") as f:
        lines = f.readlines()
    in_kpoints = False
    for line in lines:
        line_lower = line.lower().strip()
        if "begin kpoints" in line_lower:
            in_kpoints = True
            continue
        if "end kpoints" in line_lower:
            break
        if in_kpoints:
            parts = line.split()
            if len(parts) >= 3:
                try:
                    kpoints.append([float(parts[0]), float(parts[1]), float(parts[2])])
                except ValueError:
                    continue
    if not kpoints:
        raise ValueError(
            f"No k-points found in {filename}. "
            "Check that the file is a valid wannier90.nnkp."
        )
    return kpoints
