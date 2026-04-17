"""
Input/Output utilities for ABACUS and Wannier90 files.
"""

from typing import List, Dict

class Wannier90Input:
    """Helper class to construct wannier90.win file."""
    def __init__(self, **kwargs):
        self.params = kwargs

    def write(self, filename, structure: Dict):
        with open(filename, 'w') as f:
            f.write(f"num_wann = {self.params['num_wann']}\n")
            f.write(f"num_bands = {self.params['num_bands']}\n")
            if self.params.get('dis_num_iter'):
                f.write(f"dis_num_iter = {self.params['dis_num_iter']}\n")
            else:
                f.write("dis_num_iter = 200\n")
            
            f.write("! outer window\n")
            f.write(f"dis_win_min = {self.params['dis_win_min']}\n")
            f.write(f"dis_win_max = {self.params['dis_win_max']}\n")
            f.write("! inner window\n")
            f.write(f"dis_froz_min = {self.params['dis_froz_min']}\n")
            f.write(f"dis_froz_max = {self.params['dis_froz_max']}\n")
            
            f.write("write_hr = .true.\n")
            f.write(f"spinors = {'.true.' if self.params.get('spinors', True) else '.false.'}\n")
            
            f.write("begin projections\n")
            for proj in self.params['projections']:
                f.write(f"{proj}\n")
            f.write("end projections\n")
            
            f.write("begin unit_cell_cart\n")
            f.write("Ang\n")
            for vec in structure['lattice']:
                f.write(f"{vec[0]:15.10f} {vec[1]:15.10f} {vec[2]:15.10f}\n")
            f.write("end unit_cell_cart\n")
            
            f.write("begin atoms_frac\n")
            for atom in structure['atoms']:
                pos = atom['pos']
                f.write(f"{atom['name']:5s} {pos[0]:15.10f} {pos[1]:15.10f} {pos[2]:15.10f}\n")
            f.write("end atoms_frac\n")
            
            if self.params.get('kpath'):
                f.write("bands_plot = true\n")
                f.write(f"bands_num_points {self.params.get('bands_num_points', 101)}\n")
                f.write("begin kpoint_path\n")
                for kp in self.params['kpath']:
                    f.write(f"{kp['start_label']} {kp['start_pos'][0]} {kp['start_pos'][1]} {kp['start_pos'][2]} ")
                    f.write(f"{kp['end_label']} {kp['end_pos'][0]} {kp['end_pos'][1]} {kp['end_pos'][2]}\n")
                f.write("end kpoint_path\n")

            mp = self.params['mp_grid']
            f.write(f"mp_grid : {mp[0]} {mp[1]} {mp[2]}\n")

class AbacusInput:
    """Helper class to construct ABACUS INPUT file."""
    def __init__(self, **kwargs):
        self.params = kwargs

    def write(self, filename):
        with open(filename, 'w') as f:
            f.write("INPUT_PARAMETERS\n")
            f.write(f"calculation {self.params['calculation']}\n")
            f.write(f"ecutwfc {self.params['ecutwfc']}\n")
            f.write(f"nbands {self.params['nbands']}\n")
            f.write(f"nspin {self.params['nspin']}\n")
            f.write(f"lspinorb {self.params['lspinorb']}\n")
            f.write(f"noncolin {self.params['noncolin']}\n")
            
            f.write(f"scf_nmax {self.params['scf_nmax']}\n")
            f.write(f"scf_thr {self.params['scf_thr']}\n")
            f.write(f"init_chg {self.params['init_chg']}\n")
            f.write(f"symmetry {self.params['symmetry']}\n")
            
            f.write(f"towannier90 {self.params['towannier90']}\n")
            f.write(f"wannier_method {self.params['wannier_method']}\n")
            f.write(f"nnkpfile {self.params['nnkpfile']}\n")

def parse_nnkp(filename) -> List[List[float]]:
    """
    Parse wannier90.nnkp to extract k-points.
    """
    kpoints = []
    with open(filename, 'r') as f:
        lines = f.readlines()
        
    in_kpoints = False
    for line in lines:
        line_lower = line.lower()
        if "begin kpoints" in line_lower:
            in_kpoints = True
            continue
        if "end kpoints" in line_lower:
            break
        if in_kpoints:
            parts = line.split()
            if len(parts) >= 4:
                kpoints.append([float(parts[0]), float(parts[1]), float(parts[2])])
    return kpoints
