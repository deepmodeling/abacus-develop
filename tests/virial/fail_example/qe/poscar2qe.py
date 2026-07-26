import numpy as np
def convert_cartesian_to_fractional(cartesian_coords, lattice_vectors):
    lattice_matrix = np.array(lattice_vectors).T
    fractional_coords = np.linalg.solve(lattice_matrix, cartesian_coords.T).T
    return fractional_coords


def _reciprocal_box(box) :
    rbox = np.linalg.inv(box)
    rbox = rbox.T
    return rbox

def make_kspacing_kpoints(box, kspacing) :
    if type(kspacing) is not list:
        kspacing = [kspacing, kspacing, kspacing]
    box = np.array(box)
    rbox = _reciprocal_box(box)
    kpoints = [max(1,(np.ceil(2 * np.pi * np.linalg.norm(ii) / ks).astype(int))) for ii,ks in zip(rbox,kspacing)]
    return kpoints

def convert_poscar_to_qe_scfinput(poscar_file, pseudopotentials, atomic_masses):
    with open(poscar_file, 'r') as file:
        poscar_lines = file.readlines()
    
    # Extract lattice vectors
    lattice_vectors = []
    for i in range(2, 5):
        lattice_vectors.append([float(x) for x in poscar_lines[i].split()])
        
    
    # Extract atomic species and positions
    species = poscar_lines[5].split()
    num_atoms = [int(x) for x in poscar_lines[6].split()]
    
    total_atoms = sum(num_atoms)
    # Check the coordinate type and extract atomic positions
    coord_type = poscar_lines[7].strip()
    atomic_positions = []
    for i in range(8, 8 + total_atoms):
        atomic_positions.append([float(x) for x in poscar_lines[i].split()[0:3]])

    if coord_type.lower().startswith('cartesian'):
        atomic_positions = convert_cartesian_to_fractional(np.array(atomic_positions), lattice_vectors)
    
    
    # Mapping for atomic species and their respective positions
    species_positions = []
    index = 0
    for i, s in enumerate(species):
        for _ in range(num_atoms[i]):
            species_positions.append((s, atomic_positions[index]))
            index += 1
    
    # Generate ATOMIC_SPECIES content
    atomic_species_content = "ATOMIC_SPECIES\n"
    for s in species:
        if s in pseudopotentials and s in atomic_masses:
            atomic_species_content += f" {s}  {atomic_masses[s]:.6f}  {pseudopotentials[s]}\n"
        else:
            raise ValueError(f"Pseudopotential or atomic mass for {s} not found.")
    
    # Generate QE input file content
    qe_input_content = f"""
&control
    calculation = 'vc-relax'
    restart_mode = 'from_scratch',
    prefix = 'openlam'
    outdir = './tmp'
    pseudo_dir = '../pp_all'
    etot_conv_thr = 1e-5
    forc_conv_thr = 1e-3
    nstep = 100
/
&SYSTEM
 ibrav=0,   
 nosym=.true.,                     
 nat={total_atoms},                         
 ntyp={len(species)},                        
 occupations = 'smearing',       
 smearing = 'gauss',                
 degauss = 0.01,                   
 ecutwfc = 100,                   
 lspinorb = .false.,                
 noncolin = .false.,                
/
&ELECTRONS
 conv_thr = 1.0d-9,                  
 mixing_beta = 0.2,               
 electron_maxstep = 200,          
/
&electrons
    mixing_beta = 0.2
    conv_thr = 1.0d-8
/
&ions
/
&cell
    press = 0
/
{atomic_species_content}CELL_PARAMETERS {{angstrom}}        
 {lattice_vectors[0][0]:.15f}  {lattice_vectors[0][1]:.15f}  {lattice_vectors[0][2]:.15f}
 {lattice_vectors[1][0]:.15f}  {lattice_vectors[1][1]:.15f}  {lattice_vectors[1][2]:.15f}
 {lattice_vectors[2][0]:.15f}  {lattice_vectors[2][1]:.15f}  {lattice_vectors[2][2]:.15f}
ATOMIC_POSITIONS (crystal)       
"""
    
    for s, pos in species_positions:
        qe_input_content += " {0}  {1:.15f}  {2:.15f}  {3:.15f}\n".format(s, pos[0], pos[1], pos[2])
    
    kpt=make_kspacing_kpoints(lattice_vectors,0.15)
    qe_input_content += "K_POINTS {automatic} \n "            
    qe_input_content += f"{kpt[0]} {kpt[1]} {kpt[2]}  0 0 0\n " 
    
    return qe_input_content

# Define the pseudopotentials dictionary
pseudopotentials = {
'H':'H_ONCV_PBE-1.0.upf',
'He':'He_ONCV_PBE_FR-1.0.upf',
'Li':'Li_pd_04_s-high.UPF',
'Be':'Be_pd_04_s-high.UPF',
'B':'B.PD03.PBE.UPF',
'C':'C.PD04.PBE.UPF',
'N':'N_ONCV_PBE-1.0.upf',
'O':'O.PD04.PBE.UPF',
'F':'F.PD03.PBE.UPF',
'Ne':'Ne_ONCV_PBE-1.0.upf',
'Na':'Na-sp.PD04.PBE.UPF',
'Mg':'Mg.PD04.PBE.UPF',
'Al':'Al.PD04.PBE.UPF',
'Si':'Si.PD04.PBE.UPF',
'P':'P-sp.PD04.PBE.UPF',
'S':'S.PD03.PBE.UPF',
'Cl':'Cl.PD04.PBE.UPF',
'Ar':'Ar_ONCV_PBE-1.0.upf',
'K':'K-sp.PD04.PBE.UPF',
'Ca':'Ca_pd_04_sp.UPF',
'Sc':'Sc-sp.PD04.PBE.UPF',
'Ti':'Ti_pd_04_sp.UPF',
'V':'V-sp.PD04.PBE.UPF',
'Cr':'Cr.Rappe.PBE.UPF',
'Mn':'Mn-sp.PD04.PBE.UPF',
'Fe':'Fe.Rappe.PBE.UPF',
'Co':'Co.Rappe.PBE.UPF',
'Ni':'Ni-sp.PD04.PBE.UPF',
'Cu':'Cu_sg15_1.2_.UPF',
'Zn':'Zn_pd_04_.UPF',
'Ga':'Ga_pd_04_d.UPF',
'Ge':'Ge_pd_04_d.UPF',
'As':'As-d.PD04.PBE.UPF',
'Se':'Se.upf',
'Br':'Br.PD03.PBE.UPF',
'Kr':'Kr_ONCV_PBE-1.2.upf',
'Rb':'Rb-sp.PD04.PBE.UPF',
'Sr':'Sr-sp.PD04.PBE.UPF',
'Y':'Y-sp.PD04.PBE.UPF',
'Zr':'Zr-sp.PD04.PBE.UPF',
'Nb':'Nb-sp.PD04.PBE.UPF',
'Mo':'Mo-sp.PD04.PBE.UPF',
'Tc':'Tc_ONCV_PBE-1.0.upf',
'Ru':'Ru-sp.PD04.PBE.UPF',
'Rh':'Rh-sp.PD04.PBE.UPF',
'Pd':'Pd-sp.PD04.PBE.UPF',
'Ag':'Ag.Rappe.PBE.UPF',
'Cd':'Cd.Rappe.PBE.UPF',
'In':'In_pd_04_d.UPF',
'Sn':'Sn-d.PD04.PBE.UPF',
'Sb':'Sb.PD03.PBE.UPF',
'Te':'Te.PD04.PBE.UPF',
'I':'I.upf',
'Xe':'Xe.upf',
'Cs':'Cs.upf',
'Ba':'Ba.upf',
'La':'La_pd_04_sp.UPF',
'Ce':'Ce_pd_04_3+_f--core.UPF',
'Pr':'Pr3+_f--core.PD04.PBE.UPF',
'Nd':'Nd3+_f--core.PD04.PBE.UPF',
'Pm':'Pm-sp.PD04.PBE.UPF',
'Sm':'Sm3+_f--core.PD04.PBE.UPF',
'Eu':'Eu-sp.PD04.PBE.UPF',
'Gd':'Gd_pd_04_3+_f--core.UPF',
'Tb':'Tb_pd_04_3+_f--core.UPF',
'Dy':'Dy_pd_04_3+_f--core-icmod1.UPF',
'Ho':'Ho3+_f--core.PD04.PBE.UPF',
'Er':'Er3+_f--core.PD04.PBE.UPF',
'Tm':'Tm_pd_04_3+_f--core-icmod1.UPF',
'Yb':'Yb-sp.PD04.PBE.UPF',
'Lu':'Lu3+_f--core.PD04.PBE.UPF',
'Hf':'Hf_pd_04_sp.UPF',
'Ta':'Ta_pd_04_sp.UPF',
'W':'W.Rappe.PBE.UPF',
'Re':'Re-sp.PD04.PBE.UPF',
'Os':'Os.Rappe.PBE.UPF',
'Ir':'Ir.Rappe.PBE.UPF',
'Pt':'Pt_pd_04_sp.UPF',
'Au':'Au.Rappe.PBE.UPF',
'Hg':'Hg_ONCV_PBE-1.2.upf',
'Tl':'Tl_ONCV_PBE-1.0.upf',
'Pb':'Pb_pd_03_.UPF',
'Bi':'Bi-spd-high.PD04.PBE.UPF',
'Ac':'Ac-5spd.upf',
'Th':'Th-5spdf.upf',
'Pa':'Pa-5spdf.upf',
'U':'U-5spdf.upf',
'Np':'Np-5spdf.upf',
'Pu':'Pu-5spdf.upf'
}

# Define the atomic masses dictionary
atomic_masses = {
'H': 1.00794, 'He': 4.002602, 'Li': 6.94, 'Be': 9.0122, 'B': 10.81,
'C': 12.011, 'N': 14.007, 'O': 15.999, 'F': 18.998, 'Ne': 20.180,
'Na': 22.990, 'Mg': 24.305, 'Al': 26.982, 'Si': 28.085, 'P': 30.974,
'S': 32.06, 'Cl': 35.45, 'Ar': 39.948, 'K': 39.098, 'Ca': 40.078,
'Sc': 44.956, 'Ti': 47.867, 'V': 50.942, 'Cr': 51.996, 'Mn': 54.938,
'Fe': 55.845, 'Co': 58.933, 'Ni': 58.693, 'Cu': 63.546, 'Zn': 65.38,
'Ga': 69.723, 'Ge': 72.630, 'As': 74.922, 'Se': 78.971, 'Br': 79.904,
'Kr': 83.798, 'Rb': 85.468, 'Sr': 87.62, 'Y': 88.906, 'Zr': 91.224,
'Nb': 92.906, 'Mo': 95.95, 'Tc': 98, 'Ru': 101.07, 'Rh': 102.91,
'Pd': 106.42, 'Ag': 107.87, 'Cd': 112.41, 'In': 114.82, 'Sn': 118.71,
'Sb': 121.76, 'Te': 127.60, 'I': 126.90, 'Xe': 131.29, 'Cs': 132.91,
'Ba': 137.33, 'La': 138.91, 'Ce': 140.12, 'Pr': 140.91, 'Nd': 144.24,
'Pm': 145, 'Sm': 150.36, 'Eu': 151.96, 'Gd': 157.25, 'Tb': 158.93,
'Dy': 162.50, 'Ho': 164.93, 'Er': 167.26, 'Tm': 168.93, 'Yb': 173.05,
'Lu': 174.97, 'Hf': 178.49, 'Ta': 180.95, 'W': 183.84, 'Re': 186.21,
'Os': 190.23, 'Ir': 192.22, 'Pt': 195.08, 'Au': 196.97, 'Hg': 200.59,
'Tl': 204.38, 'Pb': 207.2, 'Bi': 208.98, 'Ac': 227, 'Th': 232.04, 'Pa': 231.04,
'U': 238.03, 'Np': 237, 'Pu': 244
}


# Example usage
poscar_file = 'POSCAR'
qe_input_content = convert_poscar_to_qe_scfinput(poscar_file, pseudopotentials, atomic_masses)
with open('scf.in', 'w') as qe_file:
    qe_file.write(qe_input_content)
    

print("Quantum Espresso input file 'scf.in' generated successfully.")

