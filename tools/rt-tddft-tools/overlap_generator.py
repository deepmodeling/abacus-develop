from pyatb.easy_use import stru_analyzer
from pyabacus import ModuleNAO as nao
from pyabacus import ModuleBase as base
import numpy as np
from collections import Counter
import pyatb
from pyatb.easy_use import stru_analyzer
import re
Bohr=0.5291772105638411
class overlap_R:
    def __init__(self, orb_file_dir, atoms0, atomst):
        # get orbital information from the structure
        orb_file_list = list(atoms0.info['basis'].values())
        orb_file_list = [orb_file_dir + orbfile for orbfile in orb_file_list]
        orb_file_num = len(orb_file_list)

        symbols = atoms0.get_chemical_symbols()
        element_counts = Counter(symbols)
        unique_elements = set(symbols)
        assert orb_file_num == len(unique_elements), f"The number of orbitals{orb_file_num}does not match the number of unique elements{len(unique_elements)}!"

        # Build orbital collection
        self.orb = nao.RadialCollection()
        self.orb.build(orb_file_num, orb_file_list, 'o')

        # Standardize the orbital grid
        rmax = self.orb.rcut_max * 2.0
        dr = 0.01
        nr = int(rmax/dr) + 1
        self.orb.set_uniform_grid(True, nr, rmax, 'i', True)

        # Print basic orbital information
        ntype = self.orb.ntype
        lmax = self.orb.lmax

        # Initialize the integrator
        self.S_intor = nao.TwoCenterIntegrator()
        self.S_intor.tabulate(self.orb, self.orb, 'S', nr, rmax)

        # Retrieve orbital index sorting from the structure
        # Determine the correspondence between the basis set indices and orbital information
        self.lattice_vector = atoms0.get_cell()[:] / Bohr 
        self.atom_positions_c = atoms0.get_positions() / Bohr 
        self.atom_positions_ct = atomst.get_positions() / Bohr 
        self.cal_R_direct_coor()

        self.iw2it = dict()
        self.iw2positions_c = dict()
        self.iw2positions_ct = dict()
        self.iw2iL = dict()
        self.iw2iN = dict()
        self.iw2im = dict()
        count = 0
        count_atom = 0
        for it, element in enumerate(unique_elements):
            for ia in range(element_counts[element]):
                for iL in range(self.orb.lmax_(it)+1):
                    for iN in range(self.orb.nzeta(it, iL)):
                        for im in range(2*iL+1):
                            self.iw2it[count] = it
                            self.iw2positions_c[count] = self.atom_positions_c[count_atom]
                            self.iw2positions_ct[count] = self.atom_positions_ct[count_atom]
                            self.iw2iL[count] = iL
                            self.iw2iN[count] = iN
                            if im%2 == 0:
                                self.iw2im[count] = ((im+1) // 2) * -1
                            else:
                                self.iw2im[count] = ((im+1) // 2) * 1
                            
                            count = count + 1
                count_atom = count_atom + 1

    def cal_R_direct_coor(self):
        rcut=self.orb.rcut_max*np.ones(len(self.atom_positions_c),dtype=float)
        print(rcut)
        Ncell = NeighbourCell(self.lattice_vector, self.atom_positions_c, self.atom_positions_ct, rcut)
        self.R_direct_coor = np.array(Ncell.check_interaction_neighbours())
    
    def cal_overlap(self):
        basis_num = len(self.iw2it)
        self.SR = np.zeros((self.R_direct_coor.shape[0], basis_num, basis_num), dtype=float)
        for iR, R_coor in enumerate(self.R_direct_coor):
            for row in range(basis_num):
                for col in range(basis_num):
                    result = self.__cal_S(row, col, R_coor)
                    self.SR[iR, row, col] = result[0]
        return self.SR
    
    def __cal_S(self, iw1, iw2, dR):
        dR = dR @ self.lattice_vector
        dtau = self.iw2positions_c[iw2] - self.iw2positions_ct[iw1]
        R = dR + dtau

        result = self.S_intor.calculate(
            self.iw2it[iw1], 
            self.iw2iL[iw1],
            self.iw2iN[iw1],
            self.iw2im[iw1],
            self.iw2it[iw2], 
            self.iw2iL[iw2],
            self.iw2iN[iw2],
            self.iw2im[iw2],
            R, True)
        
        return result
    def cal_Sk(self,kpoint_direct_coor):
        basis_num = len(self.iw2it)
        Sk = np.zeros([basis_num, basis_num], dtype=complex)
        for iR in range(self.R_direct_coor.shape[0]):
            arg = np.inner(kpoint_direct_coor, self.R_direct_coor[iR]) * 2 * np.pi
            phase = complex(np.cos(arg), np.sin(arg))
            Sk = Sk + phase * self.SR[iR]
        return Sk
def read_kpoints(file='.'):
    file=file+'/kpoints'
    with open(file,'r') as f:
        read_flag=0
        for line in f:
            if('nkstot now' in line):
                nks_tot=int(line.split()[3])
                kpoints=np.zeros([nks_tot,3],dtype=float)
            if(read_flag==1):
                tmp=line.split()
                kpoints[ik]=np.array([float(tmp[1]),float(tmp[2]),float(tmp[3])])
                ik+=1
                if(ik==nks_tot):
                    break
            if('WEIGHT' in line):
                read_flag=1
                ik=0
    return kpoints

class NeighbourCell:
    def __init__(self, lattice_vectors, atom_positions0, atom_positionst, cutoff_radii):
        """
        Initialize lattice vectors, atom positions, and cutoff radii for each atom.
        - lattice_vectors: Lattice vectors
        - atom_positions: Atom positions
        - cutoff_radii: Cutoff radii for each atom
        """
        if len(lattice_vectors) != 3 or len(atom_positions0) == 0 or len(cutoff_radii) != len(atom_positions0):
            raise ValueError("Invalid lattice vectors, atom positions or cutoff radii")
        
        self.lattice_vectors = np.array(lattice_vectors)
        self.atom_positions0 = np.array(atom_positions0)
        self.atom_positionst = np.array(atom_positionst)
        self.cutoff_radii = np.array(cutoff_radii)
        self.natoms = len(atom_positions0)

    def calculate_distance_to_neighbour(self, atom_index1, atom_index2, translation_vector):
        """
        Calculate the distance between any two atoms
        """
        if atom_index1 < 0 or atom_index1 >= self.natoms or \
           atom_index2 < 0 or atom_index2 >= self.natoms:
            raise IndexError("Atom index out of range")
        
        displacement = self.atom_positionst[atom_index2] - self.atom_positions0[atom_index1] + (self.lattice_vectors@np.array(translation_vector).T).T
        
         # Return the Euclidean distance
        return np.linalg.norm(displacement)

    def check_interaction_neighbours(self):
        """
        Check for interactions between atoms in the unit cell and atoms in neighbouring unit cells.
        atom_index: Index of the atom to check
        translation_vector: The translation vector for the neighbouring unit cell
        Returns a list of translation vectors for all neighbouring cells with interactions.
        """
        interactions = []  # List to store translation vectors of neighbouring cells with interactions
        # Loop over translation vectors
        for t in range(-5, 5):  
            for u in range(-5, 5):
                for v in range(-5, 5):
                    # Create a translation vector
                    neighbour_translation = np.array([t, u, v]) 
                    # Loop through each atom in the unit cell and check interactions with atoms in neighbouring cells
                    if_find = False
                    for i in range(self.natoms):
                        for j in range(self.natoms):
                            distance = self.calculate_distance_to_neighbour(i, j, neighbour_translation)
                            # if overlap exists, add the translation vector to the list
                            if distance <= (self.cutoff_radii[i] + self.cutoff_radii[j]):
                                interactions.append(neighbour_translation)
                                if_find = True
                                break 
                        if if_find:
                            break
        return interactions

def save_upper_triangle_as_matrix(matrix, filename):
    # Open the file for writing
    with open(filename, 'w') as f:
        f.write(f'{matrix.shape[0]} ')
        # Iterate over each row
        for i in range(matrix.shape[0]):
            # Extract the upper triangular part of the row
            row = matrix[i, i:] 
            row_str = ' '.join(f'({x.real:.8e},{x.imag:.8e})' for x in row)
            f.write(row_str + '\n')

def overlap_gen(stepref: int, klist, steps, stru_dir='./OUT.ABACUS/STRU', orb_file_dir='',kpoits_dir='./OUT.ABACUS'):
    """Generate overlap files for moved atoms
    :params stepref: the ground state step
    :params klist: the list for needed kpoints
    :params steps: the list for needed steps
    :params stru_dir: the directory for STRU files
    :params orb_file_dir: the directory for orb files, same as 'orbital_dir' in the INPUT
    :params kpoits_dir: the directory for kpoints files
    """
    #read in ref STRU
    with open(stru_dir+'/STRU_MD_'+str(stepref), 'r') as fd:
        atoms0 = stru_analyzer.read_abacus_stru(fd, verbose=True)
    #read in kpoints
    kpt=read_kpoints(kpoits_dir)
    #loop over steps
    for n, nstep in enumerate(steps):
        #read in time t STRU
        with open(stru_dir+'/STRU_MD_'+str(nstep), 'r') as fd:
            atomst = stru_analyzer.read_abacus_stru(fd, verbose=True)
        SR_tmp = overlap_R(orb_file_dir, atoms0, atomst)
        SR_tmp.cal_overlap()
        #loop over kpoints needed
        for ik in klist:
            Sk_ik = SR_tmp.cal_Sk(kpt[ik])
            save_upper_triangle_as_matrix(Sk_ik, f'{nstep}_data-{ik}-S')


#warning : pyatb and pyabacus are required, install them first
if __name__ == "__main__":
    #the kpoints you need, check kpoints file to get the index
    klist=[0]
    #the steps you need, check STRU_MD file to get the index
    start_step = 0
    end_step = 20
    out_interval = 5
    steps = range(start_step, end_step, out_interval)
    #the ground state step
    stepref=0
    #generate overlap files for moved atoms
    overlap_gen(stepref, klist, steps)
