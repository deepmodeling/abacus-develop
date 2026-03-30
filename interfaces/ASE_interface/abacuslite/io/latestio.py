'''implements the parser for the latest version of ABACUS'''
import re
import shutil
import unittest
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
from ase.atoms import Atoms
from ase.calculators.singlepoint import (
    SinglePointKPoint,
    SinglePointDFTCalculator
)
from ase.units import GPa
from ase.stress import full_3x3_to_voigt_6_stress

# some output formats are not updated,
# for these cases, we import from the legacyio module
from abacuslite.io.legacyio import (
    read_kpoints_from_running_log,
    read_energies_from_running_log,
    read_traj_from_md_dump,
    read_magmom_from_running_log,
    is_invalid_arr,
)

def read_esolver_type_from_running_log(src: str | Path | List[str]) \
    -> str:
    '''
    read the esolver type from the ABACUS running log file.

    Parameters
    ----------
    fn : str
        The path to the ABACUS running log file.
    
    Returns
    -------
    str
        The esolver type used in the ABACUS calculation.
    '''
    if isinstance(src, (str, Path)):
        with open(src) as f:
            raw = f.readlines()
    else: # assume the src is the return of the readlines()
        raw = src
    # with open(fn) as f:
    #     raw = f.readlines()
    raw = [l.strip() for l in raw]
    raw = [l for l in raw if l] # remove empty lines

    # search for the line with information like:
    # "#ENERGY SOLVER# ksdft_lcao"
    lines = [
        re.match(r'#ENERGY SOLVER#\s+(\S+)', l)
        for l in raw
    ]
    eslvtyp = [m.group(1) for m in lines if m is not None]
    assert len(set(eslvtyp)) == 1, \
           f'Inconsistent esolver type: {set(eslvtyp)}'
    return eslvtyp[0]

def read_band_from_running_log(src: str | Path | List[str]):
    '''in the latest branch, the band information is removed from running log.
    This function is only for backward compatibility. To get the band information,
    please call the "read_band_from_eig_occ" function.'''
    return None

def read_band_from_eig_occ(src: str | Path | List[str]) \
    -> List[Dict[str, np.ndarray]]:
    '''read the OUT.${suffix}/eig_occ.txt
    
    Parameters
    ----------
    src : str or Path or list of str
        The path to the ABACUS running log file or the return of the readlines() method.
    
    Returns
    -------
    list of dict
        A list of dictionaries containing the k-points and band energies.
        Each dictionary has keys 'k' and 'e', where 'k' is a list of k-point
        coordinates and 'e' is a numpy array of band energies.
    '''
    if isinstance(src, (str, Path)):
        with open(src) as f:
            raw = f.readlines()
    else: # assume the src is the return of the readlines()
        raw = src
    # with open(fn) as f:
    #     raw = f.readlines()
    raw = [l.strip() for l in raw]
    raw = [l for l in raw if l] # remove empty lines

    nframe = len([l for l in raw if l.endswith('# ionic step')])
    ekb_leading_pat = r'spin=(\d)\s+k-point=(\d+)/(\d+)\s+Cartesian=\s*' \
                    + r'(-?\d(\.\d+)?(e-\d+)?)\s+' \
                    + r'(-?\d(\.\d+)?(e-\d+)?)\s+' \
                    + r'(-?\d(\.\d+)?(e-\d+)?)\s+' \
                    + r'\(\d+\s+plane wave\)'
    # from the ekb leading line, there are nspin, nk information, also the 
    # coordinate of kpoints
    iekb = [i for i, l in enumerate(raw) if re.match(ekb_leading_pat, l)]
    assert len(iekb) > 0, f'No k-point found'
    
    k_raw = [re.match(ekb_leading_pat, raw[i]).groups() for i in iekb]
    
    # spin
    ispin = [int(g[0]) for g in k_raw]
    assert all(i in [1, 2] for i in ispin) # what about nspin 4?
    nspin = len(set(ispin))
    assert nspin in [1, 2]
    
    # k-points
    nk = [int(g[2]) for g in k_raw]
    nk = set(nk)
    assert len(nk) == 1, f'Sampling on k-points changed during the calculation: {nk}'
    nk = nk.pop()
    ik = [int(g[1]) for g in k_raw]
    assert all(i-1 in range(nk) for i in ik), f'k-point index out of range: {ik}'
    assert len(ik) == nframe * nspin * nk
    k  = np.array([[float(g[i]) for i in [3, 6, 9]] for g in k_raw])
    assert k.shape == (nframe * nspin * nk, 3)
    k = k.reshape((nframe, nspin, nk, 3))

    ekbpat  = r'\d+\s+'
    ekbpat += r'(-?\d+(\.\d+)?(e[+-]\d+)?)\s+'
    ekbpat += r'(\d+(\.\d+)?(e[+-]\d+)?)'
    iekb = [i for i, l in enumerate(raw) if re.match(ekbpat, l)]
    assert len(iekb) > 0, f'No band energy found'
    assert len(iekb) % (nspin*nk) == 0
    nb = len(iekb) // (nframe * nspin * nk)
    ekb_raw = np.array([list(map(float, raw[i].split()[1:])) for i in iekb])
    assert ekb_raw.shape == (nframe * nspin * nk * nb, 2)
    ekb_raw = ekb_raw.reshape(nframe, nspin, nk, nb, 2)
    ekb, occ = ekb_raw[:, :, :, :, 0], ekb_raw[:, :, :, :, 1]

    return [{'k': ki, 'e': eki, 'occ': occi}
            for ki, eki, occi in zip(k, ekb, occ)]

def read_traj_from_running_log(src: str | Path | List[str]) \
    -> List[Dict[str, np.ndarray|str]]:
    '''
    read the trajectory from the ABACUS running log file. 
    NOTE: in MD runs, the trajectory is not recorded in the running log file.

    Parameters
    ----------
    fn : str
        The path to the ABACUS running log file.
    
    Returns
    -------
    list of dict
        A list of dictionaries containing the coordinate system, cell, elements
        and the coordinates of the atoms. Each dictionary has keys 'coordinate',
        'cell', 'cell_unit', 'alat_in_angstrom', 'elem', and 'coords'. 
        The values are numpy arrays or strings.
        - 'coordinate': string, the coordinate system, e.g., 'Cartesian' or 'Direct'
        - 'cell': numpy array of shape (3, 3)
        - 'cell_unit': string, the unit of the cell, e.g., 'Angstrom'
        - 'alat_in_angstrom': float, the lattice constant in Angstrom
        - 'elem': list of strings, the chemical symbols of the elements
        - 'coords': numpy array of shape (natoms, 3), the coordinates of the atoms
    '''
    if isinstance(src, (str, Path)):
        with open(src) as f:
            raw = f.readlines()
    else: # assume the src is the return of the readlines()
        raw = src
    # with open(fn) as f:
    #     raw = f.readlines()
    raw = [l.strip() for l in raw]
    raw = [l for l in raw if l] # remove empty lines

    # search for the total number of atoms
    # r'TOTAL ATOM NUMBER = (\d+)'
    natoms = [re.search(r'TOTAL ATOM NUMBER\s*=\s*(\d+)', l) for l in raw]
    natoms = [int(n.group(1)) for n in natoms if n]
    assert len(set(natoms)) == 1, \
           f'Inconsistent number of atoms: {set(natoms)}'
    natoms = natoms[0]

    # search for the coordinate system
    # r'^([DIRECT|CARTESIAN]) COORDINATES$'
    coordinate = [re.match(r'^(DIRECT|CARTESIAN) COORDINATES', l) for l in raw]
    coordinate = [l.group(1).lower().capitalize() for l in coordinate if l]
    assert len(set(coordinate)) == 1, \
           f'Inconsistent coordinate system: {set(coordinate)}'
    coordinate = coordinate[0]

    # search for the cell, but first get the "a0": lattice constant
    # r'lattice constant (Angstrom) = (-?\d+(\.\d+)?)'
    a0 = [re.search(r'lattice constant \(Angstrom\)\s*=\s*(-?\d+(\.\d+)?)', l,
                    re.IGNORECASE) for l in raw]
    a0 = [float(n.group(1)) for n in a0 if n]
    assert len(set(a0)) == 1, f'Inconsistent lattice constant: {set(a0)}'
    a0 = a0[0]
    # then the cell
    # r'^Lattice vectors: \(Cartesian coordinate: in unit of a\_0\)$'
    icell = [i for i, l in enumerate(raw)
                if re.match(r'^Lattice vectors: \(Cartesian coordinate: in unit of a_0\)$', l)]
    assert len(icell) > 0, f'No cell found'
    # nframe = len(icell) # will be 1 for cell-invariant MD
    # assert nframe > 0, f'Invalid trajectory with length {nframe}')
    cell_raw = [raw[i+1:i+1+3] for i in icell]
    cell = [np.array([list(map(float, l.split())) for l in c]) * a0
            for c in cell_raw] # convert to Angstrom
    assert all(c.shape == (3, 3) for c in cell), \
           f'Unexpected shape of cell: {[c.shape for c in cell]}'

    # search for the elements and coordinates
    coord_leading_pat = r'atom\s+x\s+y\s+z\s+mag' # what about nspin 4?
    icoord = [i for i, l in enumerate(raw) if re.match(coord_leading_pat, l)]
    nframe = len(icoord)
    coord_raw = np.array([l.split() for i in icoord for l in raw[i+1:i+1+natoms]])
    elem = coord_raw[:, 0].astype(str).reshape((nframe, natoms))
    coords = coord_raw[:, 1:4].astype(float).reshape((nframe, natoms, 3))
    mag = coord_raw[:, 4].astype(float).reshape((nframe, natoms, ))

    cell = np.array(cell).reshape(-1, 3, 3)
    if cell.shape[0] == 1:
        cell = [cell[0] for _ in range(nframe)]
    assert len(cell) == nframe, \
           f'Unexpected number of cells: {len(cell)}, ' \
           f'expected {nframe}'
    return [{
        'coordinate': coordinate,
        'cell': c,
        'cell_unit': 'Angstrom',
        'alat_in_angstrom': a0,
        'elem': e,
        'coords': co
    } for c, e, co in zip(cell, elem, coords)]

def read_forces_from_running_log(src: str | Path | List[str]) \
    -> List[np.ndarray]:
    '''
    read the forces from the ABACUS running log file.

    Parameters
    ----------
    src : str or Path or list of str
        The path to the ABACUS running log file or the return of the readlines() method.
    
    Returns
    -------
    list of numpy array
        A list of numpy arrays containing the forces on each atom.
        Each array has shape (natoms, 3).
    '''
    if isinstance(src, (str, Path)):
        with open(src) as f:
            raw = f.readlines()
    else: # assume the src is the return of the readlines()
        raw = src
    # with open(fn) as f:
    #     raw = f.readlines()
    raw = [l.strip() for l in raw]
    raw = [l for l in raw if l] # remove empty lines

    # iteratively search for the forces, which led by the title `#TOTAL-FORCE (eV/Angstrom)#`
    forces, istart = [], 0
    while istart < len(raw):
        ith = None # index of the table header
        for i, l in enumerate(raw[istart:]):
            if re.match(r'#\s*TOTAL\-FORCE\s*\(eV\s*/Angstrom\)\s*#', l, re.IGNORECASE):
                ith = i
                break
        if ith is None: # no forces found
            break
        # otherwise
        # search for the first line that matches the pattern
        FORCEPAT_ = r'\s*([A-Z][a-z]?\d+)\s+(-?\d+(\.\d+)?)\s+(-?\d+(\.\d+)?)\s+(-?\d+(\.\d+)?)'
        itb = None # index of the first line of the table body
        for i, l in enumerate(raw[istart+ith+1:]):
            if re.match(FORCEPAT_, l):
                itb = i
                break
        if itb is None: # no content found
            break
        # otherwise
        jtb = None # index of the last line of the table body
        for j, l in enumerate(raw[istart+ith+1+itb:]):
            if not re.match(FORCEPAT_, l):
                jtb = j
                break
        if jtb is None: # no content found
            break
        
        # truncate the force table and append
        force_raw = raw[istart+ith+1+itb:istart+ith+1+itb+jtb]
        force = np.array([list(map(float, l.split()[1:])) for l in force_raw])
        forces.append(force)

        # update the istart
        istart += ith + itb + jtb + 1

    return forces

def read_stress_from_running_log(src: str | Path | List[str]) \
    -> List[np.ndarray]:
    '''
    read the stress from the ABACUS running log file.

    Parameters
    ----------
    src : str or Path or list of str
        The path to the ABACUS running log file or the return of the readlines() method.
    
    Returns
    -------
    list of numpy array
        A list of numpy arrays containing the stress on each atom.
        Each array has shape (3, 3).
    '''
    if isinstance(src, (str, Path)):
        with open(src) as f:
            raw = f.readlines()
    else: # assume the src is the return of the readlines()
        raw = src
    # with open(fn) as f:
    #     raw = f.readlines()
    raw = [l.strip() for l in raw]
    raw = [l for l in raw if l] # remove empty lines

    # iteratively search for the stress, which led by the title `#TOTAL-STRESS (kbar)#`
    stresses, istart = [], 0
    
    while istart < len(raw):
        ith = None # index of the table header
        for i, l in enumerate(raw[istart:]):
            if re.match(r'#\s*TOTAL\-STRESS\s*\(kbar\)\s*#', l, re.IGNORECASE):
                ith = i
                break
        if ith is None: # no stress found
            break
        # otherwise
        # search for the first line that matches the pattern
        STRESSPAT_ = r'\s*(-?\d+(\.\d+)?)\s+(-?\d+(\.\d+)?)\s+(-?\d+(\.\d+)?)'
        itb = None # index of the first line of the table body
        for i, l in enumerate(raw[istart+ith+1:]):
            if re.match(STRESSPAT_, l):
                itb = i
                break
        if itb is None: # no content found
            break
        # otherwise
        jtb = 3 # because the stress tensor would be a (3, 3)-matrix
        
        # truncate the stress table and append
        stress_raw = raw[istart+ith+1+itb:istart+ith+1+itb+jtb]
        stress = np.array([list(map(float, l.split())) for l in stress_raw]).reshape(3, 3)
        # unit: kbar -> GPa
        stresses.append(-0.1 * stress * GPa)

        # update the istart
        istart += ith + itb + jtb + 1

    return stresses

def read_abacus_out(fileobj, 
                    index=slice(None), 
                    results_required=True,
                    sort_atoms_with: Optional[List[int]] = None) -> Atoms | List[Atoms]:
    '''Reads the ABACUS output files. This function would be called by
    the AbacusTemplate.read_results() function. The detailed call stack
    is as follows:
       get_potential_energy()
    -> get_property()
    -> calculate()
    -> read_results()
    -> read_abacus_out() *here*

    To use this function as a standalone one, the fileobj should be
    the return of the open() function, which is a TextIOWrapper object:
    >>> with open(fn) as fileobj:
    ...     read_abacus_out(fileobj)

    Parameters
    ----------
    fileobj : str | Path | TextIOWrapper
        The file object to read.
    index : slice
        The index of the frames to read.
    results_required : bool
        Whether the results are required. If True, the results will be
        returned. If False, the results will not be returned. This parameter
        is not used.
    sort_atoms_with: Optional[List[int]]
        The sort order of the atoms. If not None, the atoms will be sorted
        according to the order in the list.

    Returns
    -------
    atoms : Atoms | List[Atoms]
        The atoms object, whose calculator is the `SinglePointDFTCalculator`.
    '''
    assert isinstance(fileobj, Path)
    # it is required that the running log file is passed as the Path object
    # so that the band structure file ``eig_occ.txt`` can be located as the
    # one in the same folder, otherwise the reading of elecstate will be
    # invalid and cause the failure
    with open(fileobj) as f:
        abacus_lines = f.readlines()
    
    # read the esolver type
    eslvtyp = read_esolver_type_from_running_log(abacus_lines)
    # FIXME: implement read_ksdft_esolver_out instead of read_abacus_out to
    # make flexible esolver support

    # read the structure, with the cell, elem, etc. (nframe)
    # if it is MD run, the trajectories will be in the file MD_dump, instead
    # of the running log
    trajectory = read_traj_from_running_log(abacus_lines) \
        if fileobj.name != 'running_md.log' else \
        read_traj_from_md_dump(fileobj.parent / 'MD_dump')

    # read the eigenvalues (nframe, nk, nbnd)
    elecstate = read_band_from_eig_occ(fileobj.parent / 'eig_occ.txt')
    # FIXME: remove thw following line till the eig_occ.txt is not written 
    #        in the append mode
    (fileobj.parent / 'eig_occ.txt').unlink()

    # read the atomic forces (nframe, nat, 3)
    forces = read_forces_from_running_log(abacus_lines)

    # read the stress (nframe, 3, 3), but may be None
    stress = read_stress_from_running_log(abacus_lines)

    # read all kpoints tables (but only want the first, spinless, with ibz2bz)
    k, _, _, _, _ = read_kpoints_from_running_log(abacus_lines)
    # unpack the kpoints information
    kvecd, wk, _ = k

    # FIXME: in principle, the two spin channels share the same set of
    # the kpoints, so it is not needed to have two sets of kpoints
    # and it is assumed that the sampling of kpoints won't change during
    # the simulation, which is, not exactly to be true for the NPT-MD
    # runs.

    _, energies   = read_energies_from_running_log(abacus_lines)
    # only keep the SCF converged energies
    energies = [edct for edct in energies if 'E_KS(sigma->0)' in edct]
    
    # read the magmom
    magmom = read_magmom_from_running_log(abacus_lines)

    # roughly check the integrity of data from their length consistency
    assert len(trajectory) == len(energies), \
        f'Inconsistent length: {len(trajectory)} != {len(energies)}'
    assert len(trajectory) == len(elecstate), \
        f'Inconsistent length: {len(trajectory)} != {len(elecstate)}'
    if len(forces) == 0:
        forces = [None] * len(trajectory)
    assert len(trajectory) == len(forces), \
        f'Inconsistent length: {len(trajectory)} != {len(forces)}'
    if len(stress) == 0:
        stress = [None] * len(trajectory)
    assert len(trajectory) == len(stress), \
        f'Inconsistent length: {len(trajectory)} != {len(stress)}'
    if len(magmom) == 0:
        magmom = [np.zeros(shape=(len(trajectory[0]['elem'])))] * len(trajectory)

    # loop over the frame...
    images, ind = [], sort_atoms_with
    for frame, estat, mag, frs, strs, ener in zip(
        trajectory, elecstate, magmom, forces, stress, energies):
        # for each frame, a structure can be defined
        ind = ind or list(range(len(frame['elem'])))
        atoms = Atoms(symbols=np.array(frame['elem'])[ind].tolist(), 
                      positions=frame['coords'][ind], 
                      cell=frame['cell'],
                      magmoms=mag[ind])
        # from result, a calculator can be assembled
        # however, sometimes the force and stress is not calculated
        # in this case, we set them to None
        frs  = None if is_invalid_arr(frs)  else frs[ind]
        strs = None if is_invalid_arr(strs) else full_3x3_to_voigt_6_stress(strs)
        calc = SinglePointDFTCalculator(atoms=atoms, energy=ener['E_KohnSham'],
                                        free_energy=ener['E_KohnSham'],
                                        forces=frs, stress=strs,
                                        magmoms=mag, efermi=ener['E_Fermi'],
                                        ibzkpts=kvecd, dipole=None)
        # import the eigenvalues and occupations kpoint-by-kpoint
        calc.kpts = []
        for ispn, (ekb, occ) in enumerate(zip(estat['e'], estat['occ'])): # loop over the spin
            calc.kpts += [SinglePointKPoint(weight=wk[ik],
                                            s=ispn,
                                            k=kvecd[ik],
                                            eps_n=ekb[ik,:],
                                            f_n=occ[ik,:])
                          for ik in range(len(kvecd))]
        # attach the calculator to the atoms
        atoms.calc = calc
        images.append(atoms)

    return images[index]

class TestLatestIO(unittest.TestCase):
    '''
    test the validities of functions for the I/O of the latest version of ABACUS
    '''
    here = Path(__file__).parent
    testfiles = here / 'testfiles'

    def test_read_esolver_type_from_running_log(self):
        self.assertEqual(
            read_esolver_type_from_running_log(
                self.testfiles / 'lcao-symm0-nspin2-multik-cellrelax'),
            'ksdft_lcao'
        )
        self.assertEqual(
            read_esolver_type_from_running_log(
                self.testfiles / 'lcao-symm0-nspin2-multik-relax'),
            'ksdft_lcao'
        )
        self.assertEqual(
            read_esolver_type_from_running_log(
                self.testfiles / 'lcao-symm1-nspin1-multik-scf'),
            'ksdft_lcao'
        )
        self.assertEqual(
            read_esolver_type_from_running_log(
                self.testfiles / 'pw-symm0-nspin4-gamma-md'),
            'ksdft_pw'
        )

    def test_read_band_from_running_log(self):
        self.assertIsNone(
            read_band_from_running_log(self.testfiles / 'lcao-symm0-nspin2-multik-cellrelax')
        )
    
    def test_read_band_from_eig_occ(self):
        fn = self.testfiles / 'nspin4-gamma-eigocc'
        elecstate = read_band_from_eig_occ(fn)
        # 2 frames, 1 spin channel (nspin 4), 1 kpoint, 35 bands
        self.assertEqual(len(elecstate), 2)
        for es in elecstate: # loop over frames
            self.assertIn('k', es)
            self.assertEqual(es['k'].shape, (1, 1, 3)) # ispin, ik, 3
            self.assertIn('e', es)
            self.assertEqual(es['e'].shape, (1, 1, 35)) # ispin, ik, nbnd
            self.assertIn('occ', es)
            self.assertEqual(es['occ'].shape, (1, 1, 35)) # ispin, ik, nbnd

    def test_read_traj_from_running_log(self):
        traj = read_traj_from_running_log(self.testfiles / 'lcao-symm0-nspin2-multik-cellrelax')
        self.assertEqual(len(traj), 3)

        reference = [
            np.array([[4.17203 , 2.086015, 2.086015],
                      [2.086015, 4.17203 , 2.086015],
                      [2.086015, 2.086015, 4.17203 ]]),
            np.array([[4.17550947, 2.05376104, 2.05376104],
                      [2.05811246, 4.12078078, 2.09979522],
                      [2.05811246, 2.09979522, 4.12078078]]),
            np.array([[4.17667347, 2.04296382, 2.04296382],
                      [2.04876712, 4.10361705, 2.10440948],
                      [2.04876712, 2.10440948, 4.10361705]])
        ]

        taud = np.array([[0.  , 1.  , 1.  ],
                         [0.5 , 0.5 , 0.5 ],
                         [0.25, 0.25, 0.25],
                         [0.75, 0.75, 0.75]])
        for t, c in zip(traj, reference):
            self.assertEqual(t['coordinate'], 'Direct')
            self.assertEqual(t['cell_unit'], 'Angstrom'),
            self.assertEqual(t['alat_in_angstrom'], 4.17203)
            self.assertTrue(all(e1 == e2 for e1, e2 in zip(t['elem'], np.array(['Ni1', 'Ni2', 'O', 'O']))))
            self.assertTrue(np.allclose(t['cell'], c))
            self.assertTrue(np.allclose(t['coords'] % 1, taud % 1))

    def test_read_forces_from_running_log(self):
        forces = read_forces_from_running_log(self.testfiles / 'lcao-symm0-nspin2-multik-relax')
        self.assertEqual(len(forces), 2) # two frames
        reference = [
            np.array([
                [-1.2994798821, -0.7502719805, 0.0000000000],
                [ 1.2994798821,  0.7502719805, 0.0000000000],
            ]),
            np.array([
                [-1.6814741322, -0.9708025265, 0.0000000000],
                [ 1.6814741322,  0.9708025265, 0.0000000000],
            ]),
        ]
        for f, ref in zip(forces, reference):
            self.assertTrue(np.allclose(f, ref))

    def test_read_stress_from_running_log(self):
        stresses = read_stress_from_running_log(
            self.testfiles / 'lcao-symm0-nspin2-multik-cellrelax')
        self.assertEqual(len(stresses), 2) # two frames
        reference = [
            np.array([
               [ 28.0257165781, -25.0814687477, -25.0814687477],
               [-25.0814687477, -57.2309475134,  52.8147857580],
               [-25.0814687477,  52.8147857580, -57.2309475134],
            ]),
            np.array([
                [ 10.5133147463,-28.8166870580, -28.8166870559],
                [-28.8166870580,  0.1640200010,   7.6016822662],
                [-28.8166870559,  7.6016822662,   0.1640199977],
            ]),
        ]
        for s, ref in zip(stresses, reference):
            self.assertTrue(np.allclose(s, -0.1 * GPa * ref))

    def test_read_abacus_out_string(self):
        fn = self.testfiles / 'pw-symm0-nspin4-gamma-md'
        # the case that does not give the instance of Path
        with self.assertRaises(AssertionError):
            read_abacus_out(str(fn))

    def test_read_pw_symm0_nspin4_gamma_md(self):
        # make files ready
        shutil.copy(self.testfiles / 'pw-symm0-nspin4-gamma-md', 
                    self.testfiles / 'running_md.log')
        shutil.copy(self.testfiles / 'nspin4-gamma-eigocc',
                    self.testfiles / 'eig_occ.txt')
        shutil.copy(self.testfiles / 'nspin4-gamma-mddump',
                    self.testfiles / 'MD_dump')

        res = read_abacus_out(self.testfiles / 'running_md.log')
        self.assertIsNotNone(res)
        self.assertEqual(len(res), 2) # two frames

        for atoms in res:
            self.assertIsInstance(atoms, Atoms)
            self.assertIsInstance(atoms.calc, SinglePointDFTCalculator)
            self.assertGreater(len(atoms.calc.kpts), 0)
            # Gamma point calculation
            for k in atoms.calc.kpts:
                self.assertIsInstance(k, SinglePointKPoint)

        # remove the files
        (self.testfiles / 'running_md.log').unlink()
        # (self.testfiles / 'eig_occ.txt').unlink()
        (self.testfiles / 'MD_dump').unlink()

if __name__ == '__main__':
    unittest.main()