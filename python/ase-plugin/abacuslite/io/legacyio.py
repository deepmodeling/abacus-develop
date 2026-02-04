import re
import shutil
import unittest
from io import TextIOWrapper
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
from ase.atoms import Atoms
from ase.calculators.singlepoint import (
    SinglePointKPoint,
    SinglePointDFTCalculator
)
from ase.units import Ry, eV, GPa, bar
from ase.constraints import full_3x3_to_voigt_6_stress
# from ase.utils import reader

__all__ = ['read_kpoints_from_running_log', 'read_istate',
           'read_band_from_running_log', 'read_traj_from_running_log',
           'read_traj_from_md_dump', 'read_forces_from_running_log',
           'read_stress_from_running_log', 'read_energies_from_running_log',
           'read_abacus_out']

def parse_kpoints_table(raw: List[str]) \
    -> Tuple[np.ndarray, np.ndarray, List[int] | None]:
    '''parse the lines that are the kpoints table
    
    Parameters
    ----------
    raw : list of str
        The lines of the kpoints table.
        
    Returns
    -------
    tuple
        A tuple containing the kpoints coordinates, weights and ibz2bz indexing,
        or without the indexing.
    '''
    # the pattern with the ibz2bz indexing
    W_PAT_ =  r'\s*(\d+)'
    W_PAT_ += r'\s+(-?\d\.\d+)\s+(-?\d\.\d+)\s+(-?\d\.\d+)'
    W_PAT_ += r'\s+(\d\.\d+)\s+(\d+)'
    # the pattern without the ibz2bz indexing
    WOPAT_ =  r'\s*(\d+)'
    WOPAT_ += r'\s+(-?\d\.\d+)\s+(-?\d\.\d+)\s+(-?\d\.\d+)'
    WOPAT_ += r'\s+(\d\.\d+)'
    # for each line, must match either the first or the second

    with_ibz = all(re.match(W_PAT_, l) for l in raw)
    assert with_ibz or all(re.match(WOPAT_, l) for l in raw), \
           f'Unexpected format of kpoints table: {raw}'
    if with_ibz:
        data = [re.match(W_PAT_, l).groups() for l in raw]
        k = np.array([list(map(float, ki[1:4])) for ki in data])
        w = np.array([float(ki[4]) for ki in data])
        ibz2bz = [int(ki[5]) for ki in data]
        return k, w, ibz2bz
    # otherwise, it must be the case without the ibz2bz
    data = [re.match(WOPAT_, l).groups() for l in raw]
    k = np.array([list(map(float, ki[1:4])) for ki in data])
    w = np.array([float(ki[4]) for ki in data])
    return k, w, None

def read_kpoints_from_running_log(src: str | Path | List[str]) \
    -> Tuple[Tuple[np.ndarray, np.ndarray, List[int]],
             Tuple[np.ndarray, np.ndarray],
             List[Tuple[np.ndarray, np.ndarray]],
             List[Tuple[np.ndarray, np.ndarray]],
             List[Tuple[np.ndarray, np.ndarray]]]:
    '''
    read the kpoint coordinates and weights from the running log file. Up to now,
    the ABACUS version 3.10.1 LTS, the ABACUS will print 4 tables of kpoints at 
    first, including the one with the ibz2bz indexing, one of spinless direct
    coordinates, then two spin cartesian and direct coordinates of kpoints. If
    the cell would change during the calculation, then after each cell-relax step,
    there will be 3 tables, including direct, cartesian and then direct (again,
    but i don't know why!) kpoint coordinates.

    Parameters
    ----------
    fn : str
        The path to the ABACUS running log file.
    
    Returns
    -------
    tuple
        A tuple contains the first table (with ibz2bz indexing), the second table
        (spinless direct coordinates), the trajectories of the third (spin-direct),
        fourth (spin-cartesian) and fifth (spin-direct) tables.
    '''
    if isinstance(src, (str, Path)):
        with open(src) as f:
            raw = f.readlines()
    else: # assume the src is the return of the readlines()
        raw = src
    # with open(fn) as f:
    #     raw = f.readlines()
    raw = [l.strip() for l in raw if l.strip()]
    raw = [l for l in raw if l] # remove empty lines

    # the pattern of table header
    THPAT_  = r'\s*(IBZ|KPOINTS)'
    THPAT_ += r'\s+(DIRECT|CARTESIAN)_X\s+(DIRECT|CARTESIAN)_Y\s+(DIRECT|CARTESIAN)_Z'
    THPAT_ += r'\s+WEIGHT(\s+ibz2bz)?'
    # the pattern of table content, optionally with the ibz2bz data
    TBPAT_  = r'\s*\d+'
    TBPAT_ += r'\s+(-?\d(\.\d+)?)\s+(-?\d(\.\d+)?)\s+(-?\d(\.\d+)?)'
    TBPAT_ += r'\s+(\d(\.\d+)?)(\s+\d+)?'

    # search for the kpoints table
    tables, istart = [], 0
    while istart < len(raw):
        # search for the table header...
        ith = None
        for i, l in enumerate(raw[istart:]):
            if re.match(THPAT_, l): # find a table header
                ith = i
                break
        if ith is None: # no more table
            break
        itb = ith + 1 # we assume there is no more header before the data
        assert re.match(TBPAT_, raw[istart+itb]) # really?

        # search for the end of the table...
        jtb = None
        for i, l in enumerate(raw[istart+itb:]):
            if not re.match(TBPAT_, l): # the next line is not a table content
                jtb = i
                break
        if jtb is None: # no more table
            break
        # parse the table
        ktab_raw = raw[istart+itb:istart+itb+jtb]
        tables.append(parse_kpoints_table(ktab_raw))

        # update the starting point to search for the next table
        istart += itb + jtb + 1
    
    # sort tables...
    kibz = tables.pop(0)
    # get the number of kpoints, then will use in reshaping the kpoints trajectories
    nk, _ = kibz[0].shape
    kspnls = (tables[0][0], tables[0][1]) # spinless direct, drop the None in [2]
    kd1traj = [(t[0].reshape(-1, nk, 3), t[1].reshape(-1, nk)) # drop the [2] because it is None
               for i, t in enumerate(tables) if i % 3 == 0]
    kctraj  = [(t[0].reshape(-1, nk, 3), t[1].reshape(-1, nk)) # drop the [2] because it is None
               for i, t in enumerate(tables) if i % 3 == 1]
    kd2traj = [(t[0].reshape(-1, nk, 3), t[1].reshape(-1, nk)) # drop the [2] because it is None
               for i, t in enumerate(tables) if i % 3 == 2]
    return kibz, kspnls, kd1traj, kctraj, kd2traj

def read_istate(src: str | Path | List[str]) \
    -> Dict[str, np.ndarray]:
    '''
    read the energy levels and occupations from istate.info file

    Parameters
    ----------
    fn : str
        The path to the istate.info file.
    
    Returns
    -------
    dict
        A dictionary containing the band energies, occupations, and k-points.
        The keys are 'ekb', 'occ', and 'k'. The values are numpy arrays with the
        following shapes:
        - 'ekb': (nspin, nk, nband)
        - 'occ': (nspin, nk, nband)
        - 'k': List of tuples containing the coordinates of k-points.
    '''
    if isinstance(src, (str, Path)):
        with open(src) as f:
            raw = f.readlines()
    else: # assume the src is the return of the readlines()
        raw = src
    # with open(fn) as f:
    #     raw = f.readlines()
    raw = [l.strip() for l in raw]
    raw = [l for l in raw if l]

    data  = [l for l in raw if l]
    ititl = [i for i, l in enumerate(data) if 
            re.match(r'BAND\s+Energy\(ev\)\s+Occupation\s+Kpoint\s+=', l)]
    idata = set(range(len(data))) - set(ititl)
    title = [data[i] for i in ititl]
    data =  [data[i] for i in idata]

    nk = len(title)
    nband = len(data) // nk
    ekb_raw = np.array([list(map(float, l.split())) for l in data])
    ekb = ekb_raw[:, 1::2].T.reshape(-1, nk, nband)
    nspin, _, _ = ekb.shape
    assert nspin in [1, 2], f'Unexpected number of spins: {nspin}'

    occ = ekb_raw[:, 2::2]
    occ = occ.T.reshape(-1, nk, nband)
    assert occ.shape == ekb.shape, \
           f'Unexpected shape of occupation numbers: {occ.shape}, expected {ekb.shape}'
    
    findk = lambda l: re.search(r'\((-?\d(\.\d+)?\s+-?\d(\.\d+)?\s+-?\d(\.\d+)?)\)', l)
    return {
        'ekb': ekb, 
        'occ': occ,
        'k': [tuple(map(float, findk(l).group(1).split())) for l in title]
    }

def read_band_from_running_log(src: str | Path | List[str]) \
    -> List[Dict[str, np.ndarray]]:
    '''
    read the band structure from the ABACUS running log file. This would be
    helpful for MD case that dftio does not support yet.

    Parameters
    ----------
    fn : str
        The path to the ABACUS running log file.
    
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
    
    # search for the number of spin channels, bands and kpoints first
    # r'nspin\s+=\s+(\d+)'
    nspin = [re.search(r'nspin\s+=\s+(\d+)', l) for l in raw]
    nspin = [int(n.group(1)) for n in nspin if n]
    assert len(set(nspin)) == 1, \
           f'Inconsistent number of spins: {set(nspin)}'
    nspin = nspin[0]
    # r'NBANDS\s+=\s+(\d+)'
    nbnd = [re.search(r'NBANDS\s+=\s+(\d+)', l) for l in raw]
    nbnd = [int(n.group(1)) for n in nbnd if n]
    assert len(set(nbnd)) == 1, \
           f'Inconsistent number of bands: {set(nbnd)}'
    nbnd = nbnd[0]
    # r'nkstot\snow\s=\s(\d+)'
    # (there may be not the nkstot_ibz for symmetry 0 case)
    # symmetry 1
    nkibz = [re.search(r'nkstot_ibz\s+=\s+(\d+)', l) for l in raw]
    nkibz = [int(n.group(1)) for n in nkibz if n]
    # symmetry 0, -1
    nktot = [re.search(r'nkstot\s+=\s+(\d+)', l) for l in raw]
    nktot = [int(n.group(1)) for n in nktot if n]
    nk = nktot if len(nkibz) == 0 else nkibz
    assert len(set(nk)) == 1, \
           f'Inconsistent number of k-points: {set(nk)}'
    nk = nk[0]

    # extract the band, first find the leading line
    # 333/511 kpoint (Cartesian) = 0.13022 0.052050 -6.5974e-10 (190 pws)
    ekb_leading_pat = r'\d+/\d+\s+kpoint\s+\(Cartesian\)\s+=\s+' \
                    + r'(-?\d(\.\d+)?(e-\d+)?)\s+' \
                    + r'(-?\d(\.\d+)?(e-\d+)?)\s+' \
                    + r'(-?\d(\.\d+)?(e-\d+)?)\s+' \
                    + r'\(\d+\s+pws\)'
    iekb = [i for i, l in enumerate(raw) if re.match(ekb_leading_pat, l)]
    assert len(iekb) > 0, f'No k-point found'
    assert len(iekb) % (nspin*nk) == 0, \
           f'Inconsistent number of k-points: {len(iekb)} vs {nk}'
    k_raw = [re.match(ekb_leading_pat, raw[i]).groups() for i in iekb]
    k = np.array([list(map(float, ki[::3])) for ki in k_raw])
    assert k.shape == (len(iekb), 3), \
           f'Unexpected shape of k-points: {k.shape}, expected ({len(iekb)}, 3)'

    nframe = len(iekb) // (nspin*nk) # number of frames, for MD or relax tasks
    ekb_raw = [l for i in iekb for l in raw[i+1:i+1+nbnd]]
    # each line should in the format of
    # r'\d+\s+(-?\d(\.\d+)?)\s+(\d(\.\d+)?)
    # Changelog: there are cases that the band energies and occupations are in scientific
    #            notation, e.g., 1.0e-01, 1.0e+00, so the regular expression should be
    #            r'\d+\s+(-?\d+(\.\d+)?(e[+-]\d+)?)\s+(\d+(\.\d+)?(e[+-]\d+)?)' instead
    ekbpat  = r'\d+\s+'
    ekbpat += r'(-?\d+(\.\d+)?(e[+-]\d+)?)\s+'
    ekbpat += r'(\d+(\.\d+)?(e[+-]\d+)?)'
    assert all(re.match(ekbpat, l) for l in ekb_raw), \
           'Unexpected format of band energies: \n' + '\n'.join(ekb_raw)
    # ekb in the second column, occ in the third column
    ekb_raw = np.array([list(map(float, l.split())) for l in ekb_raw])
    assert ekb_raw.shape == (nframe * nspin * nk * nbnd, 3), \
           f'Unexpected shape of band energies: {ekb_raw.shape}. ' \
           f'Expected ({nframe * nspin * nk * nbnd}, 3), in which ' \
           f'nframe={nframe}, nspin={nspin}, nk={nk}, nbnd={nbnd}'
    ekb = ekb_raw[:, 1].reshape(nframe, nspin, nk, nbnd)
    occ = ekb_raw[:, 2].reshape(nframe, nspin, nk, nbnd)
    # reshape the k-points to (nframe, nspin, nk, 3)
    k = k.reshape(nframe, nspin, nk, 3)

    return [{'k': ki, 'e': eki, 'occ': occi} for ki, eki, occi in zip(k, ekb, occ)]

def read_traj_from_running_log(src: str | Path | List[str]) \
    -> List[Dict[str, np.ndarray|str]]:
    '''
    read the trajectory from the ABACUS running log file. This would be
    helpful for MD case that dftio does not support yet.

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
    # nframe = len(icell) # will be 1 for NVT MD
    # assert nframe > 0, f'Invalid trajectory with length {nframe}')
    cell_raw = [raw[i+1:i+1+3] for i in icell]
    cell = [np.array([list(map(float, l.split())) for l in c]) * a0
            for c in cell_raw] # convert to Angstrom
    assert all(c.shape == (3, 3) for c in cell), \
           f'Unexpected shape of cell: {[c.shape for c in cell]}'
    
    # search for the elements and coordinates
    # r'^tau[c|d]_([A-Z][a-z]?)\d+\s+(-?\d+(\.\d+)?)\s+(-?\d+(\.\d+)?)\s+(-?\d+(\.\d+)?)*'
    m_tau = [re.match(r'^tau[c|d]_([A-Z][a-z]?)\d+'
                        r'\s+(-?\d+(\.\d+)?)\s+(-?\d+(\.\d+)?)\s+(-?\d+(\.\d+)?)', l)
                for l in raw]
    m_tau = [m for m in m_tau if m]
    assert len(m_tau) > 0, f'No atoms found'
    nframe = len(m_tau) // natoms
    assert nframe > 0, f'Invalid trajectory with length {nframe}'
    elem = [m.group(1) for m in m_tau]
    coords = np.array([[float(m.group(2)), float(m.group(4)), float(m.group(6))]
                        for m in m_tau]).reshape(-1, natoms, 3)
    assert coords.shape == (nframe, natoms, 3), \
           f'Unexpected shape of coordinates: {coords.shape}, ' \
           f'expected ({nframe}, {natoms}, 3)'
    assert len(elem) == natoms * nframe, \
           f'Unexpected number of elements: {len(elem)}, ' \
           f'expected {natoms * nframe}'
    elem = [elem[i:i+natoms] for i in range(0, len(elem), natoms)]

    # final: for volume-constant run, the cell information will only be printed for once
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

def read_traj_from_md_dump(src: str | Path | List[str]) \
    -> List[Dict[str, np.ndarray]]:
    '''
    read the trajectory from the ABACUS MD dump file

    Parameters
    ----------
    fn : str
        The path to the ABACUS MD dump file.
    
    Returns
    -------
    list of dict
        A list of dictionaries containing the coordinate system, cell, elements
        and the coordinates of the atoms. Each dictionary has keys 'coordinate',
        'cell', 'alat_in_angstrom', 'elem', and 'coords'.
        The values are numpy arrays or strings.
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

    # search for the lattice constant
    # r'LATTICE_CONSTANT:\s+(-?\d+(\.\d+)?) Angstrom'
    a0 = [re.search(r'LATTICE_CONSTANT:\s+(-?\d+(\.\d+)?) Angstrom', l) for l in raw]
    a0 = [float(n.group(1)) for n in a0 if n]
    assert len(set(a0)) == 1, \
                f'Inconsistent lattice constant: {set(a0)}'
    a0 = a0[0]
    # search for the cell
    # r'^LATTICE_VECTORS$'
    icell = [i for i, l in enumerate(raw) if re.match(r'^LATTICE_VECTORS$', l)]
    assert len(icell) > 0, f'No cell found in file'
    cell_raw = [raw[i+1:i+1+3] for i in icell]
    cell = [np.array([list(map(float, l.split())) for l in c])
            for c in cell_raw]
    assert all(c.shape == (3, 3) for c in cell), \
                f'Unexpected shape of cell: {[c.shape for c in cell]}'
    nframe = len(cell)

    # search for the elements and coordinates
    # r'^\d+\s+([A-Z][a-z]?)\s+(-?\d+(\.\d+)?)\s+(-?\d+(\.\d+)?)\s+(-?\d+(\.\d+)?)*'
    m_tau = [re.match(r'^\d+\s+([A-Z][a-z]?)\s+'
                        r'(-?\d+(\.\d+)?)\s+(-?\d+(\.\d+)?)\s+(-?\d+(\.\d+)?)', l)
                for l in raw]
    m_tau = [m for m in m_tau if m]
    assert len(m_tau) > 0, f'No atoms found'
    natoms = len(m_tau) // nframe
    assert natoms > 0, f'Invalid trajectory with length {nframe}.'
    elem = [m.group(1) for m in m_tau]
    coords = np.array([[float(m.group(2)), float(m.group(4)), float(m.group(6))]
                        for m in m_tau]).reshape(nframe, natoms, 3)
    assert coords.shape == (nframe, natoms, 3), \
           f'Unexpected shape of coordinates: {coords.shape}, ' \
           f'expected ({nframe}, {natoms}, 3)'
    assert len(elem) == natoms * nframe, \
           f'Unexpected number of elements: {len(elem)}, ' \
           f'expected {natoms * nframe}'
    elem = [elem[i:i+natoms] for i in range(0, len(elem), natoms)]
    assert len(elem) == nframe, \
           f'Unexpected number of elements: {len(elem)}, ' \
           f'expected {nframe}'
    
    return [{
        'coordinate': 'Cartesian',
        'cell': c,
        'alat_in_angstrom': a0,
        'elem': e,
        'coords': co
    } for c, e, co in zip(cell, elem, coords)]

def read_forces_from_running_log(src: str | Path | List[str]) \
    -> List[np.ndarray]:
    ''''''
    if isinstance(src, (str, Path)):
        with open(src) as f:
            raw = f.readlines()
    else: # assume the src is the return of the readlines()
        raw = src
    # with open(fn) as f:
    #     raw = f.readlines()
    raw = [l.strip() for l in raw]
    raw = [l for l in raw if l] # remove empty lines

    # iteratively search for the forces, which led by the title `TOTAL-FORCE`
    forces, istart = [], 0
    
    while istart < len(raw):
        ith = None # index of the table header
        for i, l in enumerate(raw[istart:]):
            if re.match(r'\s*TOTAL\-FORCE\s*\(eV\s*/Angstrom\)', l, re.IGNORECASE):
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
    ''''''
    if isinstance(src, (str, Path)):
        with open(src) as f:
            raw = f.readlines()
    else: # assume the src is the return of the readlines()
        raw = src
    # with open(fn) as f:
    #     raw = f.readlines()
    raw = [l.strip() for l in raw]
    raw = [l for l in raw if l] # remove empty lines

    # iteratively search for the stress, which led by the title `TOTAL-STRESS`
    stresses, istart = [], 0
    
    while istart < len(raw):
        ith = None # index of the table header
        for i, l in enumerate(raw[istart:]):
            if re.match(r'\s*TOTAL\-STRESS\s*\(KBAR\)', l, re.IGNORECASE):
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

def read_energies_from_running_log(src: str | Path | List[str]) \
    -> List[Dict[str, float]]:
    ''''''
    if isinstance(src, (str, Path)):
        with open(src) as f:
            raw = f.readlines()
    else: # assume the src is the return of the readlines()
        raw = src
    # with open(fn) as f:
    #     raw = f.readlines()
    raw = [l.strip() for l in raw]
    raw = [l for l in raw if l] # remove empty lines

    energies_ry, energies_ev, istart = [], [], 0
    
    while istart < len(raw):
        ith = None # index of the table header
        for i, l in enumerate(raw[istart:]):
            if re.match(r'\s*ENERGY\s+Rydberg\s+eV', l, re.IGNORECASE):
                ith = i
                break
        if ith is None: # no energies found
            break
        # otherwise
        # search for the first line that matches the pattern
        ENERGYPAT_ = r'\s*E_(\S+)\s+(-?\d+(\.\d+)?)\s+(-?\d+(\.\d+)?)'
        itb = None # index of the first line of the table body
        for i, l in enumerate(raw[istart+ith+1:]):
            if re.match(ENERGYPAT_, l):
                itb = i
                break
        if itb is None: # no content found
            break
        # otherwise
        jtb = None
        for j, l in enumerate(raw[istart+ith+1+itb:]):
            if not re.match(ENERGYPAT_, l):
                jtb = j
                break
        if jtb is None: # no content found
            break
        
        # truncate the energy table and append
        tb_raw = raw[istart+ith+1+itb:istart+ith+1+itb+jtb]
        # first item is the name of the energy component
        # the second and third items are the ry and ev values
        name, e_ry, e_ev = zip(*[l.split() for l in tb_raw])
        energies_ry.append(dict(zip(name, list(map(float, e_ry)))))
        energies_ev.append(dict(zip(name, list(map(float, e_ev)))))

        # update the istart
        istart += ith + itb + jtb + 1

    return energies_ry, energies_ev

def read_magmom_from_running_log(src: str | Path | List[str]) \
    -> List[np.ndarray]:
    '''
    Read the magnetic momentum from ABACUS running log. Note:
    this function is for reading the output from nspin=4 case,
    and returns the list of np.ndarray which has dimension of
    (nat, 3), where nat is the number of atoms, and 3 stands for
    the x, y and z components. For nspin=1/2, the magnetic moment 
    is not printed in the log.
    There will be problem if the output directly imported to the
    ASE SinglePointDFTCalculator, because the magmom required by
    ASE-side is the (N,) array, so another operation `np.linalg.norm
    (magmom, axis=1)` should be called to postprocess the magmom.

    Parameters
    ----------
    src : str | Path | List[str]
        The source of the running log. If it is a string, it will be
        treated as the path to the running log. If it is a list of
        strings, it will be treated as the return of the readlines()
        function.

    Returns
    -------
    magmom : List[np.ndarray]
        The trajectories of the magnetic moment of each atom.
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
    
    magmom, istart = [], 0
    while istart < len(raw):
        ith = None # index of the table header
        for i, l in enumerate(raw[istart:]):
            if re.match(r'\s*Total\sMagnetism\s\(uB\)\s+x\s+y\s+z\s*', l, 
                        re.IGNORECASE):
                ith = i
                break
        if ith is None: # no magmom found
            break
        # otherwise
        # search for the first line that matches the pattern
        MAGMOMPAT_  = r'\s*([A-Z][a-z]?\d+)'
        MAGMOMPAT_ += r'\s+([-+]?\d+\.\d+)'
        MAGMOMPAT_ += r'\s+([-+]?\d+\.\d+)'
        MAGMOMPAT_ += r'\s+([-+]?\d+\.\d+)\s*'
        itb = None # index of the first line of the table body
        for i, l in enumerate(raw[istart+ith+1:]):
            if re.match(MAGMOMPAT_, l):
                itb = i
                break
        if itb is None: # no content found
            break
        # otherwise
        jtb = None
        for j, l in enumerate(raw[istart+ith+1+itb:]):
            if not re.match(MAGMOMPAT_, l):
                jtb = j
                break
        if jtb is None: # no content found
            break
        
        # truncate the magmom table and append
        tb_raw = raw[istart+ith+1+itb:istart+ith+1+itb+jtb]
        # first item is the name of the magmom component
        # the second to fourth items are the x, y, z components
        _, mx, my, mz = zip(*[l.split() for l in tb_raw])
        magmom.append(np.array([list(map(float, (mx, my, mz))) 
                                for mx, my, mz in zip(mx, my, mz)]))
        
        # update the istart
        istart += ith + itb + jtb + 1

    return magmom

def is_invalid_arr(arr) -> bool:
    '''Check if the array is invalid, including the cases of None,
    empty array, and array with NaN values.

    Parameters
    ----------
    arr : np.ndarray
        The array to check.

    Returns
    -------
    is_invalid : bool
        Whether the array is invalid.
    '''
    if arr is None:
        return True
    if isinstance(arr, list):
        if len(arr) == 0:
            return True
    if isinstance(arr, np.ndarray):
        if len(arr) == 0:
            return True
        if np.isnan(arr).any():
            return True
    return False

# @reader
def read_abacus_out(fileobj, index=slice(None), results_required=True) \
    -> Atoms | List[Atoms]:
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

    Returns
    -------
    atoms : Atoms | List[Atoms]
        The atoms object, whose calculator is the `SinglePointDFTCalculator`.
    '''
    if isinstance(fileobj, (str, Path)):
        with open(fileobj) as f:
            abacus_lines = f.readlines()
    else: # from the `with open(fn) as fileobj:` context
        assert isinstance(fileobj, TextIOWrapper)
        abacus_lines = fileobj.readlines()
    
    # read the structure, with the cell, elem, etc. (nframe)
    trajectory = read_traj_from_running_log(abacus_lines)
    # read the eigenvalues (nframe, nk, nbnd)
    elecstate = read_band_from_running_log(abacus_lines)
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
    images = []
    for frame, estat, mag, frs, strs, ener in zip(
        trajectory, elecstate, magmom, forces, stress, energies):
        # for each frame, a structure can be defined
        atoms = Atoms(symbols=frame['elem'], positions=frame['coords'], cell=frame['cell'])
        # from result, a calculator can be assembled
        # however, sometimes the force and stress is not calculated
        # in this case, we set them to None
        frs  = None if is_invalid_arr(frs)  else frs
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

class TestCalculatorAbacusLegacyIO(unittest.TestCase):

    def setUp(self):
        here = Path(__file__).parent
        self.testfiles = here.parent.parent / 'testfiles'

    def test_read_istate(self):
        fn = self.testfiles / 'istate_nspin1.info'
        data = read_istate(fn)
        self.assertIn('ekb', data)
        self.assertIn('occ', data)
        self.assertIn('k', data)
        ekb, occ, k = data.values()
        self.assertEqual(len(k), 29)
        self.assertTrue(all(isinstance(ki, tuple) and len(ki) == 3 for ki in k))
        self.assertEqual(ekb.shape, (1, 29, 20))
        self.assertEqual(occ.shape, (1, 29, 20))

    def test_read_band_from_running_log(self):
        # nspin1
        fn = self.testfiles / 'running_scf_log_nspin1'
        data = read_band_from_running_log(fn)
        self.assertIsInstance(data, list)
        self.assertGreater(len(data), 0)
        self.assertTrue(all('k' in d and 'e' in d and 'occ' in d for d in data))
        nspin, nk, nband = data[0]['e'].shape
        self.assertEqual(nspin, 1)
        self.assertEqual(nband, 20)
        self.assertEqual(nk, 29)
        for d in data: # for each frame
            self.assertTrue(d['k'].shape == (nspin, nk, 3))
            self.assertTrue(d['e'].shape == (nspin, nk, nband))
            self.assertTrue(d['occ'].shape == (nspin, nk, nband))

        # nspin2
        fn = self.testfiles / 'running_scf_log_nspin2'
        data = read_band_from_running_log(fn)
        self.assertIsInstance(data, list)
        self.assertGreater(len(data), 0)
        self.assertTrue(all('k' in d and 'e' in d and 'occ' in d for d in data))
        nspin, nk, nband = data[0]['e'].shape
        self.assertEqual(nspin, 2)
        self.assertEqual(nband, 20)
        self.assertEqual(nk, 29)
        for d in data: # for each frame
            self.assertTrue(d['k'].shape == (nspin, nk, 3))
            self.assertTrue(d['e'].shape == (nspin, nk, nband))
            self.assertTrue(d['occ'].shape == (nspin, nk, nband))
        
        # nspin 2, cell-relax (multi-frames)
        fn = self.testfiles / 'running_cell_relax_log_nspin2'
        data = read_band_from_running_log(fn)
        self.assertIsInstance(data, list)
        self.assertEqual(len(data), 3)
        self.assertTrue(all('k' in d and 'e' in d and 'occ' in d for d in data))
        nspin, nk, nband = data[0]['e'].shape
        self.assertEqual(nspin, 2)
        self.assertEqual(nband, 20)
        self.assertEqual(nk, 29)
        for d in data: # for each frame
            self.assertTrue(d['k'].shape == (nspin, nk, 3))
            self.assertTrue(d['e'].shape == (nspin, nk, nband))
            self.assertTrue(d['occ'].shape == (nspin, nk, nband))

        # nspin 2, MD (multi-frames)
        fn = self.testfiles / 'running_md_log_nspin2'
        data = read_band_from_running_log(fn)
        self.assertIsInstance(data, list)
        self.assertEqual(len(data), 11) # from 0 to 10
        self.assertTrue(all('k' in d and 'e' in d and 'occ' in d for d in data))
        nspin, nk, nband = data[0]['e'].shape
        self.assertEqual(nspin, 2)
        self.assertEqual(nband, 20)
        self.assertEqual(nk, 260) # due to nosym = .true.
        for d in data: # for each frame
            self.assertTrue(d['k'].shape == (nspin, nk, 3))
            self.assertTrue(d['e'].shape == (nspin, nk, nband))
            self.assertTrue(d['occ'].shape == (nspin, nk, nband))

    def test_read_traj_from_running_log(self):
        fn = self.testfiles / 'running_scf_log_nspin1'
        data = read_traj_from_running_log(fn)
        self.assertIsInstance(data, list)
        self.assertEqual(len(data), 1)
        self.assertTrue(all('coordinate' in d and 'cell' in d and 
                            'elem' in d and 'coords' in d for d in data))
        for d in data:
            self.assertIn(d['coordinate'], ['Cartesian', 'Direct'])
            self.assertEqual(d['cell'].shape, (3, 3))
            self.assertIsInstance(d['elem'], list)
            self.assertEqual(len(d['elem']), 2)

        # cell-relax
        fn = self.testfiles / 'running_cell_relax_log_nspin2'
        data = read_traj_from_running_log(fn)
        self.assertIsInstance(data, list)
        self.assertEqual(len(data), 3)
        self.assertTrue(all('coordinate' in d and 'cell' in d and 
                            'elem' in d and 'coords' in d for d in data))
        for d in data:
            self.assertIn(d['coordinate'], ['Cartesian', 'Direct'])
            self.assertEqual(d['cell'].shape, (3, 3))
            self.assertIsInstance(d['elem'], list)
            self.assertEqual(len(d['elem']), 2)
            self.assertEqual(d['coords'].shape, (2, 3))

    def test_read_traj_from_md_dump(self):
        fn = self.testfiles / 'MD_dump'
        data = read_traj_from_md_dump(fn)
        self.assertIsInstance(data, list)
        self.assertEqual(len(data), 11)
        self.assertTrue(all('coordinate' in d and 'cell' in d and 
                            'alat_in_angstrom' in d and 'elem' in d and 
                            'coords' in d for d in data))
        for d in data:
            self.assertEqual(d['coordinate'], 'Cartesian')
            self.assertEqual(d['cell'].shape, (3, 3))
            self.assertEqual(d['alat_in_angstrom'], 0.999615353000)
            self.assertIsInstance(d['elem'], list)
            self.assertEqual(len(d['elem']), 2)
            self.assertEqual(d['coords'].shape, (2, 3))        

    def test_read_forces_from_running_log(self):
        fn = self.testfiles / 'running_md_log_nspin2'
        forces = read_forces_from_running_log(fn)
        self.assertEqual(len(forces), 11) # 11 frames
        self.assertTrue(all(f.shape == (2, 3) for f in forces)) # 2 atoms, 3 components

    def test_read_stress_from_running_log(self):
        fn = self.testfiles / 'running_cell_relax_log_nspin2'
        stress = read_stress_from_running_log(fn)
        self.assertEqual(len(stress), 3) # 3 cell-relax steps
        self.assertTrue(all(s.shape == (3, 3) for s in stress)) # 3x3 matrix

        reference = [
            np.array([[ 3.6005261437, -0.0000000000,  0.0000000000], 
                      [ 0.0000000000,  3.6005261437,  0.0000000000],
                      [-0.0000000000,  0.0000000000,  3.6005261437]]),
            np.array([[-2.6583368803,  0.0000000000,  0.0000000000],
                      [ 0.0000000000, -2.6583368803,  0.0000000000],
                      [ 0.0000000000,  0.0000000000, -2.6583368803]]),
            np.array([[-0.0138062589,  0.0000000000,  0.0000000000],
                      [-0.0000000000, -0.0138062589, -0.0000000000],
                      [ 0.0000000000,  0.0000000000, -0.0138062589]])
        ]
        for s, sref in zip(stress, reference):
            self.assertTrue(np.allclose(s, sref))

    def test_read_energies_from_running_log(self):
        fn = self.testfiles / 'running_scf_log_nspin2'
        energies_ry, energies_ev = read_energies_from_running_log(fn)
        self.assertIsInstance(energies_ry, list)
        self.assertIsInstance(energies_ev, list)
        self.assertEqual(len(energies_ry), 13) # 19 SCF steps
        self.assertEqual(len(energies_ev), 13)
        self.assertTrue(all(isinstance(e, dict) for e in energies_ry)) # dict of energies
        self.assertTrue(all(isinstance(e, dict) for e in energies_ev))
        for e_ev, e_ry in zip(energies_ev, energies_ry):
            # take energies of one SCF step, there are still many energy terms
            self.assertEqual(len(e_ev), len(e_ry))
            for k, ei_ev in e_ev.items():
                self.assertIn(k, e_ry)
                self.assertAlmostEqual(ei_ev, e_ry[k] * Ry / eV, delta=1e-3)

    def test_read_kpoints_from_running_log(self):
        fn = self.testfiles / 'running_cell_relax_log_nspin2'
        kpoints = read_kpoints_from_running_log(fn)
        self.assertIsInstance(kpoints, tuple)
        self.assertEqual(len(kpoints), 5)
        # thus we unpack
        kibz, kdspnls, kd1traj, kctraj, kd2traj = kpoints
        
        # kibz
        self.assertIsInstance(kibz, tuple)
        self.assertEqual(len(kibz), 3)
        # thus we unpack
        kibz, wk, ibz2bz = kibz
        self.assertIsInstance(kibz, np.ndarray)
        self.assertEqual(kibz.shape, (29, 3))
        self.assertIsInstance(wk, np.ndarray)
        self.assertEqual(wk.shape, (29,))
        self.assertAlmostEqual(wk.sum(), 1.0, delta=1e-3)
        self.assertIsInstance(ibz2bz, list)
        self.assertEqual(len(ibz2bz), 29)

        # kdspnls
        self.assertIsInstance(kdspnls, tuple)
        self.assertEqual(len(kdspnls), 2)
        # thus we unpack
        kdspnls, wk = kdspnls
        self.assertIsInstance(kdspnls, np.ndarray)
        self.assertEqual(kdspnls.shape, (29, 3))
        self.assertIsInstance(wk, np.ndarray)
        self.assertEqual(wk.shape, (29,))
        self.assertAlmostEqual(wk.sum(), 1.0, delta=1e-3)

        # kd1traj
        self.assertIsInstance(kd1traj, list)
        self.assertEqual(len(kd1traj), 3) # 3 cell-relax steps
        for i, kd1 in enumerate(kd1traj):
            # because the i=0 corresponds to the spinless case,
            # it has 29 instead of 58 kpoints like the others
            nk = 29
            nspin = 1 if i == 0 else 2
            wktot = 1 if i == 0 else 2
            self.assertIsInstance(kd1, tuple)
            self.assertEqual(len(kd1), 2)
            # thus we unpack
            kd1, wk = kd1
            self.assertIsInstance(kd1, np.ndarray)
            self.assertEqual(kd1.shape, (nspin, nk, 3))
            self.assertIsInstance(wk, np.ndarray)
            self.assertEqual(wk.shape, (nspin, nk))
            self.assertAlmostEqual(wk.sum(), wktot, delta=1e-3)

        # kctraj and kd2traj
        for ktraj in [kctraj, kd2traj]:
            self.assertIsInstance(ktraj, list)
            self.assertEqual(len(ktraj), 3) # 3 cell-relax steps
            for i, k in enumerate(ktraj):
                self.assertIsInstance(k, tuple)
                self.assertEqual(len(k), 2)
                # thus we unpack
                k, wk = k
                self.assertIsInstance(k, np.ndarray)
                self.assertEqual(k.shape, (2, 29, 3))
                self.assertIsInstance(wk, np.ndarray)
                self.assertEqual(wk.shape, (2, 29))
                self.assertAlmostEqual(wk.sum(), 2.0, delta=1e-3)

    def test_read_magmom_from_running_log(self):
        fn = self.testfiles / 'running_relax_log_nspin4'
        magmom = read_magmom_from_running_log(fn)
        self.assertIsInstance(magmom, list)
        self.assertEqual(len(magmom), 5) # 5 relax steps

        reference = np.array([[ 1.8225516369,  3.1480265966,  0.0000000000], 
                              [ 1.8147074392, -3.1525791944,  0.0000000000], 
                              [-3.6376745424,  0.0045350365,  0.0000000000], 
                              [ 1.8225516457,  3.1480266013,  0.0000000000], 
                              [-3.6376745306,  0.0045350385,  0.0000000000], 
                              [ 1.8147074473, -3.1525791936,  0.0000000000], 
                              [ 0.0000244792, -0.0000033218,  0.0000000000], 
                              [ 0.0000244819, -0.0000033210,  0.0000000000]])
        self.assertTrue(np.allclose(magmom[0], reference))
        magmom_norm = np.linalg.norm(magmom[0], axis=1)
        reference = np.array([3.6375494392,
                              3.6375704346,
                              3.6376773692,
                              3.6375494476,
                              3.6376773575,
                              3.6375704379,
                              0.0000247036,
                              0.0000247061])
        self.assertTrue(np.allclose(magmom_norm, reference))
        reference = np.array([[ 1.8404923963,  3.1109885127,  0.0000000000], 
                              [ 1.7740335195, -3.1492677583,  0.0000000000], 
                              [-3.6144130623,  0.0376416164,  0.0000000000], 
                              [ 1.8405035635,  3.1110434653,  0.0000000000], 
                              [-3.6143459346,  0.0376369850,  0.0000000000], 
                              [ 1.7740167753, -3.1492941455,  0.0000000000], 
                              [-0.0000086899,  0.0000432242,  0.0000000000], 
                              [-0.0000028719,  0.0000370664,  0.0000000000]])
        self.assertTrue(np.allclose(magmom[-1], reference))
        magmom_norm = np.linalg.norm(magmom[-1], axis=1)
        reference = np.array([3.6146454580,
                              3.6145653047,
                              3.6146090627,
                              3.6146984397,
                              3.6145418904,
                              3.6145800771,
                              0.0000440891,
                              0.0000371775])
        self.assertTrue(np.allclose(magmom_norm, reference))

    @unittest.skip('This functionality to test now works smoothly. Because it is only a container'
                   ', no need to test in CI/CD')
    def test_read_abacus_out(self):
        fn = self.testfiles / 'running_scf_log_nspin1'
        with open(fn) as f:
            for atom in read_abacus_out(f):
                for ik, k in enumerate(atom.calc.kpts):
                    print(f'KPOINT-{ik}th')
                    print(k.weight)
                    print(k.s)
                    print(k.k)
                    print(k.eps_n)
                    print(k.f_n)
                    print('----')



if __name__ == '__main__':
    unittest.main()
