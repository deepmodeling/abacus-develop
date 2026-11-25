'''implements the parser for the latest version of ABACUS'''
import re
import shutil
import unittest
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
from ase.atoms import Atoms
from ase.io import read as read_cif_atoms
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
    find_final_info_with_iter_header
)


def read_traj_from_md_stru(directory: Path,
                           running_log_lines: Optional[List[str]] = None) -> List[Dict[str, np.ndarray]]:
    """Read MD frames from appended, numbered, or final STRU output."""
    periodic = directory / 'STRU'
    if not periodic.is_file():
        periodic = directory / 'STRU' / 'STRU'
    numbered = sorted((path for path in (directory / 'STRUC').glob('STRU[0-9]*')
                       if re.fullmatch(r'STRU[0-9]+', path.name)),
                      key=lambda path: int(path.name[4:]))
    sources = [periodic] if periodic.is_file() else numbered
    final = directory / 'STRU_FINAL'
    if not final.is_file():
        final = directory / 'STRUC' / 'STRU_FINAL'
    if final.is_file():
        sources.append(final)
    if running_log_lines is not None:
        announced = []
        for line in running_log_lines:
            match = re.search(r'Write file:\s+(.+?)\s*$', line)
            if match is None:
                continue
            path = Path(match.group(1))
            if path.name in ('STRU', 'STRU_FINAL') or re.fullmatch(r'STRU[0-9]+', path.name):
                if path.parent.name in ('STRUC', 'STRU'):
                    path = directory / path.parent.name / path.name
                else:
                    path = directory / path.name
                if path.is_file() and path not in announced:
                    announced.append(path)
        if announced:
            sources = announced
    frames = {}
    for source in sources:
        step = None
        lines = []
        for raw in source.read_text().splitlines() + ['# MD step END']:
            if raw.startswith('# MD step '):
                if step is not None:
                    data = [line.split('#', 1)[0].strip() for line in lines]
                    data = [line for line in data if line]
                    lat0 = float(data[data.index('LATTICE_CONSTANT') + 1])
                    lattice_index = data.index('LATTICE_VECTORS') + 1
                    cell = np.array([list(map(float, row.split()))
                                     for row in data[lattice_index:lattice_index + 3]])
                    cell *= lat0 * 0.529177
                    index = data.index('ATOMIC_POSITIONS') + 1
                    coordinate = data[index]
                    assert coordinate == 'Cartesian_angstrom', coordinate
                    index += 1
                    labels = []
                    coords = []
                    forces = []
                    while index < len(data):
                        label = data[index]
                        count = int(data[index + 2])
                        index += 3
                        for atom in data[index:index + count]:
                            fields = atom.split()
                            labels.append(label)
                            coords.append([float(value) for value in fields[:3]])
                            if 'f' in fields:
                                force_index = fields.index('f')
                                forces.append([float(value) for value in fields[force_index + 1:force_index + 4]])
                        index += count
                    frames[step] = {'step': step, 'coordinate': coordinate,
                                    'cell': cell, 'elem': labels, 'coords': np.array(coords),
                                    'forces': np.array(forces) if len(forces) == len(labels) else None}
                if raw != '# MD step END':
                    step = int(raw.split()[3])
                    lines = []
                else:
                    step = None
            elif step is not None:
                lines.append(raw)
    return [frames[step] for step in sorted(frames)]


def read_traj_from_md_cif(directory: Path,
                          running_log_lines: Optional[List[str]] = None) -> List[Dict[str, np.ndarray]]:
    """Read single-frame MD CIF files, keeping the cumulative step markers."""
    numbered = sorted((path for path in (directory / 'STRUC').glob('STRU[0-9]*.cif')
                       if re.fullmatch(r'STRU[0-9]+\.cif', path.name)),
                      key=lambda path: int(path.stem[4:]))
    periodic = directory / 'STRU.cif'
    sources = [periodic] if periodic.is_file() else numbered
    final = directory / 'STRU_FINAL.cif'
    if not final.is_file():
        final = directory / 'STRUC' / 'STRU_FINAL.cif'
    if final.is_file():
        sources.append(final)
    if running_log_lines is not None:
        announced = []
        for line in running_log_lines:
            match = re.search(r'Write file:\s+(.+?)\s*$', line)
            if match is None:
                continue
            path = Path(match.group(1))
            if path.name in ('STRU.cif', 'STRU_FINAL.cif') or re.fullmatch(r'STRU[0-9]+\.cif', path.name):
                path = directory / 'STRUC' / path.name if path.parent.name == 'STRUC' else directory / path.name
                if path.is_file() and path not in announced:
                    announced.append(path)
        if announced:
            sources = announced
    frames = {}
    for source in sources:
        marker = re.match(r'# MD step (\d+)', source.read_text().splitlines()[0])
        if marker is None:
            continue
        atoms = read_cif_atoms(source)
        step = int(marker.group(1))
        frames[step] = {'step': step, 'coordinate': 'Cartesian_angstrom',
                        'cell': np.array(atoms.cell), 'elem': atoms.get_chemical_symbols(),
                        'coords': atoms.get_positions(), 'forces': None}
    return [frames[step] for step in sorted(frames)]


def find_eig_occ_files(directory: Path,
                       running_log_lines: Optional[List[str]] = None) -> List[Path]:
    '''Discover eigenvalue/occupation files in the output directory.

    ABACUS writes eigenvalue data in three layouts depending on
    ``out_app_flag``:

    - mode 0: one fixed ``eig_occ.txt`` holding only the latest written frame;
    - mode 1: one fixed ``eig_occ.txt`` holding every written frame;
    - mode 2: one ``eig_occg{N}.txt`` per written ionic step.

    When ``running_log_lines`` is given, only the files announced by the
    current run are returned, so stale files left in a reused output
    directory by an earlier run are ignored.  An announced file that is
    missing is dropped instead of falling back to an unannounced file: the
    log is authoritative for the layout of the current run, and an old
    ``eig_occ.txt`` must not be reported as the current results.  If the log
    carries no such announcement (older builds or hand-written fixtures),
    the filename based discovery below is used instead.

    Parameters
    ----------
    directory : Path
        The ``OUT.${suffix}`` directory to search.
    running_log_lines : list of str, optional
        Lines of the running log of the current run.

    Returns
    -------
    list of Path
        Sorted list of eigenvalue files.  Empty when none are found.
    '''
    if running_log_lines is not None:
        announced: List[Path] = []
        for line in running_log_lines:
            match = re.search(r'Write\s+Eigenvalues and occupations to file:\s*(\S+)',
                              line, re.IGNORECASE)
            if match is None:
                continue
            candidate = directory / Path(match.group(1)).name
            if candidate not in announced:
                announced.append(candidate)
        if announced:
            def _step_number(path: Path) -> int:
                match = re.search(r'g(\d+)\.txt$', path.name)
                return int(match.group(1)) if match else 0
            # The announcement set decides the layout. A missing announced
            # file is reported as missing (by its absence) rather than being
            # replaced by a stale unannounced file from an earlier run.
            existing = [path for path in announced if path.exists()]
            return sorted(existing, key=_step_number)

    single = directory / 'eig_occ.txt'
    if single.exists():
        return [single]
    numbered = sorted(directory.glob('eig_occg*.txt'),
                      key=lambda p: int(re.search(r'g(\d+)\.txt$', p.name).group(1)))
    return numbered


def read_band_from_eig_occ_files(files: List[Path]) -> List[Dict[str, np.ndarray]]:
    '''Read eigenvalue/occupation data from one or more eig_occ files.

    When multiple numbered files are given (mode 2), each file contributes
    one frame.  When a single file is given (mode 0/1), the frames are
    parsed from its internal ``# Ionic_Step`` markers.

    Parameters
    ----------
    files : list of Path
        Files returned by :func:`find_eig_occ_files`.

    Returns
    -------
    list of dict
        A list of dictionaries containing the k-points and band energies.
        Each dictionary has keys 'k', 'e', and 'occ'.
    '''
    if not files:
        return []
    if len(files) == 1:
        return read_band_from_eig_occ(files[0])

    # Mode 2: concatenate the per-step files and parse them together.
    raw: List[str] = []
    for fn in files:
        with open(fn) as f:
            raw.extend(f.readlines())
    return read_band_from_eig_occ(raw)

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

    # The current writer emits the one-based step after the marker
    # ("# Ionic_Step N"), while older outputs wrote it before
    # ("N     # Ionic_Step", and "N     # ionic step" before the keyword was
    # unified).  All three spellings are accepted case-insensitively.
    step_header_pat = r'^(?:\s*#\s*Ionic[ _]Step\s+(\d+)\s*' \
                    + r'|\s*(\d+)\s*#\s*Ionic[ _]Step\s*)$'
    step_matches = [re.match(step_header_pat, l, re.IGNORECASE) for l in raw]
    steps = [int(m.group(1) or m.group(2)) for m in step_matches if m]
    nframe = len(steps)
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

    return [{'step': si, 'k': ki, 'e': eki, 'occ': occi}
            for si, ki, eki, occi in zip(steps, k, ekb, occ)]


def align_eig_frames(elecstate: List[Dict[str, np.ndarray]],
                     step_ids: List[int]) \
        -> List[Optional[Dict[str, np.ndarray]]]:
    '''Align written eigenvalue frames with the ionic steps of a trajectory.

    ABACUS may write fewer eigenvalue frames than there are ionic steps, for
    example when ``out_freq_ion > 1`` skips steps or when ``out_app_flag=0``
    keeps only the final snapshot.  Frames are matched by their recorded
    ``Ionic_Step`` ID; a step without a written frame yields ``None`` instead
    of borrowing the spectrum of another step.

    Parameters
    ----------
    elecstate : list of dict
        Frames returned by :func:`read_band_from_eig_occ_files`.
    step_ids : list of int
        One-based ionic step of every trajectory frame, in order.

    Returns
    -------
    list of dict or None
        One entry per element of ``step_ids``; ``None`` where no eigenvalue
        frame was written.
    '''
    by_step = {frame['step']: frame for frame in elecstate if 'step' in frame}
    return [by_step.get(step) for step in step_ids]

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


def read_md_virial_by_step(lines: List[str], input_info: Path,
                           esolver_type: str) -> Dict[int, np.ndarray]:
    """Reconstruct the former electronic stress from MD virial output."""
    pressure = np.zeros(3)
    if input_info.is_file() and not esolver_type.startswith(('ofdft', 'sdft')):
        for line in input_info.read_text().splitlines():
            match = re.match(r'\s*press([123])\s+([-+\d.eE]+)', line)
            if match:
                pressure[int(match.group(1)) - 1] = float(match.group(2))
    stresses = {}
    step = None
    for index, line in enumerate(lines):
        marker = re.match(r'\s*MD STEP:\s*(\d+)', line)
        if marker:
            step = int(marker.group(1))
        if step is not None and re.match(r'\s*VIRIAL \(kbar\):', line):
            rows = [list(map(float, row.split())) for row in lines[index + 1:index + 4]]
            if len(rows) != 3 or any(len(row) != 3 for row in rows):
                raise ValueError(f'Malformed MD virial at step {step}')
            tensor = np.array(rows)
            tensor[np.diag_indices(3)] += pressure
            stresses[step] = -0.1 * tensor * GPa
    return stresses


def md_log_step_offset(lines: List[str]) -> Optional[int]:
    """Map local solver ion-move IDs to cumulative MD structure step IDs."""
    global_step = None
    for line in lines:
        marker = re.match(r'\s*STEP OF MOLECULAR DYNAMICS:\s*(\d+)', line)
        if marker:
            global_step = int(marker.group(1))
        local = re.search(r'-->\s+#ION MOVE#\s+(\d+)', line)
        if local and global_step is not None:
            return global_step - int(local.group(1))
    return None

def read_iter_header_from_running_log(src: str | Path | List[str]) \
    -> List[Tuple[int, int]]:
    '''
    read the "iteration header" from the running log, useful for determining
    which is the final iteration. The "iteration header" is defined as:
    ```
    --> #ION MOVE#         1  #ELEC ITER#         3
    ```
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

    HEADER_PAT = r'-->\s+#ION MOVE#\s+(\d+)\s+#ELEC ITER#\s+(\d+)'
    res = [re.findall(HEADER_PAT, l) for l in raw]

    return [(int(pack[0][0]), int(pack[0][1])) for pack in res if pack]

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
    if fileobj.name == 'running_md.log':
        restart_markers = [index for index, line in enumerate(abacus_lines)
                           if re.search(r'Read file: .*[/\\]STRU(?:_restart|_MD_\d+|\d+)?\s*$', line)]
        if restart_markers:
            abacus_lines = abacus_lines[restart_markers[-1]:]
    
    # read the esolver type
    eslvtyp = read_esolver_type_from_running_log(abacus_lines)
    # FIXME: implement read_ksdft_esolver_out instead of read_abacus_out to
    # make flexible esolver support

    # New MD output uses STRU frames. Keep the legacy reader for existing
    # calculation directories written by older ABACUS versions.
    announced_structure = any(re.search(r'Write file: .*[/\\]STRU(?:_FINAL|[0-9]+)?(?:\.cif)?\s*$', line)
                              for line in abacus_lines)
    has_stru_files = ((fileobj.parent / 'STRU').is_file()
                      or (fileobj.parent / 'STRU' / 'STRU').is_file()
                      or (fileobj.parent / 'STRU_FINAL').is_file()
                      or (fileobj.parent / 'STRUC' / 'STRU_FINAL').is_file()
                      or (fileobj.parent / 'STRU.cif').is_file()
                      or (fileobj.parent / 'STRU_FINAL.cif').is_file()
                      or (fileobj.parent / 'STRUC' / 'STRU_FINAL.cif').is_file()
                      or any((fileobj.parent / 'STRUC').glob('STRU[0-9]*')))
    new_md_output = (fileobj.name == 'running_md.log' and
                     (announced_structure or (has_stru_files and not (fileobj.parent / 'MD_dump').is_file())))
    if new_md_output:
        cif_output = any(re.search(r'Write file: .*[/\\]STRU(?:_FINAL|[0-9]+)?\.cif\s*$', line)
                         for line in abacus_lines)
        if cif_output or (fileobj.parent / 'STRU_FINAL.cif').is_file():
            trajectory = read_traj_from_md_cif(fileobj.parent, abacus_lines)
        else:
            trajectory = read_traj_from_md_stru(fileobj.parent, abacus_lines)
    elif fileobj.name == 'running_md.log':
        trajectory = read_traj_from_md_dump(fileobj.parent / 'MD_dump')
    else:
        trajectory = read_traj_from_running_log(abacus_lines)

    # read the eigenvalues (nframe, nk, nbnd) written by the current run
    eig_files = find_eig_occ_files(fileobj.parent, abacus_lines)
    elecstate = read_band_from_eig_occ_files(eig_files)

    # Associate eigenvalue frames with ionic steps.  The running log carries
    # one "#ION MOVE#" header per ionic step; if it is unavailable or does not
    # match the trajectory, fall back to the one-based step order.  Steps
    # without a written eigenvalue frame stay unavailable instead of reusing
    # another step's spectrum.
    iter_header = read_iter_header_from_running_log(abacus_lines)
    step_ids: List[int] = []
    for ion_move, _ in iter_header:
        if ion_move not in step_ids:
            step_ids.append(ion_move)
    if new_md_output:
        step_offset = md_log_step_offset(abacus_lines)
        step_ids = [frame['step'] for frame in trajectory]
        if step_offset is None:
            local_steps = [ion_move for ion_move, _ in iter_header]
            step_offset = max(step_ids) - max(local_steps) if step_ids and local_steps else 0
        for frame in elecstate:
            if 'step' in frame:
                frame['step'] += step_offset
    elif len(step_ids) != len(trajectory):
        step_ids = list(range(1, len(trajectory) + 1))
    elecstate = align_eig_frames(elecstate, step_ids)

    # read the atomic forces (nframe, nat, 3)
    forces = read_forces_from_running_log(abacus_lines)

    # read the stress (nframe, 3, 3), but may be None
    stress = read_stress_from_running_log(abacus_lines) if not new_md_output else []

    # read all kpoints tables (but only want the first, spinless, with ibz2bz)
    k, _, _, _, _ = read_kpoints_from_running_log(abacus_lines)
    # unpack the kpoints information
    kvecd, wk, _ = k

    # FIXME: in principle, the two spin channels share the same set of
    # the kpoints, so it is not needed to have two sets of kpoints
    # and it is assumed that the sampling of kpoints won't change during
    # the simulation, which is, not exactly to be true for the NPT-MD
    # runs.

    # only keep the energies of the final iteration for each ion step
    energies = find_final_info_with_iter_header(
        read_energies_from_running_log(abacus_lines)[1],
        iter_header
    )
    magmom = read_magmom_from_running_log(abacus_lines)

    if new_md_output:
        log_steps = []
        for ion_move, _ in iter_header:
            global_step = ion_move + step_offset
            if global_step not in log_steps:
                log_steps.append(global_step)
        if len(log_steps) != len(energies):
            log_steps = list(range(step_offset + 1, step_offset + len(energies) + 1))
        selected = [log_steps.index(step) for step in step_ids]
        energies = [energies[i] for i in selected]
        if len(forces) == len(log_steps):
            forces = [forces[i] for i in selected]
        stress_by_step = read_md_virial_by_step(abacus_lines, fileobj.parent / 'INPUT.info', eslvtyp)
        stress = [stress_by_step.get(step) for step in step_ids]
        if len(magmom) == len(log_steps):
            magmom = [magmom[i] for i in selected]
        if all(frame['forces'] is not None for frame in trajectory):
            forces = [frame['forces'] for frame in trajectory]

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
                                        magmoms=mag[ind], efermi=ener['E_Fermi'],
                                        ibzkpts=kvecd, dipole=None)
        # import the eigenvalues and occupations kpoint-by-kpoint; a step
        # without a written eigenvalue frame keeps an empty k-point list
        calc.kpts = []
        if estat is not None:
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

    def test_read_md_stru_frames_and_force(self):
        import tempfile

        def frame(step, position, force):
            return (f'# MD step {step}\n\nATOMIC_SPECIES\nH 1.0 auto auto\n'
                    'LATTICE_CONSTANT\n1.8897270000\nLATTICE_VECTORS\n'
                    '4 0 0\n0 4 0\n0 0 4\nATOMIC_POSITIONS\nCartesian_angstrom\n'
                    f'H\n0.0\n1\n{position} 0 0 m 1 1 1 v 0.1 0 0 f {force} 0 0\n')

        with tempfile.TemporaryDirectory() as scratch:
            directory = Path(scratch)
            (directory / 'STRU').write_text(frame(1, 1.0, 2.0) + frame(2, 1.5, 3.0))
            (directory / 'STRU_FINAL').write_text(frame(2, 1.5, 3.0))
            trajectory = read_traj_from_md_stru(directory)
            self.assertEqual([item['step'] for item in trajectory], [1, 2])
            self.assertTrue(np.allclose(trajectory[1]['coords'], [[1.5, 0, 0]]))
            self.assertTrue(np.allclose(trajectory[1]['forces'], [[3.0, 0, 0]]))
            self.assertTrue(np.allclose(trajectory[1]['cell'], np.eye(3) * 4.0, atol=1e-5))

            (directory / 'STRU').unlink()
            numbered = directory / 'STRUC'
            numbered.mkdir()
            (numbered / 'STRU1').write_text(frame(1, 1.0, 2.0))
            (numbered / 'STRU3').write_text(frame(3, 1.8, 4.0))
            (numbered / 'STRU2.cif').write_text('stale CIF output')
            (directory / 'STRU_FINAL').unlink()
            (numbered / 'STRU_FINAL').write_text(frame(4, 2.0, 5.0))
            self.assertEqual([item['step'] for item in read_traj_from_md_stru(directory)], [1, 3, 4])
            (directory / 'STRU').write_text(frame(99, 9.0, 9.0))
            announcements = [f'Write file: {numbered / name}'
                             for name in ('STRU1', 'STRU3', 'STRU_FINAL')]
            self.assertEqual([item['step'] for item in read_traj_from_md_stru(directory, announcements)],
                             [1, 3, 4])
            (directory / 'STRU').unlink()
            (numbered / 'STRU1').unlink()
            (numbered / 'STRU3').unlink()
            self.assertEqual([item['step'] for item in read_traj_from_md_stru(directory)], [4])

    def test_read_abacus_out_from_new_md_stru(self):
        import tempfile

        with tempfile.TemporaryDirectory() as scratch:
            directory = Path(scratch)
            shutil.copy(self.testfiles / 'pw-symm0-nspin4-gamma-md', directory / 'running_md.log')
            shutil.copy(self.testfiles / 'nspin4-gamma-eigocc', directory / 'eig_occ.txt')
            header = ('ATOMIC_SPECIES\nAs 1 auto auto\nGa 1 auto auto\n'
                      'LATTICE_CONSTANT\n1.8897270000\nLATTICE_VECTORS\n'
                      '4 0 0\n0 4 0\n0 0 4\nATOMIC_POSITIONS\nCartesian_angstrom\n')
            def frame(step, x, force):
                return (f'# MD step {step}\n' + header
                        + f'As\n0.0\n1\n1 1 1 m 1 1 1 v 0 0 0 f {force} 0 0\n'
                        + f'Ga\n0.0\n1\n{x} 0 0 m 1 1 1 v 0 0 0 f {-force} 0 0\n')
            (directory / 'STRU').write_text(frame(1, 0.0, 1.0) + frame(2, 0.1, 2.0))
            (directory / 'STRU_FINAL').write_text(frame(2, 0.1, 2.0))
            images = read_abacus_out(directory / 'running_md.log')
            (directory / 'STRU').unlink()
            final_only = read_abacus_out(directory / 'running_md.log')
        self.assertEqual(len(images), 2)
        self.assertTrue(np.allclose(images[-1].get_forces(), [[2.0, 0, 0], [-2.0, 0, 0]]))
        self.assertEqual(len(final_only), 1)
        self.assertTrue(np.allclose(final_only[0].get_forces(), [[2.0, 0, 0], [-2.0, 0, 0]]))

    def test_restarted_md_stru_aligns_energy_and_virial(self):
        import tempfile

        with tempfile.TemporaryDirectory() as scratch:
            directory = Path(scratch)
            log = (self.testfiles / 'pw-symm0-nspin4-gamma-md').read_text()
            log = log.replace('STEP OF MOLECULAR DYNAMICS: 1',
                              'STEP OF MOLECULAR DYNAMICS: 6\n MD STEP: 6\n'
                              ' VIRIAL (kbar):\n 1 0 0\n 0 2 0\n 0 0 3')
            log = log.replace('STEP OF MOLECULAR DYNAMICS: 2',
                              'STEP OF MOLECULAR DYNAMICS: 7\n MD STEP: 7\n'
                              ' VIRIAL (kbar):\n 4 0 0\n 0 5 0\n 0 0 6')
            (directory / 'running_md.log').write_text(log)
            (directory / 'INPUT.info').write_text('press1 10\npress2 20\npress3 30\n')
            shutil.copy(self.testfiles / 'nspin4-gamma-eigocc', directory / 'eig_occ.txt')
            header = ('ATOMIC_SPECIES\nAs 1 auto auto\nGa 1 auto auto\n'
                      'LATTICE_CONSTANT\n1.8897270000\nLATTICE_VECTORS\n'
                      '4 0 0\n0 4 0\n0 0 4\nATOMIC_POSITIONS\nCartesian_angstrom\n')
            def frame(step, force):
                return (f'# MD step {step}\n' + header
                        + f'As\n0.0\n1\n1 1 1 m 1 1 1 v 0 0 0 f {force} 0 0\n'
                        + f'Ga\n0.0\n1\n0 0 0 m 1 1 1 v 0 0 0 f {-force} 0 0\n')
            (directory / 'STRU').write_text(frame(6, 1.0) + frame(7, 2.0))
            (directory / 'STRU_FINAL').write_text(frame(7, 2.0))
            images = read_abacus_out(directory / 'running_md.log')
        self.assertEqual(len(images), 2)
        self.assertTrue(np.allclose(images[0].get_forces(), [[1, 0, 0], [-1, 0, 0]]))
        self.assertTrue(np.allclose(images[1].get_forces(), [[2, 0, 0], [-2, 0, 0]]))
        expected = -0.1 * GPa * np.diag([14, 25, 36])
        self.assertTrue(np.allclose(images[1].get_stress(voigt=False), expected))

    def test_md_cif_and_stale_numbered_cif(self):
        import tempfile
        from ase.io import write as write_cif

        with tempfile.TemporaryDirectory() as scratch:
            directory = Path(scratch)
            numbered = directory / 'STRUC'
            numbered.mkdir()
            (numbered / 'STRU1.cif').write_text('stale CIF')
            final = directory / 'STRU_FINAL.cif'
            write_cif(final, Atoms('He', positions=[[0, 0, 0]], cell=np.eye(3) * 4,
                                   pbc=True), format='cif')
            final.write_text('# MD step 2\n' + final.read_text())
            self.assertEqual(read_traj_from_md_stru(directory), [])
            announcements = [f'Write file: {directory / "STRU_FINAL.cif"}']
            frames = read_traj_from_md_cif(directory, announcements)
            self.assertEqual([frame['step'] for frame in frames], [2])
            self.assertEqual(frames[0]['elem'], ['He'])

    def test_read_abacus_out_reorders_calculator_magmoms(self):
        import tempfile
        from unittest.mock import patch

        frame = {
            'elem': ['Na', 'Na', 'Cl'],
            'coords': np.array([[0.0, 0.0, 0.0],
                                [1.0, 0.0, 0.0],
                                [2.0, 0.0, 0.0]]),
            'cell': np.eye(3),
        }
        elecstate = [{
            'k': np.zeros((1, 1, 3)),
            'e': np.zeros((1, 1, 1)),
            'occ': np.ones((1, 1, 1)),
        }]
        energies = [{'E_KohnSham': -1.0, 'E_Fermi': 0.0}]
        kpoints = ((np.zeros((1, 3)), np.ones(1), None), None, None, None, None)

        with tempfile.TemporaryDirectory() as tmpdir:
            running_log = Path(tmpdir) / 'running_scf.log'
            running_log.write_text('')
            (Path(tmpdir) / 'eig_occ.txt').write_text('')
            with patch(__name__ + '.read_esolver_type_from_running_log', return_value='ksdft'), \
                 patch(__name__ + '.read_traj_from_running_log', return_value=[frame]), \
                 patch(__name__ + '.find_eig_occ_files', return_value=[Path(tmpdir) / 'eig_occ.txt']), \
                 patch(__name__ + '.read_band_from_eig_occ_files', return_value=elecstate), \
                 patch(__name__ + '.read_forces_from_running_log', return_value=[]), \
                 patch(__name__ + '.read_stress_from_running_log', return_value=[]), \
                 patch(__name__ + '.read_kpoints_from_running_log', return_value=kpoints), \
                 patch(__name__ + '.read_energies_from_running_log', return_value=([], [])), \
                 patch(__name__ + '.read_iter_header_from_running_log', return_value=[]), \
                 patch(__name__ + '.find_final_info_with_iter_header', return_value=energies), \
                 patch(__name__ + '.read_magmom_from_running_log', return_value=[np.array([10.0, 20.0, 30.0])]):
                atoms = read_abacus_out(running_log, sort_atoms_with=[0, 2, 1])[0]

        self.assertEqual(atoms.get_chemical_symbols(), ['Na', 'Cl', 'Na'])
        self.assertTrue(np.allclose(atoms.calc.results['magmoms'], [10.0, 30.0, 20.0]))

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
        self.assertEqual([es['step'] for es in elecstate], [1, 2])
        for es in elecstate: # loop over frames
            self.assertIn('k', es)
            self.assertEqual(es['k'].shape, (1, 1, 3)) # ispin, ik, 3
            self.assertIn('e', es)
            self.assertEqual(es['e'].shape, (1, 1, 35)) # ispin, ik, nbnd
            self.assertIn('occ', es)
            self.assertEqual(es['occ'].shape, (1, 1, 35)) # ispin, ik, nbnd

    def test_read_band_from_eig_occ_producer_header(self):
        # A single snapshot from the current writer: the step follows the
        # marker and two provenance lines precede it.
        text = (
            '# ABACUS version: 3.10.0\n'
            '# Written at 2026-11-09 15:18:58\n'
            '# Ionic_Step 3\n'
            ' Electronic state energy (eV) and occupations\n'
            ' Spin number 1\n'
            ' spin=1 k-point=1/1 Cartesian=0.0000000 0.0000000 0.0000000 (7 plane wave)\n'
            ' 1 -8.734940678442698 1.000000000000000\n'
            ' 2 7.011672764486192 0.000000000000000\n'
        )
        elecstate = read_band_from_eig_occ(text.splitlines(keepends=True))
        self.assertEqual(len(elecstate), 1)
        self.assertEqual(elecstate[0]['step'], 3)
        self.assertEqual(elecstate[0]['e'].shape, (1, 1, 2))

    def test_read_band_from_eig_occ_legacy_headers(self):
        # Older spellings put the step before the marker; the keyword may be
        # lower case.  Both must still be parsed into the same frame layout.
        body = (
            ' spin=1 k-point=1/1 Cartesian=0.0000000 0.0000000 0.0000000 (7 plane wave)\n'
            ' 1 -8.734940678442698 1.000000000000000\n'
            ' 2 7.011672764486192 0.000000000000000\n'
        )
        for header in ['1     # Ionic_Step', '1     # ionic step']:
            elecstate = read_band_from_eig_occ((header + '\n' + body).splitlines(keepends=True))
            self.assertEqual(len(elecstate), 1)
            self.assertEqual(elecstate[0]['step'], 1)
            self.assertEqual(elecstate[0]['e'].shape, (1, 1, 2))

    def test_find_eig_occ_files(self):
        import tempfile

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            # no files
            self.assertEqual(find_eig_occ_files(tmp), [])
            # single file (mode 0/1)
            (tmp / 'eig_occ.txt').write_text('')
            self.assertEqual(find_eig_occ_files(tmp), [tmp / 'eig_occ.txt'])
            # numbered files (mode 2) take precedence when single is absent
            (tmp / 'eig_occ.txt').unlink()
            (tmp / 'eig_occg2.txt').write_text('')
            (tmp / 'eig_occg1.txt').write_text('')
            (tmp / 'eig_occg3.txt').write_text('')
            self.assertEqual(find_eig_occ_files(tmp),
                             [tmp / 'eig_occg1.txt',
                              tmp / 'eig_occg2.txt',
                              tmp / 'eig_occg3.txt'])

    def test_find_eig_occ_files_prefers_current_run(self):
        import tempfile

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            # A stale fixed file from an earlier run coexists with the files
            # written by the current mode-2 run.
            (tmp / 'eig_occ.txt').write_text('stale\n')
            (tmp / 'eig_occg1.txt').write_text('')
            (tmp / 'eig_occg2.txt').write_text('')
            log = [
                ' Write Eigenvalues and occupations to file: OUT.ABACUS/eig_occg1.txt\n',
                ' Write Eigenvalues and occupations to file: OUT.ABACUS/eig_occg2.txt\n',
            ]
            self.assertEqual(find_eig_occ_files(tmp, log),
                             [tmp / 'eig_occg1.txt', tmp / 'eig_occg2.txt'])
            # Without a log announcement the legacy discovery still applies.
            self.assertEqual(find_eig_occ_files(tmp), [tmp / 'eig_occ.txt'])

    def test_align_eig_frames_keeps_steps_distinct(self):
        frames = [
            {'step': 1, 'k': None, 'e': None, 'occ': None},
            {'step': 3, 'k': None, 'e': None, 'occ': None},
        ]
        aligned = align_eig_frames(frames, [1, 2, 3])
        self.assertIs(aligned[0], frames[0])
        self.assertIsNone(aligned[1])
        self.assertIs(aligned[2], frames[1])
        # mode 0 keeps only the final snapshot; earlier steps stay unavailable
        final_only = align_eig_frames([{'step': 3, 'k': None, 'e': None, 'occ': None}],
                                      [1, 2, 3])
        self.assertEqual([f is None for f in final_only], [True, True, False])

    def test_read_band_from_eig_occ_files(self):
        fn = self.testfiles / 'nspin4-gamma-eigocc'
        # single file path
        elecstate = read_band_from_eig_occ_files([fn])
        self.assertEqual(len(elecstate), 2)
        # empty list
        self.assertEqual(read_band_from_eig_occ_files([]), [])
        # multiple files: split the single fixture into two step files
        import tempfile
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            with open(fn) as f:
                lines = f.readlines()
            # find the frame boundary
            step_idx = [i for i, l in enumerate(lines) if '# Ionic_Step' in l]
            self.assertEqual(len(step_idx), 2)
            file1 = tmp / 'eig_occg1.txt'
            file2 = tmp / 'eig_occg2.txt'
            with open(file1, 'w') as f:
                f.writelines(lines[:step_idx[1]])
            with open(file2, 'w') as f:
                f.writelines(lines[step_idx[1]:])
            elecstate = read_band_from_eig_occ_files([file1, file2])
            self.assertEqual(len(elecstate), 2)

    def test_read_band_from_eig_occ_accepts_legacy_step_marker(self):
        # The reader promises to accept both the current "# Ionic_Step N"
        # marker and the legacy "N     # ionic step" spelling. Run the same
        # fixture through both and compare the parsed frames.
        fn = self.testfiles / 'nspin4-gamma-eigocc'
        current = read_band_from_eig_occ_files([fn])

        import tempfile
        with tempfile.TemporaryDirectory() as tmpdir:
            legacy = Path(tmpdir) / 'eig_occ.txt'
            lines = []
            for line in fn.read_text(encoding='utf-8').splitlines(keepends=True):
                match = re.match(r'^\s*#\s*Ionic[ _]Step\s+(\d+)\s*$', line,
                                 re.IGNORECASE)
                if match:
                    lines.append('{}     # ionic step\n'.format(int(match.group(1))))
                else:
                    lines.append(line)
            legacy.write_text(''.join(lines), encoding='utf-8')
            parsed = read_band_from_eig_occ_files([legacy])

        self.assertEqual(len(parsed), len(current))
        for frame, reference in zip(parsed, current):
            self.assertEqual(frame['step'], reference['step'])
            np.testing.assert_allclose(frame['k'], reference['k'])
            np.testing.assert_allclose(frame['e'], reference['e'])
            np.testing.assert_allclose(frame['occ'], reference['occ'])

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
        (self.testfiles / 'eig_occ.txt').unlink()
        (self.testfiles / 'MD_dump').unlink()

    def test_read_iter_header_from_running_log(self):
        header = read_iter_header_from_running_log(
            self.testfiles / 'lcao-symm0-nspin2-multik-cellrelax')
        self.assertTupleEqual(header[0], (1, 1))
        self.assertTupleEqual(header[1], (1, 2))
        self.assertTupleEqual(header[2], (2, 1))

if __name__ == '__main__':
    unittest.main()
