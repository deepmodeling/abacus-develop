'''this module provides some useful tools for k-sampling'''
import seekpath
import unittest
from typing import Optional, Tuple, List, Dict

import numpy as np
from ase.atoms import Atoms

def convert_kspacing_to_kpts(cell: np.ndarray,
                             kspacing: float | Tuple[float, float, float])\
     -> Tuple[int, int, int]:
    '''
    convert the kspacing to kpts, according to the given cell in Angstrom.
    This function would be helpful for those calculators that only support
    kpts, but not kspacing. A negative value in kspacing would yield nk=1
    for that direction.

    Parameters
    ----------
    cell : np.ndarray
        The cell in Angstrom, in shape (3, 3).
    kspacing : float or tuple of float
        The kspacing in Angstrom. If a tuple of float is given, it would be
        interpreted as the kspacing in each direction.
    
    Returns
    -------
    tuple of int
        The kpts in each direction.
    '''
    # get reciprocal cell vectors...
    bvec = 2*np.pi * np.linalg.solve(cell.T, np.eye(3))
    bnorm = np.linalg.norm(bvec, axis=1).tolist()

    kspacing = (kspacing, kspacing, kspacing) if isinstance(kspacing, float) else kspacing
    assert len(kspacing) == 3
    nk = [int(norm / kspac) if kspac > 0 else 1 for kspac, norm in zip(kspacing, bnorm)]
    return tuple(map(lambda x: max(1, x+1), nk))

def get_kpath(atoms: Atoms, n_interpl: int = 10) -> Tuple[np.ndarray, np.ndarray, List[str]]:
    ''''''
    from seekpath import get_path as _corefunc_seekpath
    seeked = _corefunc_seekpath(
        structure=(atoms.get_cell().tolist(), 
                   atoms.get_scaled_positions(), 
                   atoms.get_atomic_numbers()))
    
    # get the most useful informaiton
    k = np.array(list(seeked['point_coords'].values()))
    assert k.ndim == 2
    nk, dim = k.shape
    assert dim == 3

    # interpolate...
    kvec = np.vstack([np.linspace(k[i, :], k[i + 1, :], 
                                  num=n_interpl+int(i == nk - 2), 
                                  endpoint=bool(i == nk - 2)) 
                     for i in range(nk - 1)])
    
    # for easy-drawing
    kdist = np.linalg.norm(kvec[1:] - kvec[:-1], axis=1)
    kdist = np.cumsum(kdist)
    kdist = np.hstack([[0], kdist])

    kname = list(seeked['point_coords'].keys())
    kname = [[kn] + [''] * (n_interpl-1) for kn in kname[:-1]] + [[kname[-1]]]
    kname = [n for kn in kname for n in kn]

    # finally, we assert the consistency of the data
    assert len(kvec) == len(kdist) == len(kname), f'{len(kvec)}, {len(kdist)}, {len(kname)}'
    return kvec, kdist, kname

def interpolate_kpath(knodes: np.ndarray, 
                      knames: List[str],
                      n_interpl: int | List[int] = 10) -> Tuple[np.ndarray, np.ndarray, List[str]]:
    '''
    interpolate between the given k-points (nodes, as they are nodes along the path).
    The number of interpolated points is given by `n_interpl`, in which the last
    number will be ignored.
    '''
    if isinstance(n_interpl, int):
        n_interpl = [n_interpl] * (len(knodes) - 1) + [1]
    assert isinstance(knodes, np.ndarray)
    assert knodes.ndim == 2
    nk, dim = knodes.shape

    assert dim == 3
    assert len(n_interpl) == len(knames) == nk

    kvec = np.vstack([np.linspace(knodes[i, :], knodes[i + 1, :], 
                                  num=n_interpl[i]+int(i == nk - 2), 
                                  endpoint=bool(i == nk - 2)) 
                     for i in range(nk - 1)])
    
    # for easy-drawing
    kdist = np.linalg.norm(kvec[1:] - kvec[:-1], axis=1)
    kdist = np.cumsum(kdist)
    kdist = np.hstack([[0], kdist])

    # padding the knames to have the same length as kvec
    kname = [[kn] + [''] * (n_interpl[i]-1) for i, kn in enumerate(knames[:-1])] + [[knames[-1]]]
    kname = [n for kn in kname for n in kn]

    # finally, we assert the consistency of the data
    assert len(kvec) == len(kdist) == len(kname), f'{len(kvec)}, {len(kdist)}, {len(kname)}'
    return kvec, kdist, kname

def merge_ksgm(segments) -> Tuple[List[str], List[bool]]:
    '''
    SeeKpath generates the path of kpoints in the format of list of 2-element tuples.
    This is good but not quite compatible with some DFT softwares like ABACUS. This
    function will return a list of special kpoint labels accompanied with a boolean
    list that helps distinguish the breakpoint of path.
    
    Parameters
    ----------
    segments : List[Tuple[str, str]]
        The list of kpath segments. Each segment is a tuple of two strings,
        where the first string is the starting point of the segment and the
        second string is the ending point of the segment, like [(Gamma, M),
        (M, K), (K, Gamma)]

    Returns
    -------
    Tuple[List[str], List[bool]]
        The first element is the list of special kpoint labels. The second element
        is the boolean list that helps distinguish the breakpoint of path.
    '''
    klabels, is_brkpt = [segments[0][0]], [False]
    for i, (start, end) in enumerate(segments):
        if start != klabels[-1]:
            is_brkpt[-1] = True
            klabels.append(start)
            is_brkpt.append(False)
        
        klabels.append(end)
        is_brkpt.append(i == len(segments) - 1)
    
    return klabels, is_brkpt

def make_kstring(klabels: List[str],
                 is_brkpt: List[bool]) -> str:
    '''
    Make the kpoint path in the str. For path in which there are break point,
    will use the comma to distinguish

    Parameters
    ----------
    klabels : List[str]
        The list of special kpoint labels.
    is_brkpt : List[bool]
        The boolean list that helps distinguish the breakpoint of path.

    Returns
    -------
    str
        The kstring that ABACUS uses to represent the kpath.
    '''
    out = ''
    for lbl, brkpt in zip(klabels, is_brkpt):
        out += lbl
        if brkpt:
            out += ','
    return out[:-1] # remove the last comma

def make_klines(kpts, 
                is_brkpt, 
                n_interpl,
                klabels) -> List[Dict[str, np.ndarray | str | int]]:
    '''
    gather the information and output the list of dicts that organizes the
    information easy for exporting the ABAUCS KLINES file. Each dictionary
    contains the keys `coord`, `label`, `n`, where `coord` is the coordinate
    of kpoint, `label` is the label of kpoint, and `n` is the number of
    interpolation points between the next kpoint and itself (the next kpoint
    included).

    Parameters
    ----------
    kpts : np.ndarray
        The array of kpoints.
    is_brkpt : List[bool]
        The boolean list that helps distinguish the breakpoint of path.
    n_interpl : int
        The number of interpolation points between the kpoints.
    klabels : List[str]
        The list of special kpoint labels.

    Returns
    -------
    List[Dict[str, np.ndarray | str | int]]
        The list of dictionaries that organizes the information easy for
        exporting the ABAUCS KLINES file.
    '''
    def fspawnk(c: np.ndarray, n: int, label: str) -> Dict[str, np.ndarray | str | int]:
        return {'coord': c, 'label': label, 'n': n}
    
    return [fspawnk(c, 1 if is_brkpt[i] else n_interpl, klbl) 
            for i, (c, klbl) in enumerate(zip(kpts, klabels))]

def kpathgen(atoms: Atoms) -> Tuple[str, Dict[str, List[float]]]:
    '''
    Generate the k-path by default SeeK-Path flavor in the format
    compatible with ase bandpath module
    '''
    kpathseen = seekpath.get_path(
        structure=(np.array(atoms.get_cell()), 
                atoms.get_scaled_positions(), 
                atoms.get_atomic_numbers()),
        with_time_reversal=True
    )

    # convert the kpoint path to the format that is acceptable by ASE
    kpathstr, is_brkpt = merge_ksgm(kpathseen['path'])
    fklblfilter = lambda lbl: lbl if lbl != 'GAMMA' else 'G'
    kpathstr = make_kstring([fklblfilter(lbl) for lbl in kpathstr], is_brkpt)
    kpathstr = ''.join([k if k != 'GAMMA' else 'G' for k in kpathstr])
    # seekpath use 'GAMMA' as the gamma point symbol, while the ASE use 'G'
    kspecial = {k: v for k, v in kpathseen['point_coords'].items() if k != 'GAMMA'}
    kspecial['G'] = [0.0, 0.0, 0.0]
    return kpathstr, kspecial

# the following are for quick test on the silicon case
SILICON_KNODES = np.array([[0.0000000000,  0.0000000000,   0.0000000000],
                           [0.5000000000,  0.0000000000,   0.5000000000],
                           [0.6250000000,  0.2500000000,   0.6250000000],
                           [0.3750000000,  0.3750000000,   0.7500000000],
                           [0.0000000000,  0.0000000000,   0.0000000000],
                           [0.5000000000,  0.5000000000,   0.5000000000],
                           [0.5000000000,  0.2500000000,   0.7500000000],
                           [0.5000000000,  0.0000000000,   0.5000000000]])
SILICON_KNAMES = ["G", "X", "X/U", "K", "G", "L", "W", "X"]
SILICON_NINTPL = [50, 50, 1, 50, 50, 50, 50, 1] # from the DeePTB example :)

class TestKsampling(unittest.TestCase):
    
    def test_convert_kspacing(self):
        # test-reference: ABACUS 3.8.4 implementation
        from ase.geometry import cellpar_to_cell

        cell = cellpar_to_cell([4.22798145, 4.22798145, 4.22798145, 60, 60, 60])
        kpts = convert_kspacing_to_kpts(cell, 0.03*1.889725989)
        # kspacing 0.03 Bohr-1, multiply 1.889725989 to convert to Angstrom-1
        self.assertEqual(kpts, (33, 33, 33))

    def test_merge_ksgm(self):
        klables, is_brkpt = merge_ksgm([('GAMMA', 'X'), ('X', 'U'), ('K', 'GAMMA'), 
                                        ('GAMMA', 'L'), ('L', 'W'), ('W', 'X')])
        self.assertEqual(klables, ['GAMMA',  'X',  'U',   'K','GAMMA', 'L', 'W', 'X'])
        self.assertEqual(is_brkpt, [False, False, True, False,  False, False, False, True])

    def test_make_kstring(self):
        fklblfilter = lambda lbl: lbl if lbl != 'GAMMA' else 'G'
        klabels = ['GAMMA',  'X',  'U',   'K','GAMMA', 'L', 'W', 'X']
        klabels = [fklblfilter(lbl) for lbl in klabels]
        is_brkpt = [False, False, True, False,  False, False, False, True]
        self.assertEqual(make_kstring(klabels, is_brkpt), 'GXU,KGLWX')

if __name__ == '__main__':
    unittest.main()