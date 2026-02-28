'''public functions that can be used by both the LTS and the latest
version of ABACUS'''
import re
import os
import tempfile
import unittest
from pathlib import Path
from io import TextIOWrapper
from typing import Optional, Dict, List, Any

import numpy as np
from ase.atoms import Atoms
from ase.build import bulk
from ase.data import chemical_symbols, atomic_masses
ATOM_MASS = dict(zip(chemical_symbols, atomic_masses.tolist()))

def load_pseudo(pseudo_dir: str) -> Dict[str, str]:
    '''load the pseudopotential mapping from the pseudo_dir'''
    pseudo_dir = Path(pseudo_dir)
    UPFPAT_ = r'^([A-Z][a-z]?)([^a-zA-Z]+\S+)?\.(upf|UPF)$'
    # the file name is assumed to have the pattern:
    #   [element symbol][delimeter][other things...][.upf|UPF]
    return {re.match(UPFPAT_, f.name).group(1): f.name
            for f in pseudo_dir.glob('*')
            if re.match(UPFPAT_, f.name)}

def load_orbital(orbital_dir: str, efficiency: bool=True) -> Dict[str, str]:
    '''load the ABACUS numerical atomic orbital mapping from the
    orbital dir. There may be the case that there is not only one
    orbitals find for one element, in this case, the second parameter
    `efficiency` will take its effect: if true, use the "smaller" 
    one, otherwise, use the "larger" one
    
    Parameters
    ----------
    orbital_dir : str
        The directory containing the orbital basis files.
    efficiency : bool, optional
        Whether to use the "smaller" one, default is True.
    
    Returns
    -------
    Dict[str, str]
        The mapping from element symbol to orbital basis file name.
    '''
    # a small helper function
    def _basis_size(basis: List[str]) -> int:
        SPECTRUM = 'spdfghi'
        num, symb = zip(*[re.match(r'(\d+)([spdfghi])', b).groups() for b in basis])
        return sum((2*SPECTRUM.index(s) + 1) * int(n) for n, s in zip(num, symb))
    
    orbital_dir = Path(orbital_dir)
    ORBPAT_ = r'^([A-Z][a-z]?)_gga_(\d+)au_(\d+(\.\d+)?)Ry_(\w+)\.orb$'
    # a correct parse would yield:
    # elem, rcut, ecut, components = re.match(ORBPAT_, f.name)
    tosort = {}
    for f in orbital_dir.glob('*'):
        m = re.match(ORBPAT_, f.name)
        if m:
            elem, rcut, ecut, _, components = m.groups()
            components = re.findall(r'\d+[spdfghi]', components)
            tosort.setdefault(elem, []).append(
                (f.name, int(rcut), float(ecut), _basis_size(components))
            )
    # sort the basis files by their size: first by components, then by rcut, then by ecut
    ib = 0 if efficiency else -1
    return {elem: sorted(basis, key=lambda b: (b[3], b[1], b[2]))[ib][0]
            for elem, basis in tosort.items()}

def file_safe_backup(fn: Path, suffix: str = 'bak'):
    '''for the case where there are already files with the same name,
    add a suffix to the file name, like `STRU.bak.0`. If there are
    already `STRU.bak.0`, rename the elder to `STRU.bak.1` and let
    the latest one be `STRU.bak.0`.
    
    Parameters
    ----------
    fn : Path
        The path to the file to backup. Note: it must be provided
        as the Path object so that its folder is accessible by this
        function.
    suffix : str, optional
        The suffix to add to the file name. Default is 'bak'
    '''
    assert isinstance(fn, Path)
    where = fn.parent

    # get the backup files
    fbak = sorted(list(where.glob(f'{fn.name}.{suffix}.*')),
                  key=lambda p: int(p.name.split('.')[-1]))
    if fbak:
        # rename the elder by adding 1 to the suffix
        for i, f in enumerate(fbak[::-1]): # reverse order, to avoid overwrite
            j = len(fbak) - i + 1 #: STRU.bak.i -> STRU.bak.i+1
            fname = f.name.replace(f'.{j}', f'.{j+1}')
            f.rename(f.parent / fname)
    
    # backup the latest file, if there is one
    if fn.exists():
        fn.rename(fn.parent / f'{fn.name}.{suffix}.0')

def _write_stru(job_dir, stru, fname='STRU'):
    '''
    Generates a ABACUS STRU file from a STRU dict.

    Parameters
    ----------
    job_dir: str
        Directory in which the STRU file is generated.
    stru : dict
        Parameters to generate the STRU file.
        See the docstring at the beginning of this script for details.
    fname : str
        Name of the STRU file.

    '''
    with open(Path(job_dir) / fname, 'w') as f:

        #============ ATOMIC_SPECIES ============
        f.write('ATOMIC_SPECIES\n')
        width = {key + '_w' : max([len(str(s[key])) for s in stru['species']])
                 for key in ['symbol', 'mass', 'pp_file']}
        for s in stru['species']:
            f.write('{symbol:<{symbol_w}}  {mass:>{mass_w}}  {pp_file:>{pp_file_w}}'.format(**s, **width))
            if 'pp_type' in s:
                f.write(f"  {s['pp_type']}")
            f.write('\n')

        #============ NUMERICAL_ORBITAL ============
        if 'orb_file' in stru['species'][0]:
            f.write('\nNUMERICAL_ORBITAL\n')
            for s in stru['species']:
                f.write(f"{s['orb_file']}\n")
        
        #============ LATTICE_CONSTANT/PARAMETER/VECTORS ============
        f.write('\nLATTICE_CONSTANT\n')
        f.write(f"{stru['lat']['const']}\n")

        if 'vec' in stru['lat']:
            f.write('\nLATTICE_VECTORS\n')
            for v in stru['lat']['vec']:
                f.write(f'{v[0]} {v[1]} {v[2]}\n')

        if 'param' in stru['lat']:
            f.write('\nLATTICE_PARAMETER\n')
            for param in stru['lat']['param']:
                f.write(f'{param} ')
            f.write('\n')

        #============ ATOMIC_POSITIONS ============
        f.write('\nATOMIC_POSITIONS\n')
        f.write(f"{stru['coord_type']}\n")

        for s in stru['species']:
            f.write(f"\n{s['symbol']}\n")
            f.write(f"{s['mag_each']}\n")
            f.write(f"{s['natom']}\n")

            for atom in s['atom']:
                f.write(' '.join(f'{x}' for x in atom['coord']))

                for key in ['m', 'v']: # frozen atom / initial velocity
                    if key in atom:
                        f.write(f' {key}' +
                                ''.join(f' {x}' for x in atom[key]))

                if 'mag' in atom:
                    if not isinstance(atom['mag'], tuple): # collinear
                        f.write(f" mag {atom['mag']}")
                    else: # non-collinear
                        mag_coord_type, mag = atom['mag']
                        assert mag_coord_type in ['Cartesian', 'Spherical']
                        if mag_coord_type == 'Cartesian':
                            f.write(f' mag {mag[0]} {mag[1]} {mag[2]}')
                        else:
                            f.write(f' mag {mag[0]} angle1 {mag[1]} angle2 {mag[2]}')

                f.write('\n')

def write_stru(stru: Atoms,
               outdir: str,
               pp_file: Optional[Dict[str, str]],
               orb_file: Optional[Dict[str, str]] = None,
               fname: Optional[str] = 'STRU') -> str:
    '''
    write the ABACUS STRU file from an Atoms object into outdir.

    Parameters
    ----------
    stru : Atoms
        The Atoms object to write.
    outdir : str
        The directory to write the STRU file.
    pp_file : dict
        A dictionary mapping element symbols to pseudopotential file names.
        The keys are the element symbols, and the values are the paths to
        the pseudopotential files.
    orb_file : dict, optional
        A dictionary mapping element symbols to orbital file names.
        The keys are the element symbols, and the values are the paths to
        the orbital files. If None, no orbital files will be written.
    fname : str, optional
        The name of the STRU file to write. Default is 'STRU'.

    Returns
    -------
    str
        The path to the written STRU file.
    '''
    from ase.units import (
        Bohr as __BOHR__, 
        Angstrom as __ANGSTROM__
    )
    assert isinstance(stru, Atoms)
    pp_file = pp_file or {} # can be None, for those non-ESolver_KS cases

    elem = stru.get_chemical_symbols()    
    # ABACUS requires the atoms ranged species-by-species, therefore
    # we need to sort the atoms by species
    ind = np.argsort(elem)
    coords = stru.get_positions()[ind]
    elem = [elem[i] for i in ind]
    elem_uniq, nat = np.unique(elem, return_counts=True)
    stru_dict = {
        'coord_type': 'Cartesian',
        'lat': {
            'const': __ANGSTROM__ / __BOHR__, 
            'vec': np.array(stru.get_cell()).tolist()
        },
        'species': [
            {
                'symbol': e,
                'mass': ATOM_MASS[e],
                'pp_file': pp_file.get(e, ''),
                'pp_type': '',
                'natom': n,
                'mag_each': 0.0,
                'atom': [
                    {
                        'coord': coords[j].tolist(),
                        'm': [1, 1, 1],
                        'v': [0.0, 0.0, 0.0]
                    } for j in range(np.sum(nat[:i]), np.sum(nat[:i+1]))
                ]
            }
            for i, (n, e) in enumerate(zip(nat, elem_uniq))
        ]
    }
    if orb_file is not None:
        for s in stru_dict['species']:
            s['orb_file'] = orb_file[s['symbol']]
    
    _write_stru(outdir, stru_dict, fname)

    return (Path(outdir) / fname).resolve().as_posix()

def write_input(data: Dict[str, Any], fn: str) -> str:
    '''
    write the ABACUS INPUT file from a Dict[str, Any] object

    Parameters
    ----------
    data : dict
        The data to write to the INPUT file.
        The keys are the parameters and the values are the values of the parameters.
    fn : str
        The path to the ABACUS INPUT file.

    Returns
    -------
    str
        the absolute path of the written INPUT file
    '''
    with open(fn, 'w') as f:
        f.write('INPUT_PARAMETERS\n')
        for key, value in data.items():
            f.write(f'{key} {value}\n')
        f.write('\n')  # add a newline at the end

    return str(Path(fn).resolve())

def _write_kline(data: Dict[str, Any], f: TextIOWrapper):
    '''
    write the k-point file whose mode is "Line"
    '''
    f.write('K_POINTS\n')
    f.write(f"{data['nk']}\n")
    f.write('Line_Cartesian\n' if data['coordinate'].lower() == 'cartesian' 
                               else 'Line\n')
    for k, nk in zip(data['kpoints'], data['nkinterpl']):
        f.write(f"{k[0]} {k[1]} {k[2]} {nk}\n")

def _write_ksampl_mp(data: Dict[str, Any], f: TextIOWrapper):
    '''
    write the k-point file whose mode is automatically sampling with Monkhorst-Pack
    '''
    f.write('K_POINTS\n')
    f.write('0\n')
    f.write("Gamma\n" if data['gamma-centered'] else "MP\n")
    f.write(f"{data['nk'][0]} {data['nk'][1]} {data['nk'][2]} "
            f"{data['kshift'][0]} {data['kshift'][1]} {data['kshift'][2]}\n")
    
def _write_kpoint(data: Dict[str, Any], f: TextIOWrapper):
    '''
    write the k-point file whose mode is specifying kpoints one-by-one
    '''
    f.write('K_POINTS\n')
    f.write(f"{data['nk']}\n")
    f.write(f"{data['coordinate'].capitalize()}\n")
    nkinterpl = data.get('nkinterpl', [1 for _ in data['kpoints']])
    for k, nk in zip(data['kpoints'], data.get('nkinterpl', nkinterpl)):
        f.write(f"{k[0]} {k[1]} {k[2]} {nk}\n")

def read_kpoint_from_running_log(fn: str) -> Dict[str, Any]:
    '''Read the k-point information from the ABACUS running log file.'''
    raise NotImplementedError('read_kpoint_from_running_log is not implemented yet.')

def write_kpt(data: Dict[str, Any], fn: str) -> str:
    '''
    write the ABACUS KPT file from a Dict[str, Any] object

    Parameters
    ----------
    data : dict
        The data to write to the KPT file.
        The keys are the parameters and the values are the values of the parameters.
    fn : str
        The path to the ABACUS KPT file.

    Returns
    -------
    str
        the absolute path of the written KPT file
    '''
    assert 'mode' in data, 'KPT data must contain "mode" key'
    mywriter = {
         # 'line': _write_kline,
        'mp-sampling': _write_ksampl_mp,
        'point': _write_kpoint
    }
    assert data['mode'] in mywriter, \
           f'Invalid KPT mode {data["mode"]}, must be one of ' \
           f'{list(mywriter.keys())}'
    with open(fn, 'w') as f:
        mywriter[data['mode']](data, f)

    return str(Path(fn).resolve())

def _parse_coord_line(line):
    '''
    Parses a coordinate line (which may include extra parameters)
    in the ATOMIC_POSITIONS block.

    A coordinate line always contains the x, y, z coordinates of an atom,
    and may also include
        - whether an atom is frozen in MD or relaxation
        - initial velocity of an atom in MD or relaxation
        - magnetic moment of an atom

    Reference
    ---------
    https://abacus.deepmodeling.com/en/latest/advanced/input_files/stru.html

    (see section "More Key Words" for details)

    '''
    fields = line.split()
    result = { 'coord' : [float(x) for x in fields[0:3]] }

    idx = 3
    while idx < len(fields):
        if fields[idx].isdigit(): # no keyword, 0/1 -> frozen atom
            result['m'] = [int(x) for x in fields[idx:idx+3]]
            idx += 3
        elif fields[idx] == 'm': # frozen atom
            result['m'] = [int(x) for x in fields[idx+1:idx+4]]
            idx += 4
        elif fields[idx] in ['v', 'vel', 'velocity']: # initial velocity
            result['v'] = [float(x) for x in fields[idx+1:idx+4]]
            idx += 4
        elif fields[idx] in ['mag', 'magmom']:
            '''
            here we assume that frozen atom info cannot be placed after
            a collinear mag info without a keyword
            i.e., the following coordinate line
                0.0 0.0 0.0 mag 1.0 0 0 0
            is not allowed; one must explicitly specify 'm' in this case:
                0.0 0.0 0.0 mag 1.0 m 0 0 0

            '''
            if idx + 2 < len(fields) and fields[idx+2] == 'angle1':
                result['mag'] = ('Spherical',
                                 list(map(float, fields[idx+1:idx+6:2])))
                idx += 6
            elif idx + 2 < len(fields) and fields[idx+2][0].isdigit():
                result['mag'] = ('Cartesian',
                                 list(map(float, fields[idx+1:idx+4])))
                idx += 4
            else: # collinear
                result['mag'] = float(fields[idx+1])
                idx += 2
        else:
            raise ValueError('Error: unknown keyword %s'%fields[idx])

    return result

def _atomic_positions_gen(lines):
    '''
    Iteratively generates info per species from the ATOMIC_POSITIONS block.

    '''
    natom = int(lines[2])
    yield {'symbol': lines[0], 'mag_each': float(lines[1]), 'natom': natom,
           'atom': [_parse_coord_line(line) for line in lines[3:3+natom]]}
    if len(lines) > 3 + natom:
        yield from _atomic_positions_gen(lines[3+natom:])

def read_stru(fn: str) -> Dict[str, Any]:
    '''
    read the ABACUS STRU file and return a comprehensive Dict[str, Any]
    object. This function is implemented by @jinzx10 in ABACUS-CSW-NAO
    project.

    Parameters
    ----------
    fn : str
        The path to the ABACUS STRU file.

    Returns
    -------
        A dict containing the following keys-value pairs:
        'species' : list of dict
            List of atomic species. Each dict contains 'symbol', 'mass',
            'pp_file', and optionally 'pp_type'.
    '''
    block_title = ['ATOMIC_SPECIES',
                   'NUMERICAL_ORBITAL',
                   'LATTICE_CONSTANT',
                   'LATTICE_PARAMETER',
                   'LATTICE_VECTORS',
                   'ATOMIC_POSITIONS']

    def _trim(line):
        return line.split('#')[0].split('//')[0].strip(' \t\n')

    with open(fn, 'r') as f:
        lines = [_trim(line).replace('\t', ' ')
                 for line in f.readlines() if len(_trim(line)) > 0]

    # break the content into blocks
    delim = [i for i, line in enumerate(lines) if line in block_title] \
            + [len(lines)]
    blocks = {lines[delim[i]] : lines[delim[i]+1:delim[i+1]]
              for i in range(len(delim) - 1)}

    stru = {}
    #============ LATTICE_CONSTANT/PARAMETER/VECTORS ============
    stru['lat'] = {'const': float(blocks['LATTICE_CONSTANT'][0])}
    if 'LATTICE_VECTORS' in blocks:
        stru['lat']['vec'] = [[float(x) for x in line.split()]
                              for line in blocks['LATTICE_VECTORS']]
    elif 'LATTICE_PARAMETER' in blocks:
        stru['lat']['param'] = [float(x)
                                for x in blocks['LATTICE_PARAMETERS'].split()]

    #============ ATOMIC_SPECIES ============
    stru['species'] = [dict(zip(['symbol', 'mass', 'pp_file', 'pp_type'],
                                line.split()))
                       for line in blocks['ATOMIC_SPECIES']]
    for s in stru['species']:
        s['mass'] = float(s['mass'])

    #============ NUMERICAL_ORBITAL ============
    if 'NUMERICAL_ORBITAL' in blocks:
        for i, s in enumerate(stru['species']):
            s['orb_file'] = blocks['NUMERICAL_ORBITAL'][i]

    #============ ATOMIC_POSITIONS ============
    stru['coord_type'] = blocks['ATOMIC_POSITIONS'][0]
    index = { s['symbol'] : i for i, s in enumerate(stru['species']) }
    for ap in _atomic_positions_gen(blocks['ATOMIC_POSITIONS'][1:]):
        stru['species'][index[ap['symbol']]].update(ap)

    return stru

def read_input(fn: str) -> Dict[str, Any]:
    '''
    read the ABACUS INPUT file and return a comprehensive Dict[str, Any]

    Parameters
    ----------
    fn : str
        The path to the ABACUS INPUT file.
    
    Returns
    -------
    dict
    '''
    with open(fn) as f:
        raw = [l.strip() for l in f.readlines()]
    raw = [l for l in raw 
           if l and not re.match(r'^(INPUT_PARAMETERS|#|//|!)', l)]
    raw = [re.split(r'#|//|!', l)[0].strip() for l in raw]
    raw = [re.split(r'\s+', l) for l in raw]

    return dict([(l[0], ' '.join(l[1:])) for l in raw if len(l) > 0])

def _read_kline(raw: List[str]) -> Dict[str, Any]:
    '''
    read the k-point file whose mode is "Line", and return the parsed result
    '''
    assert len(raw) >= 4, \
             f'Invalid KPT file, expected at least 4 lines, got {len(raw)}'
    assert re.match(r'^\d+$', raw[1]), \
             f'Invalid KPT file, expected the second line to be an integer, ' \
             f'got {raw[1]}'
    assert int(raw[1]) == len(raw) - 3, \
             f'Invalid KPT file, expected the second line to be the number of ' \
             f'k-points, got {raw[1]}, but found {len(raw)} k-points'
    assert raw[2].lower().startswith('line'), \
             f'Invalid KPT file, expected "Line" in the third line, got {raw[2]}'
    mymatch = [re.match(r'^(-?\d+(\.\d+)?)\s+(-?\d+(\.\d+)?)\s+(-?\d+(\.\d+)?)\s+'
                        r'(\d+).*', l) 
              for l in raw[3:]]
    assert all(m for m in mymatch), \
             'Invalid KPT file, expected the k-points to be in the format ' \
             '"x y z n # comment"'
    print(raw)
    return {
        'mode': 'line',
        'coordinate': 'Cartesian' if raw[2].lower().endswith('cartesian') else 'Direct', 
        'nk': int(raw[1]),
        'kpoints': [tuple([float(x) for x in m.groups()[:3]]) for m in mymatch],
        'nkinterpl': [int(m.groups()[6]) for m in mymatch]
    }

def _read_ksampl_mp(raw: List[str]) -> Dict[str, Any]:
    '''
    read the k-point file whose mode is automatically sampling with Monkhorst-Pack
    scheme, and return the parsed result
    '''
    assert len(raw) == 4, \
             f'Invalid KPT file, expected 4 lines, got {len(raw)}'
    assert raw[2].lower() in ['gamma', 'mp'], \
             f'Invalid KPT file, expected "Gamma" or "MP" in the third line, ' \
             f'got {raw[2]}'
    # the following may raise the ValueError if unpacking is not correct
    nk1, nk2, nk3, kshift1, kshift2, kshift3 = map(int, raw[3].split())
    return {
        'mode': 'mp-sampling',
        'gamma-centered': bool(raw[2].lower() == 'gamma'),
        'nk': (nk1, nk2, nk3),
        'kshift': (kshift1, kshift2, kshift3)
    }

def _read_kpoint(raw: List[str]) -> Dict[str, Any]:
    '''
    read the k-point file whose mode is specifying kpoints one-by-one, return
    the parsed result.
    '''
    assert len(raw) >= 4, \
             f'Invalid KPT file, expected at least 4 lines, got {len(raw)}'
    assert re.match(r'^\d+$', raw[1]), \
             f'Invalid KPT file, expected the second line to be an integer, ' \
             f'got {raw[1]}'
    assert int(raw[1]) == len(raw) - 3, \
             f'Invalid KPT file, expected the second line to be the number of ' \
             f'k-points, got {raw[1]}, but found {len(raw)} k-points'
    assert raw[2].lower() in ['direct', 'cartesian'], \
             f'Invalid KPT file, expected "Direct" or "Cartesian" in the third line, ' \
             f'got {raw[2]}'
    mymatch = [re.match(r'^(-?\d+(\.\d+)?)\s+(-?\d+(\.\d+)?)\s+(-?\d+(\.\d+)?)\s*'
                        r'(\d+).*', l) 
              for l in raw[3:]]
    assert all(m for m in mymatch), \
             'Invalid KPT file, expected the k-points to be in the format ' \
             '"x y z n # comment"'

    return {
        'mode': 'point',
        'coordinate': 'Cartesian' if raw[2].lower() == 'cartesian' else 'Direct',
        'nk': int(raw[1]),
        'kpoints': [tuple(map(float, m.groups()[:3])) for m in mymatch],
        'weights': [int(m.groups()[6]) for m in mymatch]
    }

def read_kpt(fn: str) -> Dict[str, Any]:
    '''
    read the ABACUS KPT file and return a comprehensive Dict[str, Any]

    Parameters
    ----------
    fn : str
        The path to the ABACUS KPT file.
    
    Returns
    -------
    dict
        A dictionary containing the k-point information.
    '''
    with open(fn) as f:
        raw = [l.strip() for l in f.readlines()]
    raw = [l for l in raw if l]
    parser_switch = {'line': _read_kline,
                     'gamma': _read_ksampl_mp, 'mp': _read_ksampl_mp,
                     'direct': _read_kpoint, 'cartesian': _read_kpoint}
    assert raw[0] == 'K_POINTS', \
             f'Invalid KPT file {fn}, first line must be "K_POINTS"'
    # because there are modes like "Line_Cartesian", we need to
    # split the second line by '_'
    return parser_switch[raw[2].lower().split('_')[0]](raw)

class TestAbacusCalculatorIOUtil(unittest.TestCase):

    def setUp(self):
        here = Path(__file__).parent
        self.testfiles = here / 'testfiles'

    def test_input_io(self):
        fn = self.testfiles / 'input-script'
        data = read_input(fn)
        with tempfile.NamedTemporaryFile(mode='w') as f:
            fn_ = Path(f.name)
            write_input(data, fn_.as_posix())
            data_ = read_input(fn_.as_posix())
            self.assertDictEqual(data, data_)
        # will automatically delete the file after the context manager

    def test_stru_io(self):
        from ase.units import Bohr, Angstrom
        nacl = bulk('NaCl', 'rocksalt', a=5.64)
        write_stru(
            nacl, 
            outdir=self.testfiles, 
            pp_file={
                'Na': 'Na.pz-bhs.UPF',
                'Cl': 'Cl.pz-bhs.UPF'
            },
            orb_file={
                'Na': 'Na_gga_6au_100Ry_2s2p1d.orb',
                'Cl': 'Cl_gga_6au_100Ry_2s2p1d.orb',
            },
        )
        stru_ = read_stru(self.testfiles / 'STRU')
        (self.testfiles / 'STRU').unlink()
        
        self.assertIsInstance(stru_, dict)
        for k in ['lat', 'species', 'coord_type']:
            self.assertIn(k, stru_)
        for k in ['const', 'vec']:
            self.assertIn(k, stru_['lat'])
        for k in ['symbol', 'mass', 'pp_file', 'orb_file', 'mag_each', 'natom', 'atom']:
            for s in stru_['species']:
                self.assertIn(k, s)
        self.assertEqual(stru_['coord_type'], 'Cartesian')

        self.assertEqual(stru_['lat']['const'], Angstrom / Bohr)
        self.assertTrue(np.allclose(stru_['lat']['vec'], np.array(nacl.get_cell())))
        self.assertIsInstance(stru_['species'], list)

        for s in stru_['species']:
            self.assertIsInstance(s, dict)
            self.assertEqual(s['natom'], 1)
            self.assertIsInstance(s['atom'], list)
            self.assertEqual(s['natom'], len(s['atom']))
            for a in s['atom']:
                self.assertIsInstance(a, dict)
                for k in ['coord', 'm', 'v']:
                    self.assertIn(k, a)
                self.assertEqual(a['m'], [1, 1, 1])
                self.assertEqual(a['v'], [0.0, 0.0, 0.0])
                
        self.assertEqual(stru_['species'][0]['symbol'], 'Cl')
        self.assertEqual(stru_['species'][1]['symbol'], 'Na')
        self.assertEqual(stru_['species'][0]['mass'], ATOM_MASS['Cl'])
        self.assertEqual(stru_['species'][1]['mass'], ATOM_MASS['Na'])
        self.assertEqual(stru_['species'][0]['pp_file'], 'Cl.pz-bhs.UPF')
        self.assertEqual(stru_['species'][1]['pp_file'], 'Na.pz-bhs.UPF')
        self.assertEqual(stru_['species'][0]['orb_file'], 'Cl_gga_6au_100Ry_2s2p1d.orb')
        self.assertEqual(stru_['species'][1]['orb_file'], 'Na_gga_6au_100Ry_2s2p1d.orb')

    def test_kpt_io(self):
        kpt = {
            'mode': 'point',
            'coordinate': 'Cartesian',
            'nk': 1,
            'kpoints': [(0.0, 0.0, 0.0)],
            'weights': [1]
        }
        with tempfile.NamedTemporaryFile(mode='w') as f:
            fn_ = Path(f.name)
            write_kpt(kpt, fn_.as_posix())
            self.assertDictEqual(kpt, read_kpt(fn_.as_posix()))

        kpt = {
            'mode': 'mp-sampling',
            'nk': (5, 5, 5),
            'kshift': (0, 0, 0),
            'gamma-centered': True,
        }
        with tempfile.NamedTemporaryFile(mode='w') as f:
            fn_ = Path(f.name)
            write_kpt(kpt, fn_.as_posix())
            self.assertDictEqual(kpt, read_kpt(fn_.as_posix()))

        # kpt = {
        #     'mode': 'line',
        #     'nk': 6,
        #     'kpoints': [
        #         (0.0, 0.0, 0.0),
        #         (0.5, 0.5, 0.5),
        #         (0.0, 0.5, 0.5),
        #         (0.5, 0.5, 0.0),
        #         (0.5, 0.0, 0.5),
        #         (0.0, 0.0, 0.5),
        #     ],
        #     'nkinterpl': [10, 10, 10, 10, 10, 1],
        #     'coordinate': 'Cartesian',
        # }
        # with tempfile.NamedTemporaryFile(mode='w') as f:
        #     fn_ = Path(f.name)
        #     write_kpt(kpt, fn_.as_posix())
        #     self.assertDictEqual(kpt, read_kpt(fn_.as_posix()))

    @unittest.skip('')
    def test_load_pseudo(self):
        pseudo_dir = self.testfiles / 'pporb'
        pseudo_map = load_pseudo(pseudo_dir)
        self.assertEqual(pseudo_map['O'],  'O.upf')
        self.assertEqual(pseudo_map['Si'], 'Si.upf')
        self.assertEqual(pseudo_map['In'], 'In.upf')
        self.assertEqual(pseudo_map['As'], 'As.pbe-n-rrkjus_psl.0.2.UPF')
        self.assertEqual(pseudo_map['Ga'], 'Ga.pbe-dn-kjpaw_psl.1.0.0.UPF')

    @unittest.skip('')
    def test_load_orbital(self):
        orbital_dir = self.testfiles / 'pporb'
        orbital_map = load_orbital(orbital_dir)
        self.assertEqual(orbital_map['O'],  'O_gga_6au_100Ry_2s2p1d.orb')
        self.assertEqual(orbital_map['Si'], 'Si_gga_7au_100Ry_2s2p1d.orb')
        self.assertEqual(orbital_map['In'], 'In_gga_8au_100Ry_1s1p1d.orb')

        orbital_map = load_orbital(orbital_dir, efficiency=False)
        self.assertEqual(orbital_map['O'],  'O_gga_6au_100Ry_2s2p1d.orb')
        self.assertEqual(orbital_map['Si'], 'Si_gga_7au_100Ry_2s2p1d.orb')
        self.assertEqual(orbital_map['In'], 'In_gga_8au_100Ry_2s2p2d1f.orb')

if __name__ == '__main__':
    unittest.main()