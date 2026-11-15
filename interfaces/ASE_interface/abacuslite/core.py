# fmt: off

'''
Adapted from the ase-abacus implementation:
https://gitlab.com/1041176461/ase-abacus/-/blob/master/ase/calculators/abacus.py

This module defines an ASE interface to ABACUS.
Created on Fri Jun  8 16:33:38 2018

ABACUS (Atomic-orbital Based Ab-initio Computation at UStc) is an open-source 
package based on density functional theory (DFT). The package utilizes both plane 
wave and numerical atomic basis sets with the usage of pseudopotentials to describe 
the interactions between nuclear ions and valence electrons. ABACUS supports LDA, 
GGA, meta-GGA, and hybrid functionals. Apart from single-point calculations, 
the package allows geometry optimizations and ab-initio molecular dynamics with 
various ensembles. The package also provides a variety of advanced functionalities 
for simulating materials, including the DFT+U, VdW corrections, and implicit solvation
model, etc. In addition, ABACUS strives to provide a general infrastructure to 
facilitate the developments and applications of novel machine-learning-assisted 
DFT methods (DeePKS, DP-GEN, DeepH, DeePTB etc.) in molecular and material simulations.

Modified on Wed Jun 20 15:00:00 2018
@author: Shen Zhen-Xiong

Modified on Wed Jun 03 23:00:00 2022
@author: Ji Yu-yang

Refactored on 2025-12-07
@author: Huang Yi-ke
'''

import json
import os
import re
import shutil
import tempfile
import unittest
from copy import copy
from pathlib import Path
from typing import Dict, Optional, List

import numpy as np
from ase.calculators.genericfileio import (
    BaseProfile,
    CalculatorTemplate,
    GenericFileIOCalculator,
    read_stdout
)
from ase.calculators.socketio import SocketIOCalculator
from ase.atoms import Atoms
from ase.dft.kpoints import BandPath
from ase.io import read

from abacuslite.io.generalio import (
    file_safe_backup,
    read_input,
    read_stru,
    read_kpt,
    species_group_indices,
    write_input,
    write_stru,
    write_kpt
)

__LEGACYIO__ = True
SOCKET_METADATA_SCHEMA = 'abacus.socket.properties.v1'


def _schema_compatible_properties(results, natoms):
    """Project calculator results through ASE's installed Properties schema."""
    from ase.calculators.calculator import Properties, all_outputs

    exported = {}
    for name, value in results.items():
        if name not in all_outputs:
            continue
        shape = np.shape(value)
        if name == 'magmom' and shape != ():
            continue
        if name == 'magmoms' and shape != (natoms,):
            continue
        exported[name] = value
    return Properties(exported)


def _uses_legacy_io(version: str) -> bool:
    """Resolve a format without mutating another calculator's reader."""
    match = re.match(
        r'^v(\d+)\.(\d+)\.(\d+)(\.\d+|-(alpha|beta|rc)\.?\d+)?(?:\+[A-Za-z0-9.-]+)?$',
        version,
    )
    assert match, f'Invalid ABACUS version: {version}'
    assert int(match.group(1)) >= 3, 'ABACUS v2.x is not supported'
    minor = int(match.group(2))
    return not (minor >= 11 or (minor == 9 and match.group(4) is not None))


def switch_io_backend_version(version: str) -> bool:
    """Set the compatibility default for subsequently created templates.

    Existing calculators retain the format selected from their own profile.
    Explicit io.legacyio/io.latestio functions are unaffected by this helper.
    """
    global __LEGACYIO__
    __LEGACYIO__ = _uses_legacy_io(version)
    return __LEGACYIO__


class AbacusProfile(BaseProfile):
    '''AbacusProfile for interacting the ASE with ABACUS that installed in
    the practical system'''
    configvars = {'pseudo_dir', 'orbital_dir'}

    def __init__(self, 
                 command: str, 
                 pseudo_dir: Optional[str | Path] = None, 
                 orbital_dir: Optional[str | Path] = None, 
                 omp_num_threads: Optional[int] = None,
                 **kwargs):
        '''Initialize ABACUS profile.
        
        Parameters
        ----------
        command : str
            The command to run ABACUS. NOTE: there may be the case for some
            sophisticated ABACUS user they call ABACUS with command like
            `OMP_NUM_THREADS=1 mpirun -np X abacus`. Here please do not set
            the number of omp threads in `command`, instead, use `nomp=1`.
        pseudo_dir : str or Path, optional
            The directory containing pseudopotential files. If ``None`` (the
            default), ``ABACUS_PP_PATH`` is used when set; an explicit value,
            including an empty string, is preserved.
        orbital_dir : str or Path, optional
            The directory containing orbital basis files. This is only necessary
            for an ABACUS-LCAO calculation. If ``None`` (the default),
            ``ABACUS_ORB_PATH`` is used when set; an explicit value, including
            an empty string, is preserved.
        omp_num_threads : int, optional
            The number of omp threads to use.
        '''
        assert isinstance(command, str)
        # further validation on the command will be in the __init__ of
        # the base class
        super().__init__(command, **kwargs)
        # These two community conventions are the only environment fallbacks
        # in this package; explicit arguments always take precedence. Only an
        # omitted (None) argument falls back, so an explicit empty string is
        # preserved instead of being replaced by the environment.
        self.pseudo_dir = (pseudo_dir if pseudo_dir is not None
                           else (os.environ.get("ABACUS_PP_PATH") or None))
        self.orbital_dir = (orbital_dir if orbital_dir is not None
                            else (os.environ.get("ABACUS_ORB_PATH") or None))

        if omp_num_threads is not None:
            # set the number of omp threads for the present process
            assert isinstance(omp_num_threads, int)
            os.environ['OMP_NUM_THREADS'] = str(omp_num_threads)

    @staticmethod
    def parse_version(stdout) -> str:
        # MPI launchers may add informational lines before ABACUS output.
        match = re.search(r'ABACUS version (\S+)', stdout or '')
        if match is None:
            raise RuntimeError(
                'Could not parse ABACUS version from command output. '
                'Expected a line like "ABACUS version vX.Y.Z".'
            )
        return match.group(1)

    def get_calculator_command(self, inputfile) -> List[str]:
        # because ABACUS run in the folder where there are INPUT files, so the
        # additional inputfile argument is not used.
        return []

    def socketio_argv_inet(self, port: Optional[int] = None) -> List[str]:
        port = 31415 if port is None else port
        return [
            'env',
            f'ABACUS_SOCKET_ADDRESS=localhost:{port}',
            *self._split_command,
        ]

    def socketio_argv_unix(self, socket: str) -> List[str]:
        return [
            'env',
            f'ABACUS_SOCKET_ADDRESS=/tmp/ipi_{socket}:UNIX',
            *self._split_command,
        ]

    def version(self) -> str:
        '''get the abacus version information'''
        cmd_ = [*self._split_command, '--version']
        return AbacusProfile.parse_version(read_stdout(cmd_))

def _validate_atomorder(atomorder, natoms, stru_path):
    """Validate the species-to-ASE permutation before handing it to a parser."""
    if atomorder is None:
        return list(range(natoms))

    try:
        values = list(atomorder)
    except (TypeError, ValueError) as exc:
        raise RuntimeError(
            f"Cannot parse scf atomorder (STRU: {stru_path}; expected "
            f"{natoms} integer indices; fail-closed)"
        ) from exc

    if (
        len(values) != natoms
        or not all(
            isinstance(value, (int, np.integer))
            and not isinstance(value, (bool, np.bool_))
            for value in values
        )
        or sorted(values) != list(range(natoms))
    ):
        raise RuntimeError(
            f"Invalid scf atomorder (STRU: {stru_path}; received {values!r}, "
            f"expected {natoms} unique indices from 0 to {max(natoms - 1, 0)}; "
            f"fail-closed)"
        )
    return [int(value) for value in values]


def _invert_cell(cell, context):
    """Return a cell inverse or raise a diagnostic validation error."""
    try:
        return np.linalg.solve(cell, np.eye(3, dtype=float))
    except np.linalg.LinAlgError as exc:
        raise ValueError(f"singular cell in {context}") from exc


def _read_stru_reference(directory, stru_file, atomorder, log_path=None):
    """Normalize STRU coordinates, cell, symbols, and mapping for SCF matching."""
    from ase.units import Bohr
    from abacuslite.io.generalio import read_stru

    stru_path = Path(directory) / (stru_file or 'STRU')
    try:
        stru = read_stru(stru_path)
        species = stru['species']
        lat = stru['lat']
        coord_type = str(stru['coord_type']).lower()
        raw_symbols = []
        raw_coords = []
        for species_entry in species:
            symbol = species_entry['symbol']
            atoms = species_entry['atom']
            declared_natoms = species_entry['natom']
            if declared_natoms != len(atoms):
                raise ValueError(
                    f"STRU species {symbol!r} declared natom={declared_natoms}, "
                    f"but contains {len(atoms)} atom records"
                )
            for atom in atoms:
                raw_symbols.append(str(symbol))
                raw_coords.append(atom['coord'])

        natoms = len(raw_symbols)
        mapping = _validate_atomorder(atomorder, natoms, stru_path)
        coordinates = np.asarray(raw_coords, dtype=float)
        if coordinates.shape != (natoms, 3):
            raise ValueError(
                f"STRU positions shape {coordinates.shape}, expected ({natoms}, 3)"
            )
        if not np.isfinite(coordinates).all():
            raise ValueError("non-finite STRU positions")

        lattice_vectors = np.asarray(lat['vec'], dtype=float)
        if lattice_vectors.shape != (3, 3):
            raise ValueError(
                f"STRU cell shape {lattice_vectors.shape}, expected (3, 3)"
            )
        lattice_constant = float(lat['const'])
        cell = lattice_vectors * lattice_constant * Bohr
        if not np.isfinite(cell).all():
            raise ValueError("non-finite STRU cell")
        _invert_cell(cell, f"STRU {stru_path}")

        if coord_type.startswith('d'):
            species_positions = coordinates @ cell
        elif coord_type.startswith('c'):
            species_positions = coordinates * lattice_constant * Bohr
        else:
            raise ValueError(f"unsupported STRU coordinate type {stru['coord_type']!r}")
        if not np.isfinite(species_positions).all():
            raise ValueError("non-finite normalized STRU positions")
    except RuntimeError:
        raise
    except (
        AssertionError,
        KeyError,
        OSError,
        TypeError,
        ValueError,
        IndexError,
        np.linalg.LinAlgError,
    ) as exc:
        suffix = f"；log: {log_path}" if log_path is not None else ""
        raise RuntimeError(
            f"Invalid scf STRU identity data (STRU: {stru_path}{suffix}; "
            f"{exc}; fail-closed)"
        ) from exc

    return {
        'path': stru_path,
        'positions': np.asarray(species_positions, dtype=float)[mapping].copy(),
        'cell': np.asarray(cell, dtype=float).copy(),
        'symbols': [raw_symbols[index] for index in mapping],
        'atomorder': mapping,
    }


def _stru_positions_in_ase_order(directory, stru_file, atomorder):
    """Return normalized STRU positions for compatibility with existing callers."""
    return _read_stru_reference(directory, stru_file, atomorder)['positions']


def _frame_coordinate_metadata(log_path, *, legacyio=None):
    """Return the frame coordinate convention and lattice constant from a log.

    The selected reader returns fractional coordinates for DIRECT logs and
    Cartesian coordinates for CARTESIAN logs. Use the coordinate header
    actually printed in the log, not the STRU ``coord_type``: ``write_stru``
    may write a Cartesian STRU while ABACUS prints DIRECT coordinates. Both
    readers expose the same ``coordinate`` and ``alat_in_angstrom`` fields
    through ``read_traj_from_running_log``.
    """
    use_legacy = __LEGACYIO__ if legacyio is None else legacyio
    if use_legacy:
        from abacuslite.io.legacyio import read_traj_from_running_log
    else:
        from abacuslite.io.latestio import read_traj_from_running_log
    traj = read_traj_from_running_log(log_path)
    direct = str(traj[0]['coordinate']).lower().startswith('d')
    scale = float(traj[0]['alat_in_angstrom'])
    if not np.isfinite(scale) or scale <= 0:
        raise ValueError(f"Invalid lattice constant in Angstrom: {scale!r}")
    return direct, scale


def _frame_coordinate_is_direct(log_path, *, legacyio=None):
    """Return the coordinate convention exposed by the selected log reader."""
    return _frame_coordinate_metadata(log_path, legacyio=legacyio)[0]


def _validate_raw_scf_trajectory(log_path, expected, *, legacyio=None):
    """Check raw trajectory cardinality before parser atom reordering."""
    use_legacy = __LEGACYIO__ if legacyio is None else legacyio
    if use_legacy:
        from abacuslite.io.legacyio import read_traj_from_running_log
    else:
        from abacuslite.io.latestio import read_traj_from_running_log

    try:
        trajectory = read_traj_from_running_log(log_path)
        if not trajectory:
            raise ValueError("no raw trajectory frames")
        expected_shape = expected['positions'].shape
        expected_natoms = len(expected['symbols'])
        for frame_index, frame in enumerate(trajectory):
            coords = np.asarray(frame['coords'], dtype=float)
            if coords.shape != expected_shape:
                mismatch = (
                    f"atom count {coords.shape[0]}, expected {expected_shape[0]}; "
                    if coords.ndim >= 1 and coords.shape[0] != expected_shape[0]
                    else ""
                )
                raise ValueError(
                    f"frame {frame_index} {mismatch}positions shape {coords.shape}, "
                    f"expected {expected_shape}"
                )
            if not np.isfinite(coords).all():
                raise ValueError(f"frame {frame_index} positions are non-finite")
            elements = np.asarray(frame['elem'])
            if elements.ndim != 1 or len(elements) != expected_natoms:
                raise ValueError(
                    f"frame {frame_index} atom count {len(elements)}, "
                    f"expected {expected_natoms}"
                )
            cell = np.asarray(frame['cell'], dtype=float)
            if cell.shape != (3, 3):
                raise ValueError(
                    f"frame {frame_index} cell shape {cell.shape}, expected (3, 3)"
                )
            if not np.isfinite(cell).all():
                raise ValueError(f"frame {frame_index} cell is non-finite")
            _invert_cell(cell, f"frame {frame_index}")
    except RuntimeError:
        raise
    except (AssertionError, KeyError, OSError, TypeError, ValueError, IndexError) as exc:
        raise RuntimeError(
            f"Invalid raw scf trajectory (log: {log_path}; {exc}; fail-closed)"
        ) from exc
    return trajectory


def _normalize_scf_frame(frame, frame_index, expected, frame_is_direct, atol,
                         cartesian_scale=1.0):
    """Normalize one parser frame and return match diagnostics.

    A finite, well-shaped frame with a different position or cell is a valid
    candidate mismatch. Malformed data is represented as ``ValueError`` so the
    caller can validate every frame before selecting any one of them.
    """
    context = f"frame {frame_index}"
    try:
        symbols = list(frame.get_chemical_symbols())
    except (AttributeError, TypeError, ValueError) as exc:
        raise ValueError(f"{context} symbols are unavailable") from exc
    if len(symbols) != len(expected['symbols']):
        raise ValueError(
            f"{context} atom count {len(symbols)}, expected {len(expected['symbols'])}"
        )
    if symbols != expected['symbols']:
        raise ValueError(
            f"{context} symbols {symbols!r} do not match STRU symbols "
            f"{expected['symbols']!r}"
        )

    try:
        positions = np.asarray(frame.positions, dtype=float)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{context} positions cannot be converted to float") from exc
    if positions.shape != expected['positions'].shape:
        raise ValueError(
            f"{context} positions shape {positions.shape}, expected "
            f"{expected['positions'].shape}"
        )
    if not np.isfinite(positions).all():
        raise ValueError(f"{context} positions are non-finite")

    try:
        cell = np.asarray(frame.cell, dtype=float)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{context} cell cannot be converted to float") from exc
    if cell.shape != (3, 3):
        raise ValueError(f"{context} cell shape {cell.shape}, expected (3, 3)")
    if not np.isfinite(cell).all():
        raise ValueError(f"{context} cell is non-finite")
    cell_inverse = _invert_cell(cell, context)

    # ABACUS prints Cartesian tau in units of alat, while cell is already in Å.
    frame_cart = positions @ cell if frame_is_direct else positions * cartesian_scale
    if not np.isfinite(frame_cart).all():
        raise ValueError(f"{context} normalized Cartesian positions are non-finite")
    cell_delta = cell - expected['cell']
    cell_diff = float(np.max(np.abs(cell_delta)))
    same_cell = bool(np.allclose(cell, expected['cell'], atol=atol, rtol=0.0))
    if not same_cell:
        return {
            'frame': frame,
            'frame_index': frame_index,
            'frame_cart': frame_cart,
            'cell_diff': cell_diff,
            'periodic_residual': None,
            'matches': False,
        }

    delta_frac = (frame_cart - expected['positions']) @ cell_inverse
    periodic_cart = (delta_frac - np.rint(delta_frac)) @ cell
    if not np.isfinite(periodic_cart).all():
        raise ValueError(f"{context} periodic residual is non-finite")
    periodic_residual = float(np.max(np.abs(periodic_cart)))
    return {
        'frame': frame,
        'frame_index': frame_index,
        'frame_cart': frame_cart,
        'cell_diff': cell_diff,
        'periodic_residual': periodic_residual,
        'matches': periodic_residual <= atol,
    }


def _select_scf_frame_for_structure(
    frames, log_path, directory, atomorder, atol=1e-4, stru_file='STRU', _reference=None,
    *, legacyio=None
):
    """Return the latest periodic-equivalent frame; malformed SCF data fails closed.

    Parameters
    ----------
    frames : list of Atoms
        Frames returned by ``read_abacus_out``, reordered into ASE atom order.
    log_path : Path
        Running log path, used in fail-closed diagnostics.
    directory : Path
        Calculation directory containing the STRU written by ``write_input``.
    atomorder : list of int
        Permutation from species order to ASE atom order.
    stru_file : str or Path
        The STRU filename recorded by ``write_input``; defaults to ``STRU``.
    atol : float
        Absolute Angstrom bound for cell components and periodic Cartesian
        residuals, with no relative tolerance. Explicit overrides are retained.
        The ``1e-4`` default accommodates ordinary running-log rounding:
        with a lattice constant near 1 Angstrom, six-significant-digit
        vectors in ``[10, 100)`` can lose up to ``5e-5`` Angstrom per
        component. A ``1e-5`` bound rejected valid recorded LTS output.
        This fixed absolute bound is not a universal printing-error budget:
        larger cells, independently rounded lattice constants and skewed
        cells can accumulate more error. It does not change an SCF or
        optimization convergence threshold.

    Raises
    ------
    RuntimeError
        If no frame matches the current STRU, or a frame is malformed.
    """
    try:
        nframes = len(frames)
    except TypeError as exc:
        raise RuntimeError(
            f"Invalid scf frame collection (log: {log_path}; cannot read "
            f"frame count; fail-closed)"
        ) from exc
    if nframes == 0:
        raise RuntimeError(
            f"No scf frames in log {log_path} (0 frames; atol={atol}; "
            f"fail-closed)"
        )
    if not np.isfinite(atol) or atol < 0:
        raise RuntimeError(f"Invalid scf frame-matching atol: {atol!r} "
                           f"(fail-closed)")

    expected = _reference or _read_stru_reference(
        directory, stru_file, atomorder, log_path=log_path
    )
    frame_label = 'frame' if nframes == 1 else 'frames'
    try:
        frame_is_direct, cartesian_scale = _frame_coordinate_metadata(
            log_path, legacyio=legacyio
        )
    except Exception as exc:
        raise RuntimeError(
            f"Cannot parse scf log coordinate header: {log_path} "
            f"({exc}; fail-closed)"
        ) from exc

    normalized = []
    invalid = []
    for frame_index, frame in enumerate(frames):
        try:
            normalized.append(
                _normalize_scf_frame(
                    frame, frame_index, expected, frame_is_direct, atol,
                    cartesian_scale=cartesian_scale
                )
            )
        except (AttributeError, TypeError, ValueError, IndexError) as exc:
            invalid.append(f"frame {frame_index}: {exc}")

    if invalid:
        details = "; ".join(invalid)
        raise RuntimeError(
            f"Malformed scf frames in log {log_path} ({nframes} {frame_label}; "
            f"{details}; atol={atol}; fail-closed)"
        )

    for info in reversed(normalized):
        if info['matches']:
            return info['frame']

    last = normalized[-1]
    raw_diff = float(
        np.max(np.abs(last['frame_cart'] - expected['positions']))
    )
    periodic = (
        f"{last['periodic_residual']:.6g} Å"
        if last['periodic_residual'] is not None
        else "n/a (cell mismatch)"
    )
    raise RuntimeError(
        f"No frame in scf log matches the current structure: {log_path} "
        f"({nframes} {frame_label}; last-frame raw Cartesian max difference "
        f"{raw_diff:.6g} Angstrom; periodic residual={periodic}; "
        f"maximum cell difference={last['cell_diff']:.6g} Angstrom; "
        f"atol={atol}; fail-closed)"
    )


class AbacusTemplate(CalculatorTemplate):
    
    implemented_properties = [
        'energy', 'forces', 'stress', 'free_energy', 'magmom'
    ]
    _label = 'abacus'

    def __init__(self):
        self._legacyio = __LEGACYIO__
        super().__init__(
            'abacus',
            self.implemented_properties
        )
        self.non_convergence_ok = False
        # the redirect stdout and stderr
        self.inputname  = 'INPUT' # hard-coded
        self.outputname = f'{self._label}.out'
        self.errorname  = f'{self._label}.err'

        # fix: inconsistent atoms order may induce bugs, here a list
        # is kept to swap the order of atoms
        self.atomorder  = None

    '''because it may be not one-to-one mapping between the property
    desired to calculate and the keywords used in the calculation,
    in the following a series of functions for mapping the property
    calculation to the keywords settings are implemented'''
    @staticmethod
    def get_energy_keywords(self) -> Dict[str, str]:
        return {}

    @staticmethod
    def get_forces_keywords(self) -> Dict[str, str]:
        return {'cal_force': '1'}
    
    @staticmethod
    def get_stress_keywords(self) -> Dict[str, str]:
        return {'cal_stress': '1'}

    @staticmethod
    def get_free_energy_keywords(self) -> Dict[str, str]:
        return {}

    @staticmethod
    def get_magmom_keywords(self) -> Dict[str, str]:
        return {'nspin': '2'}

    def get_property_keywords(self,
                              parameters: Dict[str, str],
                              properties: List[str]) -> Dict[str, str]:
        '''Connect the relationship between the properties calculation and
        the ABACUS keywords. May be more complicated in the future, therefore
        it is better to have a separate mapping function instead of
        implementing in some other functions.
        
        Parameters
        ----------
        parameters : dict
            The parameters used to perform the calculation.
        properties : list of str
            The list of properties to calculate
        '''
        def keyword_compare_value(value):
            if isinstance(value, bool):
                return '1' if value else '0'
            if isinstance(value, (list, tuple, set)):
                return ' '.join(str(i) for i in value)
            return str(value)

        param_cache_ = {
            key: keyword_compare_value(value)
            for key, value in parameters.items()
            if value is not None
        }

        def counter(param_new: Dict[str, str]) -> Dict[str, str]:
            info = 'desired properties or explicit parameters required contradictory keywords'
            staged = {}
            for k, v in param_new.items():
                if v is None:
                    continue
                normalized_value = keyword_compare_value(v)
                if k in param_cache_ and param_cache_[k] != normalized_value:
                    raise ValueError(f'{info}: {k}={v} (now), {param_cache_[k]} (before)')
                staged[k] = normalized_value
            param_cache_.update(staged)
            return param_new

        # update the parameters with the keywords for the properties
        for p in properties:
            assert p in self.implemented_properties
            parameters.update(counter(getattr(self, f'get_{p}_keywords')(parameters)))
        
        # from the parameters, get the file path
        self.suffix = parameters.get('suffix', 'ABACUS')
        self.calculation = parameters.get('calculation', 'scf')
        # with the above two, the running log file can be positioned.
        return parameters

    def write_input(self, 
                    profile: AbacusProfile, 
                    directory: Path | str,
                    atoms: Atoms, 
                    parameters: Dict[str, str],
                    properties: List[str]) -> None:
        '''Write the input files for the calculation. This function connects
        the calculation in ASE language (atoms, properties, assisted by the
        parameters) to the input files of ABACUS.

        Parameters
        ----------
        profile : AbacusProfile
            The profile used to perform the calculation.
        directory : Path
            The working directory to store the input files.
        atoms : Atoms
            The atoms object to perform the calculation on. Because 
        parameters: dict
            The parameters used to perform the calculation.
        properties: list of str
            The list of properties to calculate
        '''
        # directory
        directory = Path(directory)
        directory.mkdir(exist_ok=True, parents=True)

        # copy the `parameters` because later we will modify it
        parameters = parameters.copy()
        self.stru_file = parameters.get('stru_file', 'STRU') or 'STRU'

        # STRU
        _ = file_safe_backup(directory / self.stru_file)
        # abacuslite itself no longer reads abacus.json; keep a previous run's
        # file out of the way so downstream completion markers cannot see it.
        _ = file_safe_backup(directory / 'abacus.json')
        # group atoms by first-occurrence species order. Keep the reverse map so
        # that we will recover the order in function read_results()
        ind = species_group_indices(atoms.get_chemical_symbols())
        self.atomorder = sorted(range(len(atoms)), key=lambda i: ind[i]) # revmap
        # then we write
        _ = write_stru(atoms[ind], 
                       outdir=directory,
                       pp_file=parameters.get('pseudopotentials'),
                       orb_file=parameters.get('basissets'),
                       fname=self.stru_file)

        # KPT, if needed
        if 'kpts' in parameters:
            _ = file_safe_backup(directory / parameters.get('kpoint_file', 'KPT'))
            _ = write_kpt(parameters['kpts'], 
                          directory / parameters.get('kpoint_file', 'KPT'))
        # should this function be responsible for checking the integrity
        # of information provided by the user? There may be the case that
        # user provides incomplete information, such that the ABACUS cannot
        # run with parameters.

        # INPUT
        # after writing the KPT and STRU, delete them from the parameters
        _ = parameters.pop('kpts', None)

        _ = parameters.pop('pseudopotentials', None)
        parameters.update({'pseudo_dir': profile.pseudo_dir})

        _ = parameters.pop('basissets', None)
        parameters.update({'orbital_dir': profile.orbital_dir})
        # update the parameters respect to the properties desired
        parameters = self.get_property_keywords(parameters, properties)
        # postprocess on the parameters: convert the key and values
        # from any to string. For the case where the value is a 
        # array, convert to the string spaced by whitespace
        for k, v in parameters.items():
            # if the v is iterable, convert to the string spaced by whitespace
            if isinstance(v, (list, tuple, set)):
                parameters[k] = ' '.join(str(i) for i in v)
        dst = directory / self.inputname
        _ = file_safe_backup(dst)
        # remove possible key-value pairs whose value is None
        parameters = {k: v for k, v in parameters.items() if v is not None}

        # Limitation: only the ksdft esolver is supported.
        if parameters.get('esolver_type', 'ksdft') != 'ksdft':
            raise NotImplementedError(
                'ABACUS Lite only supports the ksdft esolver_type presently, '
                'which means the ABACUS should always be used as a DFT '
                'calculator. For other forcefields that ABACUS supports '
                'such as the LJ, DP, etc., please either use the ABACUS '
                'directly, or the implementation of interfaces to ASE '
                'directly.'
            )

        # write the INPUT file to the target directory
        _ = write_input(parameters, dst)

    def execute(self, 
                directory: Path | str, 
                profile: AbacusProfile):
        '''Execute the ABACUS Lite calculation.

        Parameters
        ----------
        directory : Path or str
            The working directory to store the input files.
        profile : AbacusProfile
            The profile used to perform the calculation.

        Raises
        ------
        SubprocessError
            If the ABACUS Lite calculation fails.
        '''
        from subprocess import SubprocessError
        try:
            profile.run(directory=directory, 
                        inputfile=None, 
                        outputfile=self.outputname, 
                        errorfile=self.errorname)
        except SubprocessError:
            # Surface only the tails of each log: abacus.out already contains
            # the SCF iteration table, so full stdout/stderr plus 200 lines of
            # running log would make the exception message unreadably long.
            # The run dir is preserved separately by the test harness
            # (ASE_ABACUS_KEEP_DIR) when post-mortem of the full logs is
            # needed.
            def tail(p: Path, n: int) -> str:
                if not p.exists():
                    return ''
                lines = p.read_text(errors='replace').splitlines()
                header = f'=== {p.name} (last {min(n, len(lines))} lines) ==='
                return header + '\n' + '\n'.join(lines[-n:])

            directory = Path(directory)
            message = ['ABACUS Lite calculation failed']
            for fn, n in ((self.outputname, 200),
                          (self.errorname, 200)):
                chunk = tail(directory / fn, n)
                if chunk:
                    message.append(chunk)
            # The tail of the SCF running log holds the per-step drho
            # history, usually the key evidence for issues like #7794
            # ("Can't find even an electron").
            scf_log = (directory / f'OUT.{self.suffix}'
                       / f'running_{self.calculation}.log')
            chunk = tail(scf_log, 100)
            if chunk:
                message.append(chunk)
            raise SubprocessError('\n'.join(message))

    def read_results(self, directory) -> Dict:
        '''the function that returns the desired properties in dict'''
        if self._legacyio:
            from abacuslite.io.legacyio import read_abacus_out
        else:
            from abacuslite.io.latestio import read_abacus_out

        directory = Path(directory)
        outdir = directory / f'OUT.{self.suffix}'
        log = outdir / f'running_{self.calculation}.log'
        stru_file = getattr(self, 'stru_file', 'STRU') or 'STRU'
        reference = None
        parser_atomorder = self.atomorder
        if self.calculation == 'scf':
            # Read and validate STRU/mapping before passing sort_atoms_with to
            # the parser; otherwise a malformed mapping leaks a bare IndexError.
            reference = _read_stru_reference(
                directory, stru_file, self.atomorder, log_path=log
            )
            # Validate raw parser cardinality before sort_atoms_with can slice
            # away an extra/missing atom and make a corrupt frame look valid.
            _validate_raw_scf_trajectory(log, reference, legacyio=self._legacyio)
            parser_atomorder = None if self.atomorder is None else reference['atomorder']

        # Read all frames: SCF selects the frame matching the current STRU;
        # relax and MD keep the final-frame behavior.
        try:
            frames = read_abacus_out(log, sort_atoms_with=parser_atomorder)
        except (IndexError, AssertionError, KeyError, TypeError, ValueError) as exc:
            if self.calculation != 'scf':
                raise
            raise RuntimeError(
                f"scf running log parser failed: {log} (atomorder/frame "
                f"structure error: {exc}; fail-closed)"
            ) from exc
        if not frames:
            raise RuntimeError(f"no ABACUS running-log frames in {log}")
        if self.calculation != 'scf':
            atoms: Optional[Atoms] = frames[-1]  # relax/MD keep the final frame.
        else:
            atoms = _select_scf_frame_for_structure(
                frames,
                log,
                directory,
                reference['atomorder'],
                stru_file=stru_file,
                _reference=reference,
                legacyio=self._legacyio,
            )
        assert atoms is not None

        # ``SinglePointDFTCalculator.properties()`` exposes the Fermi level
        # under ASE's ``fermi_level`` property, while the CP adapter consumes
        # the backend-neutral ``efermi`` fact.  Keep both names in the result
        # mapping so file-backed ABACUS evaluations do not lose this value at
        # the GenericFileIOCalculator boundary.  Older ASE snapshots store it
        # on ``eFermi``; tolerate that spelling as well.
        # ASE's Properties schema accepts scalar magnetic moments only,
        # although SinglePointDFTCalculator supports vector moments. Project
        # other properties through ASE and retain vectors as result metadata,
        # without changing ASE's process-global output definitions.
        singlepoint = atoms.calc
        vector_moments = {
            name: singlepoint.results[name]
            for name, shape in (('magmom', (3,)), ('magmoms', (len(atoms), 3)))
            if name in singlepoint.results
            and np.shape(singlepoint.results[name]) == shape
        }
        if vector_moments:
            singlepoint = copy(singlepoint)
            singlepoint.results = {
                name: value for name, value in singlepoint.results.items()
                if name not in vector_moments
            }
        results = dict(_schema_compatible_properties(
            singlepoint.properties(), len(atoms)
        ))
        results.update(vector_moments)
        if vector_moments:
            results.setdefault('natoms', len(atoms))
        if 'energy_without_entropy' in atoms.info:
            results['energy_without_entropy'] = atoms.info['energy_without_entropy']
        efermi = results.get("efermi", results.get("fermi_level"))
        if efermi is None:
            for attribute in ("efermi", "eFermi", "fermi_level"):
                candidate = getattr(atoms.calc, attribute, None)
                if candidate is not None:
                    efermi = candidate
                    break
        if efermi is not None:
            results["efermi"] = efermi
            results.setdefault("fermi_level", efermi)
        return results

    def load_profile(self, cfg, **kwargs):
        return AbacusProfile.from_config(cfg, self.name, **kwargs)

class Abacus(GenericFileIOCalculator):
    def __init__(self, 
                 profile=None, 
                 directory='.', 
                 **kwargs):
        '''Construct the ABACUS calculator.

        The keyword arguments (kwargs) can be one of the ASE standard
        keywords: 'xc', 'kpts' or any of ABACUS'
        native keywords.

        Parameters
        ----------
        profile: AbacusProfile
            the interface that interacts with the ABACUS executable.
        directory: str or Path
            the working directory to store the input files.
        pseudopotentials: dict
            A mapping from the element to the pseudopotential file name,
            e.g. ``{'O': 'O_ONCV_PBE-1.0.upf', 'H': 'H.upf'}``.
        baisssets: dict, optional
            A mapping from the element to the ABACUS numerical atomic 
            orbital file name. This is necessary only when it is an
            ABACUS-LCAO (Linear-Combination-of-Atomic-Orbitals) calculation
            e.g. ``{'O': 'O_gga_10au_100Ry_2s2p1d.orb', 
            'H': 'H_gga_10au_100Ry_2s1p.orb'}``.
        kpts: dict
            The k-points sampling should be given as a dict. For there
            are many modes of k-sampling supported, the content may differ
            in cases. A `mode` key should be used to specify the ksampling,
            allowed modes are: `mp-sampling`, `line` and `point`. For 
            `mp-sampling` mode, `gamma-centered`, `nk` and `kshift` should
            present. `gamma-centered` is a boolean, `nk` and `kshift` should
            be lists of three integers. ... TBD
        inp: dict
            parameters setting in INPUT of ABACUS. NOTE: if there are settings
            on the `pseudo_dir` and `orbital_dir`, these will overwrite the
            value in the profile. If you do not expect this, please only use
            the profile, because the profile stands for interfacing with the
            ASE calculator instance with the computational environment.

        **kwargs:
            Other parameters to be passed to the ABACUS calculator.
        '''
        # Fallback profile for callers that do not provide one.
        profile = AbacusProfile('abacus') if profile is None else profile

        # to be compatible with both the legacy and latest format of i/o, the
        # switch is needed. 
        template = AbacusTemplate()
        template._legacyio = _uses_legacy_io(profile.version())

        # Optional ABACUS INPUT parameters, merged with calculator keywords.
        inp = kwargs.pop('inp', {})

        super().__init__(
            template=template,
            profile=profile,
            parameters=kwargs | inp,
            directory=directory,
        )

    def write_input(self, atoms, properties=None, system_changes=None):
        if properties is None:
            properties = self.template.implemented_properties
        self.template.write_input(
            profile=self.profile,
            directory=Path(self.directory),
            atoms=atoms,
            parameters=self.parameters,
            properties=properties,
        )

    def export_properties(self):
        """Export only values accepted by ASE's versioned Properties schema."""
        natoms = len(self.atoms) if self.atoms is not None else 0
        return _schema_compatible_properties(self.results, natoms)

    @classmethod
    def restart(cls, profile=None, directory='.', **kwargs):
        '''instantiate one ABACUS calculator from an existing job directory,
        optionally overwrite some keywords'''
        directory = Path(directory)
        inp_read = read_input(directory / 'INPUT')
        inp = inp_read | kwargs.pop('inp', {})

        pporb_read = read_stru(directory / (inp.get('stru_file', 'STRU') or 'STRU'))['species']
        pseudopotentials = kwargs.get(
            'pseudopotentials',
            {pporb['symbol']: pporb['pp_file'] for pporb in pporb_read}
        )
        if 'pseudopotentials' in kwargs:
            del kwargs['pseudopotentials']
        
        basissets = kwargs.get(
            'basissets',
            {pporb['symbol']: pporb.get('orb_file') for pporb in pporb_read}
        )
        if 'basissets' in kwargs:
            del kwargs['basissets']
        if all([forb is None for forb in basissets.values()]):
            basissets = {}
        assert all([forb is not None for forb in basissets.values()])

        if 'kpts' in kwargs:
            kpts_kwargs = {'kpts': kwargs.pop('kpts')}
        else:
            kpt_file = directory / (inp.get('kpoint_file', 'KPT') or 'KPT')
            if kpt_file.is_file():
                kpts_kwargs = {'kpts': read_kpt(kpt_file)}
            else:
                gamma_only = str(inp.get('gamma_only', '')).strip().lower() in (
                    '1', 'true', 't', '.true.')
                raw_spacing = inp.get('kspacing', '')
                spacing = (raw_spacing if isinstance(raw_spacing, (list, tuple))
                           else str(raw_spacing).split())
                try:
                    implicit_kpts = gamma_only or bool(spacing) and all(
                        float(value) > 0 for value in spacing)
                except ValueError:
                    implicit_kpts = gamma_only
                if not implicit_kpts:
                    read_kpt(kpt_file)  # preserve the missing/invalid KPT error
                kpts_kwargs = {}

        return cls(profile=profile, 
                   directory=directory,
                   pseudopotentials=pseudopotentials,
                   basissets=basissets,
                   **kpts_kwargs,
                   inp=inp,
                   **kwargs)

    def fixed_density(self,
                      kpts: BandPath | Dict[str, str | int | List[float]],
                      symmetry: str = 'off', 
                      profile=None, 
                      **kwargs) -> 'Abacus':
        '''spawn a new ABACUS calculator with fixed density, based on the present
        instance. This funcionality is mostly only useful when perform the 
        non-self-consistent calculations like band structure.
        This interface is referred from the ASE document at:
        https://ase-lib.org/gettingstarted/tut04_bulk/bulk.html#band-structure
        , however, we also note that it is from the implementation of the 
        GPAW python, not the ASE official.
        To make less development burden as possible, we use the same interface
        as the GPAW python.

        Parameters
        ----------
        kpts : BandPath | Dict[str, str | int | List[float]]
            The k-point path to be calculated. Can be either a BandPath object
            or a dictionary that contains the k-point information. For the latter
            case, see tbgen/calculators/abacus/generalio.py::write_kpt for more
            details.
        symmetry : str, optional
            The symmetry mode to be used. Default is 'off'. Now only the `off`
            mode is supported.
        profile : AbacusProfile, optional
            The profile to be used. Default is None. If None, the profile of
            the present instance will be used.
        **kwargs : dict
            Other parameters to be passed to the ABACUS calculator.
        
        Returns
        -------
        Abacus
            The new ABACUS calculator instance that can perform the nscf calculation
            tasks
        '''
        # we should overwrite the 'calculation' to 'nscf', and 'init_chg' to 'file'
        assert symmetry == 'off'
        
        inp = kwargs.setdefault('inp', {})
        inp.update({'calculation': 'nscf',
                    'init_chg': 'file',
                    'symmetry': 0,
                    'out_band': 1,
                    'kspacing': 0.0,       # overwrite
                    'gamma_only': False})
        source_suffix = getattr(self, 'parameters', {}).get('suffix', 'ABACUS')
        inp.setdefault('read_file_dir', f'OUT.{source_suffix}')

        profile = self.profile if profile is None else profile

        # get the kpoint coordinates
        if isinstance(kpts, BandPath):
            kwargs['kpts'] = {
                'mode': 'point',
                'nk': len(kpts.kpts),
                'nkinterpl': np.ones(len(kpts.kpts), dtype=int).tolist(),
                'coordinate': 'direct',
                'kpoints': kpts.kpts.tolist(),
            }
        else:
            assert isinstance(kpts, dict)
            kwargs['kpts'] = kpts
        
        # return
        return Abacus.restart(profile=profile, 
                              directory=self.directory,
                              **kwargs)

    def band_structure(self, efermi=None):
        '''get the band structure from ABACUS. 
        (now not only GPAW can calculate the band structure ;) )'''
        from ase.spectrum.band_structure import get_band_structure
        return get_band_structure(calc=self, reference=efermi)

class AbacusSocketIO(SocketIOCalculator):
    """ASE socket I/O calculator that launches ABACUS as an i-PI client.

    A socket calculator owns one ABACUS process with one fixed INPUT/STRU
    setup. The i-PI protocol can update positions, but electronic-structure
    parameters such as k-points, spin, basis, pseudopotentials, and species
    require a new calculator instance. Energy, forces, and stress are
    independently controlled by ABACUS INPUT. The fixed-layout i-PI response
    uses zero padding for absent fields and an extras metadata record so
    padding is never exposed as a computed property.
    """

    def __init__(self,
                 profile=None,
                 directory='.',
                 port=None,
                 unixsocket=None,
                 timeout=None,
                 log=None,
                 **kwargs):
        inp = dict(kwargs.pop('inp', {}))
        self._property_constraints = {}
        for keyword, property_name in (('cal_force', 'forces'),
                                       ('cal_stress', 'stress')):
            if keyword in inp:
                self._property_constraints[property_name] = self._input_bool(
                    inp[keyword], keyword)
        self.implemented_properties = [
            'energy', 'free_energy', 'forces', 'stress']
        self._active_properties = None
        self._last_socket_metadata = None
        self.last_scf_converged = None
        inp = self._socket_inp(inp)
        self.abacus = Abacus(
            profile=profile,
            directory=directory,
            inp=inp,
            **kwargs,
        )
        self._reference_cell = None
        super().__init__(
            port=port,
            unixsocket=unixsocket,
            timeout=timeout,
            log=log,
            launch_client=self._launch_client,
        )

    def calculate(self, atoms=None, properties=None, system_changes=None):
        from ase.calculators.calculator import (
            PropertyNotImplementedError,
            all_changes,
        )
        from ase.stress import full_3x3_to_voigt_6_stress

        # A failed new request must not expose results or convergence from the
        # previous geometry. Publish only a completely validated response.
        self.results = {}
        self.last_scf_converged = None
        self._last_socket_metadata = None
        if system_changes is None:
            system_changes = all_changes
        if atoms is None:
            atoms = self.atoms
        if atoms is None:
            raise ValueError('AbacusSocketIO.calculate requires atoms')

        requested = self._normalize_socket_properties(properties)
        self._check_requested_properties(requested)

        bad = [change for change in system_changes
               if change not in self.supported_changes]
        if self.atoms is not None and any(bad):
            raise PropertyNotImplementedError(
                'Cannot change {} through IPI protocol. '
                'Please create new socket calculator.'
                .format(bad if len(bad) > 1 else bad[0]))

        desired = set(requested)
        desired.discard('free_energy')
        desired.add('energy')
        for property_name, enabled in self._property_constraints.items():
            if enabled:
                desired.add(property_name)
        active = set(self._active_properties or ())
        if not active:
            active.update(desired)
        elif not desired.issubset(active):
            active.update(desired)
            if self.server is not None:
                self._close_socket_session()
        self._active_properties = tuple(
            name for name in ('energy', 'forces', 'stress') if name in active)

        self._check_fixed_cell(atoms)
        order = self._socket_sort_indices(atoms)
        socket_atoms = atoms[order]
        self.atoms = atoms.copy()

        if self.server is None:
            self.server = self.launch_server()
            proc = self.launch_client(socket_atoms, list(self._active_properties),
                                      port=self._port,
                                      unixsocket=self._unixsocket)
            self.server.proc = proc

        raw_results = self.server.calculate(socket_atoms)
        if not isinstance(raw_results, dict):
            raise ValueError('ABACUS socket server returned a non-mapping result')
        results = dict(raw_results)
        metadata = self._decode_socket_metadata(results.pop('morebytes', None))
        if metadata is None:
            if set(self._active_properties) != {'energy'}:
                raise ValueError(
                    'ABACUS socket response omitted property metadata; refusing '
                    'to infer forces or stress from fixed-wire padding')
            present = {'energy'}
            converged = None
        else:
            present = set(metadata['present'])
            converged = metadata['scf_converged']

        if 'energy' not in present or 'energy' not in results:
            raise ValueError('ABACUS socket response did not provide energy')
        energy = float(results['energy'])
        if not np.isfinite(energy):
            raise ValueError('ABACUS socket energy is not finite')
        free_energy = float(results.get('free_energy', energy))
        if not np.isfinite(free_energy):
            raise ValueError('ABACUS socket free energy is not finite')
        current = {'energy': energy, 'free_energy': free_energy}

        if 'forces' in present:
            if 'forces' not in results:
                raise ValueError(
                    'ABACUS socket metadata advertises forces, but wire response omitted them')
            forces = np.asarray(results['forces'], dtype=np.float64)
            expected_shape = (len(socket_atoms), 3)
            if forces.shape != expected_shape or not np.all(np.isfinite(forces)):
                raise ValueError('ABACUS socket forces have invalid shape or values')
            current['forces'] = self._forces_to_input_order(forces, order)

        if 'stress' in present:
            virial = results.get('virial')
            if virial is None:
                raise ValueError(
                    'ABACUS socket metadata advertises stress, but wire response omitted virial')
            if self.atoms.cell.rank != 3 or not any(self.atoms.pbc):
                raise PropertyNotImplementedError(
                    'ABACUS socket stress requires a periodic rank-3 cell')
            virial = np.asarray(virial, dtype=np.float64)
            if virial.shape != (3, 3) or not np.all(np.isfinite(virial)):
                raise ValueError('ABACUS socket virial is not a finite 3x3 matrix')
            vol = float(atoms.get_volume())
            if not np.isfinite(vol) or vol <= 0.0:
                raise ValueError('ABACUS socket stress requires a positive cell volume')
            current['stress'] = -full_3x3_to_voigt_6_stress(virial) / vol

        missing = [name for name in requested if name not in current]
        if missing:
            raise PropertyNotImplementedError(
                'ABACUS socket response did not provide requested {}'.format(
                    ', '.join(missing)))
        self.results = current
        self._last_socket_metadata = metadata
        self.last_scf_converged = converged

    def _check_fixed_cell(self, atoms):
        from ase.calculators.calculator import PropertyNotImplementedError

        cell = atoms.cell.array.copy()
        if self._reference_cell is None:
            self._reference_cell = cell
            return
        max_delta = np.max(np.abs(cell - self._reference_cell))
        if max_delta > 1.0e-10:
            raise PropertyNotImplementedError(
                'AbacusSocketIO is fixed-cell only; create a new socket '
                'calculator for a changed cell, or use the normal Abacus '
                'FileIO calculator for variable-cell workflows.'
            )

    def set(self, **kwargs):
        if kwargs:
            raise ValueError(
                'AbacusSocketIO input parameters are fixed after construction; '
                'create a new AbacusSocketIO calculator to change k-points, '
                'spin, basis, pseudopotentials, species, or other INPUT/STRU '
                'settings.'
            )
        return super().set(**kwargs)

    def _check_requested_properties(self, requested):
        from ase.calculators.calculator import PropertyNotImplementedError

        constraints = self._property_constraints
        keywords = {'forces': 'cal_force', 'stress': 'cal_stress'}
        for property_name, keyword in keywords.items():
            if property_name in requested and constraints.get(property_name) is False:
                raise PropertyNotImplementedError(
                    '{}=0 disables requested {}'.format(keyword, property_name))

    @staticmethod
    def _normalize_socket_properties(properties):
        from ase.calculators.calculator import PropertyNotImplementedError

        if properties is None:
            names = ['energy']
        elif isinstance(properties, str):
            names = [properties]
        else:
            names = list(properties)
        if not names:
            names = ['energy']
        allowed = {'energy', 'free_energy', 'forces', 'stress'}
        unknown = [name for name in names if name not in allowed]
        if unknown:
            raise PropertyNotImplementedError(
                'ABACUS socket does not implement {}'.format(', '.join(unknown)))
        return tuple(dict.fromkeys(names))

    def _close_socket_session(self):
        server = getattr(self, 'server', None)
        if server is not None:
            close = getattr(server, 'close', None)
            if callable(close):
                close()
        self.server = None
        self.results = {}

    @staticmethod
    def _decode_socket_metadata(raw):
        if raw is None:
            return None
        if isinstance(raw, str):
            payload = raw.encode('utf-8')
        elif isinstance(raw, (bytes, bytearray, memoryview)):
            payload = bytes(raw)
        else:
            payload = np.asarray(raw, dtype=np.uint8).tobytes()
        if not payload:
            return None
        try:
            metadata = json.loads(payload.decode('utf-8'))
        except (UnicodeDecodeError, ValueError) as error:
            raise ValueError('ABACUS socket extras are not valid UTF-8 JSON') from error
        if not isinstance(metadata, dict):
            raise ValueError('ABACUS socket extras must be a JSON object')
        if metadata.get('schema') != SOCKET_METADATA_SCHEMA:
            raise ValueError('unsupported ABACUS socket extras schema')
        present = metadata.get('present')
        if not isinstance(present, list):
            raise ValueError('ABACUS socket extras present must be a list')
        allowed = {'energy', 'forces', 'stress'}
        if any(not isinstance(name, str) or name not in allowed for name in present):
            raise ValueError('ABACUS socket extras contain an unknown property')
        if 'energy' not in present:
            raise ValueError('ABACUS socket extras must include energy')
        scf_converged = metadata.get('scf_converged')
        if not isinstance(scf_converged, bool):
            raise ValueError('ABACUS socket extras scf_converged must be Boolean')
        return {
            'present': tuple(dict.fromkeys(present)),
            'scf_converged': scf_converged,
        }

    def _launch_client(self, atoms, properties=None, port=None, unixsocket=None):
        from subprocess import Popen

        if properties is None:
            properties = list(self._active_properties or ('energy',))
        properties = set(properties)
        properties.discard('free_energy')
        properties.add('energy')
        # The i-PI response has fixed force/virial fields, but ABACUS must be
        # told explicitly which expensive quantities to evaluate.  Keep these
        # switches synchronized with the session mask before writing INPUT.
        self.abacus.parameters['cal_force'] = int('forces' in properties)
        self.abacus.parameters['cal_stress'] = int('stress' in properties)
        properties = [name for name in ('energy', 'forces', 'stress')
                      if name in properties]

        directory = Path(self.abacus.directory)
        directory.mkdir(exist_ok=True, parents=True)

        if hasattr(self.abacus, 'write_inputfiles'):
            self.abacus.write_inputfiles(atoms, properties)
        else:
            self.abacus.write_input(atoms, properties=properties)

        if unixsocket is not None:
            argv = self.abacus.profile.socketio_argv_unix(socket=unixsocket)
        else:
            argv = self.abacus.profile.socketio_argv_inet(port=port)

        stdout = open(directory / self.abacus.template.outputname, 'w')
        stderr = open(directory / self.abacus.template.errorname, 'w')
        try:
            return Popen(argv, cwd=directory, env=os.environ,
                         stdout=stdout, stderr=stderr)
        finally:
            stdout.close()
            stderr.close()

    @staticmethod
    def _socket_inp(inp):
        inp = dict(inp)
        calculation = inp.get('calculation', 'scf')
        if calculation != 'scf':
            raise ValueError('ABACUS socket I/O requires calculation="scf"')
        for keyword in ('cal_force', 'cal_stress'):
            if keyword in inp:
                inp[keyword] = int(AbacusSocketIO._input_bool(inp[keyword], keyword))
        inp.update({
            'calculation': 'scf',
            'socket_driver': 1,
        })
        return inp

    @staticmethod
    def _input_bool(value, name):
        if isinstance(value, bool):
            return value
        if isinstance(value, int) and value in (0, 1):
            return bool(value)
        if isinstance(value, str):
            normalized = value.strip().lower()
            if normalized in ('true', '1'):
                return True
            if normalized in ('false', '0'):
                return False
        raise ValueError('{} must be one of true, false, 1, or 0'.format(name))

    @staticmethod
    def _socket_sort_indices(atoms):
        return species_group_indices(atoms.get_chemical_symbols())

    @staticmethod
    def _forces_to_input_order(forces, order):
        reordered = np.empty_like(forces)
        for sorted_index, original_index in enumerate(order):
            reordered[original_index] = forces[sorted_index]
        return reordered

class TestAbacusCalculator(unittest.TestCase):

    here = Path(__file__).parent
    pporb = here.parent.parent.parent / 'tests' / 'PP_ORB'

    def test_socketio_species_order_mapping(self):
        atoms = Atoms(symbols=['Si', 'O', 'C', 'Si', 'O', 'C'])
        order = AbacusSocketIO._socket_sort_indices(atoms)
        self.assertEqual(order, [0, 3, 1, 4, 2, 5])

        socket_forces = np.arange(18).reshape(6, 3)
        input_forces = AbacusSocketIO._forces_to_input_order(
            socket_forces, order)

        expected = np.empty_like(socket_forces)
        for sorted_index, original_index in enumerate(order):
            expected[original_index] = socket_forces[sorted_index]
        np.testing.assert_array_equal(input_forces, expected)

    def test_socketio_rejects_parameter_changes(self):
        calc = object.__new__(AbacusSocketIO)
        with self.assertRaisesRegex(ValueError, 'fixed after construction'):
            calc.set(kpts={'mode': 'mp-sampling', 'nk': [2, 2, 2]})

    def test_socketio_rejects_cell_changes(self):
        from ase.calculators.calculator import PropertyNotImplementedError

        calc = object.__new__(AbacusSocketIO)
        calc.atoms = Atoms('Si', cell=[5.0, 5.0, 5.0], pbc=True)
        calc._reference_cell = calc.atoms.cell.array.copy()

        changed = calc.atoms.copy()
        changed.cell[0, 0] = 5.1
        with self.assertRaisesRegex(PropertyNotImplementedError, 'fixed-cell'):
            calc._check_fixed_cell(changed)

    def test_socketio_input_keeps_independent_property_switches(self):
        self.assertEqual(
            AbacusSocketIO._socket_inp({'cal_force': 0, 'cal_stress': 1}),
            {'calculation': 'scf', 'socket_driver': 1,
             'cal_force': 0, 'cal_stress': 1})

    def test_socketio_boolean_parser_rejects_ambiguous_values(self):
        for value in ('yes', 'no', 'on', 'off', ''):
            with self.assertRaisesRegex(ValueError, 'cal_force'):
                AbacusSocketIO._input_bool(value, 'cal_force')

    def test_socketio_metadata_rejects_unknown_property(self):
        metadata = json.dumps({
            'schema': 'abacus.socket.properties.v1',
            'present': ['energy', 'charges'],
            'scf_converged': True,
        }).encode('utf-8')
        with self.assertRaisesRegex(ValueError, 'unknown property'):
            AbacusSocketIO._decode_socket_metadata(
                np.frombuffer(metadata, dtype=np.int8))

    def test_socketio_metadata_marks_padding_absent(self):
        metadata = json.dumps({
            'schema': 'abacus.socket.properties.v1',
            'present': ['energy'],
            'scf_converged': False,
        }).encode('utf-8')
        decoded = AbacusSocketIO._decode_socket_metadata(
            np.frombuffer(metadata, dtype=np.int8))
        self.assertEqual(decoded['present'], ('energy',))
        self.assertFalse(decoded['scf_converged'])

    def test_socketio_legacy_response_cannot_infer_force_from_padding(self):
        class LegacyServer:
            def calculate(self, atoms):
                return {
                    'energy': 1.0,
                    'forces': np.zeros((len(atoms), 3)),
                    'virial': np.zeros((3, 3)),
                    'morebytes': b'',
                }

        calc = object.__new__(AbacusSocketIO)
        calc.variable_cell = False
        calc._property_constraints = {}
        calc._active_properties = ('energy', 'forces')
        calc._reference_cell = None
        calc.atoms = None
        calc.server = LegacyServer()
        atoms = Atoms('Si')
        with self.assertRaisesRegex(ValueError, 'refusing to infer forces'):
            calc.calculate(atoms=atoms, properties=('forces',), system_changes=())

    def test_socketio_failed_response_clears_previous_frame(self):
        class BrokenServer:
            def calculate(self, atoms):
                raise EOFError('incomplete frame')

        calc = object.__new__(AbacusSocketIO)
        calc._property_constraints = {}
        calc._active_properties = ('energy',)
        calc._reference_cell = None
        calc.atoms = None
        calc.results = {'energy': 123.0}
        calc.last_scf_converged = True
        calc._last_socket_metadata = {'scf_converged': True}
        calc.server = BrokenServer()
        with self.assertRaises(EOFError):
            calc.calculate(Atoms('Si'), properties=('energy',), system_changes=())
        self.assertEqual(calc.results, {})
        self.assertIsNone(calc.last_scf_converged)
        self.assertIsNone(calc._last_socket_metadata)

    def test_socketio_requested_disabled_property_is_rejected(self):
        calc = object.__new__(AbacusSocketIO)
        calc._property_constraints = {'forces': False, 'stress': False}
        from ase.calculators.calculator import PropertyNotImplementedError
        with self.assertRaises(PropertyNotImplementedError):
            calc._check_requested_properties(('forces',))

    def test_parse_version_allows_launcher_noise(self):
        stdout = 'launcher info\nABACUS version v3.11.0-beta6\n'
        self.assertEqual(AbacusProfile.parse_version(stdout), 'v3.11.0-beta6')

    def test_parse_version_rejects_missing_version(self):
        with self.assertRaisesRegex(RuntimeError, 'ABACUS version'):
            AbacusProfile.parse_version('launcher failed before abacus started')

    def test_calculator_results(self):
        from ase.build.bulk import bulk
        silicon = bulk('Si', crystalstructure='diamond', a=5.43)
        aprof = AbacusProfile(
            command='mpirun -np 2 abacus',
            pseudo_dir=self.pporb,
            orbital_dir=self.pporb,
            omp_num_threads=1
        )
        with tempfile.TemporaryDirectory() as tmpdir:
            calculator = Abacus(aprof,
                                directory=tmpdir,
                                pseudopotentials={'Si': 'Si_ONCV_PBE-1.0.upf'},
                                basissets={'Si': 'Si_gga_6au_100Ry_2s2p1d.orb'},
                                inp={'calculation': 'scf',
                                    'basis_type': 'lcao',
                                    'ks_solver': 'genelpa',
                                    'ecutwfc': 40,
                                    'symmetry': 1,
                                    'nspin': 1,
                                    'gamma_only': True,
                                    'cal_force': 1,
                                    'cal_stress': 1})
            silicon.calc = calculator
            e = silicon.get_potential_energy()
        
        # check!
        self.assertAlmostEqual(e, -194.953050937)
        self.assertIsNotNone(calculator.results)
        self.assertIsInstance(calculator.results, dict)
        for k in ['nspins', 'nkpts', 'nbands', 'eigenvalues', 'occupations',
                  'fermi_level', 'kpoint_weights', 'ibz_kpoints', 'energy', 
                  'free_energy', 'natoms', 'forces', 'stress', 'magmoms']:
            self.assertIn(k, calculator.results)
        self.assertEqual(calculator.results['nspins'], 1)
        self.assertEqual(calculator.results['nkpts'], 1)
        self.assertEqual(calculator.results['nbands'], 14)
        self.assertEqual(calculator.results['energy'], e)
        self.assertEqual(calculator.results['free_energy'], e)
        self.assertEqual(calculator.results['natoms'], 2)
        
        for k in ['eigenvalues', 'occupations', 'ibz_kpoints', 'forces', 'stress', 'magmoms']:
            self.assertIsInstance(calculator.results[k], np.ndarray)

        self.assertEqual(calculator.results['eigenvalues'].shape, (1, 1, 14))
        ekb = [-4.82194,  7.62727,  7.62727,  7.62737, 10.2436 , 10.2436 ,
                10.2436 , 10.9884 , 16.057  , 16.057  , 23.8353 , 25.421  ,
                25.421  , 25.4212 ]
        self.assertTrue(np.allclose(calculator.results['eigenvalues'][0, 0, :], np.array(ekb)))

        self.assertEqual(calculator.results['occupations'].shape, (1, 1, 14))
        occ = [2., 2., 2., 2., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.]
        self.assertTrue(np.allclose(calculator.results['occupations'][0, 0, :], np.array(occ)))
        
        self.assertEqual(calculator.results['ibz_kpoints'].shape, (1, 3))
        self.assertTrue(np.allclose(calculator.results['ibz_kpoints'][0, :], np.array([0,0,0])))

        self.assertEqual(calculator.results['forces'].shape, (2, 3))
        self.assertTrue(np.allclose(calculator.results['forces'], np.zeros((2, 3))))

        self.assertEqual(calculator.results['stress'].shape, (6,))
        stress = [-0.19327923, -0.19327923, -0.19327923, -0.        ,  0.        ,   0.        ]
        self.assertTrue(np.allclose(calculator.results['stress'], np.array(stress)))
        
        self.assertEqual(calculator.results['magmoms'].shape, (2,))
        self.assertTrue(np.allclose(calculator.results['magmoms'], np.zeros(2)))

    def test_restart(self):
        from ase.build.bulk import bulk
        silicon = bulk('Si', crystalstructure='diamond', a=5.43)
        aprof = AbacusProfile(
            command='mpirun -np 2 abacus',
            pseudo_dir=self.pporb,
            orbital_dir=self.pporb,
            omp_num_threads=1
        )
        with tempfile.TemporaryDirectory() as tmpdir:
            calculator = Abacus(aprof,
                                directory=tmpdir,
                                pseudopotentials={'Si': 'Si_ONCV_PBE-1.0.upf'},
                                basissets={'Si': 'Si_gga_6au_100Ry_2s2p1d.orb'},
                                inp={'calculation': 'scf',
                                    'basis_type': 'lcao',
                                    'ks_solver': 'genelpa',
                                    'ecutwfc': 40,
                                    'symmetry': 1,
                                    'nspin': 1,
                                    'gamma_only': True,
                                    'cal_force': 1,
                                    'cal_stress': 1})
            silicon.calc = calculator
            e = silicon.get_potential_energy()
        
            # restart
            silicon.calc = Abacus.restart(aprof, directory=tmpdir)
            e2 = silicon.get_potential_energy()
            self.assertAlmostEqual(e2, e)

    def test_version_number_check(self):
        # not a version number
        with self.assertRaises(AssertionError):
            switch_io_backend_version('not-a-version-number')
        # too old version
        with self.assertRaises(AssertionError):
            switch_io_backend_version('v2.2.2')
        self.assertTrue(switch_io_backend_version('v3.8.4'))
        self.assertFalse(switch_io_backend_version('v3.9.0.25'))
        self.assertTrue(switch_io_backend_version('v3.10.0'))
        self.assertFalse(switch_io_backend_version('v3.11.0-beta.2'))
        self.assertFalse(switch_io_backend_version('v3.11.0'))

    def test_property_keywords_reject_conflicting_user_parameters(self):
        template = AbacusTemplate()
        with self.assertRaises(ValueError):
            template.get_property_keywords({'nspin': 1}, ['magmom'])

        parameters = template.get_property_keywords({'nspin': 2}, ['magmom'])
        self.assertEqual(str(parameters['nspin']), '2')

    def test_property_keywords_accept_equivalent_boolean_user_parameters(self):
        template = AbacusTemplate()

        parameters = template.get_property_keywords(
            {'cal_force': True, 'cal_stress': True},
            ['forces', 'stress']
        )

        self.assertEqual(str(parameters['cal_force']), '1')
        self.assertEqual(str(parameters['cal_stress']), '1')

    def test_property_keywords_reject_conflicting_boolean_user_parameters(self):
        template = AbacusTemplate()

        with self.assertRaises(ValueError):
            template.get_property_keywords({'cal_force': False}, ['forces'])

        with self.assertRaises(ValueError):
            template.get_property_keywords({'cal_stress': False}, ['stress'])

    def test_property_keywords_treat_string_values_as_scalars(self):
        template = AbacusTemplate()
        template.implemented_properties = ['probe']
        template.get_probe_keywords = lambda parameters: {'custom_switch': 'true'}

        parameters = template.get_property_keywords(
            {'custom_switch': 'true'}, ['probe']
        )

        self.assertEqual(parameters['custom_switch'], 'true')

    def test_property_keywords_compare_iterables_like_input_writer(self):
        template = AbacusTemplate()
        template.implemented_properties = ['probe']
        template.get_probe_keywords = lambda parameters: {
            'custom_vector': [1, 'true']
        }

        with self.assertRaises(ValueError):
            template.get_property_keywords(
                {'custom_vector': [True, 'true']}, ['probe']
            )

    def test_property_keywords_reject_conflicting_properties(self):
        template = AbacusTemplate()
        template.implemented_properties = ['prop_a', 'prop_b']
        template.get_prop_a_keywords = lambda parameters: {'calculation': 'scf'}
        template.get_prop_b_keywords = lambda parameters: {'calculation': 'md'}

        with self.assertRaises(ValueError):
            template.get_property_keywords({}, ['prop_a', 'prop_b'])

    def test_dipole_property_is_not_implemented(self):
        template = AbacusTemplate()
        self.assertNotIn('dipole', template.implemented_properties)
        with self.assertRaises(AssertionError):
            template.get_property_keywords({}, ['dipole'])

if __name__ == '__main__':
    unittest.main()
