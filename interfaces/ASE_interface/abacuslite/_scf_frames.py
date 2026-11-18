"""Validate and match SCF output against the input structure."""

from pathlib import Path

import numpy as np
from ase.units import Bohr

from abacuslite.io.generalio import _read_stru_geometry


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

    # A permutation needs one non-boolean integer per atom, with every
    # index present exactly once; sorting checks both range and uniqueness.
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
    stru_path = Path(directory) / (stru_file or 'STRU')
    try:
        stru = _read_stru_geometry(stru_path)
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

        if coord_type == 'direct':
            species_positions = coordinates @ cell
        elif coord_type == 'cartesian':
            species_positions = coordinates * lattice_constant * Bohr
        elif coord_type == 'cartesian_au':
            species_positions = coordinates * Bohr
        elif coord_type == 'cartesian_angstrom':
            species_positions = coordinates
        elif coord_type in (
            'cartesian_angstrom_center_xy',
            'cartesian_angstrom_center_xz',
            'cartesian_angstrom_center_yz',
            'cartesian_angstrom_center_xyz',
        ):
            # The engine masks Cartesian components of the cell center.
            axes = coord_type.rsplit('_', 1)[1]
            center = cell.sum(axis=0) / 2.0
            species_positions = coordinates + center * [axis in axes for axis in 'xyz']
        else:
            raise ValueError(f"unsupported STRU coordinate type {stru['coord_type']!r}")
        if not np.isfinite(species_positions).all():
            raise ValueError("non-finite normalized STRU positions")
    except (
        AssertionError,
        KeyError,
        OSError,
        TypeError,
        ValueError,
        IndexError,
        np.linalg.LinAlgError,
    ) as exc:
        suffix = f"; log: {log_path}" if log_path is not None else ""
        raise RuntimeError(
            f"Invalid scf STRU identity data (STRU: {stru_path}{suffix}; "
            f"{exc}; fail-closed)"
        ) from exc

    return {
        'positions': species_positions[mapping],
        'cell': cell,
        'symbols': [raw_symbols[index] for index in mapping],
        'atomorder': mapping,
    }


def _validate_raw_scf_frame(frame, frame_index, expected):
    """Validate one raw frame before atom reordering can hide extra atoms."""
    context = f"frame {frame_index}"
    scale = float(frame['alat_in_angstrom'])
    if not np.isfinite(scale) or scale <= 0:
        raise ValueError(f"{context} invalid lattice constant in Angstrom: {scale!r}")
    expected_shape = expected['positions'].shape
    expected_natoms = len(expected['symbols'])
    coords = np.asarray(frame['coords'], dtype=float)
    if coords.ndim >= 1 and coords.shape[0] != expected_natoms:
        raise ValueError(
            f"{context} atom count {coords.shape[0]}, expected {expected_natoms}"
        )
    if coords.shape != expected_shape:
        raise ValueError(
            f"{context} positions shape {coords.shape}, expected {expected_shape}"
        )
    if not np.isfinite(coords).all():
        raise ValueError(f"{context} positions are non-finite")

    elements = np.asarray(frame['elem'])
    if elements.shape != (expected_natoms,):
        raise ValueError(
            f"{context} atom symbols shape {elements.shape}, "
            f"expected ({expected_natoms},)"
        )
    cell = np.asarray(frame['cell'], dtype=float)
    if cell.shape != (3, 3):
        raise ValueError(f"{context} cell shape {cell.shape}, expected (3, 3)")
    if not np.isfinite(cell).all():
        raise ValueError(f"{context} cell is non-finite")
    _invert_cell(cell, context)


def _validate_raw_scf_trajectory(log_path, expected, read_raw_trajectory):
    """Reject malformed raw frames before atom reordering can hide them."""
    try:
        trajectory = read_raw_trajectory(log_path)
        if not trajectory:
            raise ValueError("no raw trajectory frames")
        for frame_index, frame in enumerate(trajectory):
            _validate_raw_scf_frame(frame, frame_index, expected)
    except (AssertionError, KeyError, OSError, TypeError, ValueError, IndexError) as exc:
        raise RuntimeError(
            f"Invalid raw scf trajectory (log: {log_path}; {exc}; fail-closed)"
        ) from exc


def _measure_scf_frame_match(frame, frame_index, expected, atol):
    """Validate one ASE frame and return its structure-match diagnostics.

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

    # Public readers return standard ASE Cartesian positions in Angstrom.
    cell_diff = float(np.max(np.abs(cell - expected['cell'])))
    periodic_residual = None
    if np.allclose(cell, expected['cell'], atol=atol, rtol=0.0):
        delta_frac = (positions - expected['positions']) @ cell_inverse
        periodic_cart = (delta_frac - np.rint(delta_frac)) @ cell
        if not np.isfinite(periodic_cart).all():
            raise ValueError(f"{context} periodic residual is non-finite")
        periodic_residual = float(np.max(np.abs(periodic_cart)))
    return {
        'frame_cart': positions,
        'cell_diff': cell_diff,
        'periodic_residual': periodic_residual,
        'matches': periodic_residual is not None and periodic_residual <= atol,
    }


def _select_scf_frame_for_structure(
    frames, log_path, directory, atomorder, atol=1e-4, stru_file='STRU', _reference=None
):
    """Select the latest frame matching STRU in ASE atom order.

    Validate every candidate before returning, including those after a match.
    Compare cell components and periodic Cartesian residuals with absolute
    ``atol`` in Angstrom and ``rtol=0``. The 1e-4 default accommodates ordinary
    log rounding; it is a structure identity bound, not a convergence setting.
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
    selected = None
    last = None
    invalid = []
    for frame_index, frame in enumerate(frames):
        try:
            last = _measure_scf_frame_match(
                frame, frame_index, expected, atol
            )
            if last['matches']:
                selected = frame
        except (AttributeError, TypeError, ValueError, IndexError) as exc:
            invalid.append(f"frame {frame_index}: {exc}")

    if invalid:
        details = "; ".join(invalid)
        raise RuntimeError(
            f"Malformed scf frames in log {log_path} ({nframes} {frame_label}; "
            f"{details}; atol={atol}; fail-closed)"
        )
    if selected is not None:
        return selected

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
