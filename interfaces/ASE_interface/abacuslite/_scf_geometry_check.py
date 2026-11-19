"""Check single-point SCF geometry against the written STRU."""

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


def _read_stru_reference(stru_path):
    """Return STRU geometry in file order and Cartesian Angstrom."""
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
        coordinates = np.asarray(raw_coords, dtype=float)
        if coordinates.shape != (natoms, 3):
            raise ValueError(
                f"STRU positions shape {coordinates.shape}, expected ({natoms}, 3)"
            )
        if not np.isfinite(coordinates).all():
            raise ValueError("non-finite STRU positions")

        if 'vec' not in lat:
            raise ValueError("SCF geometry checking requires explicit LATTICE_VECTORS")
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
        raise RuntimeError(
            f"Invalid SCF STRU geometry ({stru_path}; "
            f"{exc}; fail-closed)"
        ) from exc

    return {
        'positions': species_positions,
        'cell': cell,
        'lattice_vectors': lattice_vectors,
        'alat': lattice_constant * Bohr,
        # ABACUS prints the product with its own Bohr-to-Angstrom constant.
        'log_alat': lattice_constant * 0.5291770,
        'symbols': raw_symbols,
    }


def _check_scf_geometry(stru_path, log_path, atomorder, read_log_geometry):
    """Check one SCF geometry in ABACUS order and return the ASE permutation."""
    expected = _read_stru_reference(stru_path)
    mapping = _validate_atomorder(atomorder, len(expected['symbols']), stru_path)
    try:
        records = read_log_geometry(log_path)
        if len(records) != 1:
            raise ValueError(f"expected exactly one geometry record, found {len(records)}")
        _check_log_geometry(records[0], expected)
    except (AssertionError, KeyError, OSError, TypeError, ValueError, IndexError) as exc:
        raise RuntimeError(
            f"Invalid SCF geometry (log: {log_path}; STRU: {stru_path}; {exc})"
        ) from exc
    return mapping if atomorder is not None else None


def _check_log_geometry(record, expected):
    """Compare one raw geometry record before atom reordering can hide errors."""
    scale = float(record['alat_in_angstrom'])
    if not np.isfinite(scale) or scale <= 0:
        raise ValueError(f"invalid lattice constant in Angstrom: {scale!r}")
    positions = np.asarray(record['coords'], dtype=float)
    cell = np.asarray(record['cell'], dtype=float)
    natoms = len(expected['symbols'])
    if positions.ndim >= 1 and positions.shape[0] != natoms:
        raise ValueError(f"atom count {positions.shape[0]}, expected {natoms}")
    if positions.shape != (natoms, 3) or not np.isfinite(positions).all():
        raise ValueError(f"expected finite positions of shape ({natoms}, 3)")
    symbols = list(record['elem'])
    if symbols != expected['symbols']:
        raise ValueError(f"atom symbols {symbols!r} do not match STRU {expected['symbols']!r}")
    if cell.shape != (3, 3) or not np.isfinite(cell).all():
        raise ValueError("expected a finite (3, 3) cell")

    # ABACUS prints the scale and the dimensionless matrix separately.
    # A single Cartesian tolerance on their product rejects valid large cells.
    if not _matches_log_precision(scale, expected['log_alat']):
        raise ValueError("lattice constant in log differs from STRU")
    if not _matches_log_precision(cell / scale, expected['lattice_vectors']):
        raise ValueError("lattice vectors in log differ from STRU")

    cell = expected['cell']
    if record['coordinate'] == 'Direct':
        positions = positions @ cell
    elif record['coordinate'] == 'Cartesian':
        positions = positions * expected['alat']
    else:
        raise ValueError(f"Unknown coordinate system: {record['coordinate']}")
    cell_inverse = _invert_cell(cell, 'STRU')
    delta_frac = (positions - expected['positions']) @ cell_inverse
    periodic_cart = (delta_frac - np.rint(delta_frac)) @ cell
    if not np.isfinite(periodic_cart).all():
        raise ValueError("non-finite periodic position residual")
    residual = float(np.max(np.abs(periodic_cart)))
    atol = 1e-4
    if residual > atol:
        raise ValueError(f"positions differ from STRU by {residual:.6g} Angstrom "
                         f"after periodic wrapping (atol={atol})")


def _matches_log_precision(actual, expected):
    """Compare using the six significant digits of ABACUS's default log."""
    magnitude = np.abs(expected)
    nonzero_magnitude = np.where(magnitude > 0, magnitude, 1.)
    half_unit = 0.5 * 10. ** (np.floor(np.log10(nonzero_magnitude)) - 5)
    tolerance = np.where(magnitude > 0, half_unit + 1e-12 * magnitude, 1e-12)
    return bool(np.all(np.abs(actual - expected) <= tolerance))
