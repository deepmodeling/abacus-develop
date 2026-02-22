from __future__ import annotations

"""
Tests for memory lifetime safety of pyabacus bindings.

These tests verify that numpy arrays returned from C++ objects
are safe to use even after the C++ object is destroyed (because
they are copies, not views).
"""

import pytest
import numpy as np
from pyabacus import ModuleNAO as nao


def test_array_lifetime_after_object_deletion():
    """Test that arrays remain valid after C++ object deletion."""
    chi = nao.NumericalRadial()

    # Build the object
    sz = 100
    dr = 0.1
    grid = np.array([i * dr for i in range(sz)], dtype=np.float64)
    values = np.exp(-grid)

    chi.build(
        l=0,
        for_r_space=True,
        ngrid=sz,
        grid=grid,
        value=values,
        p=-1,
        izeta=0,
        symbol="Test",
        itype=0,
        init_sbt=True
    )

    # Get arrays
    rgrid = chi.rgrid
    rvalue = chi.rvalue

    # Save copies for comparison
    rgrid_copy = rgrid.copy()
    rvalue_copy = rvalue.copy()

    # Delete the C++ object
    del chi

    # Arrays should still be valid and contain correct data
    # (because they are copies, not views)
    np.testing.assert_array_equal(rgrid, rgrid_copy)
    np.testing.assert_array_equal(rvalue, rvalue_copy)


def test_array_modification_isolation():
    """Test that modifying returned arrays doesn't affect original data."""
    chi = nao.NumericalRadial()

    sz = 100
    dr = 0.1
    grid = np.array([i * dr for i in range(sz)], dtype=np.float64)
    values = np.exp(-grid)

    chi.build(
        l=0,
        for_r_space=True,
        ngrid=sz,
        grid=grid,
        value=values
    )

    # Get first array and modify it
    rgrid1 = chi.rgrid
    original_value = rgrid1[0]
    rgrid1[0] = 999.0

    # Get second array - should have original value
    rgrid2 = chi.rgrid
    assert rgrid2[0] == original_value, "Modification should not affect original data"
    assert rgrid1[0] == 999.0, "Modified array should retain modification"


def test_multiple_array_accesses():
    """Test that multiple accesses return independent copies."""
    chi = nao.NumericalRadial()

    sz = 50
    dr = 0.2
    grid = np.array([i * dr for i in range(sz)], dtype=np.float64)
    values = np.exp(-grid)

    chi.build(
        l=1,
        for_r_space=True,
        ngrid=sz,
        grid=grid,
        value=values
    )

    # Get multiple arrays
    arrays = [chi.rgrid for _ in range(5)]

    # Modify each differently
    for i, arr in enumerate(arrays):
        arr[0] = float(i)

    # Verify each has its own modification
    for i, arr in enumerate(arrays):
        assert arr[0] == float(i), f"Array {i} should have value {i}"


def test_radial_collection_array_lifetime():
    """Test array lifetime for RadialCollection objects."""
    orb_dir = '../../../tests/PP_ORB/'
    file_list = [orb_dir + "C_gga_8au_100Ry_2s2p1d.orb"]

    try:
        orb = nao.RadialCollection()
        orb.build(1, file_list, 'o')

        # Get a NumericalRadial from the collection
        nr = orb(0, 0, 0)

        # Get arrays from it
        rgrid = nr.rgrid
        rvalue = nr.rvalue

        # Save copies
        rgrid_copy = rgrid.copy()
        rvalue_copy = rvalue.copy()

        # Delete the collection (which owns the NumericalRadial)
        del orb

        # Arrays should still be valid
        np.testing.assert_array_equal(rgrid, rgrid_copy)
        np.testing.assert_array_equal(rvalue, rvalue_copy)
    except (FileNotFoundError, RuntimeError):
        pytest.skip("Orbital files not available for testing")


def test_empty_array_handling():
    """Test handling of arrays when k-space is not initialized."""
    chi = nao.NumericalRadial()

    sz = 50
    dr = 0.2
    grid = np.array([i * dr for i in range(sz)], dtype=np.float64)
    values = np.exp(-grid)

    chi.build(
        l=0,
        for_r_space=True,
        ngrid=sz,
        grid=grid,
        value=values,
        init_sbt=False  # Don't initialize SBT, so k-space won't be set
    )

    # r-space should be valid
    assert chi.nr == sz
    rgrid = chi.rgrid
    assert len(rgrid) == sz

    # k-space should be empty (nk == 0)
    assert chi.nk == 0
    kgrid = chi.kgrid
    assert len(kgrid) == 0
