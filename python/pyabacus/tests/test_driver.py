"""
Unit tests for PyABACUS driver module.

This module contains tests for the abacus() function and related classes.
"""

import pytest
import numpy as np
from pathlib import Path
from unittest.mock import Mock, patch, MagicMock
import tempfile
import os


class TestCalculationResult:
    """Tests for CalculationResult dataclass."""

    def test_default_values(self):
        """Test default initialization."""
        from pyabacus.driver.runner import CalculationResult

        result = CalculationResult()
        assert result.converged is False
        assert result.niter == 0
        assert result.etot == 0.0
        assert result.forces is None
        assert result.stress is None

    def test_energy_conversion(self):
        """Test energy unit conversion."""
        from pyabacus.driver.runner import CalculationResult

        result = CalculationResult(etot=-10.0)
        # etot is now in eV, etot_ev is the same value for compatibility
        assert result.etot_ev == pytest.approx(-10.0)

    def test_energies_dict(self):
        """Test energies dictionary property."""
        from pyabacus.driver.runner import CalculationResult

        result = CalculationResult(
            etot=-10.0,
            eband=-5.0,
            hartree_energy=2.0,
            etxc=-3.0,
            ewald_energy=-4.0,
        )

        energies = result.energies
        assert 'etot' in energies
        assert energies['etot'] == -10.0  # All in eV now
        assert energies['eband'] == -5.0

    def test_forces_conversion(self):
        """Test force unit conversion."""
        from pyabacus.driver.runner import CalculationResult

        # Forces are now stored in eV/Angstrom directly
        forces = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
        result = CalculationResult(forces=forces)

        # forces_ev_ang is now the same as forces (for compatibility)
        forces_ev_ang = result.forces_ev_ang
        assert forces_ev_ang is not None
        assert forces_ev_ang[0, 0] == pytest.approx(1.0)

    def test_forces_none(self):
        """Test forces_ev_ang when forces is None."""
        from pyabacus.driver.runner import CalculationResult

        result = CalculationResult()
        assert result.forces_ev_ang is None
        assert result.has_forces is False

    def test_has_forces(self):
        """Test has_forces property."""
        from pyabacus.driver.runner import CalculationResult

        result = CalculationResult()
        assert result.has_forces is False

        result.forces = np.array([[1.0, 0.0, 0.0]])
        assert result.has_forces is True

    def test_has_stress(self):
        """Test has_stress property."""
        from pyabacus.driver.runner import CalculationResult

        result = CalculationResult()
        assert result.has_stress is False

        result.stress = np.eye(3)
        assert result.has_stress is True

    def test_summary(self):
        """Test summary string generation."""
        from pyabacus.driver.runner import CalculationResult

        result = CalculationResult(
            converged=True,
            niter=10,
            drho=1e-8,
            etot=-10.0,
            nat=2,
            ntype=1,
            nbands=8,
            nks=4,
        )

        summary = result.summary()
        assert "Converged: Yes" in summary
        assert "SCF iterations: 10" in summary
        assert "Total energy:" in summary
        assert "Atoms: 2" in summary

    def test_repr(self):
        """Test string representation."""
        from pyabacus.driver.runner import CalculationResult

        result = CalculationResult(converged=True, etot=-10.0)
        repr_str = repr(result)
        assert "converged=True" in repr_str
        assert "etot=" in repr_str


class TestAbacusFunction:
    """Tests for the abacus() function."""

    def test_file_not_found_error_nonexistent_dir(self):
        """Test FileNotFoundError for missing input directory."""
        from pyabacus.driver.runner import abacus

        with pytest.raises(FileNotFoundError, match="Input directory not found"):
            abacus("./nonexistent/")

    def test_file_not_found_error(self):
        """Test FileNotFoundError for missing input directory."""
        from pyabacus.driver.runner import abacus

        with pytest.raises(FileNotFoundError, match="Input directory not found"):
            abacus("/nonexistent/path/that/does/not/exist/")

    def test_input_file_not_found(self):
        """Test FileNotFoundError for missing INPUT file."""
        from pyabacus.driver.runner import abacus

        with tempfile.TemporaryDirectory() as tmpdir:
            mock_driver_module = MagicMock()
            with patch.dict('sys.modules', {'pyabacus.driver._driver_pack': mock_driver_module}):
                with pytest.raises(FileNotFoundError, match="INPUT file not found"):
                    abacus(tmpdir)

    def test_default_input_dir(self):
        """Test that default input_dir is current directory."""
        from pyabacus.driver.runner import abacus

        # Create a mock driver
        mock_driver = MagicMock()
        mock_result = MagicMock()
        mock_result.converged = True
        mock_result.niter = 10
        mock_result.drho = 1e-8
        mock_result.etot = -10.0
        mock_result.eband = -5.0
        mock_result.hartree_energy = 2.0
        mock_result.etxc = -3.0
        mock_result.ewald_energy = -4.0
        mock_result.demet = 0.0
        mock_result.exx = 0.0
        mock_result.evdw = 0.0
        mock_result.fermi_energy = 5.0
        mock_result.bandgap = 1.0
        mock_result.nat = 2
        mock_result.ntype = 1
        mock_result.nbands = 8
        mock_result.nks = 4
        mock_result.has_forces = False
        mock_result.has_stress = False
        mock_driver.run.return_value = mock_result

        mock_driver_class = MagicMock(return_value=mock_driver)
        mock_module = MagicMock()
        mock_module.PyDriver = mock_driver_class

        with tempfile.TemporaryDirectory() as tmpdir:
            # Create INPUT file
            input_file = Path(tmpdir) / "INPUT"
            input_file.write_text("INPUT_PARAMETERS\ncalculation scf\n")

            with patch.dict('sys.modules', {'pyabacus.driver._driver_pack': mock_module}):
                # Change to temp directory and run
                old_cwd = os.getcwd()
                try:
                    os.chdir(tmpdir)
                    result = abacus()
                    assert result.converged is True
                finally:
                    os.chdir(old_cwd)


class TestRunScf:
    """Tests for run_scf convenience function."""

    def test_run_scf_calls_abacus(self):
        """Test that run_scf calls abacus with correct arguments."""
        from pyabacus.driver.runner import run_scf, abacus

        with patch('pyabacus.driver.runner.abacus') as mock_abacus:
            mock_result = MagicMock()
            mock_abacus.return_value = mock_result

            result = run_scf("./test_dir/", verbosity=0)

            mock_abacus.assert_called_once_with("./test_dir/", verbosity=0)


class TestRunRelax:
    """Tests for run_relax convenience function."""

    def test_run_relax_enables_forces(self):
        """Test that run_relax enables force calculation by default."""
        from pyabacus.driver.runner import run_relax

        with patch('pyabacus.driver.runner.abacus') as mock_abacus:
            mock_result = MagicMock()
            mock_abacus.return_value = mock_result

            result = run_relax("./test_dir/")

            # Check that calculate_force=True was passed
            call_kwargs = mock_abacus.call_args[1]
            assert call_kwargs.get('calculate_force', False) is True


class TestDriverModule:
    """Tests for the driver module initialization."""

    def test_module_exports(self):
        """Test that module exports expected symbols."""
        from pyabacus import driver

        assert hasattr(driver, 'abacus')
        assert hasattr(driver, 'CalculationResult')

    def test_pyabacus_exports_abacus(self):
        """Test that pyabacus exports abacus function."""
        import pyabacus

        # The abacus function should be accessible
        assert 'abacus' in pyabacus.__all__
        assert 'CalculationResult' in pyabacus.__all__


class TestIntegration:
    """Integration tests (require actual ABACUS installation)."""

    @pytest.mark.skip(reason="Requires ABACUS installation and test data")
    def test_si_scf(self):
        """Test SCF calculation on Si."""
        from pyabacus import abacus

        result = abacus("./tests/integrate/101_PW_15_f_pseudopots/")

        assert result.converged
        assert result.etot < 0  # Energy should be negative
        assert result.nat > 0

    @pytest.mark.skip(reason="Requires ABACUS installation and test data")
    def test_si_forces(self):
        """Test force calculation on Si."""
        from pyabacus import abacus

        result = abacus(
            "./tests/integrate/101_PW_15_f_pseudopots/",
            calculate_force=True,
        )

        assert result.converged
        assert result.has_forces
        assert result.forces.shape == (result.nat, 3)


class TestCalculationResultOutputTracking:
    """Tests for CalculationResult output file tracking."""

    def test_default_output_values(self):
        """Test default output tracking values."""
        from pyabacus.driver.runner import CalculationResult

        result = CalculationResult()
        assert result.output_dir == ""
        assert result.log_file == ""
        assert result.output_files == {}
        assert result.has_output_dir is False

    def test_output_dir_set(self):
        """Test output directory tracking."""
        from pyabacus.driver.runner import CalculationResult

        with tempfile.TemporaryDirectory() as tmpdir:
            result = CalculationResult(output_dir=tmpdir)
            assert result.output_dir == tmpdir
            assert result.has_output_dir is True

    def test_output_dir_nonexistent(self):
        """Test has_output_dir with non-existent directory."""
        from pyabacus.driver.runner import CalculationResult

        result = CalculationResult(output_dir="/nonexistent/path")
        assert result.output_dir == "/nonexistent/path"
        assert result.has_output_dir is False

    def test_output_files_dict(self):
        """Test output files dictionary."""
        from pyabacus.driver.runner import CalculationResult

        output_files = {
            "running_scf.log": "/path/to/OUT.ABACUS/running_scf.log",
            "BANDS_1.dat": "/path/to/OUT.ABACUS/BANDS_1.dat",
        }
        result = CalculationResult(output_files=output_files)

        assert len(result.output_files) == 2
        assert "running_scf.log" in result.output_files
        assert result.output_files["running_scf.log"] == "/path/to/OUT.ABACUS/running_scf.log"

    def test_get_output_file(self):
        """Test get_output_file method."""
        from pyabacus.driver.runner import CalculationResult

        output_files = {
            "running_scf.log": "/path/to/OUT.ABACUS/running_scf.log",
        }
        result = CalculationResult(output_files=output_files)

        assert result.get_output_file("running_scf.log") == "/path/to/OUT.ABACUS/running_scf.log"
        assert result.get_output_file("nonexistent.dat") is None

    def test_list_output_files(self):
        """Test list_output_files method."""
        from pyabacus.driver.runner import CalculationResult

        output_files = {
            "running_scf.log": "/path/to/OUT.ABACUS/running_scf.log",
            "BANDS_1.dat": "/path/to/OUT.ABACUS/BANDS_1.dat",
            "CHARGE.cube": "/path/to/OUT.ABACUS/CHARGE.cube",
        }
        result = CalculationResult(output_files=output_files)

        file_list = result.list_output_files()
        assert len(file_list) == 3
        assert "running_scf.log" in file_list
        assert "BANDS_1.dat" in file_list
        assert "CHARGE.cube" in file_list

    def test_summary_includes_output_info(self):
        """Test that summary includes output tracking info."""
        from pyabacus.driver.runner import CalculationResult

        with tempfile.TemporaryDirectory() as tmpdir:
            # Create a mock log file
            log_file = os.path.join(tmpdir, "running_scf.log")
            with open(log_file, 'w') as f:
                f.write("test log")

            output_files = {
                "running_scf.log": log_file,
                "BANDS_1.dat": os.path.join(tmpdir, "BANDS_1.dat"),
            }

            result = CalculationResult(
                converged=True,
                etot=-10.0,
                output_dir=tmpdir,
                log_file=log_file,
                output_files=output_files,
            )

            summary = result.summary()
            assert "Output:" in summary
            assert "Directory:" in summary
            assert tmpdir in summary
            assert "Log file: running_scf.log" in summary
            assert "2 output files" in summary

    def test_summary_without_output_info(self):
        """Test that summary works without output tracking info."""
        from pyabacus.driver.runner import CalculationResult

        result = CalculationResult(converged=True, etot=-10.0)
        summary = result.summary()

        # Output section is always present, but shows N/A when not set
        assert "Output:" in summary
        assert "Directory: N/A" in summary
        assert "Log file: N/A" in summary
        assert "0 output files" in summary


class TestCollectOutputFiles:
    """Tests for _collect_output_files helper function."""

    def test_collect_output_files(self):
        """Test collecting output files from directory."""
        from pyabacus.driver.runner import _collect_output_files

        with tempfile.TemporaryDirectory() as tmpdir:
            # Create some test files
            files = ["running_scf.log", "BANDS_1.dat", "CHARGE.cube"]
            for f in files:
                Path(tmpdir, f).touch()

            # Create a subdirectory (should be ignored)
            os.makedirs(os.path.join(tmpdir, "subdir"))

            output_files = _collect_output_files(tmpdir)

            assert len(output_files) == 3
            for f in files:
                assert f in output_files
                assert output_files[f] == os.path.join(tmpdir, f)

    def test_collect_output_files_empty_dir(self):
        """Test collecting from empty directory."""
        from pyabacus.driver.runner import _collect_output_files

        with tempfile.TemporaryDirectory() as tmpdir:
            output_files = _collect_output_files(tmpdir)
            assert output_files == {}

    def test_collect_output_files_nonexistent(self):
        """Test collecting from non-existent directory."""
        from pyabacus.driver.runner import _collect_output_files

        output_files = _collect_output_files("/nonexistent/path")
        assert output_files == {}


class TestParseForces:
    """Tests for _parse_forces_from_log helper function."""

    def test_parse_forces_from_log(self):
        """Test parsing forces from log file."""
        from pyabacus.driver.runner import _parse_forces_from_log

        with tempfile.TemporaryDirectory() as tmpdir:
            log_file = os.path.join(tmpdir, "running_scf.log")
            log_content = """
 TOTAL-FORCE (eV/Angstrom)
 ------------------------------------------------------------------------------------
      atom                   x              y              z
        Si1         0.00100000     0.00200000     0.00300000
        Si2        -0.00100000    -0.00200000    -0.00300000
 ------------------------------------------------------------------------------------
"""
            with open(log_file, 'w') as f:
                f.write(log_content)

            forces = _parse_forces_from_log(log_file, 2)
            assert forces is not None
            assert forces.shape == (2, 3)
            # Forces are now in eV/Angstrom directly (no conversion)
            assert abs(forces[0, 0] - 0.001) < 1e-10
            assert abs(forces[0, 1] - 0.002) < 1e-10
            assert abs(forces[1, 0] - (-0.001)) < 1e-10

    def test_parse_forces_no_forces(self):
        """Test parsing when no forces in log."""
        from pyabacus.driver.runner import _parse_forces_from_log

        with tempfile.TemporaryDirectory() as tmpdir:
            log_file = os.path.join(tmpdir, "running_scf.log")
            with open(log_file, 'w') as f:
                f.write("No forces here\n")

            forces = _parse_forces_from_log(log_file, 2)
            assert forces is None

    def test_parse_forces_nonexistent_file(self):
        """Test parsing from non-existent file."""
        from pyabacus.driver.runner import _parse_forces_from_log

        forces = _parse_forces_from_log("/nonexistent/file.log", 2)
        assert forces is None


class TestParseStress:
    """Tests for _parse_stress_from_log helper function."""

    def test_parse_stress_from_log(self):
        """Test parsing stress from log file."""
        from pyabacus.driver.runner import _parse_stress_from_log

        with tempfile.TemporaryDirectory() as tmpdir:
            log_file = os.path.join(tmpdir, "running_scf.log")
            log_content = """
 TOTAL-STRESS (KBAR)
 ------------------------------------------------------------------------------------
      1.23456789      0.00000000      0.00000000
      0.00000000      1.23456789      0.00000000
      0.00000000      0.00000000      1.23456789
 ------------------------------------------------------------------------------------
"""
            with open(log_file, 'w') as f:
                f.write(log_content)

            stress = _parse_stress_from_log(log_file)
            assert stress is not None
            assert stress.shape == (3, 3)
            assert abs(stress[0, 0] - 1.23456789) < 1e-6
            assert abs(stress[1, 1] - 1.23456789) < 1e-6

    def test_parse_stress_no_stress(self):
        """Test parsing when no stress in log."""
        from pyabacus.driver.runner import _parse_stress_from_log

        with tempfile.TemporaryDirectory() as tmpdir:
            log_file = os.path.join(tmpdir, "running_scf.log")
            with open(log_file, 'w') as f:
                f.write("No stress here\n")

            stress = _parse_stress_from_log(log_file)
            assert stress is None


class TestParallelParameters:
    """Tests for nprocs and nthreads parameters."""

    def test_abacus_accepts_parallel_params(self):
        """Test that abacus() accepts nprocs and nthreads parameters."""
        from pyabacus.driver.runner import abacus
        import inspect

        sig = inspect.signature(abacus)
        params = sig.parameters

        assert 'nprocs' in params
        assert 'nthreads' in params
        assert params['nprocs'].default == 1
        assert params['nthreads'].default == 1


# Fixtures for test data
@pytest.fixture
def mock_cpp_result():
    """Create a mock C++ CalculationResult."""
    result = MagicMock()
    result.converged = True
    result.niter = 15
    result.drho = 1e-9
    result.etot = -15.5
    result.eband = -8.0
    result.hartree_energy = 3.5
    result.etxc = -4.2
    result.ewald_energy = -6.8
    result.demet = 0.0
    result.exx = 0.0
    result.evdw = 0.0
    result.fermi_energy = 6.5
    result.bandgap = 1.2
    result.nat = 2
    result.ntype = 1
    result.nbands = 10
    result.nks = 8
    result.has_forces = True
    result.has_stress = False
    result.forces = np.array([[0.01, 0.02, 0.03], [-0.01, -0.02, -0.03]])
    return result


@pytest.fixture
def temp_input_dir():
    """Create a temporary directory with INPUT file."""
    with tempfile.TemporaryDirectory() as tmpdir:
        input_file = Path(tmpdir) / "INPUT"
        input_file.write_text("""INPUT_PARAMETERS
calculation scf
basis_type pw
ecutwfc 50
scf_thr 1e-8
""")
        yield tmpdir


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
