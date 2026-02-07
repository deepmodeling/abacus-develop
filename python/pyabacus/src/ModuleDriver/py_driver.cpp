/**
 * @file py_driver.cpp
 * @brief Implementation of PyDriver and pybind11 bindings
 *
 * This file implements the PyDriver class that wraps the complete ABACUS
 * calculation workflow for Python access.
 */

#include "py_driver.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

// ABACUS headers
#include "source_main/driver.h"
#include "source_cell/unitcell.h"
#include "source_cell/check_atomic_stru.h"
#include "source_esolver/esolver.h"
#include "source_io/read_input.h"
#include "source_io/input_conv.h"
#include "source_io/module_parameter/parameter.h"
#include "source_base/global_variable.h"
#include "source_base/global_file.h"
#include "source_base/timer.h"
#include "source_base/memory.h"
#include "source_base/matrix.h"
#include "source_relax/relax_driver.h"

#include <fstream>
#include <sstream>
#include <iostream>
#include <ctime>
#include <stdexcept>
#include <filesystem>

namespace py = pybind11;
namespace fs = std::filesystem;

namespace py_driver
{

/**
 * @brief RAII class for managing global state
 *
 * Saves and restores global state to allow multiple calculations
 * in the same Python session.
 */
class GlobalStateGuard
{
public:
    GlobalStateGuard()
    {
        // Save current state
        saved_my_rank_ = GlobalV::MY_RANK;
        saved_nproc_ = GlobalV::NPROC;
    }

    ~GlobalStateGuard()
    {
        // Restore state
        GlobalV::MY_RANK = saved_my_rank_;
        GlobalV::NPROC = saved_nproc_;
    }

private:
    int saved_my_rank_ = 0;
    int saved_nproc_ = 1;
};

/**
 * @brief Implementation class for PyDriver (PIMPL pattern)
 */
class PyDriver::Impl
{
public:
    Impl() = default;
    ~Impl() { cleanup(); }

    void cleanup()
    {
        if (p_esolver_)
        {
            delete p_esolver_;
            p_esolver_ = nullptr;
        }
        ucell_.reset();
    }

    // ESolver instance
    ModuleESolver::ESolver* p_esolver_ = nullptr;

    // UnitCell instance
    std::unique_ptr<UnitCell> ucell_;

    // Output stream for running log
    std::ofstream ofs_running_;
    std::ofstream ofs_warning_;

    // Original working directory
    std::string original_cwd_;

    // Null stream for silent mode
    std::ofstream null_stream_;

    // Store original stream buffers
    std::streambuf* orig_running_buf_ = nullptr;
    std::streambuf* orig_warning_buf_ = nullptr;
};

PyDriver::PyDriver() : impl_(std::make_unique<Impl>())
{
}

PyDriver::~PyDriver()
{
    cleanup_context();
}

void PyDriver::initialize_context()
{
    // Set up for serial mode (no MPI in Python context)
    PARAM.set_pal_param(0, 1, 1);  // rank=0, nproc=1, nthread=1
    GlobalV::MY_RANK = 0;
    GlobalV::NPROC = 1;

    initialized_ = true;
}

void PyDriver::cleanup_context()
{
    if (impl_)
    {
        impl_->cleanup();

        // Restore original stream buffers
        // Note: We use static_cast to std::ostream& because std::ofstream::rdbuf()
        // doesn't accept arguments, but std::ostream::rdbuf(streambuf*) does
        if (impl_->orig_running_buf_)
        {
            static_cast<std::ostream&>(GlobalV::ofs_running).rdbuf(impl_->orig_running_buf_);
            impl_->orig_running_buf_ = nullptr;
        }
        if (impl_->orig_warning_buf_)
        {
            static_cast<std::ostream&>(GlobalV::ofs_warning).rdbuf(impl_->orig_warning_buf_);
            impl_->orig_warning_buf_ = nullptr;
        }

        // Close output streams
        if (impl_->ofs_running_.is_open())
        {
            impl_->ofs_running_.close();
        }
        if (impl_->ofs_warning_.is_open())
        {
            impl_->ofs_warning_.close();
        }
        if (impl_->null_stream_.is_open())
        {
            impl_->null_stream_.close();
        }

        // Restore working directory if changed
        if (!impl_->original_cwd_.empty())
        {
            try
            {
                fs::current_path(impl_->original_cwd_);
            }
            catch (...)
            {
                // Ignore errors
            }
            impl_->original_cwd_.clear();
        }
    }

    initialized_ = false;
}

void PyDriver::setup_output(const std::string& output_dir, int verbosity)
{
    std::string out_dir = output_dir.empty() ? "OUT.PYABACUS" : output_dir;

    // Create output directory
    fs::create_directories(out_dir);

    // Save original stream buffers
    impl_->orig_running_buf_ = GlobalV::ofs_running.rdbuf();
    impl_->orig_warning_buf_ = GlobalV::ofs_warning.rdbuf();

    // Open log files based on verbosity
    if (verbosity >= 1)
    {
        std::string running_log = out_dir + "/running.log";
        impl_->ofs_running_.open(running_log);
        if (impl_->ofs_running_.is_open())
        {
            static_cast<std::ostream&>(GlobalV::ofs_running).rdbuf(impl_->ofs_running_.rdbuf());
        }
    }
    else
    {
        // Silent mode - redirect to null
        impl_->null_stream_.open("/dev/null");
        if (impl_->null_stream_.is_open())
        {
            static_cast<std::ostream&>(GlobalV::ofs_running).rdbuf(impl_->null_stream_.rdbuf());
        }
    }

    std::string warning_log = out_dir + "/warning.log";
    impl_->ofs_warning_.open(warning_log);
    if (impl_->ofs_warning_.is_open())
    {
        static_cast<std::ostream&>(GlobalV::ofs_warning).rdbuf(impl_->ofs_warning_.rdbuf());
    }
}

void PyDriver::read_input(
    const std::string& input_dir,
    const std::string& input_file,
    const std::string& stru_file,
    const std::string& kpt_file,
    const std::string& pseudo_dir,
    const std::string& orbital_dir,
    const std::string& output_dir)
{
    // Save original working directory
    impl_->original_cwd_ = fs::current_path().string();

    // Determine input file path
    std::string input_path;
    if (!input_file.empty())
    {
        input_path = fs::absolute(input_file).string();
    }
    else
    {
        input_path = (fs::absolute(input_dir) / "INPUT").string();
    }

    // Check if input file exists
    if (!fs::exists(input_path))
    {
        throw std::runtime_error("INPUT file not found: " + input_path);
    }

    // Change to input directory for relative paths
    std::string work_dir = input_dir;
    if (work_dir.empty())
    {
        work_dir = fs::path(input_path).parent_path().string();
    }
    if (!work_dir.empty() && work_dir != ".")
    {
        fs::current_path(work_dir);
    }

    // Read INPUT file
    // Note: ReadInput will set PARAM.globalv.global_in_card internally
    ModuleIO::ReadInput reader(0);  // rank 0
    std::string input_filename = fs::path(input_path).filename().string();
    reader.read_parameters(PARAM, input_filename);

    // Create output directory
    reader.create_directory(PARAM);

    // Convert input parameters to internal format
    Input_Conv::Convert();
}

CalculationResult PyDriver::collect_results(bool calculate_force, bool calculate_stress)
{
    CalculationResult result;

    if (!impl_->p_esolver_ || !impl_->ucell_)
    {
        return result;
    }

    // Get convergence info
    result.converged = impl_->p_esolver_->conv_esolver;

    // Get energy
    result.etot = impl_->p_esolver_->cal_energy();

    // Get system info from UnitCell
    result.nat = impl_->ucell_->nat;
    result.ntype = impl_->ucell_->ntype;

    // Calculate forces if requested
    if (calculate_force)
    {
        ModuleBase::matrix force(result.nat, 3);
        impl_->p_esolver_->cal_force(*impl_->ucell_, force);

        // Convert to numpy array
        std::vector<ssize_t> shape = {static_cast<ssize_t>(result.nat), 3};
        result.forces = py::array_t<double>(shape);
        auto buf = result.forces.request();
        double* ptr = static_cast<double*>(buf.ptr);

        for (int i = 0; i < result.nat; ++i)
        {
            for (int j = 0; j < 3; ++j)
            {
                ptr[i * 3 + j] = force(i, j);
            }
        }
        result.has_forces = true;
    }

    // Calculate stress if requested
    if (calculate_stress)
    {
        ModuleBase::matrix stress(3, 3);
        impl_->p_esolver_->cal_stress(*impl_->ucell_, stress);

        // Convert to numpy array
        std::vector<ssize_t> shape = {3, 3};
        result.stress = py::array_t<double>(shape);
        auto buf = result.stress.request();
        double* ptr = static_cast<double*>(buf.ptr);

        for (int i = 0; i < 3; ++i)
        {
            for (int j = 0; j < 3; ++j)
            {
                ptr[i * 3 + j] = stress(i, j);
            }
        }
        result.has_stress = true;
    }

    // Collect output file tracking information
    result.output_dir = PARAM.sys.global_out_dir;

    // Find the log file
    if (!result.output_dir.empty() && fs::exists(result.output_dir))
    {
        // Look for running_*.log files
        std::vector<std::string> log_patterns = {
            "running_scf.log",
            "running_relax.log",
            "running_cell-relax.log",
            "running_nscf.log",
            "running_md.log"
        };

        for (const auto& log_name : log_patterns)
        {
            std::string log_path = result.output_dir + "/" + log_name;
            if (fs::exists(log_path))
            {
                result.log_file = log_path;
                break;
            }
        }

        // Iterate directory to populate output_files map
        try
        {
            for (const auto& entry : fs::directory_iterator(result.output_dir))
            {
                if (entry.is_regular_file())
                {
                    std::string filename = entry.path().filename().string();
                    std::string full_path = entry.path().string();
                    result.output_files[filename] = full_path;
                }
            }
        }
        catch (const std::exception& e)
        {
            // Ignore errors during directory iteration
        }
    }

    return result;
}

CalculationResult PyDriver::run(
    const std::string& input_dir,
    const std::string& input_file,
    const std::string& stru_file,
    const std::string& kpt_file,
    const std::string& pseudo_dir,
    const std::string& orbital_dir,
    const std::string& output_dir,
    bool calculate_force,
    bool calculate_stress,
    int verbosity)
{
    // Use RAII guard for global state
    GlobalStateGuard state_guard;

    // Clean up any previous calculation
    cleanup_context();

    // Initialize context
    initialize_context();

    // Setup output
    setup_output(output_dir, verbosity);

    // Start timer
    ModuleBase::timer::start();

    try
    {
        // Read input files
        read_input(input_dir, input_file, stru_file, kpt_file,
                   pseudo_dir, orbital_dir, output_dir);

        // Create UnitCell
        impl_->ucell_ = std::make_unique<UnitCell>();
        impl_->ucell_->setup(
            PARAM.inp.latname,
            PARAM.inp.ntype,
            PARAM.inp.lmaxmax,
            PARAM.inp.init_vel,
            PARAM.inp.fixed_axes
        );

        // Read structure
        impl_->ucell_->setup_cell(PARAM.globalv.global_in_stru, GlobalV::ofs_running);

        // Check atomic structure
        unitcell::check_atomic_stru(*impl_->ucell_, PARAM.inp.min_dist_coef);

        // Initialize ESolver
        impl_->p_esolver_ = ModuleESolver::init_esolver(PARAM.inp, *impl_->ucell_);

        // Run before_all_runners
        impl_->p_esolver_->before_all_runners(*impl_->ucell_, PARAM.inp);

        // Run calculation based on calculation type
        const std::string& cal = PARAM.inp.calculation;

        if (cal == "scf" || cal == "relax" || cal == "cell-relax" || cal == "nscf")
        {
            Relax_Driver rl_driver;
            rl_driver.relax_driver(impl_->p_esolver_, *impl_->ucell_, PARAM.inp);
        }
        else if (cal == "get_s")
        {
            impl_->p_esolver_->runner(*impl_->ucell_, 0);
        }
        else
        {
            throw std::runtime_error("Unsupported calculation type: " + cal);
        }

        // Collect results
        last_result_ = collect_results(calculate_force, calculate_stress);

        // Run after_all_runners
        impl_->p_esolver_->after_all_runners(*impl_->ucell_);
    }
    catch (const std::exception& e)
    {
        // Stop timer on error
        ModuleBase::timer::finish(GlobalV::ofs_running);
        // Clean up on error
        cleanup_context();
        throw;
    }

    // Stop timer
    ModuleBase::timer::finish(GlobalV::ofs_running);

    // Print memory usage
    ModuleBase::Memory::print_all(GlobalV::ofs_running);

    return last_result_;
}

} // namespace py_driver

// ============================================================================
// Pybind11 Module Definition
// ============================================================================

PYBIND11_MODULE(_driver_pack, m)
{
    m.doc() = R"pbdoc(
        PyABACUS Driver Module
        ----------------------

        This module provides Python bindings for running complete ABACUS
        DFT calculations.

        Main Classes
        ------------
        PyDriver : Main driver class for running calculations
        CalculationResult : Container for calculation results

        Example
        -------
        >>> from pyabacus.driver import PyDriver
        >>> driver = PyDriver()
        >>> result = driver.run("./Si_scf/")
        >>> print(f"Energy: {result.etot_eV()} eV")
        >>> print(result.summary())
    )pbdoc";

    // Bind CalculationResult
    py::class_<py_driver::CalculationResult>(m, "CalculationResult",
        R"pbdoc(
        Container for DFT calculation results.

        Attributes
        ----------
        converged : bool
            Whether SCF converged
        niter : int
            Number of SCF iterations
        drho : float
            Final charge density difference
        etot : float
            Total energy in Rydberg
        forces : numpy.ndarray
            Forces on atoms (nat, 3) in Ry/Bohr
        stress : numpy.ndarray
            Stress tensor (3, 3) in kbar
        )pbdoc")
        .def(py::init<>())
        .def_readonly("converged", &py_driver::CalculationResult::converged,
            "Whether SCF converged")
        .def_readonly("niter", &py_driver::CalculationResult::niter,
            "Number of SCF iterations")
        .def_readonly("drho", &py_driver::CalculationResult::drho,
            "Final charge density difference")
        .def_readonly("etot", &py_driver::CalculationResult::etot,
            "Total energy (Ry)")
        .def_readonly("eband", &py_driver::CalculationResult::eband,
            "Band energy (Ry)")
        .def_readonly("hartree_energy", &py_driver::CalculationResult::hartree_energy,
            "Hartree energy (Ry)")
        .def_readonly("etxc", &py_driver::CalculationResult::etxc,
            "Exchange-correlation energy (Ry)")
        .def_readonly("ewald_energy", &py_driver::CalculationResult::ewald_energy,
            "Ewald energy (Ry)")
        .def_readonly("demet", &py_driver::CalculationResult::demet,
            "-TS term for metals (Ry)")
        .def_readonly("exx", &py_driver::CalculationResult::exx,
            "Exact exchange energy (Ry)")
        .def_readonly("evdw", &py_driver::CalculationResult::evdw,
            "van der Waals energy (Ry)")
        .def_readonly("forces", &py_driver::CalculationResult::forces,
            "Forces on atoms (nat, 3) in Ry/Bohr")
        .def_readonly("has_forces", &py_driver::CalculationResult::has_forces,
            "Whether forces are available")
        .def_readonly("stress", &py_driver::CalculationResult::stress,
            "Stress tensor (3, 3) in kbar")
        .def_readonly("has_stress", &py_driver::CalculationResult::has_stress,
            "Whether stress is available")
        .def_readonly("fermi_energy", &py_driver::CalculationResult::fermi_energy,
            "Fermi energy (eV)")
        .def_readonly("bandgap", &py_driver::CalculationResult::bandgap,
            "Band gap (eV)")
        .def_readonly("nat", &py_driver::CalculationResult::nat,
            "Number of atoms")
        .def_readonly("ntype", &py_driver::CalculationResult::ntype,
            "Number of atom types")
        .def_readonly("nbands", &py_driver::CalculationResult::nbands,
            "Number of bands")
        .def_readonly("nks", &py_driver::CalculationResult::nks,
            "Number of k-points")
        .def_readonly("output_dir", &py_driver::CalculationResult::output_dir,
            "Path to output directory (OUT.$suffix)")
        .def_readonly("log_file", &py_driver::CalculationResult::log_file,
            "Path to the main log file")
        .def_readonly("output_files", &py_driver::CalculationResult::output_files,
            "Dictionary of output files (filename -> full path)")
        .def("etot_eV", &py_driver::CalculationResult::etot_eV,
            "Get total energy in eV")
        .def("get_energies", &py_driver::CalculationResult::get_energies,
            "Get all energies as a dictionary")
        .def("get_forces_eV_Ang", &py_driver::CalculationResult::get_forces_eV_Ang,
            "Get forces in eV/Angstrom")
        .def("summary", &py_driver::CalculationResult::summary,
            "Get a summary string of the calculation result")
        .def("__repr__", [](const py_driver::CalculationResult& r) {
            std::ostringstream ss;
            ss << "<CalculationResult converged=" << (r.converged ? "True" : "False")
               << " etot=" << r.etot << " Ry>";
            return ss.str();
        });

    // Bind PyDriver
    py::class_<py_driver::PyDriver>(m, "PyDriver",
        R"pbdoc(
        Python wrapper for ABACUS Driver.

        This class provides a Python interface for running complete ABACUS
        DFT calculations.

        Example
        -------
        >>> driver = PyDriver()
        >>> result = driver.run(
        ...     input_dir="./Si_scf/",
        ...     calculate_force=True,
        ...     calculate_stress=True
        ... )
        >>> print(f"Energy: {result.etot_eV()} eV")
        >>> print(f"Converged: {result.converged}")
        )pbdoc")
        .def(py::init<>())
        .def("run", &py_driver::PyDriver::run,
            R"pbdoc(
            Run a complete DFT calculation.

            Parameters
            ----------
            input_dir : str, optional
                Directory containing INPUT, STRU, KPT files (default: ".")
            input_file : str, optional
                Explicit path to INPUT file
            stru_file : str, optional
                Explicit path to STRU file
            kpt_file : str, optional
                Explicit path to KPT file
            pseudo_dir : str, optional
                Directory containing pseudopotentials
            orbital_dir : str, optional
                Directory containing orbital files
            output_dir : str, optional
                Directory for output files
            calculate_force : bool, optional
                Whether to calculate forces (default: True)
            calculate_stress : bool, optional
                Whether to calculate stress (default: False)
            verbosity : int, optional
                Output verbosity level (0=silent, 1=normal, 2=verbose)

            Returns
            -------
            CalculationResult
                Container with all calculation results
            )pbdoc",
            py::arg("input_dir") = ".",
            py::arg("input_file") = "",
            py::arg("stru_file") = "",
            py::arg("kpt_file") = "",
            py::arg("pseudo_dir") = "",
            py::arg("orbital_dir") = "",
            py::arg("output_dir") = "",
            py::arg("calculate_force") = true,
            py::arg("calculate_stress") = false,
            py::arg("verbosity") = 1)
        .def("is_ready", &py_driver::PyDriver::is_ready,
            "Check if the driver is ready for calculation")
        .def("get_last_result", &py_driver::PyDriver::get_last_result,
            py::return_value_policy::reference_internal,
            "Get the last calculation result");
}
