/**
 * @file py_driver.hpp
 * @brief Python bindings for ABACUS Driver - complete DFT calculation workflow
 *
 * This file provides the PyDriver class that wraps the complete ABACUS
 * calculation workflow, enabling Python to run full DFT calculations.
 */

#ifndef PY_DRIVER_HPP
#define PY_DRIVER_HPP

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <string>
#include <iomanip>
#include <memory>
#include <vector>
#include <map>
#include <sstream>

namespace py = pybind11;

namespace py_driver
{

/**
 * @brief Result container for calculation results
 *
 * Stores all results from a DFT calculation including energies,
 * forces, stress, and convergence information.
 */
struct CalculationResult
{
    // Convergence info
    bool converged = false;
    int niter = 0;
    double drho = 0.0;

    // Energies (in Rydberg)
    double etot = 0.0;
    double eband = 0.0;
    double hartree_energy = 0.0;
    double etxc = 0.0;
    double ewald_energy = 0.0;
    double demet = 0.0;
    double exx = 0.0;
    double evdw = 0.0;

    // Forces (nat x 3, in Ry/Bohr)
    py::array_t<double> forces;
    bool has_forces = false;

    // Stress (3 x 3, in kbar)
    py::array_t<double> stress;
    bool has_stress = false;

    // Electronic structure info
    double fermi_energy = 0.0;  // in eV
    double bandgap = 0.0;       // in eV
    int nat = 0;
    int ntype = 0;
    int nbands = 0;
    int nks = 0;

    // Output file tracking
    std::string output_dir = "";   // Path to OUT.$suffix folder
    std::string log_file = "";     // Path to the main log file
    std::map<std::string, std::string> output_files;  // filename -> full path

    // Unit conversion constants
    static constexpr double Ry_to_eV = 13.605693122994;
    static constexpr double Bohr_to_Ang = 0.529177249;

    // Convenience methods
    double etot_eV() const { return etot * Ry_to_eV; }

    py::dict get_energies() const
    {
        py::dict result;
        result["etot"] = etot;
        result["etot_eV"] = etot * Ry_to_eV;
        result["eband"] = eband;
        result["hartree_energy"] = hartree_energy;
        result["etxc"] = etxc;
        result["ewald_energy"] = ewald_energy;
        result["demet"] = demet;
        result["exx"] = exx;
        result["evdw"] = evdw;
        return result;
    }

    py::array_t<double> get_forces_eV_Ang() const
    {
        if (!has_forces)
        {
            throw std::runtime_error("Forces not available. Set calculate_force=True.");
        }
        // Convert from Ry/Bohr to eV/Ang
        auto buf = forces.request();
        auto result = py::array_t<double>(buf.shape);
        auto result_buf = result.request();
        double* src = static_cast<double*>(buf.ptr);
        double* dst = static_cast<double*>(result_buf.ptr);
        double factor = Ry_to_eV / Bohr_to_Ang;
        for (ssize_t i = 0; i < buf.size; ++i)
        {
            dst[i] = src[i] * factor;
        }
        return result;
    }

    std::string summary() const
    {
        std::ostringstream ss;
        ss << "=== ABACUS Calculation Result ===\n";
        ss << "Converged: " << (converged ? "Yes" : "No") << "\n";
        ss << "SCF iterations: " << niter << "\n";
        ss << "Final drho: " << std::scientific << drho << "\n";
        ss << "\nEnergies:\n";
        ss << std::fixed << std::setprecision(8);
        ss << "  Total energy: " << etot << " Ry (" << etot * Ry_to_eV << " eV)\n";
        ss << "  Band energy:  " << eband << " Ry\n";
        ss << "  Hartree:      " << hartree_energy << " Ry\n";
        ss << "  XC energy:    " << etxc << " Ry\n";
        ss << "  Ewald:        " << ewald_energy << " Ry\n";
        if (has_forces)
        {
            ss << "\nForces: calculated (" << nat << " atoms)\n";
        }
        if (has_stress)
        {
            ss << "Stress: calculated\n";
        }
        ss << "\nSystem info:\n";
        ss << "  Atoms: " << nat << ", Types: " << ntype << "\n";
        ss << "  Bands: " << nbands << ", K-points: " << nks << "\n";
        if (fermi_energy != 0.0)
        {
            ss << "  Fermi energy: " << fermi_energy << " eV\n";
        }
        if (bandgap > 0.0)
        {
            ss << "  Band gap: " << bandgap << " eV\n";
        }
        // Output file tracking
        if (!output_dir.empty())
        {
            ss << "\nOutput:\n";
            ss << "  Directory: " << output_dir << "\n";
            if (!log_file.empty())
            {
                // Extract just the filename from the path
                size_t pos = log_file.find_last_of("/\\");
                std::string log_filename = (pos != std::string::npos) ? log_file.substr(pos + 1) : log_file;
                ss << "  Log file: " << log_filename << "\n";
            }
            if (!output_files.empty())
            {
                ss << "  Files: " << output_files.size() << " output files\n";
            }
        }
        return ss.str();
    }
};

/**
 * @brief Python wrapper for ABACUS Driver
 *
 * This class provides a Python interface for running complete ABACUS
 * DFT calculations. It handles:
 * - Input file reading (INPUT, STRU, KPT)
 * - ESolver initialization and execution
 * - Result collection (energy, forces, stress)
 * - Global state management
 */
class PyDriver
{
public:
    PyDriver();
    ~PyDriver();

    // Disable copy
    PyDriver(const PyDriver&) = delete;
    PyDriver& operator=(const PyDriver&) = delete;

    /**
     * @brief Run a complete DFT calculation
     *
     * @param input_dir Directory containing INPUT, STRU, KPT files
     * @param input_file Optional: explicit path to INPUT file
     * @param stru_file Optional: explicit path to STRU file
     * @param kpt_file Optional: explicit path to KPT file
     * @param pseudo_dir Optional: directory containing pseudopotentials
     * @param orbital_dir Optional: directory containing orbital files
     * @param output_dir Optional: directory for output files
     * @param calculate_force Whether to calculate forces
     * @param calculate_stress Whether to calculate stress
     * @param verbosity Output verbosity level (0=silent, 1=normal, 2=verbose)
     * @return CalculationResult containing all results
     */
    CalculationResult run(
        const std::string& input_dir = ".",
        const std::string& input_file = "",
        const std::string& stru_file = "",
        const std::string& kpt_file = "",
        const std::string& pseudo_dir = "",
        const std::string& orbital_dir = "",
        const std::string& output_dir = "",
        bool calculate_force = true,
        bool calculate_stress = false,
        int verbosity = 1
    );

    /**
     * @brief Check if the driver is ready for calculation
     */
    bool is_ready() const { return initialized_; }

    /**
     * @brief Get the last calculation result
     */
    const CalculationResult& get_last_result() const { return last_result_; }

private:
    class Impl;
    std::unique_ptr<Impl> impl_;

    bool initialized_ = false;
    CalculationResult last_result_;

    // Internal methods
    void initialize_context();
    void cleanup_context();
    void read_input(const std::string& input_dir,
                    const std::string& input_file,
                    const std::string& stru_file,
                    const std::string& kpt_file,
                    const std::string& pseudo_dir,
                    const std::string& orbital_dir,
                    const std::string& output_dir);
    void setup_output(const std::string& output_dir, int verbosity);
    CalculationResult collect_results(bool calculate_force, bool calculate_stress);
};

} // namespace py_driver

#endif // PY_DRIVER_HPP
