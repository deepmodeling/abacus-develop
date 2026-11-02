#include "source_io/module_current/td_spin_pw_io.h"

#include "source_base/tool_quit.h"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits.h>
#include <sstream>
#include <unistd.h>

namespace ModuleIO
{
namespace
{
std::ofstream open_series(const std::string& path, int step, const std::string& header, std::ostream& log)
{
    const std::ios::openmode mode = step == 0 ? std::ios::out : std::ios::app;
    std::ofstream output(path.c_str(), mode);
    if (!output)
    {
        ModuleBase::WARNING_QUIT("write_td_spin_pw", "Cannot open " + path);
    }
    if (step == 0)
    {
        output << header << '\n';
        char cwd[PATH_MAX];
        std::string full_path = path;
        if (!path.empty() && path[0] != '/' && getcwd(cwd, sizeof(cwd)) != nullptr)
        {
            full_path = std::string(cwd) + "/" + path;
        }
        log << "Write file: " << full_path << '\n';
    }
    return output;
}

void finish_series(std::ofstream* output)
{
    output->flush();
    if (!output->good())
    {
        ModuleBase::WARNING_QUIT("write_td_spin_pw", "Failed to write magnetic time series.");
    }
}

void write_endpoint(std::ostream& output, int step, double time_fs)
{
    // Internal step zero is the initial state at time zero, displayed as step one.
    output << std::right << std::setw(6) << step + 1 << "  " << std::fixed << std::setprecision(6) << std::setw(10) << time_fs;
    output << std::scientific << std::setprecision(12);
}

void write_atoms(std::ofstream* output, int step, double time, const std::vector<double>& values)
{
    for (size_t atom = 0; atom < values.size() / 3; ++atom)
    {
        write_endpoint(*output, step, time);
        *output << "  " << std::setw(6) << atom + 1;
        for (int a = 0; a < 3; ++a)
        {
            *output << "  " << std::setw(19) << values[3 * atom + a];
        }
        *output << '\n';
    }
    finish_series(output);
}
} // namespace

void write_td_spin_pw(const std::string& directory,
                      int step,
                      double time_fs,
                      const std::vector<std::string>& orbitals,
                      const std::string& mag_partition,
                      const std::string& orbital_dir,
                      double onsite_radius,
                      bool calc_torque,
                      bool spin_current,
                      const pw::TDSpinResultPW& spin,
                      const pw::TDVelResultPW& current,
                      std::ostream& log)
{
    if (mag_partition != "none")
    {
        std::string header = "# method=" + mag_partition + "; total=density_integral; M=ABACUS_Pauli_moment(mu_B); no_electron_minus_sign";
        if (mag_partition == "onsite")
        {
            std::ostringstream radius;
            radius << std::scientific << std::setprecision(12) << onsite_radius;
            header += "\n# onsite_radius(Bohr)=" + radius.str() + "; projector=q(k+G+A); atom_phase_cancels";
            for (int it = 0; it < orbitals.size(); ++it)
            {
                header += "\n# orbital=" + orbital_dir + orbitals[it];
            }
        }
        std::ostringstream total_columns;
        total_columns << std::left << std::setw(6) << "#step" << "  " << std::right << std::setw(10) << "time_fs";
        for (const char* label: {"Mx", "My", "Mz", "M_norm", "M_abs_integral", "N_elec"})
        {
            total_columns << "  " << std::setw(19) << label;
        }
        std::ostringstream atom_columns;
        atom_columns << std::left << std::setw(6) << "#step" << "  " << std::right << std::setw(10) << "time_fs"
                     << "  " << std::setw(6) << "atom";
        for (const char* label: {"Mx", "My", "Mz"})
        {
            atom_columns << "  " << std::setw(19) << label;
        }
        const std::string total_header = header + '\n' + total_columns.str();
        const std::string atom_header = header + '\n' + atom_columns.str();
        const std::string total_path = directory + "mag_tot.txt";
        const std::string atom_path = directory + "mag_atom.txt";
        std::ofstream total = open_series(total_path, step, total_header, log);
        std::ofstream atoms = open_series(atom_path, step, atom_header, log);
        write_endpoint(total, step, time_fs);
        double norm = 0.0;
        for (double value: spin.mag)
        {
            total << "  " << std::setw(19) << value;
            norm += value * value;
        }
        total << "  " << std::setw(19) << std::sqrt(norm) << "  " << std::setw(19) << spin.mag_abs << "  " << std::setw(19)
              << spin.electrons << '\n';
        finish_series(&total);
        write_atoms(&atoms, step, time_fs, spin.atom_mag);
    }
    if (calc_torque)
    {
        std::ostringstream total_columns;
        total_columns << std::left << std::setw(6) << "#step" << "  " << std::right << std::setw(10) << "time_fs";
        for (const char* label: {"tau_xc_x", "tau_xc_y", "tau_xc_z", "tau_soc_x", "tau_soc_y", "tau_soc_z"})
        {
            total_columns << "  " << std::setw(19) << label;
        }
        std::ostringstream atom_columns;
        atom_columns << std::left << std::setw(6) << "#step" << "  " << std::right << std::setw(10) << "time_fs"
                     << "  " << std::setw(6) << "atom";
        for (const char* label: {"tau_soc_x", "tau_soc_y", "tau_soc_z"})
        {
            atom_columns << "  " << std::setw(19) << label;
        }
        const std::string total_header = "# torque=dS/dt; S=M/2; unit=hbar/Hartree_atomic_time\n" + total_columns.str();
        const std::string atom_header = "# SOC_pseudopotential_atom_decomposition; unit=hbar/Hartree_atomic_time\n"
                                        "# not_the_derivative_of_atomic_moments\n"
                                        + atom_columns.str();
        const std::string total_path = directory + "torque_tot.txt";
        const std::string atom_path = directory + "torque_atom.txt";
        std::ofstream total = open_series(total_path, step, total_header, log);
        std::ofstream atoms = open_series(atom_path, step, atom_header, log);
        write_endpoint(total, step, time_fs);
        for (double value: spin.xc)
        {
            total << "  " << std::setw(19) << value;
        }
        for (double value: spin.soc)
        {
            total << "  " << std::setw(19) << value;
        }
        total << '\n';
        finish_series(&total);
        write_atoms(&atoms, step, time_fs, spin.atom_soc);
    }
    if (spin_current)
    {
        std::ostringstream columns;
        columns << std::left << std::setw(6) << "#step" << "  " << std::right << std::setw(10) << "time_fs";
        for (const char* label: {"Jxx", "Jxy", "Jxz", "Jyx", "Jyy", "Jyz", "Jzx", "Jzy", "Jzz"})
        {
            columns << "  " << std::setw(19) << label;
        }
        const std::string header = "# S=M/2; unit=hbar/(Bohr^2*Hartree_atomic_time); flow_index_first\n" + columns.str();
        const std::string path = directory + "spincur_tot.txt";
        std::ofstream output = open_series(path, step, header, log);
        write_endpoint(output, step, time_fs);
        for (double value: current.spin)
        {
            output << "  " << std::setw(19) << value;
        }
        output << '\n';
        finish_series(&output);
    }
}
} // namespace ModuleIO
