#include "source_estate/scf_diagnostics.h"

#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <stdexcept>

namespace
{
void write_matrix(std::ostream& out, const std::vector<double>& matrix)
{
    if (matrix.empty())
    {
        out << "null";
        return;
    }
    out << '[';
    for (std::size_t i = 0; i < matrix.size(); ++i)
    {
        if (i != 0)
        {
            out << ',';
        }
        out << matrix[i];
    }
    out << ']';
}
}

bool ScfDiagnostics::scf_diagnostics_enabled()
{
    const char* value = std::getenv("ABACUS_SCF_DIAGNOSTICS");
    return value != nullptr && std::string(value) == "1";
}

void ScfDiagnostics::write_scf_diagnostics(const ScfDiagnosticStep& step,
                                     const std::vector<ScfOccupationBlock>& blocks,
                                     std::ostream& running_log)
{
    if (step.rank != 0)
    {
        return;
    }
    const std::string path = step.output_dir + "scf_diagnostics.jsonl";
    const bool first = step.ionic_step == 1 && step.iteration == 1;
    const std::ios::openmode mode = first ? std::ios::trunc : std::ios::app;
    std::ofstream out(path.c_str(), std::ios::out | mode);
    if (!out)
    {
        throw std::runtime_error("Cannot write SCF diagnostics: " + path);
    }
    if (first)
    {
        running_log << "Write SCF diagnostic records to file: " << path << '\n';
    }
    out << std::scientific << std::setprecision(17);
    out << "{\"schema\":3,\"basis_type\":\"" << step.basis_type << "\",\"ionic_step\":" << step.ionic_step
        << ",\"iteration\":" << step.iteration
        << ",\"threshold\":" << step.threshold << ",\"drho\":" << step.drho
        << ",\"charge_residual_L1_per_electron\":" << step.charge_residual
        << ",\"magnetic_residual_L1_per_electron\":" << step.magnetic_residual
        << ",\"scf_thr_type\":" << step.scf_thr_type;
    out << ",\"charge_residual_L2_normalized\":";
    if (step.has_l2_channels) { out << step.charge_residual_l2; } else { out << "null"; }
    out << ",\"magnetic_residual_L2_normalized\":";
    if (step.has_l2_channels) { out << step.magnetic_residual_l2; } else { out << "null"; }
    out
        << ",\"scf_wall_seconds_before_occupation_observer\":" << step.wall_seconds
        << ",\"mixing_restart_count\":" << step.restart_count
        << ",\"mixing_restart_scheduled_step\":" << step.restart_scheduled_step
        << ",\"mixing_restart_executed\":" << step.restart_executed
        << ",\"scf_converged\":" << step.converged
        << ",\"logged_E_gap_k_eV\":" << step.gap_eV
        << ",\"fundamental_gap_eV\":null"
        << ",\"total_magnetism_muB_cell\":" << step.total_magnetism
        << ",\"absolute_magnetism_muB_cell\":" << step.absolute_magnetism
        << ",\"total_energy_eV\":" << step.total_energy_eV
        << ",\"source_band_correction_eV\":" << step.source_band_correction_eV
        << ",\"mixing_type\":\"" << step.mixing_type << "\",\"mixing_beta\":" << step.beta
        << ",\"mixing_beta_mag\":" << step.beta_mag << ",\"mixing_dftu\":" << step.mixing_dftu
        << ",\"mixing_dmr\":" << step.mixing_dmr << ",\"mixing_restart\":" << step.mixing_restart
        << ",\"blocks\":[";
    for (std::size_t i = 0; i < blocks.size(); ++i)
    {
        const ScfOccupationBlock& b = blocks[i];
        if (i != 0)
        {
            out << ',';
        }
        out << "{\"atom\":" << b.atom << ",\"spin\":" << b.spin << ",\"l\":" << b.l
            << ",\"dimension\":" << 2 * b.l + 1 << ",\"U_eff_eV\":" << b.u_eV
            << ",\"n_in\":";
        write_matrix(out, b.input);
        out << ",\"n_raw_out\":";
        write_matrix(out, b.raw_output);
        // This observer only supports onsite collinear runs without occupation
        // mixing. Do not invent a post-mixing occupation when none exists.
        out << ",\"n_after_occupation_mixing\":null}";
    }
    out << "]}\n";
    out.flush();
    if (!out)
    {
        throw std::runtime_error("Incomplete SCF diagnostic write: " + path);
    }
}
