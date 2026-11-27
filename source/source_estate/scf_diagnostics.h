#ifndef DFTU_SCF_DIAGNOSTICS_H
#define DFTU_SCF_DIAGNOSTICS_H

#include <ostream>
#include <string>
#include <vector>

namespace ScfDiagnostics
{
// Observational copies only. These objects never feed the Hamiltonian or mixer.
struct ScfOccupationBlock
{
    int atom = 0;
    int spin = 0;
    int l = 0;
    double u_eV = 0.0;
    std::vector<double> input;
    std::vector<double> raw_output;
};

struct ScfDiagnosticStep
{
    std::string output_dir;
    std::string basis_type = "lcao";
    std::string mixing_type;
    int rank = 0;
    int ionic_step = 1;
    int iteration = 1;
    int mixing_dftu = 0;
    bool mixing_dmr = false;
    double mixing_restart = 0.0;
    double beta = 0.0;
    double beta_mag = 0.0;
    double threshold = 0.0;
    double drho = 0.0;
    int scf_thr_type = 2;
    bool has_l2_channels = false;
    double charge_residual_l2 = 0.0;
    double magnetic_residual_l2 = 0.0;
    double charge_residual = 0.0;
    double magnetic_residual = 0.0;
    double wall_seconds = 0.0;
    int restart_count = 0;
    int restart_scheduled_step = 0;
    bool restart_executed = false;
    bool converged = false;
    double gap_eV = 0.0;
    double total_magnetism = 0.0;
    double absolute_magnetism = 0.0;
    double total_energy_eV = 0.0;
    double source_band_correction_eV = 0.0;
};

// Opt-in observer for this diagnostic build; no INPUT/SCF defaults are changed.
bool scf_diagnostics_enabled();
void write_scf_diagnostics(const ScfDiagnosticStep& step,
                           const std::vector<ScfOccupationBlock>& blocks,
                           std::ostream& running_log);
} // namespace ScfDiagnostics
#endif
