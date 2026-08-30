#ifndef MODULE_RI_RPA_ABFS_PREORTHOGONALIZATION_H
#define MODULE_RI_RPA_ABFS_PREORTHOGONALIZATION_H

#include "source_basis/module_ao/ORB_atomic_lm.h"
#include "source_io/module_parameter/input_parameter.h"

#include <string>
#include <vector>

namespace ModuleRI
{

using RpaAbfsOrbitalSet = std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>>;

struct RpaAbfsPreorthChannelReport
{
    int atom_type = -1;
    int angular_momentum = -1;
    int input_count = 0;
    int output_count = 0;
    std::vector<int> rejected_indices;
    double minimum_residual_norm = 1.0;
    double maximum_identity_error = 0.0;
};

struct RpaAbfsPreorthReport
{
    std::string mode = "none";
    double threshold = 1.0e-2;
    int input_expanded_size = 0;
    int output_expanded_size = 0;
    std::vector<RpaAbfsPreorthChannelReport> channels;
};

std::vector<std::vector<double>> compute_onsite_coulomb_metric(const std::vector<Numerical_Orbital_Lm>& channel);

RpaAbfsPreorthReport apply_rpa_abfs_preorthogonalization(RpaAbfsOrbitalSet& abfs,
                                                         const std::string& mode,
                                                         double threshold,
                                                         bool force_flag);

RpaAbfsPreorthReport finalize_rpa_abfs_from_input(RpaAbfsOrbitalSet& abfs, const Input_para& input, bool force_flag);

std::string format_rpa_abfs_preorth_report(const RpaAbfsPreorthReport& report);

} // namespace ModuleRI

#endif
