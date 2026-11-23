#ifndef MODULE_RI_RPA_ABFS_PREORTH_H
#define MODULE_RI_RPA_ABFS_PREORTH_H

#include "source_basis/module_ao/orb_atomic_lm.h"

#include <string>
#include <vector>

// Only a const reference is used in the declarations below; the full
// Input_para definition is included by the implementation.
struct Input_para;

namespace ModuleRI
{

using AbfsOrbitalSet = std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>>;

struct PreorthChannel
{
    int atom_type = -1;
    int angular_momentum = -1;
    int input_count = 0;
    int output_count = 0;
    std::vector<int> rejected_indices;
    double min_residual_norm = 1.0;
    double max_identity_error = 0.0;
};

struct PreorthReport
{
    std::string mode = "none";
    double threshold = 1.0e-2;
    int input_size = 0;
    int output_size = 0;
    std::vector<PreorthChannel> channels;
};

std::vector<std::vector<double>> coulomb_metric(const std::vector<Numerical_Orbital_Lm>& channel);

PreorthReport apply_preorth(AbfsOrbitalSet& abfs, const std::string& mode, double threshold, bool force_flag);

PreorthReport preorth_from_input(AbfsOrbitalSet& abfs, const Input_para& input, bool force_flag);

std::string preorth_report(const PreorthReport& report);

} // namespace ModuleRI

#endif
