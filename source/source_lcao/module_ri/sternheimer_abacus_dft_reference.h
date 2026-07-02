#ifndef STERNHEIMER_ABACUS_DFT_REFERENCE_H
#define STERNHEIMER_ABACUS_DFT_REFERENCE_H

#include <vector>

namespace elecstate
{
class ElecState;
}

namespace ModuleRI
{

std::vector<double> copy_sternheimer_dft_eigenvalues(const elecstate::ElecState& elec_state,
                                                     int k_index,
                                                     int num_bands);

std::vector<double> copy_sternheimer_dft_occupations(const elecstate::ElecState& elec_state,
                                                     int k_index,
                                                     int num_bands);

} // namespace ModuleRI

#endif
