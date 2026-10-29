#ifndef ABACUS_LR_SPECTRUM_WEIGHTS_H
#define ABACUS_LR_SPECTRUM_WEIGHTS_H

#include <vector>

namespace LR_Util
{
/// KS weights include a factor two only for nspin=1. Optical singlet
/// amplitudes account for that degeneracy separately through sqrt(2).
inline double spatial_k_weight(const double ks_weight, const int nspin_global)
{
    return nspin_global == 1 ? 0.5 * ks_weight : ks_weight;
}

/// Convert KS weights at the LR/spectrum boundary without changing the shared k list.
/// BSE already owns spatial weights and must pass them directly to LR_Spectrum.
inline std::vector<double> spatial_k_weights(const std::vector<double>& ks_weights, const int nspin_global)
{
    std::vector<double> weights(ks_weights);
    for (double& weight : weights)
    {
        weight = spatial_k_weight(weight, nspin_global);
    }
    return weights;
}

}

#endif
