#ifndef ABACUS_LR_SPECTRUM_WEIGHTS_H
#define ABACUS_LR_SPECTRUM_WEIGHTS_H

#include <vector>
#include <complex>

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

/// One production velocity-gauge transition, before the common i/omega and
/// singlet factors. Spatial weights are indexed by response spin and k point.
template <typename T>
inline std::complex<double> weighted_velocity_transition(
    const std::complex<double>& velocity, const T& x, const T& y,
    const double gap, const bool use_ks_gap, const std::vector<double>& weights,
    const int spin, const int ik, const int nk)
{
    const double weight = weights[spin * nk + ik] * nk;
    const std::complex<double> excitation = velocity * x;
    const std::complex<double> deexcitation = std::conj(velocity) * y;
    if (use_ks_gap) { return (excitation - deexcitation) * (weight / gap); }
    return (excitation + deexcitation) * weight;
}

}

#endif
