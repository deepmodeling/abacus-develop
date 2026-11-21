#ifndef SCCS_LOWPASS_H
#define SCCS_LOWPASS_H

#include <vector>

namespace ModulePW
{
class PW_Basis;
}

namespace ModuleSccs
{
struct CavityParameters;
struct CavityDerivatives;
struct SccsResponse;

// Reciprocal-space filter on the local G vectors; empty when disabled. As
// Environ 3.1.1 core_fft_lowpass.f90 (filter_lowpass): f(G) =
// erfc(p1 G^2 / Gcut^2 - p2) / 2, applied before every FFT derivative.
std::vector<double> switch_filter(const CavityParameters& cavity, const ModulePW::PW_Basis& basis);

// Exact derivative of the filtered, discrete sqrt-CG reaction energy with
// respect to the boundary s; collective over the PW pool.
void lowpass_bnd_pot(const std::vector<double>& charge,
                     const CavityParameters& cavity,
                     const CavityDerivatives& derivatives,
                     const ModulePW::PW_Basis& basis,
                     double tpiba,
                     SccsResponse& response);
} // namespace ModuleSccs

#endif
