#ifndef DFTU_FULL_DC_H
#define DFTU_FULL_DC_H

#include <complex>
#include <vector>

namespace DFTU_LCAO
{
// One per-spin Hermitian shell, row-major. alpha=0 is AMF, alpha=1 is FLL.
// Returns the interaction energy; potential is dE/dn with dE = Re Tr(V dn).
// band_correction = E - Re Tr(V n), using the SAME occupation as V.
// This is the Full operator's existing energy convention, not a correction
// evaluated with a different (post-diagonalization) occupation.
template <typename T>
double evaluate_full_dc(const std::vector<T>& occupation,
                        int dimension,
                        double u_eff,
                        double alpha,
                        std::vector<T>& potential,
                        double& band_correction);

extern template double evaluate_full_dc<double>(
    const std::vector<double>&, int, double, double, std::vector<double>&, double&);
extern template double evaluate_full_dc<std::complex<double>>(
    const std::vector<std::complex<double>>&, int, double, double,
    std::vector<std::complex<double>>&, double&);
} // namespace DFTU_LCAO

#endif
