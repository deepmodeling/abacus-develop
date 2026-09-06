#ifndef DFTU_HAMILT_H
#define DFTU_HAMILT_H

#include <complex>
#include <vector>

class Plus_U;
class Parallel_Orbitals;
class UnitCell;

#ifdef __LCAO
namespace DFTU_LCAO {

/// @brief Compute the LCAO-basis U-term effective potential matrix (complex).
/// Calls DFTU_LCAO::pot_onsite_complex plus the S-projection GEMM.
void pot_uterm_complex(Plus_U& dftu,
                       const UnitCell& ucell,
                       const Parallel_Orbitals* pv,
                       const int ik,
                       std::complex<double>* pot_uterm,
                       const std::vector<int>& isk,
                       const std::complex<double>* sk,
                       const int npol);

/// @brief Compute the LCAO-basis U-term effective potential matrix (real).
/// Calls DFTU_LCAO::pot_onsite_real plus the S-projection GEMM.
void pot_uterm_real(Plus_U& dftu,
                    const UnitCell& ucell,
                    const Parallel_Orbitals* pv,
                    const int ik,
                    double* pot_uterm,
                    const std::vector<int>& isk,
                    const double* sk,
                    const int npol);

/// @brief Accumulate the DFT+U term into the real-space HR (double).
/// Wraps pot_onsite_real plus the (pot_onsite*SR + SR*pot_onsite)/2 GEMM pair.
void pot_uterm_HR_real(const Plus_U& dftu,
                       const UnitCell& ucell,
                       const Parallel_Orbitals* pv,
                       const int ispin,
                       double* SR,
                       double* HR,
                       const int npol);

/// @brief Accumulate the DFT+U term into the real-space HR (complex).
/// Wraps pot_onsite_complex plus the (pot_onsite*SR + SR*pot_onsite)/2 GEMM pair.
void pot_uterm_HR_complex(const Plus_U& dftu,
                          const UnitCell& ucell,
                          const Parallel_Orbitals* pv,
                          const int ispin,
                          std::complex<double>* SR,
                          std::complex<double>* HR,
                          const int npol);

} // namespace DFTU_LCAO
#endif

#endif
