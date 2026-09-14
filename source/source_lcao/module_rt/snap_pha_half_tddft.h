#ifndef SNAP_PHA_HALF_TDDFT_H
#define SNAP_PHA_HALF_TDDFT_H

#include <complex>
#include <vector>

class LCAO_Orbitals;

namespace ModuleBase
{
template <class T> class Vector3;
}

namespace module_rt
{

struct SnapIntegrationOptions;

/**
 * @brief Compute RT-TDDFT velocity-gauge DeePKS alpha-projector overlaps.
 *
 * DeePKS alpha projectors currently use the global descriptor basis
 * orb.Alpha[0]. Each radial channel is passed to the shared projector
 * integrator through getPsi_r(), following the required r * alpha(r)
 * convention.
 */
void snap_phialpha_half_tddft(const LCAO_Orbitals& orb,
                              std::vector<std::vector<std::complex<double>>>& nlm,
                              const ModuleBase::Vector3<double>& R1,
                              const int& T1,
                              const int& L1,
                              const int& m1,
                              const int& N1,
                              const ModuleBase::Vector3<double>& R0,
                              const ModuleBase::Vector3<double>& A,
                              const bool& calc_r);

/**
 * @brief Compute DeePKS alpha-projector overlaps with explicit quadrature.
 */
void snap_phialpha_half_tddft(const LCAO_Orbitals& orb,
                              std::vector<std::vector<std::complex<double>>>& nlm,
                              const ModuleBase::Vector3<double>& R1,
                              const int& T1,
                              const int& L1,
                              const int& m1,
                              const int& N1,
                              const ModuleBase::Vector3<double>& R0,
                              const ModuleBase::Vector3<double>& A,
                              const bool& calc_r,
                              const SnapIntegrationOptions& options);

} // namespace module_rt

#endif
