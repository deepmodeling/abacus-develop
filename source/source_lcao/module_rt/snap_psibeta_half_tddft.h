#ifndef SNAP_PSIBETA_HALF_TDDFT
#define SNAP_PSIBETA_HALF_TDDFT

#include "source_base/vector3.h"
#include "source_basis/module_ao/ORB_read.h"
#include "source_cell/setup_nonlocal.h"
#include "source_lcao/module_rt/snap_projector_half_tddft.h"

#include <complex>
#include <vector>

namespace module_rt
{
/**
 * @brief Compute RT-TDDFT velocity-gauge beta-projector overlaps.
 *
 * This overload uses the production quadrature settings.
 */
void snap_psibeta_half_tddft(const LCAO_Orbitals& orb,
                             const InfoNonlocal& infoNL_,
                             std::vector<std::vector<std::complex<double>>>& nlm,
                             const ModuleBase::Vector3<double>& R1,
                             const int& T1,
                             const int& L1,
                             const int& m1,
                             const int& N1,
                             const ModuleBase::Vector3<double>& R0, // The projector.
                             const int& T0,
                             const ModuleBase::Vector3<double>& A,
                             const bool& calc_r);

/**
 * @brief Compute RT-TDDFT velocity-gauge beta-projector overlaps.
 *
 * This overload is used by tests to select the radial and Lebedev-Laikov grids.
 */
void snap_psibeta_half_tddft(const LCAO_Orbitals& orb,
                             const InfoNonlocal& infoNL_,
                             std::vector<std::vector<std::complex<double>>>& nlm,
                             const ModuleBase::Vector3<double>& R1,
                             const int& T1,
                             const int& L1,
                             const int& m1,
                             const int& N1,
                             const ModuleBase::Vector3<double>& R0,
                             const int& T0,
                             const ModuleBase::Vector3<double>& A,
                             const bool& calc_r,
                             const SnapIntegrationOptions& options);

} // namespace module_rt

#endif
