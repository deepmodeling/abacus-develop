#ifndef SNAP_PROJECTOR_HALF_TDDFT_H
#define SNAP_PROJECTOR_HALF_TDDFT_H

#include "source_base/vector3.h"
#include "source_basis/module_ao/ORB_read.h"

#include <complex>
#include <vector>

namespace module_rt
{

/**
 * @brief Radial projector channel integrated against one LCAO orbital.
 */
struct ProjectorChannel
{
    int l = 0;
    int mesh = 0;
    double dk = 0.0;
    double rcut = 0.0;
    const double* radial_values = nullptr;
    const double* radial_grid = nullptr;
};

/**
 * @brief Numerical quadrature settings for projector snapshots.
 *
 * The default values reproduce the production RT-TDDFT path.
 */
struct SnapIntegrationOptions
{
    int radial_grid_num = 140;
    int lebedev_grid_points = 110;
};

/**
 * @brief Compute <phi|exp(-i A r)|projector> with default quadrature settings.
 */
void snap_projector_half_tddft(const LCAO_Orbitals& orb,
                               const std::vector<ProjectorChannel>& projector_channels,
                               std::vector<std::vector<std::complex<double>>>& nlm,
                               const ModuleBase::Vector3<double>& R1,
                               const int& T1,
                               const int& L1,
                               const int& m1,
                               const int& N1,
                               const ModuleBase::Vector3<double>& R0,
                               const ModuleBase::Vector3<double>& A,
                               const bool& calc_r,
                               const char* timer_name);

/**
 * @brief Compute <phi|exp(-i A r)|projector> with explicit quadrature settings.
 *
 * If calc_r is true, nlm[1..3] also store the Cartesian position moments.
 */
void snap_projector_half_tddft(const LCAO_Orbitals& orb,
                               const std::vector<ProjectorChannel>& projector_channels,
                               std::vector<std::vector<std::complex<double>>>& nlm,
                               const ModuleBase::Vector3<double>& R1,
                               const int& T1,
                               const int& L1,
                               const int& m1,
                               const int& N1,
                               const ModuleBase::Vector3<double>& R0,
                               const ModuleBase::Vector3<double>& A,
                               const bool& calc_r,
                               const SnapIntegrationOptions& options,
                               const char* timer_name);

} // namespace module_rt

#endif
