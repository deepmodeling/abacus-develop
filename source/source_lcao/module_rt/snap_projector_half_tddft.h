#ifndef SNAP_PROJECTOR_HALF_TDDFT_H
#define SNAP_PROJECTOR_HALF_TDDFT_H

#include "source_base/vector3.h"
#include "source_basis/module_ao/ORB_read.h"

#include <complex>
#include <vector>

namespace module_rt
{

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
 * @brief Radial projector channel integrated against one LCAO orbital.
 *
 * radial_times_r stores r * p_l(r), where p_l(r) is the radial projector.
 * The radial part of the volume integral is therefore evaluated as
 * (r * p_l(r)) * r dr = p_l(r) * r^2 dr.
 */
struct ProjectorChannel
{
    int l = 0;
    int mesh = 0;
    double dk = 0.0;
    double rcut = 0.0;
    const double* radial_times_r = nullptr;
    const double* radial_grid = nullptr;
};

/**
 * @brief Compute projector overlaps with default quadrature settings.
 *
 * The shared integral is
 * I_m(A) = <phi_{T1,L1,m1,N1}(r-R1)|exp(-i A.r)|p_{l,m}(r-R0)>.
 *
 * The phase A is given in Cartesian coordinates. The returned nlm[0] stores
 * I_m(A) for all projector magnetic components. If calc_r is true, nlm[1..3]
 * store
 * R_a,m(A) = <phi|r_a exp(-i A.r)|p_{l,m}>.
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
 * @brief Compute projector overlaps with explicit quadrature settings.
 *
 * The ProjectorChannel radial convention is always r * p_l(r). Callers that
 * own different physical projectors are responsible for passing that form.
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
