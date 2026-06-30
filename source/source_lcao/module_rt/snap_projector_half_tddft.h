#ifndef SNAP_PROJECTOR_HALF_TDDFT_H
#define SNAP_PROJECTOR_HALF_TDDFT_H

#include "source_base/vector3.h"
#include "source_basis/module_ao/ORB_read.h"

#include <complex>
#include <vector>

namespace module_rt
{

struct ProjectorChannel
{
    int l = 0;
    int mesh = 0;
    double dk = 0.0;
    double rcut = 0.0;
    const double* radial_values = nullptr;
    const double* radial_grid = nullptr;
};

struct SnapIntegrationOptions
{
    int radial_grid_num = 140;
    int lebedev_grid_points = 110;
};

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
