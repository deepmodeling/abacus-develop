#ifndef ABACUS_NCGGA_SF_GPU_H
#define ABACUS_NCGGA_SF_GPU_H

#include "source_base/matrix.h"
#include "source_base/vector3.h"
#include "source_hamilt/module_xc/xc_ncgga_radial.h"

#include <vector>

namespace ModuleXC
{
#if defined(__CUDA)
/// Evaluate the local LDA and PBE-GGA terms of the regularized NCGGA-SF path
/// on the device. The radial spin map is supplied by the caller; the function
/// returns the four potential channels, the two gradient flux fields, and the
/// energy sum. Returns false when no device evaluation was performed.
bool ncgga_sf_local_gpu(int nrxx,
                        const NcggaSpinMapPoint* spin_map,
                        const double* spin_density_up,
                        const double* spin_density_down,
                        const ModuleBase::Vector3<double>* gradient_up,
                        const ModuleBase::Vector3<double>* gradient_down,
                        double e2,
                        bool is_gga,
                        ModuleBase::matrix& potential,
                        std::vector<ModuleBase::Vector3<double>>& flux_up,
                        std::vector<ModuleBase::Vector3<double>>& flux_down,
                        double& etxc);
#endif
} // namespace ModuleXC

#endif
