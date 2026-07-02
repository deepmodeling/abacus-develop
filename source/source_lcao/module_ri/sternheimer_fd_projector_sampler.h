#ifndef STERNHEIMER_FD_PROJECTOR_SAMPLER_H
#define STERNHEIMER_FD_PROJECTOR_SAMPLER_H

#include "source_base/matrix.h"
#include "source_base/vector3.h"
#include "source_lcao/module_ri/sternheimer_fd_hamiltonian.h"
#include "source_lcao/module_ri/sternheimer_fd_nonlocal_projector.h"

#include <complex>
#include <vector>

namespace ModuleRI
{

struct SternheimerFDRadialProjectorSet
{
    using Complex = std::complex<double>;
    using Matrix = SternheimerFDNonlocalProjector::Matrix;

    std::vector<double> radial_grid;
    std::vector<std::vector<double>> beta_radials;
    std::vector<int> angular_momenta;
    Matrix d_radial;
};

SternheimerFDRadialProjectorSet
make_sternheimer_fd_radial_projector_set_from_abacus_matrices(const std::vector<double>& radial_grid,
                                                              const ModuleBase::matrix& beta_radials,
                                                              const std::vector<int>& angular_momenta,
                                                              const ModuleBase::matrix& d_radial);

SternheimerFDNonlocalProjector::ProjectorBlock
sample_sternheimer_fd_projector_block(const SternheimerFDRadialProjectorSet& radial_set,
                                      const SternheimerFDHamiltonian::Grid& grid,
                                      const ModuleBase::Vector3<double>& atom_position);

} // namespace ModuleRI

#endif
