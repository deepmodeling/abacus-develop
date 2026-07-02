#include "source_lcao/module_ri/sternheimer_abacus_fd_adapter.h"

#include "source_base/matrix3.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace
{

double dot_product(const double ax, const double ay, const double az, const double bx, const double by, const double bz)
{
    return ax * bx + ay * by + az * bz;
}

double vector_norm(const double x, const double y, const double z)
{
    return std::sqrt(dot_product(x, y, z, x, y, z));
}

bool is_orthogonal_pair(const double dot, const double norm_a, const double norm_b, const double tolerance)
{
    return std::abs(dot) <= tolerance * std::max(1.0, norm_a * norm_b);
}

} // namespace

namespace ModuleRI
{

SternheimerABACUSFDGridData make_sternheimer_fd_grid_from_lattice(const int nx,
                                                                  const int ny,
                                                                  const int nz,
                                                                  const int nrxx,
                                                                  const double lat0,
                                                                  const ModuleBase::Matrix3& latvec,
                                                                  const double orthogonality_tolerance)
{
    if (nx <= 0 || ny <= 0 || nz <= 0)
    {
        throw std::invalid_argument("Sternheimer ABACUS FD grid requires positive grid dimensions.");
    }
    if (lat0 <= 0.0)
    {
        throw std::invalid_argument("Sternheimer ABACUS FD grid requires a positive lattice constant.");
    }

    const int nxyz = nx * ny * nz;
    if (nrxx != nxyz)
    {
        throw std::invalid_argument(
            "Sternheimer ABACUS FD dense prototype requires an undistributed full real-space grid.");
    }

    const double a1_norm = vector_norm(latvec.e11, latvec.e12, latvec.e13);
    const double a2_norm = vector_norm(latvec.e21, latvec.e22, latvec.e23);
    const double a3_norm = vector_norm(latvec.e31, latvec.e32, latvec.e33);
    if (a1_norm <= 0.0 || a2_norm <= 0.0 || a3_norm <= 0.0)
    {
        throw std::invalid_argument("Sternheimer ABACUS FD grid requires nonzero lattice vectors.");
    }

    const double a12 = dot_product(latvec.e11, latvec.e12, latvec.e13, latvec.e21, latvec.e22, latvec.e23);
    const double a13 = dot_product(latvec.e11, latvec.e12, latvec.e13, latvec.e31, latvec.e32, latvec.e33);
    const double a23 = dot_product(latvec.e21, latvec.e22, latvec.e23, latvec.e31, latvec.e32, latvec.e33);
    if (!is_orthogonal_pair(a12, a1_norm, a2_norm, orthogonality_tolerance)
        || !is_orthogonal_pair(a13, a1_norm, a3_norm, orthogonality_tolerance)
        || !is_orthogonal_pair(a23, a2_norm, a3_norm, orthogonality_tolerance))
    {
        throw std::invalid_argument(
            "Sternheimer ABACUS FD stencil currently supports only orthogonal real-space lattices.");
    }

    SternheimerABACUSFDGridData grid_data;
    grid_data.grid.nx = nx;
    grid_data.grid.ny = ny;
    grid_data.grid.nz = nz;
    grid_data.grid.hx = lat0 * a1_norm / nx;
    grid_data.grid.hy = lat0 * a2_norm / ny;
    grid_data.grid.hz = lat0 * a3_norm / nz;
    grid_data.grid.periodic = true;
    grid_data.volume_element = grid_data.grid.hx * grid_data.grid.hy * grid_data.grid.hz;
    return grid_data;
}

SternheimerFDHamiltonian make_sternheimer_fd_hamiltonian_from_local_potential(
    const SternheimerABACUSFDGridData& grid_data,
    std::vector<double> local_potential,
    const double kinetic_prefactor,
    std::shared_ptr<const SternheimerFDNonlocalProjector> nonlocal_projector)
{
    return SternheimerFDHamiltonian(grid_data.grid,
                                    std::move(local_potential),
                                    kinetic_prefactor,
                                    std::move(nonlocal_projector));
}

} // namespace ModuleRI
