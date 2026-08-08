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
    (void)orthogonality_tolerance;
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

    const double determinant
        = latvec.e11 * (latvec.e22 * latvec.e33 - latvec.e23 * latvec.e32)
          - latvec.e12 * (latvec.e21 * latvec.e33 - latvec.e23 * latvec.e31)
          + latvec.e13 * (latvec.e21 * latvec.e32 - latvec.e22 * latvec.e31);
    if (std::abs(determinant) <= 1.0e-14 * a1_norm * a2_norm * a3_norm)
    {
        throw std::invalid_argument("Sternheimer ABACUS FD grid requires a nonsingular lattice.");
    }

    SternheimerABACUSFDGridData grid_data;
    grid_data.grid.nx = nx;
    grid_data.grid.ny = ny;
    grid_data.grid.nz = nz;
    grid_data.grid.hx = lat0 * a1_norm / nx;
    grid_data.grid.hy = lat0 * a2_norm / ny;
    grid_data.grid.hz = lat0 * a3_norm / nz;
    grid_data.grid.periodic = true;
    grid_data.grid.lattice_vectors = {{{lat0 * latvec.e11, lat0 * latvec.e12, lat0 * latvec.e13},
                                       {lat0 * latvec.e21, lat0 * latvec.e22, lat0 * latvec.e23},
                                       {lat0 * latvec.e31, lat0 * latvec.e32, lat0 * latvec.e33}}};
    grid_data.volume_element = std::abs(determinant) * lat0 * lat0 * lat0 / static_cast<double>(nxyz);
    return grid_data;
}

std::vector<double> embed_sternheimer_local_z_slab(const std::vector<double>& local_values,
                                                   const int nxy,
                                                   const int nz,
                                                   const int nplane,
                                                   const int startz)
{
    if (nxy <= 0 || nz <= 0 || nplane < 0 || startz < 0 || startz + nplane > nz
        || local_values.size() != static_cast<std::size_t>(nxy) * static_cast<std::size_t>(nplane))
    {
        throw std::invalid_argument("Invalid Sternheimer local z-slab dimensions.");
    }

    std::vector<double> full_values(static_cast<std::size_t>(nxy) * static_cast<std::size_t>(nz), 0.0);
    for (int ixy = 0; ixy != nxy; ++ixy)
    {
        for (int iz = 0; iz != nplane; ++iz)
        {
            full_values[static_cast<std::size_t>(ixy) * static_cast<std::size_t>(nz)
                        + static_cast<std::size_t>(startz + iz)]
                = local_values[static_cast<std::size_t>(ixy) * static_cast<std::size_t>(nplane)
                               + static_cast<std::size_t>(iz)];
        }
    }
    return full_values;
}

SternheimerFDHamiltonian make_sternheimer_fd_hamiltonian_from_local_potential(
    const SternheimerABACUSFDGridData& grid_data,
    std::vector<double> local_potential,
    const double kinetic_prefactor,
    std::shared_ptr<const SternheimerFDNonlocalProjector> nonlocal_projector,
    const int finite_difference_order)
{
    return SternheimerFDHamiltonian(grid_data.grid,
                                    std::move(local_potential),
                                    kinetic_prefactor,
                                    std::move(nonlocal_projector),
                                    finite_difference_order);
}

} // namespace ModuleRI
