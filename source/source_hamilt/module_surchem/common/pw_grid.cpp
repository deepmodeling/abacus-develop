#include "pw_grid.h"
#include "lattice_row.h"

#include "source_base/matrix3.h"
#include "source_basis/module_pw/pw_basis.h"

#include <cmath>
#include <stdexcept>

namespace ModuleSurchem
{
std::vector<ModuleBase::Vector3<double>> pw_grid_positions(
    const ModulePW::PW_Basis& basis,
    const ModuleBase::Matrix3& lattice_vectors,
    const double lattice_scale)
{
    if (basis.nx <= 0 || basis.ny <= 0 || basis.nz <= 0 || basis.nplane <= 0
        || basis.nrxx != basis.nx * basis.ny * basis.nplane
        || !std::isfinite(lattice_scale) || lattice_scale <= 0.0)
    {
        throw std::invalid_argument("SCCS PW grid positions require initialized grid dimensions and lattice");
    }
    const ModuleBase::Vector3<double> a1 = lattice_row(lattice_vectors, 0, lattice_scale);
    const ModuleBase::Vector3<double> a2 = lattice_row(lattice_vectors, 1, lattice_scale);
    const ModuleBase::Vector3<double> a3 = lattice_row(lattice_vectors, 2, lattice_scale);
    std::vector<ModuleBase::Vector3<double>> positions(basis.nrxx);
    for (int ir = 0; ir < basis.nrxx; ++ir)
    {
        const int ix = ir / (basis.ny * basis.nplane);
        const int iy = ir / basis.nplane - ix * basis.ny;
        const int iz = ir % basis.nplane + basis.startz_current;
        // FFT samples are at integer grid nodes, not voxel centers.
        const double fx = static_cast<double>(ix) / static_cast<double>(basis.nx);
        const double fy = static_cast<double>(iy) / static_cast<double>(basis.ny);
        const double fz = static_cast<double>(iz) / static_cast<double>(basis.nz);
        const double grid_x = fx * a1.x + fy * a2.x + fz * a3.x;
        const double grid_y = fx * a1.y + fy * a2.y + fz * a3.y;
        const double grid_z = fx * a1.z + fy * a2.z + fz * a3.z;
        positions[ir] = ModuleBase::Vector3<double>(grid_x, grid_y, grid_z);
    }
    return positions;
}

} // namespace ModuleSurchem
