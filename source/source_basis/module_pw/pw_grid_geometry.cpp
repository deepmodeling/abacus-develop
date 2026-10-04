#include "pw_grid_geometry.h"

#include "pw_basis.h"
#include "source_base/matrix3.h"

#include <cmath>

namespace ModulePW
{
bool grid_positions(const PW_Basis& basis,
                     const ModuleBase::Matrix3& lattice,
                     const double lattice_scale,
                     std::vector<ModuleBase::Vector3<double>>& positions,
                     std::string& error)
{
    error.clear();
    if (basis.nx <= 0 || basis.ny <= 0 || basis.nz <= 0 || basis.nplane < 0
        || basis.nrxx != basis.nx * basis.ny * basis.nplane
        || basis.startz_current < 0 || basis.startz_current + basis.nplane > basis.nz
        || !std::isfinite(lattice_scale) || lattice_scale <= 0.0)
    {
        error = "PW grid coordinates require valid dimensions, slab offset and lattice scale";
        return false;
    }
    const ModuleBase::Vector3<double> a1(lattice.e11, lattice.e12, lattice.e13);
    const ModuleBase::Vector3<double> a2(lattice.e21, lattice.e22, lattice.e23);
    const ModuleBase::Vector3<double> a3(lattice.e31, lattice.e32, lattice.e33);
    std::vector<ModuleBase::Vector3<double>> candidate(basis.nrxx);
    for (int ir = 0; ir < basis.nrxx; ++ir)
    {
        const int ix = ir / (basis.ny * basis.nplane);
        const int iy = ir / basis.nplane - ix * basis.ny;
        const int iz = ir % basis.nplane + basis.startz_current;
        const double fx = static_cast<double>(ix) / basis.nx;
        const double fy = static_cast<double>(iy) / basis.ny;
        const double fz = static_cast<double>(iz) / basis.nz;
        const ModuleBase::Vector3<double> position = (a1 * fx + a2 * fy + a3 * fz) * lattice_scale;
        if (!std::isfinite(position.x) || !std::isfinite(position.y) || !std::isfinite(position.z))
        {
            error = "PW grid coordinates must be finite";
            return false;
        }
        candidate[ir] = position;
    }
    positions.swap(candidate);
    return true;
}
} // namespace ModulePW
