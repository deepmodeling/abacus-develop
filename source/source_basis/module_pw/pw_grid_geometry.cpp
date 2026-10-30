#include "pw_grid_geometry.h"

#include "pw_basis.h"
#include "source_base/matrix3.h"

namespace ModulePW
{
void grid_positions(const PW_Basis& basis,
                    const ModuleBase::Matrix3& lattice,
                    const double lattice_scale,
                    std::vector<ModuleBase::Vector3<double>>& positions)
{
    const ModuleBase::Vector3<double> a1(lattice.e11, lattice.e12, lattice.e13);
    const ModuleBase::Vector3<double> a2(lattice.e21, lattice.e22, lattice.e23);
    const ModuleBase::Vector3<double> a3(lattice.e31, lattice.e32, lattice.e33);
    positions.resize(basis.nrxx);
    for (int ir = 0; ir < basis.nrxx; ++ir)
    {
        const int ix = ir / (basis.ny * basis.nplane);
        const int iy = ir / basis.nplane - ix * basis.ny;
        const int iz = ir % basis.nplane + basis.startz_current;
        const double fx = static_cast<double>(ix) / basis.nx;
        const double fy = static_cast<double>(iy) / basis.ny;
        const double fz = static_cast<double>(iz) / basis.nz;
        positions[ir] = (a1 * fx + a2 * fy + a3 * fz) * lattice_scale;
    }
}
} // namespace ModulePW
