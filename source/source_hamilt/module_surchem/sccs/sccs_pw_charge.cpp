#include "sccs_pw_charge.h"

#include "source_base/constants.h"
#include "source_basis/module_pw/pw_basis.h"

#include <cmath>
#include <complex>
#include <stdexcept>

namespace ModuleSccs
{
namespace
{

ModuleBase::Vector3<double> lattice_row(const ModuleBase::Matrix3& lattice,
                                        const int row,
                                        const double scale)
{
    if (row == 0)
    {
        return ModuleBase::Vector3<double>(scale * lattice.e11,
                                           scale * lattice.e12,
                                           scale * lattice.e13);
    }
    if (row == 1)
    {
        return ModuleBase::Vector3<double>(scale * lattice.e21,
                                           scale * lattice.e22,
                                           scale * lattice.e23);
    }
    return ModuleBase::Vector3<double>(scale * lattice.e31,
                                       scale * lattice.e32,
                                       scale * lattice.e33);
}

double dot(const ModuleBase::Vector3<double>& left, const ModuleBase::Vector3<double>& right)
{
    return left.x * right.x + left.y * right.y + left.z * right.z;
}

} // namespace

double validate_cubic_cell(const ModuleBase::Matrix3& lattice_vectors,
                           const double lattice_scale,
                           const double relative_tolerance)
{
    if (!std::isfinite(lattice_scale) || lattice_scale <= 0.0
        || !std::isfinite(relative_tolerance) || relative_tolerance <= 0.0)
    {
        throw std::invalid_argument("SCCS cubic-cell validation requires positive finite scale and tolerance");
    }
    const ModuleBase::Vector3<double> a1 = lattice_row(lattice_vectors, 0, lattice_scale);
    const ModuleBase::Vector3<double> a2 = lattice_row(lattice_vectors, 1, lattice_scale);
    const ModuleBase::Vector3<double> a3 = lattice_row(lattice_vectors, 2, lattice_scale);
    const double length1 = std::sqrt(dot(a1, a1));
    const double length2 = std::sqrt(dot(a2, a2));
    const double length3 = std::sqrt(dot(a3, a3));
    if (!std::isfinite(length1) || !std::isfinite(length2) || !std::isfinite(length3)
        || length1 <= 0.0 || length2 <= 0.0 || length3 <= 0.0)
    {
        throw std::invalid_argument("SCCS PCC cell vectors must have positive finite lengths");
    }
    const double length = (length1 + length2 + length3) / 3.0;
    if (std::abs(length1 - length) > relative_tolerance * length
        || std::abs(length2 - length) > relative_tolerance * length
        || std::abs(length3 - length) > relative_tolerance * length
        || std::abs(dot(a1, a2)) > relative_tolerance * length * length
        || std::abs(dot(a1, a3)) > relative_tolerance * length * length
        || std::abs(dot(a2, a3)) > relative_tolerance * length * length)
    {
        throw std::invalid_argument("zero-dimensional SCCS PCC requires an orthogonal equal-edge cubic cell");
    }
    return length;
}

ModuleBase::Vector3<double> cell_center(const ModuleBase::Matrix3& lattice_vectors,
                                       const double lattice_scale)
{
    const ModuleBase::Vector3<double> a1 = lattice_row(lattice_vectors, 0, lattice_scale);
    const ModuleBase::Vector3<double> a2 = lattice_row(lattice_vectors, 1, lattice_scale);
    const ModuleBase::Vector3<double> a3 = lattice_row(lattice_vectors, 2, lattice_scale);
    return ModuleBase::Vector3<double>(0.5 * (a1.x + a2.x + a3.x),
                                       0.5 * (a1.y + a2.y + a3.y),
                                       0.5 * (a1.z + a2.z + a3.z));
}

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
        const double fx = (static_cast<double>(ix) + 0.5) / static_cast<double>(basis.nx);
        const double fy = (static_cast<double>(iy) + 0.5) / static_cast<double>(basis.ny);
        const double fz = (static_cast<double>(iz) + 0.5) / static_cast<double>(basis.nz);
        positions[ir] = ModuleBase::Vector3<double>(fx * a1.x + fy * a2.x + fz * a3.x,
                                                    fx * a1.y + fy * a2.y + fz * a3.y,
                                                    fx * a1.z + fy * a2.z + fz * a3.z);
    }
    return positions;
}

std::vector<double> ionic_charge_from_local_potential(
    const std::vector<double>& local_potential_rydberg,
    const double total_ionic_charge,
    const double cell_volume,
    const double tpiba,
    const ModulePW::PW_Basis& basis)
{
    if (local_potential_rydberg.size() != static_cast<std::size_t>(basis.nrxx))
    {
        throw std::invalid_argument("SCCS local potential does not match the local PW real-space grid");
    }
    if (!std::isfinite(total_ionic_charge) || total_ionic_charge < 0.0
        || !std::isfinite(cell_volume) || cell_volume <= 0.0 || !std::isfinite(tpiba)
        || tpiba <= 0.0 || basis.npw <= 0 || basis.nrxx <= 0 || basis.gg == nullptr)
    {
        throw std::invalid_argument("SCCS ionic charge reconstruction inputs are invalid");
    }

    std::vector<std::complex<double>> potential_g(basis.npw);
    std::vector<std::complex<double>> ionic_g(basis.npw);
    basis.real2recip(local_potential_rydberg.data(), potential_g.data());
    const double tpiba2 = tpiba * tpiba;
    for (int ig = 0; ig < basis.npw; ++ig)
    {
        if (ig == basis.ig_gge0 || basis.gg[ig] == 0.0)
        {
            ionic_g[ig] = total_ionic_charge / cell_volume;
        }
        else
        {
            const double coulomb_rydberg
                = ModuleBase::e2 * ModuleBase::FOUR_PI / (tpiba2 * basis.gg[ig]);
            ionic_g[ig] = -potential_g[ig] / coulomb_rydberg;
        }
    }

    std::vector<double> ionic_density(basis.nrxx);
    basis.recip2real(ionic_g.data(), ionic_density.data());
    return ionic_density;
}

} // namespace ModuleSccs
