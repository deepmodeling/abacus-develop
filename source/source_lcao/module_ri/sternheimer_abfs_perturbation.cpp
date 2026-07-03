#include "source_lcao/module_ri/sternheimer_abfs_perturbation.h"

#include "source_base/math_ylmreal.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace
{

int grid_size(const ModuleRI::SternheimerFDHamiltonian::Grid& grid)
{
    return grid.nx * grid.ny * grid.nz;
}

int grid_index(const ModuleRI::SternheimerFDHamiltonian::Grid& grid, const int ix, const int iy, const int iz)
{
    return (ix * grid.ny + iy) * grid.nz + iz;
}

void validate_grid(const ModuleRI::SternheimerFDHamiltonian::Grid& grid)
{
    if (grid.nx <= 0 || grid.ny <= 0 || grid.nz <= 0)
    {
        throw std::invalid_argument("Sternheimer ABFS perturbation requires positive grid dimensions.");
    }
    if (grid.hx <= 0.0 || grid.hy <= 0.0 || grid.hz <= 0.0)
    {
        throw std::invalid_argument("Sternheimer ABFS perturbation requires positive grid spacings.");
    }
}

double minimum_image_displacement(double displacement, const double length)
{
    if (length > 0.0)
    {
        displacement -= length * std::round(displacement / length);
    }
    return displacement;
}

double interpolate_radial(const std::vector<double>& radial_grid, const std::vector<double>& values, const double radius)
{
    if (radial_grid.size() != values.size() || radial_grid.size() < 2)
    {
        throw std::invalid_argument("Sternheimer ABFS radial interpolation found an invalid radial function.");
    }
    if (radius < radial_grid.front() || radius > radial_grid.back())
    {
        return 0.0;
    }
    if (radius == radial_grid.front())
    {
        return values.front();
    }

    const auto upper = std::upper_bound(radial_grid.begin(), radial_grid.end(), radius);
    if (upper == radial_grid.end())
    {
        return values.back();
    }
    const int hi = static_cast<int>(upper - radial_grid.begin());
    const int lo = hi - 1;
    const double width = radial_grid[hi] - radial_grid[lo];
    const double t = (radius - radial_grid[lo]) / width;
    return (1.0 - t) * values[lo] + t * values[hi];
}

void evaluate_real_spherical_harmonics(const int lmax,
                                       const double x,
                                       const double y,
                                       const double z,
                                       std::vector<double>& ylm)
{
    ylm.assign((lmax + 1) * (lmax + 1), 0.0);
    if (lmax == 0 || x * x + y * y + z * z > 1.0e-28)
    {
        ModuleBase::YlmReal::rlylm(lmax, x, y, z, ylm.data());
    }
    else
    {
        ylm[0] = 0.28209479177387814347;
    }
}

void validate_radial(const ModuleRI::SternheimerRadialPerturbation& radial)
{
    if (radial.angular_momentum < 0)
    {
        throw std::invalid_argument("Sternheimer ABFS perturbation found negative angular momentum.");
    }
    if (radial.radial_grid.size() != radial.radial_values.size() || radial.radial_grid.size() < 2)
    {
        throw std::invalid_argument("Sternheimer ABFS perturbation found inconsistent radial data.");
    }
    for (std::size_t ir = 1; ir != radial.radial_grid.size(); ++ir)
    {
        if (radial.radial_grid[ir] <= radial.radial_grid[ir - 1])
        {
            throw std::invalid_argument("Sternheimer ABFS perturbation radial grid is not strictly increasing.");
        }
    }
}

} // namespace

namespace ModuleRI
{

std::vector<std::vector<SternheimerRadialPerturbation>> make_sternheimer_radial_perturbations_from_orbitals(
    const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>>& orbitals)
{
    std::vector<std::vector<SternheimerRadialPerturbation>> radials_by_type(orbitals.size());
    for (std::size_t type = 0; type != orbitals.size(); ++type)
    {
        for (std::size_t l = 0; l != orbitals[type].size(); ++l)
        {
            for (std::size_t n = 0; n != orbitals[type][l].size(); ++n)
            {
                const Numerical_Orbital_Lm& orbital = orbitals[type][l][n];
                SternheimerRadialPerturbation radial;
                radial.type_index = static_cast<int>(type);
                radial.angular_momentum = orbital.getL();
                radial.radial_index = orbital.getChi();
                radial.label = orbital.getLabel();
                radial.radial_grid = orbital.get_r_radial();
                radial.radial_values = orbital.get_psi();
                validate_radial(radial);
                radials_by_type[type].push_back(std::move(radial));
            }
        }
    }
    return radials_by_type;
}

std::vector<SternheimerABFGridChannel> sample_sternheimer_abf_grid_channels(
    const std::vector<std::vector<SternheimerRadialPerturbation>>& radials_by_type,
    const std::vector<int>& atom_types,
    const std::vector<ModuleBase::Vector3<double>>& atom_positions,
    const SternheimerFDHamiltonian::Grid& grid,
    const int max_channels)
{
    validate_grid(grid);
    if (atom_types.size() != atom_positions.size())
    {
        throw std::invalid_argument("Sternheimer ABFS perturbation atom type/position count mismatch.");
    }

    std::vector<SternheimerABFGridChannel> channels;
    const int size = grid_size(grid);
    const double lx = grid.nx * grid.hx;
    const double ly = grid.ny * grid.hy;
    const double lz = grid.nz * grid.hz;
    int channel_index = 0;

    for (std::size_t iat = 0; iat != atom_types.size(); ++iat)
    {
        const int type = atom_types[iat];
        if (type < 0 || type >= static_cast<int>(radials_by_type.size()))
        {
            throw std::invalid_argument("Sternheimer ABFS perturbation atom type is out of range.");
        }
        int atom_local_index = 0;
        for (const SternheimerRadialPerturbation& radial: radials_by_type[type])
        {
            validate_radial(radial);
            std::vector<double> ylm;
            evaluate_real_spherical_harmonics(radial.angular_momentum, 0.0, 0.0, 0.0, ylm);
            for (int m_index = 0; m_index != 2 * radial.angular_momentum + 1; ++m_index)
            {
                if (max_channels > 0 && static_cast<int>(channels.size()) >= max_channels)
                {
                    return channels;
                }

                SternheimerABFGridChannel channel;
                channel.channel_index = channel_index++;
                channel.atom_index = static_cast<int>(iat);
                channel.atom_local_index = atom_local_index++;
                channel.type_index = type;
                channel.angular_momentum = radial.angular_momentum;
                channel.radial_index = radial.radial_index;
                channel.magnetic_index = m_index;
                channel.potential_r.assign(size, 0.0);

                const int ylm_index = radial.angular_momentum * radial.angular_momentum + m_index;
                for (int iz = 0; iz != grid.nz; ++iz)
                {
                    for (int iy = 0; iy != grid.ny; ++iy)
                    {
                        for (int ix = 0; ix != grid.nx; ++ix)
                        {
                            double dx = ix * grid.hx - atom_positions[iat].x;
                            double dy = iy * grid.hy - atom_positions[iat].y;
                            double dz = iz * grid.hz - atom_positions[iat].z;
                            if (grid.periodic)
                            {
                                dx = minimum_image_displacement(dx, lx);
                                dy = minimum_image_displacement(dy, ly);
                                dz = minimum_image_displacement(dz, lz);
                            }
                            const double radius = std::sqrt(dx * dx + dy * dy + dz * dz);
                            evaluate_real_spherical_harmonics(radial.angular_momentum, dx, dy, dz, ylm);
                            const double value
                                = interpolate_radial(radial.radial_grid, radial.radial_values, radius) * ylm[ylm_index];
                            const int ir = grid_index(grid, ix, iy, iz);
                            channel.potential_r[ir] = value;
                            channel.max_abs = std::max(channel.max_abs, std::abs(value));
                        }
                    }
                }
                channels.push_back(std::move(channel));
            }
        }
    }
    return channels;
}

} // namespace ModuleRI
