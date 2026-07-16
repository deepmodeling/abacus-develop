#include "source_lcao/module_ri/sternheimer_abfs_perturbation.h"

#include "source_base/math_ylmreal.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
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

struct PeriodicImageRange
{
    int first = 0;
    int last = 0;
};

PeriodicImageRange periodic_image_range(const double displacement,
                                        const double length,
                                        const double cutoff,
                                        const bool periodic)
{
    if (!periodic)
    {
        return {0, 0};
    }
    return {static_cast<int>(std::ceil((displacement - cutoff) / length)),
            static_cast<int>(std::floor((displacement + cutoff) / length))};
}

void validate_qpoint(const ModuleRI::SternheimerReducedKPoint& qpoint, const bool periodic)
{
    bool has_nonzero_qpoint = false;
    for (const double coordinate: qpoint)
    {
        if (!std::isfinite(coordinate))
        {
            throw std::invalid_argument("Sternheimer ABFS perturbation requires a finite reduced q point.");
        }
        has_nonzero_qpoint = has_nonzero_qpoint || coordinate != 0.0;
    }
    if (!periodic && has_nonzero_qpoint)
    {
        throw std::invalid_argument("Sternheimer ABFS nonperiodic grids cannot use a Bloch q point.");
    }
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

std::vector<SternheimerABFBlochGridChannel> sample_sternheimer_abf_bloch_grid_channels(
    const std::vector<std::vector<SternheimerRadialPerturbation>>& radials_by_type,
    const std::vector<int>& atom_types,
    const std::vector<ModuleBase::Vector3<double>>& atom_positions,
    const SternheimerFDHamiltonian::Grid& grid,
    const SternheimerReducedKPoint& qpoint,
    const int max_channels)
{
    validate_grid(grid);
    validate_qpoint(qpoint, grid.periodic);
    if (atom_types.size() != atom_positions.size())
    {
        throw std::invalid_argument("Sternheimer ABFS perturbation atom type/position count mismatch.");
    }

    std::vector<SternheimerABFBlochGridChannel> channels;
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

                SternheimerABFBlochGridChannel channel;
                channel.channel_index = channel_index++;
                channel.atom_index = static_cast<int>(iat);
                channel.atom_local_index = atom_local_index++;
                channel.type_index = type;
                channel.angular_momentum = radial.angular_momentum;
                channel.radial_index = radial.radial_index;
                channel.magnetic_index = m_index;
                channel.label = radial.label;
                channel.potential_r.assign(size, std::complex<double>(0.0, 0.0));

                const int ylm_index = radial.angular_momentum * radial.angular_momentum + m_index;
                const double cutoff = radial.radial_grid.back();
                for (int iz = 0; iz != grid.nz; ++iz)
                {
                    for (int iy = 0; iy != grid.ny; ++iy)
                    {
                        for (int ix = 0; ix != grid.nx; ++ix)
                        {
                            const double displacement_x = ix * grid.hx - atom_positions[iat].x;
                            const double displacement_y = iy * grid.hy - atom_positions[iat].y;
                            const double displacement_z = iz * grid.hz - atom_positions[iat].z;
                            const PeriodicImageRange image_x
                                = periodic_image_range(displacement_x, lx, cutoff, grid.periodic);
                            const PeriodicImageRange image_y
                                = periodic_image_range(displacement_y, ly, cutoff, grid.periodic);
                            const PeriodicImageRange image_z
                                = periodic_image_range(displacement_z, lz, cutoff, grid.periodic);
                            const int ir = grid_index(grid, ix, iy, iz);

                            for (int rz = image_z.first; rz <= image_z.last; ++rz)
                            {
                                for (int ry = image_y.first; ry <= image_y.last; ++ry)
                                {
                                    for (int rx = image_x.first; rx <= image_x.last; ++rx)
                                    {
                                        const double dx = displacement_x - rx * lx;
                                        const double dy = displacement_y - ry * ly;
                                        const double dz = displacement_z - rz * lz;
                                        const double radius = std::sqrt(dx * dx + dy * dy + dz * dz);
                                        if (radius > cutoff)
                                        {
                                            continue;
                                        }
                                        evaluate_real_spherical_harmonics(
                                            radial.angular_momentum, dx, dy, dz, ylm);
                                        const double radial_angular_value
                                            = interpolate_radial(
                                                  radial.radial_grid, radial.radial_values, radius)
                                              * ylm[ylm_index];
                                        const std::complex<double> phase
                                            = sternheimer_bloch_phase(
                                                qpoint, std::array<int, 3>{rx, ry, rz});
                                        channel.potential_r[ir] += phase * radial_angular_value;
                                    }
                                }
                            }
                            channel.max_abs = std::max(channel.max_abs, std::abs(channel.potential_r[ir]));
                        }
                    }
                }
                channels.push_back(std::move(channel));
            }
        }
    }
    return channels;
}

std::vector<SternheimerABFGridChannel> sample_sternheimer_abf_grid_channels(
    const std::vector<std::vector<SternheimerRadialPerturbation>>& radials_by_type,
    const std::vector<int>& atom_types,
    const std::vector<ModuleBase::Vector3<double>>& atom_positions,
    const SternheimerFDHamiltonian::Grid& grid,
    const int max_channels)
{
    const std::vector<SternheimerABFBlochGridChannel> bloch_channels
        = sample_sternheimer_abf_bloch_grid_channels(
            radials_by_type, atom_types, atom_positions, grid, {0.0, 0.0, 0.0}, max_channels);

    std::vector<SternheimerABFGridChannel> channels;
    channels.reserve(bloch_channels.size());
    for (const SternheimerABFBlochGridChannel& bloch: bloch_channels)
    {
        SternheimerABFGridChannel channel;
        channel.channel_index = bloch.channel_index;
        channel.atom_index = bloch.atom_index;
        channel.atom_local_index = bloch.atom_local_index;
        channel.type_index = bloch.type_index;
        channel.angular_momentum = bloch.angular_momentum;
        channel.radial_index = bloch.radial_index;
        channel.magnetic_index = bloch.magnetic_index;
        channel.label = bloch.label;
        channel.max_abs = bloch.max_abs;
        channel.potential_r.reserve(bloch.potential_r.size());
        for (const std::complex<double>& value: bloch.potential_r)
        {
            if (std::abs(value.imag()) > 1.0e-13 * std::max(1.0, channel.max_abs))
            {
                throw std::runtime_error("Sternheimer ABFS Gamma potential acquired an unexpected imaginary part.");
            }
            channel.potential_r.push_back(value.real());
        }
        channels.push_back(std::move(channel));
    }
    return channels;
}

} // namespace ModuleRI
