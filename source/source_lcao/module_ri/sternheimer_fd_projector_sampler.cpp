#include "source_lcao/module_ri/sternheimer_fd_projector_sampler.h"

#include "source_base/math_ylmreal.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <stdexcept>

namespace
{

struct ProjectorChannel
{
    int radial_index = 0;
    int angular_momentum = 0;
    int magnetic_index = 0;
    int ylm_index = 0;
};

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
        throw std::invalid_argument("Sternheimer projector sampler requires positive grid dimensions.");
    }
    if (grid.hx <= 0.0 || grid.hy <= 0.0 || grid.hz <= 0.0)
    {
        throw std::invalid_argument("Sternheimer projector sampler requires positive grid spacings.");
    }
    bool has_nonzero_kpoint = false;
    for (const double coordinate : grid.kpoint)
    {
        if (!std::isfinite(coordinate))
        {
            throw std::invalid_argument("Sternheimer projector sampler requires finite reduced k-point coordinates.");
        }
        has_nonzero_kpoint = has_nonzero_kpoint || coordinate != 0.0;
    }
    if (!grid.periodic && has_nonzero_kpoint)
    {
        throw std::invalid_argument("Sternheimer projector sampler nonperiodic grids cannot use a Bloch k-point.");
    }
}

void validate_radial_set(const ModuleRI::SternheimerFDRadialProjectorSet& radial_set)
{
    const int nr = static_cast<int>(radial_set.radial_grid.size());
    const int nproj = static_cast<int>(radial_set.beta_radials.size());
    if (nr < 2)
    {
        throw std::invalid_argument("Sternheimer projector sampler requires at least two radial grid points.");
    }
    if (nproj == 0)
    {
        throw std::invalid_argument("Sternheimer projector sampler requires at least one radial projector.");
    }
    if (static_cast<int>(radial_set.angular_momenta.size()) != nproj)
    {
        throw std::invalid_argument("Sternheimer projector sampler angular-momentum count mismatch.");
    }
    if (static_cast<int>(radial_set.d_radial.size()) != nproj)
    {
        throw std::invalid_argument("Sternheimer projector sampler D matrix row count mismatch.");
    }
    if (radial_set.radial_grid.front() < 0.0)
    {
        throw std::invalid_argument("Sternheimer projector sampler radial grid must be non-negative.");
    }
    for (int ir = 1; ir != nr; ++ir)
    {
        if (radial_set.radial_grid[ir] <= radial_set.radial_grid[ir - 1])
        {
            throw std::invalid_argument("Sternheimer projector sampler radial grid must be strictly increasing.");
        }
    }
    for (int ip = 0; ip != nproj; ++ip)
    {
        if (radial_set.angular_momenta[ip] < 0)
        {
            throw std::invalid_argument("Sternheimer projector sampler angular momenta must be non-negative.");
        }
        if (static_cast<int>(radial_set.beta_radials[ip].size()) != nr)
        {
            throw std::invalid_argument("Sternheimer projector sampler radial beta size mismatch.");
        }
        if (static_cast<int>(radial_set.d_radial[ip].size()) != nproj)
        {
            throw std::invalid_argument("Sternheimer projector sampler D matrix must be square.");
        }
    }
}

std::vector<ProjectorChannel> build_channels(const std::vector<int>& angular_momenta)
{
    std::vector<ProjectorChannel> channels;
    for (int ip = 0; ip != static_cast<int>(angular_momenta.size()); ++ip)
    {
        const int l = angular_momenta[ip];
        for (int m_index = 0; m_index != 2 * l + 1; ++m_index)
        {
            ProjectorChannel channel;
            channel.radial_index = ip;
            channel.angular_momentum = l;
            channel.magnetic_index = m_index;
            channel.ylm_index = l * l + m_index;
            channels.push_back(channel);
        }
    }
    return channels;
}

struct PeriodicImageRange
{
    int first = 0;
    int last = 0;
};

std::array<PeriodicImageRange, 3> periodic_image_ranges(
    const std::array<double, 3>& displacement,
    const ModuleRI::SternheimerFDHamiltonian::Grid& grid,
    const double cutoff)
{
    std::array<PeriodicImageRange, 3> ranges{};
    if (!grid.periodic)
    {
        return ranges;
    }
    const ModuleRI::SternheimerFDLatticeVectors dual = ModuleRI::sternheimer_fd_grid_dual_vectors(grid);
    for (int direction = 0; direction != 3; ++direction)
    {
        double reduced_displacement = 0.0;
        double dual_norm_squared = 0.0;
        for (int component = 0; component != 3; ++component)
        {
            reduced_displacement += dual[direction][component] * displacement[component];
            dual_norm_squared += dual[direction][component] * dual[direction][component];
        }
        const double reduced_cutoff = std::sqrt(dual_norm_squared) * cutoff;
        ranges[direction].first = static_cast<int>(std::ceil(reduced_displacement - reduced_cutoff));
        ranges[direction].last = static_cast<int>(std::floor(reduced_displacement + reduced_cutoff));
    }
    return ranges;
}

double interpolate_radial(const std::vector<double>& radial_grid, const std::vector<double>& values, const double radius)
{
    if (radius < radial_grid.front() || radius > radial_grid.back())
    {
        return 0.0;
    }
    if (radius == radial_grid.front())
    {
        return values.front();
    }

    const std::vector<double>::const_iterator upper = std::upper_bound(radial_grid.begin(), radial_grid.end(), radius);
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

double abacus_betar_to_real_space_beta(const std::vector<double>& radial_grid,
                                       const ModuleBase::matrix& beta_radials,
                                       const std::vector<int>& angular_momenta,
                                       const int ip,
                                       const int ir)
{
    const double radius = radial_grid[ir];
    if (radius > 0.0)
    {
        return beta_radials(ip, ir) / radius;
    }

    if (angular_momenta[ip] > 0)
    {
        return 0.0;
    }

    for (int jr = ir + 1; jr != static_cast<int>(radial_grid.size()); ++jr)
    {
        if (radial_grid[jr] > 0.0)
        {
            return beta_radials(ip, jr) / radial_grid[jr];
        }
    }
    return 0.0;
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

} // namespace

namespace ModuleRI
{

SternheimerFDRadialProjectorSet
make_sternheimer_fd_radial_projector_set_from_abacus_matrices(const std::vector<double>& radial_grid,
                                                              const ModuleBase::matrix& beta_radials,
                                                              const std::vector<int>& angular_momenta,
                                                              const ModuleBase::matrix& d_radial)
{
    const int nproj = static_cast<int>(angular_momenta.size());
    const int nr = static_cast<int>(radial_grid.size());
    if (nproj <= 0)
    {
        throw std::invalid_argument("Sternheimer ABACUS projector set requires at least one projector.");
    }
    if (beta_radials.nr != nproj || beta_radials.nc != nr)
    {
        throw std::invalid_argument("Sternheimer ABACUS projector set expects betar(nbeta, mesh).");
    }
    if (d_radial.nr != nproj || d_radial.nc != nproj)
    {
        throw std::invalid_argument("Sternheimer ABACUS projector set expects dion(nbeta, nbeta).");
    }

    SternheimerFDRadialProjectorSet radial_set;
    radial_set.radial_grid = radial_grid;
    radial_set.angular_momenta = angular_momenta;
    radial_set.beta_radials.assign(nproj, std::vector<double>(nr, 0.0));
    radial_set.d_radial.assign(nproj, SternheimerFDNonlocalProjector::Vector(nproj, {0.0, 0.0}));

    for (int ip = 0; ip != nproj; ++ip)
    {
        for (int ir = 0; ir != nr; ++ir)
        {
            radial_set.beta_radials[ip][ir]
                = abacus_betar_to_real_space_beta(radial_grid, beta_radials, angular_momenta, ip, ir);
        }
        for (int jp = 0; jp != nproj; ++jp)
        {
            radial_set.d_radial[ip][jp] = d_radial(ip, jp);
        }
    }
    return radial_set;
}

SternheimerFDNonlocalProjector::ProjectorBlock
sample_sternheimer_fd_projector_block(const SternheimerFDRadialProjectorSet& radial_set,
                                      const SternheimerFDHamiltonian::Grid& grid,
                                      const ModuleBase::Vector3<double>& atom_position)
{
    validate_grid(grid);
    validate_radial_set(radial_set);

    const std::vector<ProjectorChannel> channels = build_channels(radial_set.angular_momenta);
    const int nchannel = static_cast<int>(channels.size());
    const int size = grid_size(grid);

    SternheimerFDNonlocalProjector::ProjectorBlock block;
    block.projectors.assign(nchannel, SternheimerFDNonlocalProjector::Vector(size, {0.0, 0.0}));
    block.d_matrix.assign(nchannel, SternheimerFDNonlocalProjector::Vector(nchannel, {0.0, 0.0}));

    int lmax = 0;
    for (const int l : radial_set.angular_momenta)
    {
        lmax = std::max(lmax, l);
    }

    std::vector<double> ylm;
    const double cutoff = radial_set.radial_grid.back();
    for (int iz = 0; iz != grid.nz; ++iz)
    {
        for (int iy = 0; iy != grid.ny; ++iy)
        {
            for (int ix = 0; ix != grid.nx; ++ix)
            {
                const std::array<double, 3> position
                    = sternheimer_fd_grid_cartesian_position(grid, ix, iy, iz);
                const std::array<double, 3> displacement{
                    position[0] - atom_position.x,
                    position[1] - atom_position.y,
                    position[2] - atom_position.z};
                const std::array<PeriodicImageRange, 3> image_ranges
                    = periodic_image_ranges(displacement, grid, cutoff);
                const int ir = grid_index(grid, ix, iy, iz);

                for (int rz = image_ranges[2].first; rz <= image_ranges[2].last; ++rz)
                {
                    for (int ry = image_ranges[1].first; ry <= image_ranges[1].last; ++ry)
                    {
                        for (int rx = image_ranges[0].first; rx <= image_ranges[0].last; ++rx)
                        {
                            const std::array<int, 3> image{rx, ry, rz};
                            const std::array<double, 3> translation
                                = sternheimer_fd_grid_lattice_translation(grid, image);
                            const double dx = displacement[0] - translation[0];
                            const double dy = displacement[1] - translation[1];
                            const double dz = displacement[2] - translation[2];
                            const double radius = std::sqrt(dx * dx + dy * dy + dz * dz);
                            if (radius > cutoff)
                            {
                                continue;
                            }

                            evaluate_real_spherical_harmonics(lmax, dx, dy, dz, ylm);
                            const std::complex<double> phase = sternheimer_bloch_phase(grid.kpoint, image);
                            for (int ich = 0; ich != nchannel; ++ich)
                            {
                                const ProjectorChannel& channel = channels[ich];
                                const double radial_value
                                    = interpolate_radial(radial_set.radial_grid,
                                                         radial_set.beta_radials[channel.radial_index],
                                                         radius);
                                block.projectors[ich][ir] += phase * radial_value * ylm[channel.ylm_index];
                            }
                        }
                    }
                }
            }
        }
    }

    for (int ich = 0; ich != nchannel; ++ich)
    {
        const ProjectorChannel& left = channels[ich];
        for (int jch = 0; jch != nchannel; ++jch)
        {
            const ProjectorChannel& right = channels[jch];
            if (left.angular_momentum == right.angular_momentum && left.magnetic_index == right.magnetic_index)
            {
                block.d_matrix[ich][jch] = radial_set.d_radial[left.radial_index][right.radial_index];
            }
        }
    }

    return block;
}

} // namespace ModuleRI
