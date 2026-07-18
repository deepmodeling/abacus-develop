#include "source_lcao/module_ri/sternheimer_abfs_perturbation.h"

#include "source_base/constants.h"
#include "source_base/math_ylmreal.h"

#include <fftw3.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <new>
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

int signed_fft_index(const int index, const int count)
{
    return index < (count + 1) / 2 ? index : index - count;
}

bool is_gamma_qpoint(const ModuleRI::SternheimerReducedKPoint& qpoint)
{
    constexpr double tolerance = 1.0e-14;
    return std::all_of(qpoint.begin(), qpoint.end(), [](const double value) {
        return std::abs(value) <= tolerance;
    });
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
                            const std::array<double, 3> position
                                = sternheimer_fd_grid_cartesian_position(grid, ix, iy, iz);
                            const std::array<double, 3> displacement{
                                position[0] - atom_positions[iat].x,
                                position[1] - atom_positions[iat].y,
                                position[2] - atom_positions[iat].z};
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
                                        evaluate_real_spherical_harmonics(
                                            radial.angular_momentum, dx, dy, dz, ylm);
                                        const double radial_angular_value
                                            = interpolate_radial(
                                                  radial.radial_grid, radial.radial_values, radius)
                                              * ylm[ylm_index];
                                        const std::complex<double> phase
                                            = sternheimer_bloch_phase(qpoint, image);
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

std::vector<SternheimerABFBlochGridChannel> solve_sternheimer_abf_periodic_full_coulomb(
    const std::vector<SternheimerABFBlochGridChannel>& density_channels,
    const SternheimerFDHamiltonian::Grid& grid,
    const SternheimerReducedKPoint& qpoint)
{
    validate_grid(grid);
    validate_qpoint(qpoint, grid.periodic);
    if (!grid.periodic)
    {
        throw std::invalid_argument("Sternheimer periodic Poisson solve requires a periodic grid.");
    }
    if (is_gamma_qpoint(qpoint))
    {
        throw std::invalid_argument(
            "Sternheimer periodic full-Coulomb body solve requires a nonzero q point.");
    }

    const int size = grid_size(grid);
    for (const SternheimerABFBlochGridChannel& density: density_channels)
    {
        if (density.potential_r.size() != static_cast<std::size_t>(size))
        {
            throw std::invalid_argument(
                "Sternheimer periodic Poisson density size does not match the grid.");
        }
    }
    if (density_channels.empty())
    {
        return {};
    }

    fftw_complex* buffer = fftw_alloc_complex(static_cast<std::size_t>(size));
    if (buffer == nullptr)
    {
        throw std::bad_alloc();
    }
    fftw_plan forward = fftw_plan_dft_3d(
        grid.nx, grid.ny, grid.nz, buffer, buffer, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_plan backward = fftw_plan_dft_3d(
        grid.nx, grid.ny, grid.nz, buffer, buffer, FFTW_BACKWARD, FFTW_ESTIMATE);
    if (forward == nullptr || backward == nullptr)
    {
        if (forward != nullptr)
        {
            fftw_destroy_plan(forward);
        }
        if (backward != nullptr)
        {
            fftw_destroy_plan(backward);
        }
        fftw_free(buffer);
        throw std::runtime_error("Failed to initialize the Sternheimer periodic Poisson FFT.");
    }

    const SternheimerFDLatticeVectors dual = sternheimer_fd_grid_dual_vectors(grid);
    std::vector<SternheimerABFBlochGridChannel> potentials;
    potentials.reserve(density_channels.size());
    try
    {
        for (const SternheimerABFBlochGridChannel& density: density_channels)
        {
            for (int ix = 0; ix != grid.nx; ++ix)
            {
                for (int iy = 0; iy != grid.ny; ++iy)
                {
                    for (int iz = 0; iz != grid.nz; ++iz)
                    {
                        const int ir = grid_index(grid, ix, iy, iz);
                        const double phase_angle = -ModuleBase::TWO_PI
                                                   * (qpoint[0] * static_cast<double>(ix) / grid.nx
                                                      + qpoint[1] * static_cast<double>(iy) / grid.ny
                                                      + qpoint[2] * static_cast<double>(iz) / grid.nz);
                        const std::complex<double> periodic_value
                            = std::exp(std::complex<double>(0.0, phase_angle))
                              * density.potential_r[static_cast<std::size_t>(ir)];
                        buffer[ir][0] = periodic_value.real();
                        buffer[ir][1] = periodic_value.imag();
                    }
                }
            }
            fftw_execute(forward);

            for (int ix = 0; ix != grid.nx; ++ix)
            {
                for (int iy = 0; iy != grid.ny; ++iy)
                {
                    for (int iz = 0; iz != grid.nz; ++iz)
                    {
                        const std::array<double, 3> reduced_wavevector{
                            static_cast<double>(signed_fft_index(ix, grid.nx)) + qpoint[0],
                            static_cast<double>(signed_fft_index(iy, grid.ny)) + qpoint[1],
                            static_cast<double>(signed_fft_index(iz, grid.nz)) + qpoint[2]};
                        std::array<double, 3> wavevector{};
                        for (int direction = 0; direction != 3; ++direction)
                        {
                            for (int component = 0; component != 3; ++component)
                            {
                                wavevector[component]
                                    += ModuleBase::TWO_PI * reduced_wavevector[direction] * dual[direction][component];
                            }
                        }
                        const double wavevector_squared
                            = wavevector[0] * wavevector[0]
                              + wavevector[1] * wavevector[1]
                              + wavevector[2] * wavevector[2];
                        if (wavevector_squared <= 1.0e-28)
                        {
                            throw std::runtime_error(
                                "Sternheimer periodic Poisson solve encountered a zero G+q vector.");
                        }
                        const double factor
                            = ModuleBase::FOUR_PI / (wavevector_squared * static_cast<double>(size));
                        const int ig = grid_index(grid, ix, iy, iz);
                        buffer[ig][0] *= factor;
                        buffer[ig][1] *= factor;
                    }
                }
            }
            fftw_execute(backward);

            SternheimerABFBlochGridChannel potential = density;
            potential.max_abs = 0.0;
            for (int ix = 0; ix != grid.nx; ++ix)
            {
                for (int iy = 0; iy != grid.ny; ++iy)
                {
                    for (int iz = 0; iz != grid.nz; ++iz)
                    {
                        const int ir = grid_index(grid, ix, iy, iz);
                        const double phase_angle = ModuleBase::TWO_PI
                                                   * (qpoint[0] * static_cast<double>(ix) / grid.nx
                                                      + qpoint[1] * static_cast<double>(iy) / grid.ny
                                                      + qpoint[2] * static_cast<double>(iz) / grid.nz);
                        const std::complex<double> periodic_potential(buffer[ir][0], buffer[ir][1]);
                        potential.potential_r[static_cast<std::size_t>(ir)]
                            = std::exp(std::complex<double>(0.0, phase_angle)) * periodic_potential;
                        potential.max_abs = std::max(
                            potential.max_abs,
                            std::abs(potential.potential_r[static_cast<std::size_t>(ir)]));
                    }
                }
            }
            potentials.push_back(std::move(potential));
        }
    }
    catch (...)
    {
        fftw_destroy_plan(forward);
        fftw_destroy_plan(backward);
        fftw_free(buffer);
        throw;
    }
    fftw_destroy_plan(forward);
    fftw_destroy_plan(backward);
    fftw_free(buffer);
    return potentials;
}

} // namespace ModuleRI
