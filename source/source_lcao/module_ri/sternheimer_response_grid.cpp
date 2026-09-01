#include "source_lcao/module_ri/sternheimer_response_grid.h"

#include "source_basis/module_pw/pw_basis.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <map>
#include <stdexcept>

#ifdef __MPI
#include <mpi.h>
#endif

namespace ModuleRI
{
namespace
{

using ReciprocalIndex = std::array<int, 3>;

std::unique_ptr<ModulePW::PW_Basis> make_serial_basis_with_dimensions(const ModulePW::PW_Basis& reference,
                                                                      const double ecutwfc,
                                                                      const int fft_mode)
{
    auto basis = std::make_unique<ModulePW::PW_Basis>("cpu", "double");
#ifdef __MPI
    basis->initmpi(1, 0, MPI_COMM_SELF);
#endif
    basis->initgrids(reference.lat0, reference.latvec, reference.nx, reference.ny, reference.nz);
    basis->initparameters(false, 4.0 * ecutwfc);
    basis->fft_bundle.initfftmode(fft_mode);
    basis->setuptransform();
    basis->collect_local_pw();
    if (basis->nrxx != basis->nxyz)
    {
        throw std::runtime_error("Serial Sternheimer fine FFT basis does not own the complete real-space grid.");
    }
    return basis;
}

std::unique_ptr<ModulePW::PW_Basis> make_serial_basis_with_cutoff(const ModulePW::PW_Basis& reference,
                                                                  const double ecutwfc,
                                                                  const int fft_mode)
{
    auto basis = std::make_unique<ModulePW::PW_Basis>("cpu", "double");
#ifdef __MPI
    basis->initmpi(1, 0, MPI_COMM_SELF);
#endif
    basis->initgrids(reference.lat0, reference.latvec, 4.0 * ecutwfc);
    basis->initparameters(false, 4.0 * ecutwfc);
    basis->fft_bundle.initfftmode(fft_mode);
    basis->setuptransform();
    basis->collect_local_pw();
    if (basis->nrxx != basis->nxyz)
    {
        throw std::runtime_error("Serial Sternheimer response FFT basis does not own the complete real-space grid.");
    }
    return basis;
}

ReciprocalIndex reciprocal_index(const ModuleBase::Vector3<double>& vector)
{
    ReciprocalIndex index = {0, 0, 0};
    const std::array<double, 3> components = {vector.x, vector.y, vector.z};
    for (std::size_t direction = 0; direction != components.size(); ++direction)
    {
        const double rounded = std::round(components[direction]);
        if (std::abs(components[direction] - rounded) > 1.0e-10)
        {
            throw std::runtime_error("Sternheimer FFT basis returned a noninteger reciprocal-grid index.");
        }
        index[direction] = static_cast<int>(rounded);
    }
    return index;
}

} // namespace

SternheimerResponseGrid::SternheimerResponseGrid() = default;
SternheimerResponseGrid::~SternheimerResponseGrid() = default;
SternheimerResponseGrid::SternheimerResponseGrid(SternheimerResponseGrid&&) noexcept = default;
SternheimerResponseGrid& SternheimerResponseGrid::operator=(SternheimerResponseGrid&&) noexcept = default;

bool sternheimer_uses_independent_response_grid(const double response_ecutwfc, const double pbe_ecutwfc)
{
    if (pbe_ecutwfc <= 0.0)
    {
        throw std::invalid_argument("PBE ecutwfc must be positive before selecting a Sternheimer response grid.");
    }
    if (response_ecutwfc < 0.0)
    {
        throw std::invalid_argument("Sternheimer response ecutwfc must be zero or positive.");
    }
    const double tolerance = 1.0e-12 * std::max(1.0, pbe_ecutwfc);
    if (response_ecutwfc > pbe_ecutwfc + tolerance)
    {
        throw std::invalid_argument("Sternheimer response ecutwfc must not exceed the PBE ecutwfc.");
    }
    return response_ecutwfc > 0.0;
}

SternheimerResponseGrid make_sternheimer_response_grid(const ModulePW::PW_Basis& pbe_basis,
                                                       const double pbe_ecutwfc,
                                                       const double response_ecutwfc,
                                                       const int fft_mode)
{
    SternheimerResponseGrid result;
    result.basis = &pbe_basis;
    result.requested_ecutwfc = response_ecutwfc;
    if (!sternheimer_uses_independent_response_grid(response_ecutwfc, pbe_ecutwfc))
    {
        return result;
    }

    std::unique_ptr<ModulePW::PW_Basis> response_basis
        = make_serial_basis_with_cutoff(pbe_basis, response_ecutwfc, fft_mode);
    if (response_basis->nx == pbe_basis.nx && response_basis->ny == pbe_basis.ny && response_basis->nz == pbe_basis.nz)
    {
        return result;
    }

    result.serial_fine_basis = make_serial_basis_with_dimensions(pbe_basis, pbe_ecutwfc, fft_mode);
    result.serial_response_basis = std::move(response_basis);
    result.basis = result.serial_response_basis.get();
    result.independent = true;
    return result;
}

std::vector<double> restrict_sternheimer_real_field(const ModulePW::PW_Basis& fine_basis,
                                                    const ModulePW::PW_Basis& coarse_basis,
                                                    const std::vector<double>& fine_values)
{
    if (&fine_basis == &coarse_basis)
    {
        if (fine_values.size() != static_cast<std::size_t>(fine_basis.nxyz))
        {
            throw std::invalid_argument("Sternheimer real field size does not match the identity FFT grid.");
        }
        return fine_values;
    }
    if (fine_basis.nrxx != fine_basis.nxyz || coarse_basis.nrxx != coarse_basis.nxyz)
    {
        throw std::invalid_argument("Sternheimer spectral restriction requires serial complete FFT grids.");
    }
    if (fine_values.size() != static_cast<std::size_t>(fine_basis.nxyz))
    {
        throw std::invalid_argument("Sternheimer fine real field size does not match its FFT grid.");
    }

    std::vector<std::complex<double>> fine_coefficients(static_cast<std::size_t>(fine_basis.npw));
    fine_basis.real2recip(fine_values.data(), fine_coefficients.data());

    std::map<ReciprocalIndex, std::complex<double>> coefficients_by_index;
    for (int ig = 0; ig != fine_basis.npw; ++ig)
    {
        coefficients_by_index.emplace(reciprocal_index(fine_basis.gdirect[ig]), fine_coefficients[ig]);
    }

    std::vector<std::complex<double>> coarse_coefficients(static_cast<std::size_t>(coarse_basis.npw),
                                                          std::complex<double>(0.0, 0.0));
    for (int ig = 0; ig != coarse_basis.npw; ++ig)
    {
        const auto found = coefficients_by_index.find(reciprocal_index(coarse_basis.gdirect[ig]));
        if (found != coefficients_by_index.end())
        {
            coarse_coefficients[ig] = found->second;
        }
    }

    std::vector<double> coarse_values(static_cast<std::size_t>(coarse_basis.nxyz), 0.0);
    coarse_basis.recip2real(coarse_coefficients.data(), coarse_values.data());
    return coarse_values;
}

} // namespace ModuleRI
