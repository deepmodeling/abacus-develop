#include "source_lcao/module_ri/sternheimer_response_grid.h"

#include "source_basis/module_pw/pw_basis.h"

#include <fftw3.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <map>
#include <new>
#include <stdexcept>

#ifdef __MPI
#include <mpi.h>
#endif

namespace ModuleRI
{
namespace
{

using ReciprocalIndex = std::array<int, 3>;

class FFTWBuffer
{
  public:
    FFTWBuffer(const int nx, const int ny, const int nz, const int direction)
        : size_(static_cast<std::size_t>(nx) * static_cast<std::size_t>(ny) * static_cast<std::size_t>(nz))
    {
        data_ = fftw_alloc_complex(size_);
        if (data_ == nullptr)
        {
            throw std::bad_alloc();
        }
        for (std::size_t index = 0; index != size_; ++index)
        {
            data_[index][0] = 0.0;
            data_[index][1] = 0.0;
        }
#ifdef _OPENMP
#pragma omp critical(sternheimer_fftw_plan_management)
#endif
        {
            plan_ = fftw_plan_dft_3d(nx, ny, nz, data_, data_, direction, FFTW_ESTIMATE);
        }
        if (plan_ == nullptr)
        {
            fftw_free(data_);
            data_ = nullptr;
            throw std::runtime_error("Failed to initialize a Sternheimer response-grid FFT.");
        }
    }

    ~FFTWBuffer()
    {
#ifdef _OPENMP
#pragma omp critical(sternheimer_fftw_plan_management)
#endif
        {
            if (plan_ != nullptr)
            {
                fftw_destroy_plan(plan_);
            }
        }
        fftw_free(data_);
    }

    FFTWBuffer(const FFTWBuffer&) = delete;
    FFTWBuffer& operator=(const FFTWBuffer&) = delete;

    fftw_complex* data()
    {
        return data_;
    }

    const fftw_complex* data() const
    {
        return data_;
    }

    void execute() const
    {
        fftw_execute(plan_);
    }

  private:
    std::size_t size_ = 0;
    fftw_complex* data_ = nullptr;
    fftw_plan plan_ = nullptr;
};

std::size_t fft_index(const int ix, const int iy, const int iz, const int ny, const int nz)
{
    return static_cast<std::size_t>((ix * ny + iy) * nz + iz);
}

int signed_frequency(const int index, const int dimension)
{
    return index <= dimension / 2 ? index : index - dimension;
}

int frequency_index(const int frequency, const int dimension)
{
    return frequency >= 0 ? frequency : frequency + dimension;
}

struct FrequencySelection
{
    std::array<int, 2> values = {0, 0};
    int size = 1;
};

FrequencySelection fine_frequencies_for_coarse_index(const int coarse_index,
                                                      const int coarse_dimension,
                                                      const int fine_dimension)
{
    const int frequency = signed_frequency(coarse_index, coarse_dimension);
    if (coarse_dimension % 2 == 0 && coarse_index == coarse_dimension / 2
        && fine_dimension > coarse_dimension)
    {
        return {{frequency, -frequency}, 2};
    }
    return {{frequency, 0}, 1};
}

std::unique_ptr<ModulePW::PW_Basis> make_serial_basis_with_dimensions(
    const ModulePW::PW_Basis& reference,
    const double ecutwfc,
    const std::array<int, 3>& dimensions,
    const int fft_mode)
{
    auto basis = std::make_unique<ModulePW::PW_Basis>("cpu", "double");
#ifdef __MPI
    basis->initmpi(1, 0, MPI_COMM_SELF);
#endif
    basis->initgrids(reference.lat0, reference.latvec, dimensions[0], dimensions[1], dimensions[2]);
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

bool explicit_dimensions_requested(const std::array<int, 3>& dimensions)
{
    const int positive = static_cast<int>(std::count_if(dimensions.begin(), dimensions.end(), [](const int value) {
        return value > 0;
    }));
    if (std::any_of(dimensions.begin(), dimensions.end(), [](const int value) { return value < 0; }))
    {
        throw std::invalid_argument("Explicit Sternheimer response-grid dimensions must be zero or positive.");
    }
    if (positive != 0 && positive != 3)
    {
        throw std::invalid_argument("Explicit Sternheimer response-grid dimensions must be all zero or all positive.");
    }
    return positive == 3;
}

void validate_explicit_dimensions(const ModulePW::PW_Basis& pbe_basis,
                                  const std::array<int, 3>& dimensions,
                                  const int fd_order)
{
    if (fd_order != 2 && fd_order != 4 && fd_order != 6 && fd_order != 8)
    {
        throw std::invalid_argument("Sternheimer finite-difference order must be 2, 4, 6, or 8.");
    }
    const std::array<int, 3> pbe_dimensions = {pbe_basis.nx, pbe_basis.ny, pbe_basis.nz};
    const int minimum_dimension = fd_order + 1;
    for (std::size_t direction = 0; direction != dimensions.size(); ++direction)
    {
        if (dimensions[direction] < minimum_dimension)
        {
            throw std::invalid_argument("Each explicit Sternheimer response-grid dimension must exceed twice the "
                                        "finite-difference stencil radius.");
        }
        if (dimensions[direction] > pbe_dimensions[direction])
        {
            throw std::invalid_argument("Explicit Sternheimer response-grid dimensions must not exceed the PBE grid.");
        }
    }
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

const char* sternheimer_response_grid_source_name(const SternheimerResponseGridSource source)
{
    switch (source)
    {
        case SternheimerResponseGridSource::Pbe:
            return "pbe";
        case SternheimerResponseGridSource::Cutoff:
            return "cutoff";
        case SternheimerResponseGridSource::Explicit:
            return "explicit";
    }
    throw std::invalid_argument("Unknown Sternheimer response-grid source.");
}

SternheimerResponseGrid make_sternheimer_response_grid(const ModulePW::PW_Basis& pbe_basis,
                                                       const double pbe_ecutwfc,
                                                       const double response_ecutwfc,
                                                       const int fft_mode,
                                                       const std::array<int, 3> response_dimensions,
                                                       const int fd_order)
{
    SternheimerResponseGrid result;
    result.basis = &pbe_basis;
    result.requested_ecutwfc = response_ecutwfc;
    result.requested_dimensions = response_dimensions;
    const bool uses_cutoff = sternheimer_uses_independent_response_grid(response_ecutwfc, pbe_ecutwfc);
    const bool uses_explicit_dimensions = explicit_dimensions_requested(response_dimensions);
    if (uses_cutoff && uses_explicit_dimensions)
    {
        throw std::invalid_argument(
            "sternheimer_response_ecutwfc and explicit response-grid dimensions are mutually exclusive.");
    }
    if (!uses_cutoff && !uses_explicit_dimensions)
    {
        return result;
    }

    std::unique_ptr<ModulePW::PW_Basis> response_basis;
    if (uses_explicit_dimensions)
    {
        validate_explicit_dimensions(pbe_basis, response_dimensions, fd_order);
        response_basis
            = make_serial_basis_with_dimensions(pbe_basis, pbe_ecutwfc, response_dimensions, fft_mode);
        result.source = SternheimerResponseGridSource::Explicit;
    }
    else
    {
        response_basis = make_serial_basis_with_cutoff(pbe_basis, response_ecutwfc, fft_mode);
        result.source = SternheimerResponseGridSource::Cutoff;
    }
    if (response_basis->nx == pbe_basis.nx && response_basis->ny == pbe_basis.ny && response_basis->nz == pbe_basis.nz)
    {
        result.source = SternheimerResponseGridSource::Pbe;
        return result;
    }

    result.serial_fine_basis = make_serial_basis_with_dimensions(
        pbe_basis, pbe_ecutwfc, {pbe_basis.nx, pbe_basis.ny, pbe_basis.nz}, fft_mode);
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

std::vector<double> restrict_sternheimer_real_field_rectangular(const ModulePW::PW_Basis& fine_basis,
                                                                const ModulePW::PW_Basis& coarse_basis,
                                                                const std::vector<double>& fine_values)
{
    if (fine_basis.nrxx != fine_basis.nxyz || coarse_basis.nrxx != coarse_basis.nxyz)
    {
        throw std::invalid_argument("Rectangular Sternheimer restriction requires serial complete FFT grids.");
    }
    if (fine_values.size() != static_cast<std::size_t>(fine_basis.nxyz))
    {
        throw std::invalid_argument("Sternheimer fine real field size does not match its rectangular FFT grid.");
    }
    if (coarse_basis.nx > fine_basis.nx || coarse_basis.ny > fine_basis.ny || coarse_basis.nz > fine_basis.nz)
    {
        throw std::invalid_argument("Rectangular Sternheimer restriction does not support grid upsampling.");
    }
    if (fine_basis.nx == coarse_basis.nx && fine_basis.ny == coarse_basis.ny && fine_basis.nz == coarse_basis.nz)
    {
        return fine_values;
    }

    FFTWBuffer fine_fft(fine_basis.nx, fine_basis.ny, fine_basis.nz, FFTW_FORWARD);
    for (std::size_t index = 0; index != fine_values.size(); ++index)
    {
        fine_fft.data()[index][0] = fine_values[index];
        fine_fft.data()[index][1] = 0.0;
    }
    fine_fft.execute();

    FFTWBuffer coarse_fft(coarse_basis.nx, coarse_basis.ny, coarse_basis.nz, FFTW_BACKWARD);
    const double inverse_fine_size = 1.0 / static_cast<double>(fine_basis.nxyz);
    for (int ix = 0; ix != coarse_basis.nx; ++ix)
    {
        const auto fine_x = fine_frequencies_for_coarse_index(ix, coarse_basis.nx, fine_basis.nx);
        for (int iy = 0; iy != coarse_basis.ny; ++iy)
        {
            const auto fine_y = fine_frequencies_for_coarse_index(iy, coarse_basis.ny, fine_basis.ny);
            for (int iz = 0; iz != coarse_basis.nz; ++iz)
            {
                const auto fine_z = fine_frequencies_for_coarse_index(iz, coarse_basis.nz, fine_basis.nz);
                std::complex<double> coefficient(0.0, 0.0);
                for (int x_index = 0; x_index != fine_x.size; ++x_index)
                {
                    const int kx = fine_x.values[static_cast<std::size_t>(x_index)];
                    for (int y_index = 0; y_index != fine_y.size; ++y_index)
                    {
                        const int ky = fine_y.values[static_cast<std::size_t>(y_index)];
                        for (int z_index = 0; z_index != fine_z.size; ++z_index)
                        {
                            const int kz = fine_z.values[static_cast<std::size_t>(z_index)];
                            const std::size_t fine_index
                                = fft_index(frequency_index(kx, fine_basis.nx),
                                            frequency_index(ky, fine_basis.ny),
                                            frequency_index(kz, fine_basis.nz),
                                            fine_basis.ny,
                                            fine_basis.nz);
                            coefficient += std::complex<double>(fine_fft.data()[fine_index][0],
                                                                fine_fft.data()[fine_index][1]);
                        }
                    }
                }
                coefficient *= inverse_fine_size;
                const std::size_t coarse_index = fft_index(ix, iy, iz, coarse_basis.ny, coarse_basis.nz);
                coarse_fft.data()[coarse_index][0] = coefficient.real();
                coarse_fft.data()[coarse_index][1] = coefficient.imag();
            }
        }
    }
    coarse_fft.execute();

    std::vector<double> coarse_values(static_cast<std::size_t>(coarse_basis.nxyz), 0.0);
    for (std::size_t index = 0; index != coarse_values.size(); ++index)
    {
        coarse_values[index] = coarse_fft.data()[index][0];
    }
    return coarse_values;
}

} // namespace ModuleRI
