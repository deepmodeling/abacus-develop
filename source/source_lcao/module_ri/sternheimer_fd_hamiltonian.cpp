#include "source_lcao/module_ri/sternheimer_fd_hamiltonian.h"

#include "source_base/module_external/lapack_connector.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <utility>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace
{

template <int Radius>
std::array<double, Radius + 1> second_derivative_coefficients();

template <>
std::array<double, 2> second_derivative_coefficients<1>()
{
    return {{-2.0, 1.0}};
}

template <>
std::array<double, 3> second_derivative_coefficients<2>()
{
    return {{-5.0 / 2.0, 4.0 / 3.0, -1.0 / 12.0}};
}

template <>
std::array<double, 4> second_derivative_coefficients<3>()
{
    return {{-49.0 / 18.0, 3.0 / 2.0, -3.0 / 20.0, 1.0 / 90.0}};
}

template <>
std::array<double, 5> second_derivative_coefficients<4>()
{
    return {{-205.0 / 72.0, 8.0 / 5.0, -1.0 / 5.0, 8.0 / 315.0, -1.0 / 560.0}};
}

} // namespace

namespace ModuleRI
{

int SternheimerFDHamiltonian::Grid::size() const
{
    return nx * ny * nz;
}

SternheimerFDHamiltonian::SternheimerFDHamiltonian(
    Grid grid,
    std::vector<double> local_potential,
    const double kinetic_prefactor,
    std::shared_ptr<const SternheimerFDNonlocalProjector> nonlocal_projector,
    const int finite_difference_order)
    : grid_(grid), local_potential_(std::move(local_potential)), kinetic_prefactor_(kinetic_prefactor),
      finite_difference_order_(finite_difference_order), nonlocal_projector_(std::move(nonlocal_projector))
{
    if (grid_.nx <= 0 || grid_.ny <= 0 || grid_.nz <= 0)
    {
        throw std::invalid_argument("SternheimerFDHamiltonian requires positive grid dimensions.");
    }
    if (grid_.hx <= 0.0 || grid_.hy <= 0.0 || grid_.hz <= 0.0)
    {
        throw std::invalid_argument("SternheimerFDHamiltonian requires positive grid spacings.");
    }
    if (static_cast<int>(local_potential_.size()) != grid_.size())
    {
        throw std::invalid_argument("SternheimerFDHamiltonian local potential size does not match the grid.");
    }
    if (kinetic_prefactor_ < 0.0)
    {
        throw std::invalid_argument("SternheimerFDHamiltonian requires a non-negative kinetic prefactor.");
    }
    if (finite_difference_order_ != 2 && finite_difference_order_ != 4 && finite_difference_order_ != 6
        && finite_difference_order_ != 8)
    {
        throw std::invalid_argument("SternheimerFDHamiltonian finite-difference order must be 2, 4, 6, or 8.");
    }
    if (nonlocal_projector_ != nullptr && nonlocal_projector_->grid_size() != grid_.size())
    {
        throw std::invalid_argument("SternheimerFDHamiltonian nonlocal projector size does not match the grid.");
    }

    const int radius = finite_difference_order_ / 2;
    const auto shifted_coordinate = [this](const int shifted, const int extent) {
        if (grid_.periodic)
        {
            return (shifted % extent + extent) % extent;
        }
        return shifted >= 0 && shifted < extent ? shifted : -1;
    };
    const auto initialize_coordinates
        = [radius, shifted_coordinate](const int extent,
                                       std::array<std::vector<int>, max_stencil_radius_>& positive_coordinates,
                                       std::array<std::vector<int>, max_stencil_radius_>& negative_coordinates) {
              for (int offset = 1; offset <= radius; ++offset)
              {
                  std::vector<int>& positive = positive_coordinates[offset - 1];
                  std::vector<int>& negative = negative_coordinates[offset - 1];
                  positive.resize(extent);
                  negative.resize(extent);
                  for (int coordinate = 0; coordinate != extent; ++coordinate)
                  {
                      const int positive_shifted = coordinate + offset;
                      const int negative_shifted = coordinate - offset;
                      positive[coordinate] = shifted_coordinate(positive_shifted, extent);
                      negative[coordinate] = shifted_coordinate(negative_shifted, extent);
                  }
              }
          };

    initialize_coordinates(grid_.nx, x_positive_coordinates_, x_negative_coordinates_);
    initialize_coordinates(grid_.ny, y_positive_coordinates_, y_negative_coordinates_);
    initialize_coordinates(grid_.nz, z_positive_coordinates_, z_negative_coordinates_);
}

const SternheimerFDHamiltonian::Grid& SternheimerFDHamiltonian::grid() const
{
    return grid_;
}

const std::vector<double>& SternheimerFDHamiltonian::local_potential() const
{
    return local_potential_;
}

double SternheimerFDHamiltonian::kinetic_prefactor() const
{
    return kinetic_prefactor_;
}

int SternheimerFDHamiltonian::finite_difference_order() const
{
    return finite_difference_order_;
}

const SternheimerFDNonlocalProjector* SternheimerFDHamiltonian::nonlocal_projector() const
{
    return nonlocal_projector_.get();
}

int SternheimerFDHamiltonian::index(const int ix, const int iy, const int iz) const
{
    return (ix * grid_.ny + iy) * grid_.nz + iz;
}

void SternheimerFDHamiltonian::apply(const Vector& psi, Vector& hpsi) const
{
    apply(psi, hpsi, nullptr);
}

void SternheimerFDHamiltonian::apply(const Vector& psi, Vector& hpsi, int* threads_used) const
{
    if (static_cast<int>(psi.size()) != grid_.size())
    {
        throw std::invalid_argument("SternheimerFDHamiltonian::apply_kinetic input size does not match the grid.");
    }

    hpsi.assign(psi.size(), Complex(0.0, 0.0));

    switch (finite_difference_order_)
    {
    case 2:
        apply_local<1>(psi, hpsi, threads_used);
        break;
    case 4:
        apply_local<2>(psi, hpsi, threads_used);
        break;
    case 6:
        apply_local<3>(psi, hpsi, threads_used);
        break;
    case 8:
        apply_local<4>(psi, hpsi, threads_used);
        break;
    default:
        throw std::logic_error("Unsupported Sternheimer finite-difference order.");
    }

    if (nonlocal_projector_ != nullptr)
    {
        nonlocal_projector_->add_to(psi, hpsi);
    }
}

template <int Radius>
void SternheimerFDHamiltonian::apply_local(const Vector& psi, Vector& hpsi, int* threads_used) const
{
    const std::array<double, Radius + 1> coefficients = second_derivative_coefficients<Radius>();
    const double hx2_inv = 1.0 / (grid_.hx * grid_.hx);
    const double hy2_inv = 1.0 / (grid_.hy * grid_.hy);
    const double hz2_inv = 1.0 / (grid_.hz * grid_.hz);

#ifdef _OPENMP
#pragma omp parallel
    {
#pragma omp single
        {
            if (threads_used != nullptr)
            {
                *threads_used = omp_get_num_threads();
            }
        }

#pragma omp for collapse(2) schedule(static)
#endif
        for (int ix = 0; ix != grid_.nx; ++ix)
        {
            for (int iy = 0; iy != grid_.ny; ++iy)
            {
                for (int iz = 0; iz != grid_.nz; ++iz)
                {
                    const int center = index(ix, iy, iz);
                    const Complex psi_center = psi[center];

                    Complex laplacian = coefficients[0] * (hx2_inv + hy2_inv + hz2_inv) * psi_center;
                    for (int offset = 0; offset != Radius; ++offset)
                    {
                        const double coefficient = coefficients[offset + 1];
                        const int xp = x_positive_coordinates_[offset][ix];
                        const int xm = x_negative_coordinates_[offset][ix];
                        const int yp = y_positive_coordinates_[offset][iy];
                        const int ym = y_negative_coordinates_[offset][iy];
                        const int zp = z_positive_coordinates_[offset][iz];
                        const int zm = z_negative_coordinates_[offset][iz];

                        if (xp != -1)
                        {
                            laplacian += coefficient * hx2_inv * psi[index(xp, iy, iz)];
                        }
                        if (xm != -1)
                        {
                            laplacian += coefficient * hx2_inv * psi[index(xm, iy, iz)];
                        }
                        if (yp != -1)
                        {
                            laplacian += coefficient * hy2_inv * psi[index(ix, yp, iz)];
                        }
                        if (ym != -1)
                        {
                            laplacian += coefficient * hy2_inv * psi[index(ix, ym, iz)];
                        }
                        if (zp != -1)
                        {
                            laplacian += coefficient * hz2_inv * psi[index(ix, iy, zp)];
                        }
                        if (zm != -1)
                        {
                            laplacian += coefficient * hz2_inv * psi[index(ix, iy, zm)];
                        }
                    }

                    hpsi[center] = -kinetic_prefactor_ * laplacian + local_potential_[center] * psi_center;
                }
            }
        }
#ifdef _OPENMP
    }
#else
    if (threads_used != nullptr)
    {
        *threads_used = 1;
    }
#endif
}

SternheimerFDHamiltonian::Matrix SternheimerFDHamiltonian::dense_matrix(const int max_size) const
{
    const int size = grid_.size();
    if (size > max_size)
    {
        throw std::invalid_argument("SternheimerFDHamiltonian::dense_matrix is only intended for small debug grids.");
    }

    Matrix matrix(size, Vector(size, Complex(0.0, 0.0)));
    Vector basis(size, Complex(0.0, 0.0));
    Vector h_basis;
    for (int col = 0; col != size; ++col)
    {
        std::fill(basis.begin(), basis.end(), Complex(0.0, 0.0));
        basis[col] = Complex(1.0, 0.0);
        apply(basis, h_basis);
        for (int row = 0; row != size; ++row)
        {
            matrix[row][col] = h_basis[row];
        }
    }
    return matrix;
}

SternheimerFDHamiltonian::Eigenpairs SternheimerFDHamiltonian::diagonalize_dense(const int max_size) const
{
    const int size = grid_.size();
    if (size > max_size)
    {
        throw std::invalid_argument(
            "SternheimerFDHamiltonian::diagonalize_dense is only intended for small debug grids.");
    }

    const Matrix matrix = dense_matrix(max_size);
    std::vector<Complex> lapack_matrix(size * size, Complex(0.0, 0.0));
    for (int col = 0; col != size; ++col)
    {
        for (int row = 0; row != size; ++row)
        {
            lapack_matrix[row + size * col] = matrix[row][col];
        }
    }

    char jobz = 'V';
    char uplo = 'U';
    const int n = size;
    const int lda = size;
    const int minus_one = -1;
    int info = 0;
    std::vector<double> eigenvalues(size, 0.0);
    Complex work_query(0.0, 0.0);
    std::vector<double> rwork(std::max(1, 3 * n - 2), 0.0);

    zheev_(&jobz,
           &uplo,
           &n,
           lapack_matrix.data(),
           &lda,
           eigenvalues.data(),
           &work_query,
           &minus_one,
           rwork.data(),
           &info);
    if (info != 0)
    {
        throw std::runtime_error("SternheimerFDHamiltonian::diagonalize_dense LAPACK workspace query failed.");
    }

    const int lwork = std::max(1, static_cast<int>(std::ceil(work_query.real())));
    std::vector<Complex> work(lwork, Complex(0.0, 0.0));
    zheev_(&jobz, &uplo, &n, lapack_matrix.data(), &lda, eigenvalues.data(), work.data(), &lwork, rwork.data(), &info);
    if (info != 0)
    {
        throw std::runtime_error("SternheimerFDHamiltonian::diagonalize_dense LAPACK diagonalization failed.");
    }

    Eigenpairs eigenpairs;
    eigenpairs.eigenvalues = std::move(eigenvalues);
    eigenpairs.eigenvectors.assign(size, Vector(size, Complex(0.0, 0.0)));
    for (int band = 0; band != size; ++band)
    {
        for (int row = 0; row != size; ++row)
        {
            eigenpairs.eigenvectors[band][row] = lapack_matrix[row + size * band];
        }
    }
    return eigenpairs;
}

} // namespace ModuleRI
