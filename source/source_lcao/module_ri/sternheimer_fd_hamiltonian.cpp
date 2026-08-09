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

namespace ModuleRI
{

int SternheimerFDHamiltonian::Grid::size() const
{
    return nx * ny * nz;
}

SternheimerFDHamiltonian::SternheimerFDHamiltonian(Grid grid,
                                                   std::vector<double> local_potential,
                                                   const double kinetic_prefactor,
                                                   std::shared_ptr<const SternheimerFDNonlocalProjector>
                                                       nonlocal_projector,
                                                   const int finite_difference_order)
    : grid_(grid),
      local_potential_(std::move(local_potential)),
      kinetic_prefactor_(kinetic_prefactor),
      finite_difference_order_(finite_difference_order),
      nonlocal_projector_(std::move(nonlocal_projector))
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
    if (finite_difference_order_ != 2 && finite_difference_order_ != 4 && finite_difference_order_ != 6)
    {
        throw std::invalid_argument("SternheimerFDHamiltonian finite-difference order must be 2, 4, or 6.");
    }
    if (nonlocal_projector_ != nullptr && nonlocal_projector_->grid_size() != grid_.size())
    {
        throw std::invalid_argument("SternheimerFDHamiltonian nonlocal projector size does not match the grid.");
    }
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

SternheimerFDHamiltonian::ShiftedGridPoint SternheimerFDHamiltonian::shifted_grid_point(int ix,
                                                                                        int iy,
                                                                                        int iz) const
{
    if (grid_.periodic)
    {
        const std::array<int, 3> dimensions{grid_.nx, grid_.ny, grid_.nz};
        std::array<int, 3> coordinates{ix, iy, iz};
        std::array<int, 3> lattice_translation{};
        for (std::size_t direction = 0; direction != coordinates.size(); ++direction)
        {
            const int dimension = dimensions[direction];
            const int wrapped = (coordinates[direction] % dimension + dimension) % dimension;
            lattice_translation[direction] = (coordinates[direction] - wrapped) / dimension;
            coordinates[direction] = wrapped;
        }
        return {index(coordinates[0], coordinates[1], coordinates[2]),
                sternheimer_bloch_phase(grid_.kpoint, lattice_translation)};
    }

    if (ix < 0 || ix >= grid_.nx || iy < 0 || iy >= grid_.ny || iz < 0 || iz >= grid_.nz)
    {
        return {-1, Complex(1.0, 0.0)};
    }
    return {index(ix, iy, iz), Complex(1.0, 0.0)};
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

                    const bool fourth_order = finite_difference_order_ == 4;
                    const bool sixth_order = finite_difference_order_ == 6;
                    const double center_coefficient
                        = sixth_order ? -49.0 / 18.0 : (fourth_order ? -2.5 : -2.0);
                    Complex laplacian
                        = center_coefficient * (hx2_inv + hy2_inv + hz2_inv) * psi_center;
                    const int xp = shifted_index(ix + 1, iy, iz);
                    const int xm = shifted_index(ix - 1, iy, iz);
                    const int yp = shifted_index(ix, iy + 1, iz);
                    const int ym = shifted_index(ix, iy - 1, iz);
                    const int zp = shifted_index(ix, iy, iz + 1);
                    const int zm = shifted_index(ix, iy, iz - 1);

                    const double nearest_coefficient
                        = sixth_order ? 3.0 / 2.0 : (fourth_order ? 4.0 / 3.0 : 1.0);
                    if (xp >= 0)
                    {
                        laplacian += nearest_coefficient * hx2_inv * psi[xp];
                    }
                    if (xm >= 0)
                    {
                        laplacian += nearest_coefficient * hx2_inv * psi[xm];
                    }
                    if (yp >= 0)
                    {
                        laplacian += nearest_coefficient * hy2_inv * psi[yp];
                    }
                    if (ym >= 0)
                    {
                        laplacian += nearest_coefficient * hy2_inv * psi[ym];
                    }
                    if (zp >= 0)
                    {
                        laplacian += nearest_coefficient * hz2_inv * psi[zp];
                    }
                    if (zm >= 0)
                    {
                        laplacian += nearest_coefficient * hz2_inv * psi[zm];
                    }

                    if (fourth_order || sixth_order)
                    {
                        const double next_nearest_coefficient = sixth_order ? -3.0 / 20.0 : -1.0 / 12.0;
                        const int xpp = shifted_index(ix + 2, iy, iz);
                        const int xmm = shifted_index(ix - 2, iy, iz);
                        const int ypp = shifted_index(ix, iy + 2, iz);
                        const int ymm = shifted_index(ix, iy - 2, iz);
                        const int zpp = shifted_index(ix, iy, iz + 2);
                        const int zmm = shifted_index(ix, iy, iz - 2);
                        if (xpp >= 0)
                        {
                            laplacian += next_nearest_coefficient * hx2_inv * psi[xpp];
                        }
                        if (xmm >= 0)
                        {
                            laplacian += next_nearest_coefficient * hx2_inv * psi[xmm];
                        }
                        if (ypp >= 0)
                        {
                            laplacian += next_nearest_coefficient * hy2_inv * psi[ypp];
                        }
                        if (ymm >= 0)
                        {
                            laplacian += next_nearest_coefficient * hy2_inv * psi[ymm];
                        }
                        if (zpp >= 0)
                        {
                            laplacian += next_nearest_coefficient * hz2_inv * psi[zpp];
                        }
                        if (zmm >= 0)
                        {
                            laplacian += next_nearest_coefficient * hz2_inv * psi[zmm];
                        }
                    }

                    if (sixth_order)
                    {
                        constexpr double third_nearest_coefficient = 1.0 / 90.0;
                        const int xppp = shifted_index(ix + 3, iy, iz);
                        const int xmmm = shifted_index(ix - 3, iy, iz);
                        const int yppp = shifted_index(ix, iy + 3, iz);
                        const int ymmm = shifted_index(ix, iy - 3, iz);
                        const int zppp = shifted_index(ix, iy, iz + 3);
                        const int zmmm = shifted_index(ix, iy, iz - 3);
                        if (xppp >= 0)
                        {
                            laplacian += third_nearest_coefficient * hx2_inv * psi[xppp];
                        }
                        if (xmmm >= 0)
                        {
                            laplacian += third_nearest_coefficient * hx2_inv * psi[xmmm];
                        }
                        if (yppp >= 0)
                        {
                            laplacian += third_nearest_coefficient * hy2_inv * psi[yppp];
                        }
                        if (ymmm >= 0)
                        {
                            laplacian += third_nearest_coefficient * hy2_inv * psi[ymmm];
                        }
                        if (zppp >= 0)
                        {
                            laplacian += third_nearest_coefficient * hz2_inv * psi[zppp];
                        }
                        if (zmmm >= 0)
                        {
                            laplacian += third_nearest_coefficient * hz2_inv * psi[zmmm];
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

                    Complex laplacian = -2.0 * (laplacian_coefficients[0][0]
                                                 + laplacian_coefficients[1][1]
                                                 + laplacian_coefficients[2][2])
                                        * psi_center;
                    const auto add_shift = [&](const int dx, const int dy, const int dz, const double coefficient) {
                        const ShiftedGridPoint shifted = shifted_grid_point(ix + dx, iy + dy, iz + dz);
                        if (shifted.index >= 0)
                        {
                            laplacian += coefficient * shifted.phase * psi[shifted.index];
                        }
                    };
                    add_shift(1, 0, 0, laplacian_coefficients[0][0]);
                    add_shift(-1, 0, 0, laplacian_coefficients[0][0]);
                    add_shift(0, 1, 0, laplacian_coefficients[1][1]);
                    add_shift(0, -1, 0, laplacian_coefficients[1][1]);
                    add_shift(0, 0, 1, laplacian_coefficients[2][2]);
                    add_shift(0, 0, -1, laplacian_coefficients[2][2]);

                    const auto add_mixed = [&](const int dx1,
                                               const int dy1,
                                               const int dz1,
                                               const int dx2,
                                               const int dy2,
                                               const int dz2,
                                               const double coefficient) {
                        add_shift(dx1 + dx2, dy1 + dy2, dz1 + dz2, 0.5 * coefficient);
                        add_shift(dx1 - dx2, dy1 - dy2, dz1 - dz2, -0.5 * coefficient);
                        add_shift(-dx1 + dx2, -dy1 + dy2, -dz1 + dz2, -0.5 * coefficient);
                        add_shift(-dx1 - dx2, -dy1 - dy2, -dz1 - dz2, 0.5 * coefficient);
                    };
                    add_mixed(1, 0, 0, 0, 1, 0, laplacian_coefficients[0][1]);
                    add_mixed(1, 0, 0, 0, 0, 1, laplacian_coefficients[0][2]);
                    add_mixed(0, 1, 0, 0, 0, 1, laplacian_coefficients[1][2]);

                    output[center] = -kinetic_prefactor_ * laplacian;
                    if (include_local_potential)
                    {
                        output[center] += local_potential_[center] * psi_center;
                    }
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

void SternheimerFDHamiltonian::apply_local_potential(const Vector& psi, Vector& local_psi) const
{
    if (static_cast<int>(psi.size()) != grid_.size())
    {
        throw std::invalid_argument(
            "SternheimerFDHamiltonian::apply_local_potential input size does not match the grid.");
    }
    local_psi.resize(psi.size());
#pragma omp parallel for schedule(static)
    for (std::size_t ir = 0; ir != psi.size(); ++ir)
    {
        local_psi[ir] = local_potential_[ir] * psi[ir];
    }
}

void SternheimerFDHamiltonian::apply_nonlocal(const Vector& psi, Vector& nonlocal_psi) const
{
    if (static_cast<int>(psi.size()) != grid_.size())
    {
        throw std::invalid_argument(
            "SternheimerFDHamiltonian::apply_nonlocal input size does not match the grid.");
    }
    nonlocal_psi.assign(psi.size(), Complex(0.0, 0.0));
    if (nonlocal_projector_ != nullptr)
    {
        nonlocal_projector_->add_to(psi, nonlocal_psi);
    }
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
