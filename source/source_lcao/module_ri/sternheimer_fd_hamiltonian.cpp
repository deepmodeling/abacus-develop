#include "source_lcao/module_ri/sternheimer_fd_hamiltonian.h"

#include "source_base/module_external/lapack_connector.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <utility>

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
                                                       nonlocal_projector)
    : grid_(grid),
      local_potential_(std::move(local_potential)),
      kinetic_prefactor_(kinetic_prefactor),
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
    bool has_nonzero_twist = false;
    for (const double coordinate: grid_.kpoint)
    {
        if (!std::isfinite(coordinate))
        {
            throw std::invalid_argument("SternheimerFDHamiltonian requires finite reduced k-point coordinates.");
        }
        has_nonzero_twist = has_nonzero_twist || coordinate != 0.0;
    }
    if (!grid_.periodic && has_nonzero_twist)
    {
        throw std::invalid_argument("SternheimerFDHamiltonian nonperiodic grids cannot use a Bloch twist.");
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

const SternheimerReducedKPoint& SternheimerFDHamiltonian::kpoint() const
{
    return grid_.kpoint;
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
    if (static_cast<int>(psi.size()) != grid_.size())
    {
        throw std::invalid_argument("SternheimerFDHamiltonian::apply input size does not match the grid.");
    }

    hpsi.assign(psi.size(), Complex(0.0, 0.0));
    const double hx2_inv = 1.0 / (grid_.hx * grid_.hx);
    const double hy2_inv = 1.0 / (grid_.hy * grid_.hy);
    const double hz2_inv = 1.0 / (grid_.hz * grid_.hz);

    for (int iz = 0; iz != grid_.nz; ++iz)
    {
        for (int iy = 0; iy != grid_.ny; ++iy)
        {
            for (int ix = 0; ix != grid_.nx; ++ix)
            {
                const int center = index(ix, iy, iz);
                const Complex psi_center = psi[center];

                Complex laplacian = -2.0 * (hx2_inv + hy2_inv + hz2_inv) * psi_center;
                const ShiftedGridPoint xp = shifted_grid_point(ix + 1, iy, iz);
                const ShiftedGridPoint xm = shifted_grid_point(ix - 1, iy, iz);
                const ShiftedGridPoint yp = shifted_grid_point(ix, iy + 1, iz);
                const ShiftedGridPoint ym = shifted_grid_point(ix, iy - 1, iz);
                const ShiftedGridPoint zp = shifted_grid_point(ix, iy, iz + 1);
                const ShiftedGridPoint zm = shifted_grid_point(ix, iy, iz - 1);

                if (xp.index >= 0)
                {
                    laplacian += hx2_inv * xp.phase * psi[xp.index];
                }
                if (xm.index >= 0)
                {
                    laplacian += hx2_inv * xm.phase * psi[xm.index];
                }
                if (yp.index >= 0)
                {
                    laplacian += hy2_inv * yp.phase * psi[yp.index];
                }
                if (ym.index >= 0)
                {
                    laplacian += hy2_inv * ym.phase * psi[ym.index];
                }
                if (zp.index >= 0)
                {
                    laplacian += hz2_inv * zp.phase * psi[zp.index];
                }
                if (zm.index >= 0)
                {
                    laplacian += hz2_inv * zm.phase * psi[zm.index];
                }

                hpsi[center] = -kinetic_prefactor_ * laplacian + local_potential_[center] * psi_center;
            }
        }
    }

    if (nonlocal_projector_ != nullptr)
    {
        nonlocal_projector_->add_to(psi, hpsi);
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
