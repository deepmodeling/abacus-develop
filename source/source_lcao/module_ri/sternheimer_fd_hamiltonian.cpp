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
    initialize_operator_cache();
}

void SternheimerFDHamiltonian::initialize_operator_cache()
{
    const SternheimerFDStencilWeights weights = sternheimer_fd_stencil_weights(finite_difference_order_);
    fd_radius_ = weights.radius;
    fd_second_weights_ = weights.second;
    fd_first_weights_ = weights.first;

    const SternheimerFDLatticeVectors dual = sternheimer_fd_grid_dual_vectors(grid_);
    const std::array<double, 3> dimensions{static_cast<double>(grid_.nx),
                                           static_cast<double>(grid_.ny),
                                           static_cast<double>(grid_.nz)};
    for (int left = 0; left != 3; ++left)
    {
        for (int right = 0; right != 3; ++right)
        {
            for (int component = 0; component != 3; ++component)
            {
                laplacian_coefficients_[left][right]
                    += dual[left][component] * dual[right][component] * dimensions[left] * dimensions[right];
            }
        }
    }
    laplacian_center_coefficient_
        = fd_second_weights_[0]
          * (laplacian_coefficients_[0][0] + laplacian_coefficients_[1][1]
             + laplacian_coefficients_[2][2]);

    for (int left = 0; left != 3; ++left)
    {
        for (int right = left + 1; right != 3; ++right)
        {
            if (laplacian_coefficients_[left][right] == 0.0)
            {
                continue;
            }
            active_mixed_derivative_pairs_[static_cast<std::size_t>(active_mixed_derivative_pair_count_)]
                = {left, right};
            ++active_mixed_derivative_pair_count_;
            if (grid_.periodic && fd_radius_ > 1)
            {
                cached_first_derivative_directions_[right] = true;
            }
        }
    }

    for (const double coordinate: grid_.kpoint)
    {
        has_zero_bloch_twist_ = has_zero_bloch_twist_ && coordinate == 0.0;
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

const SternheimerReducedKPoint& SternheimerFDHamiltonian::kpoint() const
{
    return grid_.kpoint;
}

const SternheimerFDNonlocalProjector* SternheimerFDHamiltonian::nonlocal_projector() const
{
    return nonlocal_projector_.get();
}

int SternheimerFDHamiltonian::active_mixed_derivative_pair_count() const
{
    return active_mixed_derivative_pair_count_;
}

const std::array<bool, 3>& SternheimerFDHamiltonian::cached_first_derivative_directions() const
{
    return cached_first_derivative_directions_;
}

int SternheimerFDHamiltonian::index(const int ix, const int iy, const int iz) const
{
    return (ix * grid_.ny + iy) * grid_.nz + iz;
}

SternheimerFDHamiltonian::ShiftedGridPoint SternheimerFDHamiltonian::shifted_grid_point(int ix, int iy, int iz) const
{
    if (grid_.periodic)
    {
        if (ix >= 0 && ix < grid_.nx && iy >= 0 && iy < grid_.ny && iz >= 0 && iz < grid_.nz)
        {
            return {index(ix, iy, iz), Complex(1.0, 0.0)};
        }
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
        const Complex phase = has_zero_bloch_twist_ ? Complex(1.0, 0.0)
                                                    : sternheimer_bloch_phase(grid_.kpoint, lattice_translation);
        return {index(coordinates[0], coordinates[1], coordinates[2]), phase};
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
        throw std::invalid_argument("SternheimerFDHamiltonian::apply input size does not match the grid.");
    }

    apply_grid_terms(psi, hpsi, true, threads_used);
    if (nonlocal_projector_ != nullptr)
    {
        nonlocal_projector_->add_to(psi, hpsi);
    }
}

void SternheimerFDHamiltonian::apply_batch(const Matrix& psi, Matrix& hpsi) const
{
    apply_batch(psi, hpsi, nullptr);
}

void SternheimerFDHamiltonian::apply_batch(const Matrix& psi, Matrix& hpsi, int* threads_used) const
{
    apply_grid_terms_batch(psi, hpsi, true, threads_used);
    if (nonlocal_projector_ != nullptr && !psi.empty())
    {
        nonlocal_projector_->add_to_batch(psi, hpsi);
    }
}

void SternheimerFDHamiltonian::apply_kinetic(const Vector& psi, Vector& kinetic_psi) const
{
    apply_kinetic(psi, kinetic_psi, nullptr);
}

void SternheimerFDHamiltonian::apply_kinetic(const Vector& psi, Vector& kinetic_psi, int* threads_used) const
{
    apply_grid_terms(psi, kinetic_psi, false, threads_used);
}

void SternheimerFDHamiltonian::apply_grid_terms(const Vector& psi,
                                                Vector& output,
                                                const bool include_local_potential,
                                                int* threads_used) const
{
    if (static_cast<int>(psi.size()) != grid_.size())
    {
        throw std::invalid_argument("SternheimerFDHamiltonian::apply_kinetic input size does not match the grid.");
    }

    output.assign(psi.size(), Complex(0.0, 0.0));
    const bool use_separable_mixed_derivatives = grid_.periodic && fd_radius_ > 1;
    std::array<Vector, 3> first_derivatives;
    if (use_separable_mixed_derivatives)
    {
        for (int direction = 0; direction != 3; ++direction)
        {
            if (!cached_first_derivative_directions_[direction])
            {
                continue;
            }
            Vector& derivative = first_derivatives[direction];
            derivative.assign(psi.size(), Complex(0.0, 0.0));
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static)
#endif
            for (int ix = 0; ix != grid_.nx; ++ix)
            {
                for (int iy = 0; iy != grid_.ny; ++iy)
                {
                    for (int iz = 0; iz != grid_.nz; ++iz)
                    {
                        Complex value(0.0, 0.0);
                        for (int offset = 1; offset <= fd_radius_; ++offset)
                        {
                            std::array<int, 3> shift{};
                            shift[direction] = offset;
                            const ShiftedGridPoint positive
                                = shifted_grid_point(ix + shift[0], iy + shift[1], iz + shift[2]);
                            const ShiftedGridPoint negative
                                = shifted_grid_point(ix - shift[0], iy - shift[1], iz - shift[2]);
                            const double coefficient = fd_first_weights_[static_cast<std::size_t>(offset - 1)];
                            if (positive.index >= 0)
                            {
                                value += coefficient * positive.phase * psi[positive.index];
                            }
                            if (negative.index >= 0)
                            {
                                value -= coefficient * negative.phase * psi[negative.index];
                            }
                        }
                        derivative[static_cast<std::size_t>(index(ix, iy, iz))] = value;
                    }
                }
            }
        }
    }

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

                    Complex laplacian = laplacian_center_coefficient_ * psi_center;
                    const auto add_shift = [&](const int dx, const int dy, const int dz, const double coefficient) {
                        const ShiftedGridPoint shifted = shifted_grid_point(ix + dx, iy + dy, iz + dz);
                        if (shifted.index >= 0)
                        {
                            laplacian += coefficient * shifted.phase * psi[shifted.index];
                        }
                    };
                    for (int direction = 0; direction != 3; ++direction)
                    {
                        for (int offset = 1; offset <= fd_radius_; ++offset)
                        {
                            std::array<int, 3> shift{};
                            shift[direction] = offset;
                            const double coefficient = fd_second_weights_[static_cast<std::size_t>(offset)]
                                                       * laplacian_coefficients_[direction][direction];
                            add_shift(shift[0], shift[1], shift[2], coefficient);
                            add_shift(-shift[0], -shift[1], -shift[2], coefficient);
                        }
                    }

                    for (int pair_index = 0; pair_index != active_mixed_derivative_pair_count_; ++pair_index)
                    {
                        const int left = active_mixed_derivative_pairs_[static_cast<std::size_t>(pair_index)][0];
                        const int right = active_mixed_derivative_pairs_[static_cast<std::size_t>(pair_index)][1];
                        if (use_separable_mixed_derivatives)
                        {
                            Complex mixed_derivative(0.0, 0.0);
                            for (int offset = 1; offset <= fd_radius_; ++offset)
                            {
                                std::array<int, 3> shift{};
                                shift[left] = offset;
                                const ShiftedGridPoint positive
                                    = shifted_grid_point(ix + shift[0], iy + shift[1], iz + shift[2]);
                                const ShiftedGridPoint negative
                                    = shifted_grid_point(ix - shift[0], iy - shift[1], iz - shift[2]);
                                const double coefficient = fd_first_weights_[static_cast<std::size_t>(offset - 1)];
                                if (positive.index >= 0)
                                {
                                    mixed_derivative
                                        += coefficient * positive.phase * first_derivatives[right][positive.index];
                                }
                                if (negative.index >= 0)
                                {
                                    mixed_derivative
                                        -= coefficient * negative.phase * first_derivatives[right][negative.index];
                                }
                            }
                            laplacian += 2.0 * laplacian_coefficients_[left][right] * mixed_derivative;
                        }
                        else
                        {
                            for (int left_offset = 1; left_offset <= fd_radius_; ++left_offset)
                            {
                                for (int right_offset = 1; right_offset <= fd_radius_; ++right_offset)
                                {
                                    const double coefficient
                                        = 2.0 * laplacian_coefficients_[left][right]
                                          * fd_first_weights_[static_cast<std::size_t>(left_offset - 1)]
                                          * fd_first_weights_[static_cast<std::size_t>(right_offset - 1)];
                                    for (const int left_sign: {-1, 1})
                                    {
                                        for (const int right_sign: {-1, 1})
                                        {
                                            std::array<int, 3> shift{};
                                            shift[left] = left_sign * left_offset;
                                            shift[right] = right_sign * right_offset;
                                            add_shift(shift[0],
                                                      shift[1],
                                                      shift[2],
                                                      coefficient * left_sign * right_sign);
                                        }
                                    }
                                }
                            }
                        }
                    }

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

void SternheimerFDHamiltonian::apply_grid_terms_batch(const Matrix& psi,
                                                      Matrix& output,
                                                      const bool include_local_potential,
                                                      int* threads_used) const
{
    if (psi.empty())
    {
        output.clear();
        if (threads_used != nullptr)
        {
            *threads_used = 1;
        }
        return;
    }
    for (const Vector& column: psi)
    {
        if (static_cast<int>(column.size()) != grid_.size())
        {
            throw std::invalid_argument(
                "SternheimerFDHamiltonian::apply_grid_terms_batch input size does not match the grid.");
        }
    }

    const std::size_t batch_size = psi.size();
    const std::size_t grid_size = static_cast<std::size_t>(grid_.size());
    output.assign(batch_size, Vector(grid_size, Complex(0.0, 0.0)));
    const bool use_separable_mixed_derivatives = grid_.periodic && fd_radius_ > 1;
    std::array<Matrix, 3> first_derivatives;
    if (use_separable_mixed_derivatives)
    {
        for (int direction = 0; direction != 3; ++direction)
        {
            if (!cached_first_derivative_directions_[direction])
            {
                continue;
            }
            Matrix& derivatives = first_derivatives[direction];
            derivatives.assign(batch_size, Vector(grid_size, Complex(0.0, 0.0)));
#ifdef _OPENMP
#pragma omp parallel
#endif
            {
                std::vector<Complex> values(batch_size, Complex(0.0, 0.0));
#ifdef _OPENMP
#pragma omp for collapse(2) schedule(static)
#endif
                for (int ix = 0; ix != grid_.nx; ++ix)
                {
                    for (int iy = 0; iy != grid_.ny; ++iy)
                    {
                        for (int iz = 0; iz != grid_.nz; ++iz)
                        {
                            std::fill(values.begin(), values.end(), Complex(0.0, 0.0));
                            for (int offset = 1; offset <= fd_radius_; ++offset)
                            {
                                std::array<int, 3> shift{};
                                shift[direction] = offset;
                                const ShiftedGridPoint positive
                                    = shifted_grid_point(ix + shift[0], iy + shift[1], iz + shift[2]);
                                const ShiftedGridPoint negative
                                    = shifted_grid_point(ix - shift[0], iy - shift[1], iz - shift[2]);
                                const double coefficient
                                    = fd_first_weights_[static_cast<std::size_t>(offset - 1)];
                                for (std::size_t column = 0; column != batch_size; ++column)
                                {
                                    if (positive.index >= 0)
                                    {
                                        values[column]
                                            += coefficient * positive.phase * psi[column][positive.index];
                                    }
                                    if (negative.index >= 0)
                                    {
                                        values[column]
                                            -= coefficient * negative.phase * psi[column][negative.index];
                                    }
                                }
                            }
                            const std::size_t center = static_cast<std::size_t>(index(ix, iy, iz));
                            for (std::size_t column = 0; column != batch_size; ++column)
                            {
                                derivatives[column][center] = values[column];
                            }
                        }
                    }
                }
            }
        }
    }

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        std::vector<Complex> laplacians(batch_size, Complex(0.0, 0.0));
        std::vector<Complex> mixed_derivatives(batch_size, Complex(0.0, 0.0));
#ifdef _OPENMP
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
                    for (std::size_t column = 0; column != batch_size; ++column)
                    {
                        laplacians[column] = laplacian_center_coefficient_ * psi[column][center];
                    }

                    const auto add_shift
                        = [&](const int dx, const int dy, const int dz, const double coefficient) {
                              const ShiftedGridPoint shifted = shifted_grid_point(ix + dx, iy + dy, iz + dz);
                              if (shifted.index >= 0)
                              {
                                  for (std::size_t column = 0; column != batch_size; ++column)
                                  {
                                      laplacians[column]
                                          += coefficient * shifted.phase * psi[column][shifted.index];
                                  }
                              }
                          };
                    for (int direction = 0; direction != 3; ++direction)
                    {
                        for (int offset = 1; offset <= fd_radius_; ++offset)
                        {
                            std::array<int, 3> shift{};
                            shift[direction] = offset;
                            const double coefficient = fd_second_weights_[static_cast<std::size_t>(offset)]
                                                       * laplacian_coefficients_[direction][direction];
                            add_shift(shift[0], shift[1], shift[2], coefficient);
                            add_shift(-shift[0], -shift[1], -shift[2], coefficient);
                        }
                    }

                    for (int pair_index = 0; pair_index != active_mixed_derivative_pair_count_; ++pair_index)
                    {
                        const int left
                            = active_mixed_derivative_pairs_[static_cast<std::size_t>(pair_index)][0];
                        const int right
                            = active_mixed_derivative_pairs_[static_cast<std::size_t>(pair_index)][1];
                        if (use_separable_mixed_derivatives)
                        {
                            std::fill(mixed_derivatives.begin(), mixed_derivatives.end(), Complex(0.0, 0.0));
                            for (int offset = 1; offset <= fd_radius_; ++offset)
                            {
                                std::array<int, 3> shift{};
                                shift[left] = offset;
                                const ShiftedGridPoint positive
                                    = shifted_grid_point(ix + shift[0], iy + shift[1], iz + shift[2]);
                                const ShiftedGridPoint negative
                                    = shifted_grid_point(ix - shift[0], iy - shift[1], iz - shift[2]);
                                const double coefficient
                                    = fd_first_weights_[static_cast<std::size_t>(offset - 1)];
                                for (std::size_t column = 0; column != batch_size; ++column)
                                {
                                    if (positive.index >= 0)
                                    {
                                        mixed_derivatives[column]
                                            += coefficient * positive.phase
                                               * first_derivatives[right][column][positive.index];
                                    }
                                    if (negative.index >= 0)
                                    {
                                        mixed_derivatives[column]
                                            -= coefficient * negative.phase
                                               * first_derivatives[right][column][negative.index];
                                    }
                                }
                            }
                            for (std::size_t column = 0; column != batch_size; ++column)
                            {
                                laplacians[column] += 2.0 * laplacian_coefficients_[left][right]
                                                      * mixed_derivatives[column];
                            }
                        }
                        else
                        {
                            for (int left_offset = 1; left_offset <= fd_radius_; ++left_offset)
                            {
                                for (int right_offset = 1; right_offset <= fd_radius_; ++right_offset)
                                {
                                    const double coefficient
                                        = 2.0 * laplacian_coefficients_[left][right]
                                          * fd_first_weights_[static_cast<std::size_t>(left_offset - 1)]
                                          * fd_first_weights_[static_cast<std::size_t>(right_offset - 1)];
                                    for (const int left_sign: {-1, 1})
                                    {
                                        for (const int right_sign: {-1, 1})
                                        {
                                            std::array<int, 3> shift{};
                                            shift[left] = left_sign * left_offset;
                                            shift[right] = right_sign * right_offset;
                                            add_shift(shift[0],
                                                      shift[1],
                                                      shift[2],
                                                      coefficient * left_sign * right_sign);
                                        }
                                    }
                                }
                            }
                        }
                    }

                    for (std::size_t column = 0; column != batch_size; ++column)
                    {
                        output[column][center] = -kinetic_prefactor_ * laplacians[column];
                        if (include_local_potential)
                        {
                            output[column][center] += local_potential_[center] * psi[column][center];
                        }
                    }
                }
            }
        }
    }
#ifndef _OPENMP
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
        throw std::invalid_argument("SternheimerFDHamiltonian::apply_nonlocal input size does not match the grid.");
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
