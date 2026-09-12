#include "source_lcao/module_ri/sternheimer_galerkin_audit.h"

#include "source_base/module_external/blas_connector.h"
#include "source_lcao/module_ri/sternheimer_galerkin_operator.h"
#include "source_lcao/module_ri/sternheimer_grid_transfer.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace
{

using Hamiltonian = ModuleRI::SternheimerFDHamiltonian;
using Grid = Hamiltonian::Grid;
using Complex = Hamiltonian::Complex;
using Vector = Hamiltonian::Vector;
using Galerkin = ModuleRI::SternheimerGalerkinOperator;
using Transfer = ModuleRI::SternheimerGridTransfer;
constexpr std::size_t state_block_size = 8;
constexpr std::size_t grid_block_size = 8192;

std::size_t checked_add(const std::size_t left, const std::size_t right)
{
    if (right > std::numeric_limits<std::size_t>::max() - left)
    {
        throw std::length_error("Sternheimer Galerkin audit workspace size overflows size_t.");
    }
    return left + right;
}

std::size_t checked_multiply(const std::size_t left, const std::size_t right)
{
    if (left != 0 && right > std::numeric_limits<std::size_t>::max() / left)
    {
        throw std::length_error("Sternheimer Galerkin audit workspace size overflows size_t.");
    }
    return left * right;
}

std::size_t complex_bytes(const std::size_t count)
{
    return checked_multiply(count, sizeof(Complex));
}

bool same_dimensions(const Grid& coarse, const Grid& fine)
{
    // Geometry and Bloch-sector equality are validated by Transfer first.
    return coarse.nx == fine.nx && coarse.ny == fine.ny && coarse.nz == fine.nz;
}

double volume_element(const Grid& grid)
{
    const auto a = ModuleRI::sternheimer_fd_grid_lattice_vectors(grid);
    const double volume = std::abs(a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
                                 + a[0][1] * (a[1][2] * a[2][0] - a[1][0] * a[2][2])
                                 + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]));
    return volume / grid.size();
}

bool finite(const Complex value)
{
    return std::isfinite(value.real()) && std::isfinite(value.imag());
}

double norm_squared(const Vector& state, const double weight)
{
    double sum = 0.0;
    for (const Complex value : state)
    {
        sum += std::norm(value);
    }
    const double norm = weight * sum;
    if (!std::isfinite(norm))
    {
        throw std::runtime_error("Sternheimer Galerkin audit encountered a nonfinite norm.");
    }
    return norm;
}

// Grid tiles keep BLAS leading dimensions small; state columns stay in their
// original coordinates. This workspace is private to one synchronous audit.
class BlockedProducts
{
  public:
    BlockedProducts(const std::size_t count, const std::size_t points, const std::size_t block)
        : points_(points),
          left_packed_(std::min(grid_block_size, points) * count),
          right_packed_(std::min(grid_block_size, points) * block)
    {
    }

    void accumulate(const std::vector<const Complex*>& left,
                    const std::vector<const Complex*>& right,
                    const std::size_t right_count,
                    const double weight,
                    Complex* output)
    {
        for (std::size_t begin = 0; begin < points_; begin += grid_block_size)
        {
            const std::size_t count = std::min(grid_block_size, points_ - begin);
            for (std::size_t column = 0; column != left.size(); ++column)
            {
                std::copy_n(left[column] + begin, count, left_packed_.data() + count * column);
            }
            for (std::size_t column = 0; column != right_count; ++column)
            {
                std::copy_n(right[column] + begin, count, right_packed_.data() + count * column);
            }
            BlasConnector::gemm_cm('C', 'N', static_cast<int>(left.size()), static_cast<int>(right_count),
                                   static_cast<int>(count), Complex(weight, 0.0),
                                   left_packed_.data(), static_cast<int>(count),
                                   right_packed_.data(), static_cast<int>(count), Complex(1.0, 0.0),
                                   output, static_cast<int>(left.size()));
        }
    }

  private:
    std::size_t points_;
    Vector left_packed_;
    Vector right_packed_;
};

} // namespace

namespace ModuleRI
{

std::size_t sternheimer_galerkin_audit_workspace_bytes(const Hamiltonian& fine,
                                                      const Grid& coarse,
                                                      const std::size_t dimension)
{
    if (dimension == 0 || dimension > static_cast<std::size_t>(std::numeric_limits<int>::max()))
    {
        throw std::invalid_argument("Sternheimer Galerkin audit needs a nonempty BLAS-sized state set.");
    }
    // Validate grids, Bloch sector and FD8 before size arithmetic or allocation.
    const std::size_t operator_bytes = Galerkin::workspace_bytes_required(fine, coarse);
    const std::size_t transfer_bytes = Transfer::workspace_bytes_required(coarse, fine.grid());
    const std::size_t nf = static_cast<std::size_t>(fine.grid().size());
    const std::size_t nc = static_cast<std::size_t>(coarse.size());
    const bool identity = same_dimensions(coarse, fine.grid());
    const std::size_t block = std::min(state_block_size, dimension);
    const std::size_t matrix_size = checked_multiply(dimension, dimension);
    std::size_t persistent = checked_add(complex_bytes(checked_multiply(5, matrix_size)),
                                         checked_multiply(checked_multiply(3, dimension), sizeof(double)));
    persistent = checked_add(persistent, checked_multiply(checked_multiply(2, dimension), sizeof(const Complex*)));
    if (!identity)
    {
        persistent = checked_add(persistent, complex_bytes(checked_multiply(checked_add(nc, nf), dimension)));
    }

    const std::size_t projection = identity ? 0 : checked_add(transfer_bytes, complex_bytes(checked_add(nc, nf)));
    const std::size_t coarse_pack = checked_multiply(std::min(grid_block_size, nc), checked_add(dimension, block));
    std::size_t coarse_vectors = checked_add(checked_multiply(nc, block), coarse_pack);
    coarse_vectors = checked_add(coarse_vectors, checked_multiply(identity ? 1 : 2, nc));
    const std::size_t pointers = checked_multiply(block, sizeof(const Complex*));
    const std::size_t coarse_stage = checked_add(checked_add(operator_bytes, complex_bytes(coarse_vectors)), pointers);

    // Reuse the frozen operator's bound for scalar H_f transient scratch only.
    const std::size_t fine_scratch = operator_bytes - transfer_bytes - complex_bytes(checked_multiply(2, nf));
    const std::size_t fine_pack = checked_multiply(std::min(grid_block_size, nf), checked_add(dimension, block));
    std::size_t fine_vectors = checked_add(checked_multiply(nf, block), fine_pack);
    fine_vectors = checked_add(fine_vectors, checked_multiply(identity ? 1 : 2, nf));
    fine_vectors = checked_add(fine_vectors, matrix_size);
    const std::size_t fine_stage = checked_add(checked_add(fine_scratch, complex_bytes(fine_vectors)), pointers);
    return checked_add(persistent, std::max(projection, std::max(coarse_stage, fine_stage)));
}

SternheimerGalerkinAuditMatrices audit_sternheimer_galerkin_matrices(
    std::shared_ptr<const Hamiltonian> fine,
    const std::vector<const Vector*>& fine_states,
    Grid coarse,
    const std::size_t max_workspace_bytes)
{
    if (!fine)
    {
        throw std::invalid_argument("Sternheimer Galerkin audit requires a fine Hamiltonian.");
    }
    const std::size_t n = fine_states.size();
    const std::size_t required = sternheimer_galerkin_audit_workspace_bytes(*fine, coarse, n);
    if (required > max_workspace_bytes)
    {
        throw std::length_error("Sternheimer Galerkin audit numerical payload exceeds its budget.");
    }
    const std::size_t nf = static_cast<std::size_t>(fine->grid().size());
    const std::size_t nc = static_cast<std::size_t>(coarse.size());
    for (const Vector* state : fine_states)
    {
        if (state == nullptr || state->size() != nf)
        {
            throw std::invalid_argument("Sternheimer Galerkin audit requires nonnull full fine-grid state views.");
        }
        for (const Complex value : *state)
        {
            if (!finite(value))
            {
                throw std::invalid_argument("Sternheimer Galerkin audit input states must be finite.");
            }
        }
    }

    SternheimerGalerkinAuditMatrices result;
    result.dimension = n;
    result.workspace_bytes = required;
    result.overlap.resize(n * n);
    result.kinetic.resize(n * n);
    result.local_potential.resize(n * n);
    result.nonlocal.resize(n * n);
    result.hamiltonian.resize(n * n);
    result.fine_norm_squared.resize(n);
    result.projected_norm_squared.resize(n);
    result.reconstruction_relative_error.resize(n);
    const bool identity = same_dimensions(coarse, fine->grid());
    const double dvf = volume_element(fine->grid());
    const double dvc = volume_element(coarse);
    const std::size_t block = std::min(state_block_size, n);
    Vector coarse_cache;
    Vector filtered_cache;
    std::vector<const Complex*> coarse_columns(n);
    std::vector<const Complex*> filtered_columns(n);
    for (std::size_t column = 0; column != n; ++column)
    {
        result.fine_norm_squared[column] = norm_squared(*fine_states[column], dvf);
    }
    if (identity)
    {
        for (std::size_t column = 0; column != n; ++column)
        {
            coarse_columns[column] = fine_states[column]->data();
            filtered_columns[column] = fine_states[column]->data();
        }
        result.projected_norm_squared = result.fine_norm_squared;
    }
    else
    {
        coarse_cache.resize(nc * n);
        filtered_cache.resize(nf * n);
        Transfer transfer(coarse, fine->grid());
        Vector projected(nc);
        Vector reconstructed(nf);
        for (std::size_t column = 0; column != n; ++column)
        {
            transfer.restrict_adjoint(*fine_states[column], projected);
            transfer.interpolate(projected, reconstructed);
            std::copy(projected.begin(), projected.end(), coarse_cache.data() + nc * column);
            std::copy(reconstructed.begin(), reconstructed.end(), filtered_cache.data() + nf * column);
            coarse_columns[column] = coarse_cache.data() + nc * column;
            filtered_columns[column] = filtered_cache.data() + nf * column;
            result.projected_norm_squared[column] = norm_squared(projected, dvc);
            double error_squared = 0.0;
            for (std::size_t point = 0; point != nf; ++point)
            {
                error_squared += std::norm(reconstructed[point] - (*fine_states[column])[point]);
            }
            if (!std::isfinite(error_squared))
            {
                throw std::runtime_error("Sternheimer Galerkin audit encountered a nonfinite reconstruction error.");
            }
            const double fine_norm = result.fine_norm_squared[column];
            result.reconstruction_relative_error[column] = fine_norm == 0.0 ? 0.0 : std::sqrt(dvf * error_squared / fine_norm);
        }
    }

    // Assemble C* H_c,term C through the actual Galerkin component API. Total H
    // is applied separately, not manufactured by summing reported components.
    {
        Galerkin op(fine, coarse, Galerkin::workspace_bytes_required(*fine, coarse));
        BlockedProducts products(n, nc, block);
        Vector input(identity ? 0 : nc);
        Vector output(nc);
        Vector images(nc * block);
        std::vector<const Complex*> right(block);
        using Apply = void (Galerkin::*)(const Vector&, Vector&);
        const std::array<Apply, 4> apply{{&Galerkin::apply_kinetic, &Galerkin::apply_local_potential,
                                        &Galerkin::apply_nonlocal, &Galerkin::apply}};
        const std::array<Vector*, 4> matrices{{&result.kinetic, &result.local_potential,
                                             &result.nonlocal, &result.hamiltonian}};
        for (std::size_t begin = 0; begin < n; begin += block)
        {
            const std::size_t count = std::min(block, n - begin);
            for (std::size_t column = 0; column != count; ++column)
            {
                right[column] = coarse_columns[begin + column];
            }
            products.accumulate(coarse_columns, right, count, dvc, result.overlap.data() + n * begin);
            for (std::size_t term = 0; term != apply.size(); ++term)
            {
                for (std::size_t column = 0; column != count; ++column)
                {
                    const Vector* state = fine_states[begin + column];
                    if (!identity)
                    {
                        std::copy_n(coarse_columns[begin + column], nc, input.data());
                        state = &input;
                    }
                    (op.*apply[term])(*state, output);
                    std::copy(output.begin(), output.end(), images.data() + nc * column);
                    right[column] = images.data() + nc * column;
                }
                products.accumulate(coarse_columns, right, count, dvc, matrices[term]->data() + n * begin);
            }
        }
    }

    // Independent FILTERED fine contraction, not a comparison to unfiltered KS.
    // Cached J C columns are reused: this stage performs no FFTs and just one
    // direct H_f application per state, with bounded blocks of operator images.
    {
        BlockedProducts products(n, nf, block);
        Vector input(identity ? 0 : nf);
        Vector output(nf);
        Vector images(nf * block);
        Vector reference(n * n);
        std::vector<const Complex*> right(block);
        for (std::size_t begin = 0; begin < n; begin += block)
        {
            const std::size_t count = std::min(block, n - begin);
            for (std::size_t column = 0; column != count; ++column)
            {
                const Vector* state = fine_states[begin + column];
                if (!identity)
                {
                    std::copy_n(filtered_columns[begin + column], nf, input.data());
                    state = &input;
                }
                fine->apply(*state, output);
                std::copy(output.begin(), output.end(), images.data() + nf * column);
                right[column] = images.data() + nf * column;
            }
            products.accumulate(filtered_columns, right, count, dvf, reference.data() + n * begin);
        }
        for (std::size_t i = 0; i != reference.size(); ++i)
        {
            const double error = std::abs(reference[i] - result.hamiltonian[i]);
            if (!std::isfinite(error))
            {
                throw std::runtime_error("Sternheimer Galerkin audit encountered a nonfinite filtered matrix check.");
            }
            result.filtered_matrix_identity_max_abs_error = std::max(result.filtered_matrix_identity_max_abs_error, error);
        }
    }
    for (const Vector* matrix : {&result.overlap, &result.kinetic, &result.local_potential,
                                 &result.nonlocal, &result.hamiltonian})
    {
        for (const Complex value : *matrix)
        {
            if (!finite(value))
            {
                throw std::runtime_error("Sternheimer Galerkin audit encountered a nonfinite matrix element.");
            }
        }
    }
    return result;
}

} // namespace ModuleRI
