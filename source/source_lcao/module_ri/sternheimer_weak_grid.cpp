#include "source_lcao/module_ri/sternheimer_weak_grid.h"

#include "source_base/module_external/blas_connector.h"
#include "source_base/module_external/lapack_connector.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace
{
using H = ModuleRI::SternheimerFDHamiltonian;
using Operator = ModuleRI::SternheimerWeakGridOperator;
using Complex = H::Complex;

std::size_t Add(std::size_t a, std::size_t b)
{
    if (b > std::numeric_limits<std::size_t>::max() - a)
        throw std::length_error("Sternheimer weak grid size overflow.");
    return a + b;
}

std::size_t Mul(std::size_t a, std::size_t b)
{
    if (a && b > std::numeric_limits<std::size_t>::max() / a)
        throw std::length_error("Sternheimer weak grid size overflow.");
    return a * b;
}

const H& RequireFine(const std::shared_ptr<const H>& fine)
{
    if (!fine) throw std::invalid_argument("Sternheimer weak grid requires a fine Hamiltonian.");
    return *fine;
}

void Validate(const H::Vector& input, std::size_t size)
{
    if (input.size() != size) throw std::invalid_argument("Sternheimer weak grid vector dimension mismatch.");
    for (const auto& value : input)
        if (!std::isfinite(value.real()) || !std::isfinite(value.imag()))
            throw std::invalid_argument("Sternheimer weak grid requires finite values.");
}

void CheckOutput(const H::Vector& output)
{
    for (const auto& value : output)
        if (!std::isfinite(value.real()) || !std::isfinite(value.imag()))
            throw std::overflow_error("Sternheimer weak grid arithmetic overflow.");
}

double VolumeElement(const H::Grid& grid)
{
    const auto cell = ModuleRI::sternheimer_fd_grid_lattice_vectors(grid);
    const auto& a = cell[0]; const auto& b = cell[1]; const auto& c = cell[2];
    return std::abs(a[0] * (b[1] * c[2] - b[2] * c[1])
                    - a[1] * (b[0] * c[2] - b[2] * c[0])
                    + a[2] * (b[0] * c[1] - b[1] * c[0])) / grid.size();
}

std::size_t CheckWorkspace(const H& fine, const H::Grid& coarse, std::size_t budget)
{
    const auto bytes = Operator::workspace_bytes_required(fine, coarse);
    if (bytes > budget) throw std::length_error("Sternheimer weak grid workspace budget exceeded.");
    const double fine_dv = VolumeElement(fine.grid());
    const double coarse_dv = VolumeElement(coarse);
    if (!std::isfinite(fine_dv) || fine_dv <= 0.0 || !std::isfinite(coarse_dv) || coarse_dv <= 0.0)
        throw std::invalid_argument("Sternheimer weak grid requires finite positive volume elements.");
    if (!std::isfinite(fine.kinetic_prefactor()))
        throw std::invalid_argument("Sternheimer weak grid requires a finite kinetic prefactor.");
    for (double v : fine.local_potential())
        if (!std::isfinite(v)) throw std::invalid_argument("Sternheimer weak grid requires a finite potential.");
    const auto* projector = fine.nonlocal_projector();
    if (projector)
    {
        const double projector_dv = projector->volume_element();
        if (!std::isfinite(projector_dv) || projector_dv <= 0.0
            || std::abs(projector_dv - fine_dv) > 1e-12 * std::max(projector_dv, fine_dv))
            throw std::invalid_argument("Sternheimer weak grid nonlocal projector has a different or invalid volume element.");
        for (const auto& block : projector->blocks())
        {
            for (const auto& beta : block.projectors) Validate(beta, fine.grid().size());
            double scale = 1.0;
            for (const auto& row : block.d_matrix)
            {
                Validate(row, block.projectors.size());
                for (const auto value : row) scale = std::max(scale, std::abs(value));
            }
            if (!std::isfinite(scale))
                throw std::invalid_argument("Sternheimer weak grid nonlocal D magnitude overflow.");
            for (std::size_t i = 0; i < block.d_matrix.size(); ++i)
                for (std::size_t j = 0; j < block.d_matrix.size(); ++j)
                    if (std::abs(block.d_matrix[i][j] / scale - std::conj(block.d_matrix[j][i]) / scale) > 1e-12)
                        throw std::invalid_argument("Sternheimer weak grid requires a Hermitian nonlocal D matrix.");
        }
    }
    return bytes;
}
} // namespace

namespace ModuleRI
{
double orthonormalize_sternheimer_weak_states_in_place(
    std::vector<SternheimerDeltaGridFunction>& states, const double volume_element,
    const double min_metric, const std::size_t max_workspace_bytes)
{
    if (states.empty() || states.size() > static_cast<std::size_t>(std::numeric_limits<int>::max() / 4)
        || !std::isfinite(volume_element) || volume_element <= 0.0
        || !std::isfinite(min_metric) || min_metric <= 0.0 || min_metric >= 1.0)
        throw std::invalid_argument("Sternheimer weak normalization has invalid dimensions, volume or metric threshold.");
    const std::size_t grid_size = states.front().values.size();
    if (grid_size == 0)
        throw std::invalid_argument("Sternheimer weak normalization requires nonempty fields.");
    for (const auto& state : states)
    {
        Validate(state.values, grid_size);
        for (const auto& gradient : state.gradients) Validate(gradient, grid_size);
    }
    if (states.size() > grid_size)
        throw std::domain_error("Sternheimer weak normalization metric is singular: more states than grid points.");

    const int n = static_cast<int>(states.size());
    const std::size_t square = Mul(states.size(), states.size());
    const int lwork = std::max(1, 2 * n - 1);
    const int lrwork = std::max(1, 3 * n - 2);
    const std::size_t fixed_bytes = Add(Mul(Add(Mul(2, square), lwork), sizeof(Complex)),
                                        Mul(Add(states.size(), lrwork), sizeof(double)));
    const std::size_t bytes_per_row = Mul(states.size(), sizeof(Complex));
    if (Add(fixed_bytes, bytes_per_row) > max_workspace_bytes)
        throw std::length_error("Sternheimer weak normalization workspace budget exceeded.");
    const std::size_t tile_rows = std::min({grid_size, std::size_t(4096),
                                          (max_workspace_bytes - fixed_bytes) / bytes_per_row});
    H::Vector gram(square, Complex(0.0)), eigen_matrix(square), work(lwork);
    H::Vector tile(Mul(tile_rows, states.size()));
    std::vector<double> eigenvalues(n), rwork(lrwork);

    const auto assemble_gram = [&](H::Vector& output) {
        std::fill(output.begin(), output.end(), Complex(0.0));
        for (std::size_t begin = 0; begin < grid_size; begin += tile_rows)
        {
            const int rows = static_cast<int>(std::min(tile_rows, grid_size - begin));
            for (std::size_t j = 0; j < states.size(); ++j)
                std::copy_n(states[j].values.data() + begin, rows,
                            tile.data() + static_cast<std::size_t>(rows) * j);
            BlasConnector::gemm_cm('C', 'N', n, n, rows, Complex(volume_element), tile.data(), rows,
                                   tile.data(), rows, Complex(1.0), output.data(), n);
        }
        CheckOutput(output);
    };
    assemble_gram(gram);
    eigen_matrix = gram;
    const char jobz = 'N', upper = 'U';
    int info = 0;
    zheev_(&jobz, &upper, &n, eigen_matrix.data(), &n, eigenvalues.data(), work.data(), &lwork,
           rwork.data(), &info);
    if (info != 0)
        throw std::runtime_error("Sternheimer weak normalization metric eigensolver failed.");
    for (const double value : eigenvalues)
        if (!std::isfinite(value))
            throw std::overflow_error("Sternheimer weak normalization metric spectrum is nonfinite.");
    const double threshold = min_metric * std::max(1.0, eigenvalues.back());
    if (eigenvalues.front() <= threshold)
    {
        std::ostringstream message;
        message.precision(17);
        message << "Sternheimer weak normalization metric is singular or ill-conditioned: min="
                << eigenvalues.front() << ", max=" << eigenvalues.back() << ", threshold=" << threshold
                << "; no states dropped or shifts applied.";
        throw std::domain_error(message.str());
    }
    zpotrf_(&upper, &n, gram.data(), &n, &info);
    if (info > 0)
        throw std::domain_error("Sternheimer weak normalization metric is not positive definite.");
    if (info < 0)
        throw std::runtime_error("Sternheimer weak normalization Cholesky factorization failed.");
    CheckOutput(gram);

    // Every tile contains all states before any are overwritten. Right-solving
    // with upper R mixes column j only with columns <=j, preserving each prefix.
    const char right = 'R', no_transpose = 'N', nonunit = 'N';
    const Complex one(1.0);
    for (std::size_t begin = 0; begin < grid_size; begin += tile_rows)
    {
        const int rows = static_cast<int>(std::min(tile_rows, grid_size - begin));
        for (int component = -1; component < 3; ++component)
        {
            for (std::size_t j = 0; j < states.size(); ++j)
            {
                const auto& field = component < 0 ? states[j].values : states[j].gradients[component];
                std::copy_n(field.data() + begin, rows, tile.data() + static_cast<std::size_t>(rows) * j);
            }
            ztrsm_(&right, &upper, &no_transpose, &nonunit, &rows, &n, &one,
                   gram.data(), &n, tile.data(), &rows);
            CheckOutput(tile);
            for (std::size_t j = 0; j < states.size(); ++j)
            {
                auto& field = component < 0 ? states[j].values : states[j].gradients[component];
                std::copy_n(tile.data() + static_cast<std::size_t>(rows) * j, rows, field.data() + begin);
            }
        }
    }
    assemble_gram(eigen_matrix);
    double maximum_error = 0.0;
    for (std::size_t j = 0; j < states.size(); ++j)
        for (std::size_t i = 0; i < states.size(); ++i)
            maximum_error = std::max(maximum_error,
                                    std::abs(eigen_matrix[i + states.size() * j] - Complex(i == j ? 1.0 : 0.0)));
    if (!std::isfinite(maximum_error))
        throw std::overflow_error("Sternheimer weak normalization orthogonality diagnostic overflow.");
    return maximum_error;
}

SternheimerWeakGridOperator::SternheimerWeakGridOperator(std::shared_ptr<const Hamiltonian> fine,
                                                       const Grid& coarse, std::size_t budget)
    : fine_(std::move(fine)), workspace_bytes_(CheckWorkspace(RequireFine(fine_), coarse, budget)),
      transfer_(coarse, fine_->grid()), coarse_dv_(VolumeElement(coarse)),
      fine_dv_(VolumeElement(fine_->grid())), fine_input_(fine_->grid().size()),
      fine_output_(fine_->grid().size()), coarse_temporary_(coarse.size()), coarse_result_(coarse.size())
{}

const SternheimerWeakGridOperator::Grid& SternheimerWeakGridOperator::grid() const { return transfer_.coarse_grid(); }
const SternheimerWeakGridOperator::Hamiltonian& SternheimerWeakGridOperator::fine_hamiltonian() const { return *fine_; }
double SternheimerWeakGridOperator::coarse_volume_element() const { return coarse_dv_; }
double SternheimerWeakGridOperator::fine_volume_element() const { return fine_dv_; }

std::size_t SternheimerWeakGridOperator::workspace_bytes_required(const Hamiltonian& fine, const Grid& coarse)
{
    const auto transfer = SternheimerGridTransfer::gradient_workspace_bytes_required(coarse, fine.grid());
    std::size_t projectors = 0;
    if (fine.nonlocal_projector())
        for (const auto& block : fine.nonlocal_projector()->blocks())
            projectors = std::max(projectors, block.projectors.size());
    const auto vectors = Mul(2, Add(Add(fine.grid().size(), coarse.size()), projectors));
    return Add(transfer, Mul(vectors, sizeof(Complex)));
}

void SternheimerWeakGridOperator::apply(const Vector& coefficients, Vector& result)
{
    Validate(coefficients, grid().size());
    transfer_.apply_negative_laplacian(coefficients, coarse_result_, workspace_bytes_);
    for (auto& value : coarse_result_) value *= fine_->kinetic_prefactor();
    transfer_.interpolate(coefficients, fine_input_);
    fine_->apply_local_potential(fine_input_, fine_output_);
    transfer_.restrict_adjoint(fine_output_, coarse_temporary_);
    for (std::size_t i = 0; i < coarse_result_.size(); ++i) coarse_result_[i] += coarse_temporary_[i];
    fine_->apply_nonlocal(fine_input_, fine_output_);
    transfer_.restrict_adjoint(fine_output_, coarse_temporary_);
    for (std::size_t i = 0; i < coarse_result_.size(); ++i) coarse_result_[i] += coarse_temporary_[i];
    CheckOutput(coarse_result_);
    result = coarse_result_;
}

void SternheimerWeakGridOperator::lift(const Vector& coefficients, Vector& fine_values)
{
    Validate(coefficients, grid().size());
    transfer_.interpolate(coefficients, fine_values);
    for (auto& value : fine_values) value /= std::sqrt(coarse_dv_);
    CheckOutput(fine_values);
}

void SternheimerWeakGridOperator::project(const Vector& fine_values, Vector& coefficients)
{
    Validate(fine_values, fine_->grid().size());
    transfer_.restrict_adjoint(fine_values, coefficients);
    for (auto& value : coefficients) value *= std::sqrt(coarse_dv_);
    CheckOutput(coefficients);
}

SternheimerWeakGridBlocks SternheimerWeakGridOperator::assemble_blocks(const std::vector<Function>& states,
                                                                     double metric_tolerance,
                                                                     std::size_t matrix_budget)
{
    if (states.empty() || !(metric_tolerance > 0.0) || !std::isfinite(metric_tolerance))
        throw std::invalid_argument("Sternheimer weak grid requires states and a finite positive metric tolerance.");
    const auto n = states.size(), nc = std::size_t(grid().size());
    // Five square component matrices plus two rectangular retained couplings.
    const auto count = Add(Mul(5, Mul(n, n)), Mul(2, Mul(n, nc)));
    if (Mul(count, sizeof(Complex)) > matrix_budget)
        throw std::length_error("Sternheimer weak grid matrix budget exceeded.");
    for (const auto& state : states)
    {
        Validate(state.values, fine_->grid().size());
        for (const auto& gradient : state.gradients) Validate(gradient, state.values.size());
    }
    auto matrices = assemble_delta_sternheimer_grid_matrices_fast(*fine_, states, fine_dv_);
    SternheimerWeakGridBlocks blocks;
    blocks.state_count = n;
    blocks.coarse_count = nc;
    for (std::size_t j = 0; j < n; ++j)
        for (std::size_t i = 0; i < n; ++i)
        {
            const auto error = std::abs(matrices.overlap[i + n * j] - Complex(i == j ? 1.0 : 0.0));
            if (!std::isfinite(error) || error > metric_tolerance)
                throw std::invalid_argument("Sternheimer weak grid states are not fine-metric orthonormal.");
            blocks.maximum_metric_error = std::max(blocks.maximum_metric_error, error);
        }
    blocks.state_hamiltonian = std::move(matrices.hamiltonian);
    CheckOutput(blocks.state_hamiltonian);
    blocks.state_coarse_overlap.resize(n * nc);
    blocks.state_coarse_hamiltonian.resize(n * nc);
    const double normalization = std::sqrt(coarse_dv_);
    for (std::size_t i = 0; i < n; ++i)
    {
        transfer_.restrict_adjoint(states[i].values, coarse_temporary_);
        for (std::size_t j = 0; j < nc; ++j)
            blocks.state_coarse_overlap[i + n * j] = normalization * std::conj(coarse_temporary_[j]);
        // The derivative adjoint is essential: FFT differentiation of the NAO
        // values would discard the supplied analytic-gradient weak form.
        transfer_.restrict_gradient_adjoint(states[i].gradients, coarse_result_, workspace_bytes_);
        for (auto& value : coarse_result_) value *= fine_->kinetic_prefactor();
        fine_->apply_local_potential(states[i].values, fine_output_);
        transfer_.restrict_adjoint(fine_output_, coarse_temporary_);
        for (std::size_t j = 0; j < nc; ++j) coarse_result_[j] += coarse_temporary_[j];
        fine_->apply_nonlocal(states[i].values, fine_output_);
        transfer_.restrict_adjoint(fine_output_, coarse_temporary_);
        for (std::size_t j = 0; j < nc; ++j)
            blocks.state_coarse_hamiltonian[i + n * j] = normalization * std::conj(coarse_result_[j] + coarse_temporary_[j]);
    }
    CheckOutput(blocks.state_coarse_overlap);
    CheckOutput(blocks.state_coarse_hamiltonian);
    return blocks;
}
} // namespace ModuleRI
