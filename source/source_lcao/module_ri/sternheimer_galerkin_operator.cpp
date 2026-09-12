#include "source_lcao/module_ri/sternheimer_galerkin_operator.h"

#include <algorithm>
#include <limits>
#include <stdexcept>
#include <utility>

namespace
{

using Hamiltonian = ModuleRI::SternheimerFDHamiltonian;
using Galerkin = ModuleRI::SternheimerGalerkinOperator;

std::size_t CheckedAdd(const std::size_t left, const std::size_t right)
{
    if (right > std::numeric_limits<std::size_t>::max() - left)
    {
        throw std::length_error("Sternheimer Galerkin workspace size overflows size_t.");
    }
    return left + right;
}

std::size_t CheckedMultiply(const std::size_t left, const std::size_t right)
{
    if (left != 0 && right > std::numeric_limits<std::size_t>::max() / left)
    {
        throw std::length_error("Sternheimer Galerkin workspace size overflows size_t.");
    }
    return left * right;
}

const Hamiltonian& RequireFine(const std::shared_ptr<const Hamiltonian>& fine)
{
    if (!fine)
    {
        throw std::invalid_argument("Sternheimer Galerkin requires a fine Hamiltonian.");
    }
    return *fine;
}

std::size_t CheckedWorkspace(const Hamiltonian& fine, const Hamiltonian::Grid& coarse, const std::size_t budget)
{
    const std::size_t bytes = Galerkin::workspace_bytes_required(fine, coarse);
    if (bytes > budget)
    {
        throw std::length_error("Sternheimer Galerkin numerical workspace exceeds its instance budget.");
    }
    return bytes;
}

std::shared_ptr<const Hamiltonian> CopyFine(const Hamiltonian& fine,
                                          const Hamiltonian::Grid& coarse,
                                          const std::size_t budget)
{
    CheckedWorkspace(fine, coarse, budget);
    return std::make_shared<const Hamiltonian>(fine);
}

} // namespace

namespace ModuleRI
{

SternheimerGalerkinOperator::SternheimerGalerkinOperator(std::shared_ptr<const Hamiltonian> fine_hamiltonian,
                                                       const Grid& coarse_grid,
                                                       const std::size_t max_workspace_bytes)
    : fine_hamiltonian_(std::move(fine_hamiltonian)),
      workspace_bytes_(CheckedWorkspace(RequireFine(fine_hamiltonian_), coarse_grid, max_workspace_bytes)),
      transfer_(coarse_grid, fine_hamiltonian_->grid()),
      fine_input_(static_cast<std::size_t>(fine_hamiltonian_->grid().size())),
      fine_output_(fine_input_.size())
{
}

SternheimerGalerkinOperator::SternheimerGalerkinOperator(const Hamiltonian& fine_hamiltonian,
                                                       const Grid& coarse_grid,
                                                       const std::size_t max_workspace_bytes)
    : SternheimerGalerkinOperator(CopyFine(fine_hamiltonian, coarse_grid, max_workspace_bytes),
                                  coarse_grid, max_workspace_bytes)
{
}

const SternheimerGalerkinOperator::Grid& SternheimerGalerkinOperator::grid() const
{
    return transfer_.coarse_grid();
}

const SternheimerGalerkinOperator::Hamiltonian& SternheimerGalerkinOperator::fine_hamiltonian() const
{
    return *fine_hamiltonian_;
}

std::size_t SternheimerGalerkinOperator::workspace_bytes_required(const Hamiltonian& fine, const Grid& coarse)
{
    if (fine.finite_difference_order() != 8)
    {
        throw std::invalid_argument("Sternheimer Galerkin reference prototype requires the fine FD8 operator.");
    }
    // Transfer validation checks both grid sizes, periodicity, lattice and k
    // before grid.size() or any FFT/fine-vector allocation is used here.
    const std::size_t transfer_bytes = SternheimerGridTransfer::workspace_bytes_required(coarse, fine.grid());
    const std::size_t fine_size = static_cast<std::size_t>(fine.grid().size());
    const std::size_t fine_vector_bytes = CheckedMultiply(fine_size, sizeof(Complex));
    const std::size_t persistent_bytes = CheckedMultiply(2, fine_vector_bytes);
    const auto& directions = fine.cached_first_derivative_directions();
    const auto direction_count = static_cast<std::size_t>(std::count(directions.begin(), directions.end(), true));
    const std::size_t derivative_bytes = CheckedMultiply(direction_count, fine_vector_bytes);
    std::size_t max_projectors = 0;
    if (fine.nonlocal_projector() != nullptr)
    {
        for (const auto& block: fine.nonlocal_projector()->blocks())
        {
            max_projectors = std::max(max_projectors, block.projectors.size());
        }
    }
    const std::size_t projector_bytes = CheckedMultiply(CheckedMultiply(2, max_projectors), sizeof(Complex));
    // Scalar FD derivative buffers are destroyed before the nonlocal apply;
    // nonlocal blocks are processed serially, each with two coefficient vectors.
    return CheckedAdd(CheckedAdd(transfer_bytes, persistent_bytes), std::max(derivative_bytes, projector_bytes));
}

std::size_t SternheimerGalerkinOperator::workspace_bytes() const
{
    return workspace_bytes_;
}

void SternheimerGalerkinOperator::validate_input(const Vector& psi) const
{
    if (psi.size() != static_cast<std::size_t>(grid().size()))
    {
        throw std::invalid_argument("Sternheimer Galerkin input size does not match the coarse grid.");
    }
}

void SternheimerGalerkinOperator::apply_component(const Vector& psi, Vector& output, const Component component)
{
    validate_input(psi);
    transfer_.interpolate(psi, fine_input_);
    switch (component)
    {
    case Component::Total: fine_hamiltonian_->apply(fine_input_, fine_output_); break;
    case Component::Kinetic: fine_hamiltonian_->apply_kinetic(fine_input_, fine_output_); break;
    case Component::Local: fine_hamiltonian_->apply_local_potential(fine_input_, fine_output_); break;
    case Component::Nonlocal: fine_hamiltonian_->apply_nonlocal(fine_input_, fine_output_); break;
    }
    transfer_.restrict_adjoint(fine_output_, output);
}

void SternheimerGalerkinOperator::apply_component_batch(const Matrix& psi, Matrix& output, const Component component)
{
    for (const auto& vector: psi)
    {
        validate_input(vector);
    }
    if (&psi != &output)
    {
        output.resize(psi.size());
    }
    // Stream RHS vectors to bound fine-grid workspace independently of batch
    // length. The fine Hamiltonian still controls its internal OpenMP execution.
    for (std::size_t column = 0; column != psi.size(); ++column)
    {
        apply_component(psi[column], output[column], component);
    }
}

void SternheimerGalerkinOperator::apply(const Vector& psi, Vector& hpsi)
{
    apply_component(psi, hpsi, Component::Total);
}

void SternheimerGalerkinOperator::apply_kinetic(const Vector& psi, Vector& kinetic_psi)
{
    apply_component(psi, kinetic_psi, Component::Kinetic);
}

void SternheimerGalerkinOperator::apply_local_potential(const Vector& psi, Vector& local_psi)
{
    apply_component(psi, local_psi, Component::Local);
}

void SternheimerGalerkinOperator::apply_nonlocal(const Vector& psi, Vector& nonlocal_psi)
{
    apply_component(psi, nonlocal_psi, Component::Nonlocal);
}

void SternheimerGalerkinOperator::apply_batch(const Matrix& psi, Matrix& hpsi)
{
    apply_component_batch(psi, hpsi, Component::Total);
}

void SternheimerGalerkinOperator::apply_kinetic_batch(const Matrix& psi, Matrix& kinetic_psi)
{
    apply_component_batch(psi, kinetic_psi, Component::Kinetic);
}

void SternheimerGalerkinOperator::apply_local_potential_batch(const Matrix& psi, Matrix& local_psi)
{
    apply_component_batch(psi, local_psi, Component::Local);
}

void SternheimerGalerkinOperator::apply_nonlocal_batch(const Matrix& psi, Matrix& nonlocal_psi)
{
    apply_component_batch(psi, nonlocal_psi, Component::Nonlocal);
}

} // namespace ModuleRI
