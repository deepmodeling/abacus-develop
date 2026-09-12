#ifndef STERNHEIMER_GALERKIN_OPERATOR_H
#define STERNHEIMER_GALERKIN_OPERATOR_H

#include "source_lcao/module_ri/sternheimer_fd_hamiltonian.h"
#include "source_lcao/module_ri/sternheimer_grid_transfer.h"

#include <cstddef>
#include <memory>

namespace ModuleRI
{

// Reference prototype: H_c = J* H_f J, where J* is the grid-weighted adjoint.
// The supplied H_f remains FD8, including its original local potential and
// nonlocal projectors/D matrices. No point-local coarse potential is defined.
// H_f is the FullGrid FD8 reference, not the old fine weak-A/hybrid response.
// Each concurrent worker must own a separate instance; apply uses mutable
// instance-owned scratch. A shared fine Hamiltonian must remain immutable.
// Apply has no cross-instance lock; FFTW plan management belongs to transfer.
class SternheimerGalerkinOperator
{
  public:
    using Hamiltonian = SternheimerFDHamiltonian;
    using Grid = Hamiltonian::Grid;
    using Complex = Hamiltonian::Complex;
    using Vector = Hamiltonian::Vector;
    using Matrix = Hamiltonian::Matrix;

    SternheimerGalerkinOperator(std::shared_ptr<const Hamiltonian> fine_hamiltonian,
                               const Grid& coarse_grid,
                               std::size_t max_workspace_bytes = 512ULL * 1024 * 1024);
    // Copies H_f; its already-immutable nonlocal projector can remain shared.
    SternheimerGalerkinOperator(const Hamiltonian& fine_hamiltonian,
                               const Grid& coarse_grid,
                               std::size_t max_workspace_bytes = 512ULL * 1024 * 1024);

    SternheimerGalerkinOperator(const SternheimerGalerkinOperator&) = delete;
    SternheimerGalerkinOperator& operator=(const SternheimerGalerkinOperator&) = delete;

    const Grid& grid() const;
    const Hamiltonian& fine_hamiltonian() const;

    // Budget for numerical payload: transfer storage, two retained fine-grid
    // vectors and peak scalar fine-H scratch. Excludes H_f storage (also when
    // copied), caller input/output, FFTW opaque plans, allocator/container and
    // OpenMP runtime overhead. This is not a process-RSS limit. Batch outputs
    // are caller-owned; batches reuse the same single-vector workspace.
    static std::size_t workspace_bytes_required(const Hamiltonian& fine_hamiltonian, const Grid& coarse_grid);
    std::size_t workspace_bytes() const;

    // Scalar and batch input/output may alias. All input dimensions are checked
    // before output modification. Matrix batches are vector-major, as in H_f.
    void apply(const Vector& psi, Vector& hpsi);
    void apply_kinetic(const Vector& psi, Vector& kinetic_psi);
    void apply_local_potential(const Vector& psi, Vector& local_psi);
    void apply_nonlocal(const Vector& psi, Vector& nonlocal_psi);
    void apply_batch(const Matrix& psi, Matrix& hpsi);
    void apply_kinetic_batch(const Matrix& psi, Matrix& kinetic_psi);
    void apply_local_potential_batch(const Matrix& psi, Matrix& local_psi);
    void apply_nonlocal_batch(const Matrix& psi, Matrix& nonlocal_psi);

  private:
    enum class Component { Total, Kinetic, Local, Nonlocal };

    void validate_input(const Vector& psi) const;
    void apply_component(const Vector& psi, Vector& output, Component component);
    void apply_component_batch(const Matrix& psi, Matrix& output, Component component);

    std::shared_ptr<const Hamiltonian> fine_hamiltonian_;
    // Initialization order ensures budget checks precede transfer/FFT allocation.
    std::size_t workspace_bytes_;
    SternheimerGridTransfer transfer_;
    Vector fine_input_;
    Vector fine_output_;
};

} // namespace ModuleRI

#endif
