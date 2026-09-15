#ifndef STERNHEIMER_WEAK_PRECONDITIONER_H
#define STERNHEIMER_WEAK_PRECONDITIONER_H

#include "source_lcao/module_ri/sternheimer_fd_hamiltonian.h"

#include <memory>

namespace ModuleRI
{

// Analytic coarse-grid inverse kinetic shift for the weak-form complement.
// The symbol is alpha*|G+k|^2-epsilon+i*omega+regularization; it is not the
// finite-difference stencil symbol. Each worker must own one instance because
// FFT execution mutates private scratch. Input and output may alias.
class SternheimerWeakSpectralPreconditioner
{
  public:
    using Complex = SternheimerFDHamiltonian::Complex;
    using Vector = SternheimerFDHamiltonian::Vector;
    using Grid = SternheimerFDHamiltonian::Grid;

    SternheimerWeakSpectralPreconditioner(const Grid& grid,
                                          double kinetic_prefactor,
                                          double reference_eigenvalue,
                                          double omega,
                                          double regularization = 0.0);
    ~SternheimerWeakSpectralPreconditioner();

    SternheimerWeakSpectralPreconditioner(
        const SternheimerWeakSpectralPreconditioner&) = delete;
    SternheimerWeakSpectralPreconditioner& operator=(
        const SternheimerWeakSpectralPreconditioner&) = delete;

    void apply(const Vector& input, Vector& output) const;

  private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

} // namespace ModuleRI

#endif
