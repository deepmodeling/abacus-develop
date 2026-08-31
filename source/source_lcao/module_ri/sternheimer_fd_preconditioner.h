#ifndef STERNHEIMER_FD_PRECONDITIONER_H
#define STERNHEIMER_FD_PRECONDITIONER_H

#include "source_lcao/module_ri/sternheimer_fd_hamiltonian.h"

#include <memory>

namespace ModuleRI
{

class SternheimerFDSpectralPreconditioner
{
  public:
    SternheimerFDSpectralPreconditioner(const SternheimerFDHamiltonian& hamiltonian,
                                        double reference_eigenvalue,
                                        double omega,
                                        double regularization);
    ~SternheimerFDSpectralPreconditioner();

    SternheimerFDSpectralPreconditioner(const SternheimerFDSpectralPreconditioner&) = delete;
    SternheimerFDSpectralPreconditioner& operator=(const SternheimerFDSpectralPreconditioner&) = delete;

    void apply(const SternheimerFDHamiltonian::Vector& input,
               SternheimerFDHamiltonian::Vector& output) const;
    void apply_batch(const SternheimerFDHamiltonian::Matrix& input,
                     SternheimerFDHamiltonian::Matrix& output) const;
    void apply_batch(const SternheimerFDHamiltonian::Matrix& input,
                     SternheimerFDHamiltonian::Matrix& output,
                     int* threads_used) const;

    bool is_compatible(const SternheimerFDHamiltonian& hamiltonian,
                       double reference_eigenvalue,
                       double omega,
                       double regularization) const;

  private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

std::shared_ptr<SternheimerFDSpectralPreconditioner>
get_thread_local_sternheimer_fd_spectral_preconditioner(const SternheimerFDHamiltonian& hamiltonian,
                                                         double reference_eigenvalue,
                                                         double omega,
                                                         double regularization);

} // namespace ModuleRI

#endif
