/**
 * @file propagator.h
 * @brief compute propagtor to evolve the wave function
 *  This file originally belonged to file LCAO_evolve.cpp
 */
#ifndef PROPAGATOR_H
#define PROPAGATOR_H

#include "module_base/constants.h"
#include "module_base/module_container/ATen/core/tensor.h" // ct::Tensor
#include "module_basis/module_ao/parallel_orbitals.h"

#include <complex>

namespace module_tddft
{
class Propagator
{
  public:
    Propagator(const int ptype, const Parallel_Orbitals* pv, const double& dt)
    {
        this->ptype = ptype;
        this->ParaV = pv;
        this->dt = dt / ModuleBase::AU_to_FS;
    }
    ~Propagator();

#ifdef __MPI
    /**
     *  @brief compute propagator
     *
     * @param[in] nlocal number of orbitals
     * @param[in] Stmp overlap matrix
     * @param[in] Htmp H(t+dt/2) or H(t+dt)
     * @param[in] H_laststep H(t)
     * @param[in] print_matirx print internal matrix or not
     * @param[out] U_operator operator of propagator
     */
    void compute_propagator(const int nlocal,
                            const std::complex<double>* Stmp,
                            const std::complex<double>* Htmp,
                            const std::complex<double>* H_laststep,
                            std::complex<double>* U_operator,
                            const int print_matrix) const;

    template <typename Device>
    void compute_propagator_tensor(const int nlocal,
                                   const ct::Tensor& Stmp,
                                   const ct::Tensor& Htmp,
                                   const ct::Tensor& H_laststep,
                                   ct::Tensor& U_operator,
                                   const int print_matrix,
                                   const bool use_lapack) const;
#endif // __MPI

  private:
    int ptype; // type of propagator
    const Parallel_Orbitals* ParaV;
    double dt; // time step

#ifdef __MPI

    /**
     *  @brief compute propagator of method Crank-Nicolson
     *
     * @param[in] nlocal number of orbitals
     * @param[in] Stmp overlap matrix
     * @param[in] Htmp H(t+dt/2) or H(t+dt)
     * @param[in] print_matirx print internal matrix or not
     * @param[out] U_operator operator of propagator
     */
    void compute_propagator_cn2(const int nlocal,
                                const std::complex<double>* Stmp,
                                const std::complex<double>* Htmp,
                                std::complex<double>* U_operator,
                                const int print_matrix) const;

    void compute_propagator_cn2_tensor(const int nlocal,
                                       const ct::Tensor& Stmp,
                                       const ct::Tensor& Htmp,
                                       ct::Tensor& U_operator,
                                       const int print_matrix) const;

    template <typename Device>
    void compute_propagator_cn2_tensor_lapack(const int nlocal,
                                              const ct::Tensor& Stmp,
                                              const ct::Tensor& Htmp,
                                              ct::Tensor& U_operator,
                                              const int print_matrix) const;

    /**
     *  @brief compute propagator of method 4th Taylor
     *
     * @param[in] nlocal number of orbitals
     * @param[in] Stmp overlap matrix
     * @param[in] Htmp H(t+dt/2) or H(t+dt)
     * @param[in] print_matirx print internal matrix or not
     * @param[in] tag a parametre different for 4th Taylor and ETRS
     * @param[out] U_operator operator of propagator
     */
    void compute_propagator_taylor(const int nlocal,
                                   const std::complex<double>* Stmp,
                                   const std::complex<double>* Htmp,
                                   std::complex<double>* U_operator,
                                   const int print_matrix,
                                   const int tag) const;

    /**
     *  @brief compute propagator of method ETRS
     *
     * @param[in] nlocal number of orbitals
     * @param[in] Stmp overlap matrix
     * @param[in] Htmp H(t+dt/2) or H(t+dt)
     * @param[in] H_laststep H(t)
     * @param[in] print_matirx print internal matrix or not
     * @param[out] U_operator operator of propagator
     */
    void compute_propagator_etrs(const int nlocal,
                                 const std::complex<double>* Stmp,
                                 const std::complex<double>* Htmp,
                                 const std::complex<double>* H_laststep,
                                 std::complex<double>* U_operator,
                                 const int print_matrix) const;
#endif // __MPI
};
} // namespace module_tddft

#endif
