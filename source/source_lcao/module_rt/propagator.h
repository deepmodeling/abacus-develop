/**
 * @file propagator.h
 * @brief compute propagtor to evolve the wave function
 *  This file originally belonged to file LCAO_evolve.cpp
 */
#ifndef PROPAGATOR_H
#define PROPAGATOR_H

#include "source_base/constants.h"
#include "source_base/module_container/ATen/core/tensor.h" // ct::Tensor
#include "source_lcao/module_rt/kernels/cublasmp_context.h"

#include <complex>

class Parallel_Blkcyclic; // Defined in source/source_base/parallel_blkcyc.h

namespace module_rt
{
//--------------------------------- Utility function ---------------------------------//
#ifdef __MPI
inline int globalIndex(int localindex, int nblk, int nprocs, int myproc)
{
    int iblock, gIndex;
    iblock = localindex / nblk;
    gIndex = (iblock * nprocs + myproc) * nblk + localindex % nblk;
    return gIndex;
}
#endif // __MPI

//--------------------------------- Utility function ---------------------------------//

class Propagator
{
  public:
    Propagator(int ptype, const Parallel_Blkcyclic& layout, const double& dt)
        : ptype(ptype), dt(dt / ModuleBase::AU_to_FS), layout_(layout)
    {
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
                            std::ofstream& ofs_running,
                            const int print_matrix) const;

    template <typename Device>
    void compute_propagator_tensor(const int nlocal,
                                   const ct::Tensor& Stmp,
                                   const ct::Tensor& Htmp,
                                   const ct::Tensor& H_laststep,
                                   ct::Tensor& U_operator,
                                   std::ofstream& ofs_running,
                                   const int print_matrix,
                                   const bool use_lapack,
                                   CublasMpResources& cublas_res) const;
#endif // __MPI

  private:
    int ptype; // type of propagator
    double dt; // time step
    const Parallel_Blkcyclic& layout_; // non-owning block-cyclic layout

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
                                std::ofstream& ofs_running,
                                const int print_matrix) const;

    void compute_propagator_cn2_tensor(const int nlocal,
                                       const ct::Tensor& Stmp,
                                       const ct::Tensor& Htmp,
                                       ct::Tensor& U_operator,
                                       std::ofstream& ofs_running,
                                       const int print_matrix,
                                       CublasMpResources& cublas_res) const;

    template <typename Device>
    void compute_propagator_cn2_tensor_lapack(const int nlocal,
                                              const ct::Tensor& Stmp,
                                              const ct::Tensor& Htmp,
                                              ct::Tensor& U_operator,
                                              std::ofstream& ofs_running,
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
                                   std::ofstream& ofs_running,
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
                                 std::ofstream& ofs_running,
                                 const int print_matrix) const;
#endif // __MPI
};
} // namespace module_rt

#endif
