/**
 * @file propagator.h
 * @brief compute propagtor to evolve the wave function
 *  This file originally belonged to file LCAO_evolve.cpp
 */
#ifndef PROPAGATOR_H
#define PROPAGATOR_H

#include "module_basis/module_ao/parallel_orbitals.h"

namespace module_tddft
{
class Propagator
{
  public:
    Propagator(const int ptype, const Parallel_Orbitals* pv)
    {
        this->ptype = ptype;
        this->ParaV = pv;
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
#endif

  private:
    int ptype; // type of propagator
    const Parallel_Orbitals* ParaV;

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
#endif
};
} // namespace module_tddft

#endif