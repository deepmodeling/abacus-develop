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
    void compute_propagator(const int nband,
                            const int nlocal,
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

    void compute_propagator_cn2(const int nband,
                                const int nlocal,
                                const std::complex<double>* Stmp,
                                const std::complex<double>* Htmp,
                                std::complex<double>* U_operator,
                                const int print_matrix) const;

    void compute_propagator_taylor(const int nband,
                                   const int nlocal,
                                   const std::complex<double>* Stmp,
                                   const std::complex<double>* Htmp,
                                   std::complex<double>* U_operator,
                                   const int print_matrix,
                                   const int tag) const;

    void compute_propagator_etrs(const int nband,
                                 const int nlocal,
                                 const std::complex<double>* Stmp,
                                 const std::complex<double>* Htmp,
                                 const std::complex<double>* H_laststep,
                                 std::complex<double>* U_operator,
                                 const int print_matrix) const;
#endif
};
} // namespace module_tddft

#endif