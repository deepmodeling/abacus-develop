#ifndef DIGAOPEXSI_H
#define DIGAOPEXSI_H

#include "diagh.h"
#include "module_basis/module_ao/parallel_orbitals.h"
#include "module_pexsi/pexsi_solver.h"

namespace hsolver
{

template <typename T>
class DiagoPexsi : public DiagH<T>
{
  private:
    using Real = typename GetTypeReal<T>::type;

  public:
    DiagoPexsi(const Parallel_Orbitals* ParaV_in)
    {
        this->ParaV = ParaV_in;
    }
    void diag(hamilt::Hamilt<T>* phm_in, psi::Psi<T>& psi, Real* eigenvalue_in) override;
    const Parallel_Orbitals* ParaV;
    T* DM;
    double* EDM;
    double totalEnergyH;
    double totalEnergyS;
    double totalFreeEnergy;
    pexsi::PEXSI_Solver* ps;
};
} // namespace hsolver

#endif
