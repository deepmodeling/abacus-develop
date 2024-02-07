#ifndef DIGAOPEXSI_H
#define DIGAOPEXSI_H

#include <vector>
#include "diagh.h"
#include "module_base/global_variable.h"
#include "module_basis/module_ao/parallel_orbitals.h"
#include "module_pexsi/pexsi_solver.h"

namespace hsolver
{

template <typename T>
class DiagoPexsi : public DiagH<T>
{
  private:
    using Real = typename GetTypeReal<T>::type;
    std::vector<double> mu_buffer;

  public:
    DiagoPexsi(const Parallel_Orbitals* ParaV_in)
    {
        mu_buffer.resize(GlobalV::NSPIN);
        for (int i = 0; i < GlobalV::NSPIN; i++)
        {
            mu_buffer[i] = this->ps->pexsi_mu;
        }
        this->ParaV = ParaV_in;
    }
    void diag(hamilt::Hamilt<T>* phm_in, psi::Psi<T>& psi, Real* eigenvalue_in) override;
    const Parallel_Orbitals* ParaV;
    std::vector<T*> DM;
    std::vector<T*> EDM;
    double totalEnergyH;
    double totalEnergyS;
    double totalFreeEnergy;
    pexsi::PEXSI_Solver* ps;

};
} // namespace hsolver

#endif
