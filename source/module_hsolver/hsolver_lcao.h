#ifndef HSOLVERLCAO_H
#define HSOLVERLCAO_H

#include "hsolver.h"
#include "module_basis/module_ao/parallel_orbitals.h"

namespace hsolver
{

template <typename T, typename Device = base_device::DEVICE_CPU>
class HSolverLCAO
{
  public:
    HSolverLCAO(const Parallel_Orbitals* ParaV_in, std::string method_in) : ParaV(ParaV_in), method(method_in) {};

    void solve(hamilt::Hamilt<T>* pHamilt,
               psi::Psi<T>& psi,
               elecstate::ElecState* pes,
               const std::string method_in,
               const bool skip_charge);

  private:
    void hamiltSolvePsiK(hamilt::Hamilt<T>* hm, psi::Psi<T>& psi, double* eigenvalue);

    const Parallel_Orbitals* ParaV;

    void parakSolve(hamilt::Hamilt<T>* pHamilt, psi::Psi<T>& psi, elecstate::ElecState* pes, int kpar);

    bool is_first_scf = true;

    using Real = typename GetTypeReal<T>::type;
    std::vector<Real> precondition_lcao;

    DiagH<T, Device>* pdiagh = nullptr; // for single Hamiltonian matrix diagonal solver

    std::string method = "none";
};

} // namespace hsolver

#endif