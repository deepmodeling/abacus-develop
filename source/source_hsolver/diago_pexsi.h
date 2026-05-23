#ifndef DIGAOPEXSI_H
#define DIGAOPEXSI_H

#include <vector>
#include <memory>
#include "source_base/macros.h"   // GetRealType
#include "source_hamilt/hamilt.h"
#include "source_basis/module_ao/parallel_orbitals.h"
#include "module_pexsi/pexsi_solver_interface.h"

namespace hsolver
{

template <typename T>
class DiagoPexsi
{
  private:
    using Real = typename GetTypeReal<T>::type;
    static std::vector<double> mu_buffer;

  public:
    explicit DiagoPexsi(const Parallel_Orbitals* ParaV_in,
                        std::unique_ptr<pexsi::IPexsiSolver> solver_in = nullptr);
    void diag(hamilt::Hamilt<T>* phm_in, psi::Psi<T>& psi, Real* eigenvalue_in);
    const Parallel_Orbitals* ParaV = nullptr;
    std::vector<T*> DM;
    std::vector<T*> EDM;
    double totalEnergyH = 0.0;
    double totalEnergyS = 0.0;
    double totalFreeEnergy = 0.0;
    std::unique_ptr<pexsi::IPexsiSolver> ps;
    ~DiagoPexsi();

  private:
    std::vector<std::vector<T>> dm_buffer_;
    std::vector<std::vector<T>> edm_buffer_;
};
} // namespace hsolver

#endif
