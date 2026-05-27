#ifndef HSOLVERLCAO_H
#define HSOLVERLCAO_H

#include "source_estate/elecstate.h"
#include "source_hamilt/hamilt.h"
#include "source_basis/module_ao/parallel_orbitals.h"

#include "source_estate/module_charge/charge.h" // mohan add 20251024
#include "source_estate/module_dm/density_matrix.h" // mohan add 20251103

#include <vector>

#ifdef __PEXSI
#include "module_pexsi/simple_pexsi.h"
#endif

namespace hsolver
{

template <typename TK, typename Device = base_device::DEVICE_CPU>
class HSolverLCAO
{
  public:
#ifdef __PEXSI
    HSolverLCAO(const Parallel_Orbitals* ParaV_in,
                std::string method_in,
                std::vector<pexsi::PexsiComplexReuseContext>* pexsi_reuse_contexts_in = nullptr)
        : ParaV(ParaV_in),
          method(method_in),
          pexsi_reuse_contexts_ptr_(pexsi_reuse_contexts_in != nullptr ? pexsi_reuse_contexts_in
                                                                        : &owned_pexsi_reuse_contexts_) {};
#else
    HSolverLCAO(const Parallel_Orbitals* ParaV_in, std::string method_in) : ParaV(ParaV_in), method(method_in) {};
#endif

    void solve(hamilt::Hamilt<TK>* pHamilt,
               psi::Psi<TK>& psi,
               elecstate::ElecState* pes,
			   elecstate::DensityMatrix<TK, double>& dm, // mohan add 2025-11-03
			   Charge &chr, // charge density
			   const int nspin,
			   const bool skip_charge);

  private:
    void hamiltSolvePsiK(hamilt::Hamilt<TK>* hm, psi::Psi<TK>& psi, double* eigenvalue); // for kpar_lcao == 1

    void parakSolve(hamilt::Hamilt<TK>* pHamilt, psi::Psi<TK>& psi, elecstate::ElecState* pes, int kpar); // for kpar_lcao > 1

    // The solving algorithm using cusolver is different from others, so a separate function is needed
    void parakSolve_cusolver(hamilt::Hamilt<TK>* pHamilt,
                             psi::Psi<TK>& psi,
                             elecstate::ElecState* pes);

    const Parallel_Orbitals* ParaV = nullptr;
    
    const std::string method;

#ifdef __PEXSI
    std::vector<pexsi::PexsiComplexReuseContext> owned_pexsi_reuse_contexts_;
    std::vector<pexsi::PexsiComplexReuseContext>* pexsi_reuse_contexts_ptr_ = nullptr;
#endif
};

} // namespace hsolver

#endif
