#ifndef HSOLVERLCAO_H
#define HSOLVERLCAO_H

#include "source_estate/elecstate.h"
#include "source_hamilt/hamilt.h"
#include "source_basis/module_ao/parallel_orbitals.h"

#include "source_estate/module_charge/charge.h" // mohan add 20251024
#include "source_estate/module_dm/density_matrix.h" // mohan add 20251103

namespace hsolver
{

template <typename TK>
class HSolverLCAO
{
  public:
    HSolverLCAO(const Parallel_Orbitals* ParaV_in,
                const std::string method_in,
                const int kpar_lcao_in,
                const int nlocal_in,
                const int nbands_in,
                const double nelec_in,
                const bool use_gpu_in,
                const int world_nproc_in,
                const int world_rank_in)
        : ParaV(ParaV_in), method(method_in), kpar_lcao(kpar_lcao_in), nlocal(nlocal_in), nbands(nbands_in),
          nelec(nelec_in), use_gpu(use_gpu_in), world_nproc(world_nproc_in), world_rank(world_rank_in){};

    void solve(hamilt::Hamilt<TK>* pHamilt,
               psi::Psi<TK>& psi,
               elecstate::ElecState* pes,
			   elecstate::DensityMatrix<TK, double>& dm, // mohan add 2025-11-03
			   Charge &chr, // charge density
			   const int nspin,
			   const bool skip_charge);

  private:
    void hamiltSolvePsiK(hamilt::Hamilt<TK>* hm, psi::Psi<TK>& psi, double* eigenvalue); // for kpar_lcao == 1

    void parakSolve(hamilt::Hamilt<TK>* pHamilt,
                    psi::Psi<TK>& psi,
                    elecstate::ElecState* pes,
                    const int kpar,
                    const int nspin); // for kpar_lcao > 1

    // The solving algorithm using cusolver is different from others, so a separate function is needed
    void parakSolve_cusolver(hamilt::Hamilt<TK>* pHamilt,
                             psi::Psi<TK>& psi,
                             elecstate::ElecState* pes);

    const Parallel_Orbitals* ParaV = nullptr;

    const std::string method;

    const int kpar_lcao; // number of pools for LCAO diagonalization
    const int nlocal;    // global dimension of the NAO Hamiltonian
    const int nbands;    // number of bands to be solved for
    const double nelec;  // total number of electrons, only used by the pexsi branch
    const bool use_gpu;  // true if running on GPU, only used by the native-ELPA branch
    const int world_nproc;
    const int world_rank;
};

} // namespace hsolver

#endif
