#ifndef DIGAOPEXSI_H
#define DIGAOPEXSI_H

#include <vector>
#include <memory>
#include "source_base/macros.h"   // GetRealType
#include "source_base/matrix_block.h"
#include "source_basis/module_ao/parallel_orbitals.h"
#include "source_psi/psi.h"
#include "module_pexsi/pexsi_solver.h"

namespace hsolver
{

template <typename T>
class DiagoPexsi
{
  private:
    using Real = typename GetTypeReal<T>::type;
    static std::vector<double> mu_buffer;

  public:
    DiagoPexsi(const Parallel_Orbitals* ParaV_in,
               const int nspin_in,
               const int nlocal_in,
               const double nelec_in,
               const int world_nproc_in);
    void diag(ModuleBase::MatrixBlock<T>& h_mat,
              ModuleBase::MatrixBlock<T>& s_mat,
              psi::Psi<T>& psi,
              Real* eigenvalue_in);
    const Parallel_Orbitals* ParaV = nullptr;
    std::vector<T*> DM;
    std::vector<T*> EDM;
    double totalEnergyH;
    double totalEnergyS;
    double totalFreeEnergy;
    std::unique_ptr<pexsi::PEXSI_Solver> ps;
    ~DiagoPexsi();

  private:
    /// number of density matrices to keep: nspin, except that nspin == 4 is
    /// treated as a single (spinor) density matrix
    int nspin_dm = 1;
    /// global dimension of the NAO Hamiltonian
    int nlocal = 0;
    double nelec = 0.0;
    int world_nproc = 1;
};
} // namespace hsolver

#endif
