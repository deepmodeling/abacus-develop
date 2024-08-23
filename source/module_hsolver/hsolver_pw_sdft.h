#ifndef HSOLVERPW_SDFT_H
#define HSOLVERPW_SDFT_H
#include "hsolver_pw.h"
#include "module_hamilt_pw/hamilt_stodft/sto_iter.h"
namespace hsolver
{
class HSolverPW_SDFT : public HSolverPW<std::complex<double>>
{
  public:
    HSolverPW_SDFT(K_Vectors* pkv,
                   ModulePW::PW_Basis_K* wfc_basis_in,
                   wavefunc* pwf_in,
                   Stochastic_WF& stowf,
                   StoChe<double>& stoche,
                   const bool initialed_psi_in)
        : HSolverPW(wfc_basis_in, pwf_in, initialed_psi_in)
    {
        stoiter.init(pkv, wfc_basis_in, stowf, stoche);
    }

    virtual void solve(hamilt::Hamilt<std::complex<double>>* pHamilt,
                       psi::Psi<std::complex<double>>& psi,
                       elecstate::ElecState* pes,
                       ModulePW::PW_Basis_K* wfc_basis,
                       Stochastic_WF& stowf,
                       const int istep,
                       const int iter,
                       const std::string method_in,
                       const int scf_iter_in,
                       const bool need_subspace_in,
                       const int diag_iter_max_in,
                       const double pw_diag_thr_in,
                       const bool skip_charge);

    Stochastic_Iter stoiter;
};
} // namespace hsolver
#endif