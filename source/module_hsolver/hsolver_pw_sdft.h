#ifndef HSOLVERPW_SDFT_H
#define HSOLVERPW_SDFT_H
#include "hsolver_pw.h"
#include "module_hamilt_pw/hamilt_stodft/sto_iter.h"
namespace hsolver
{
    class HSolverPW_SDFT : public HSolverPW<double>
    {
        public:
          HSolverPW_SDFT(K_Vectors* pkv,
                         ModulePW::PW_Basis_K* wfc_basis_in,
                         wavefunc* pwf_in,
                         Stochastic_WF& stowf,
                         const int method_sto)
              : HSolverPW(wfc_basis_in, pwf_in)
          {
              this->classname = "HSolverPW_SDFT";
              stoiter.init(method_sto, pkv, wfc_basis_in, stowf);
        }
        virtual void solve(hamilt::Hamilt<double>* pHamilt,
                           psi::Psi<std::complex<double>>& psi,
                           elecstate::ElecState* pes,
                           ModulePW::PW_Basis_K* wfc_basis,
                           Stochastic_WF& stowf,
                           const int istep,
                           const int iter,
                           const std::string method_in,
                           const bool skip_charge) override;
        virtual double set_diagethr(const int istep, const int iter, const double drho) override;                   
        virtual double cal_hsolerror() override
        {
            return 0.0;
        }
        Stochastic_Iter stoiter;
    };
}
#endif