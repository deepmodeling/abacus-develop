#ifndef HSOLVERPW_SDFT_H
#define HSOLVERPW_SDFT_H
#include "hsolver_pw.h"
#include "src_pw/sto_iter.h"
namespace hsolver
{
    class HSolverPW_SDFT : public HSolverPW
    {
        public:
        HSolverPW_SDFT(ModulePW::PW_Basis_K* wfc_basis_in, Stochastic_WF& stowf, const int method_sto):HSolverPW(wfc_basis_in)
        {
            this->classname = "HSolverPW_SDFT";
            stoiter.init(wfc_basis_in->npwk_max, stowf.nchip, method_sto, stowf);
        }
        virtual void solve(hamilt::Hamilt* pHamilt, 
                           psi::Psi<std::complex<double>>& psi, 
                           elecstate::ElecState* pes, 
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