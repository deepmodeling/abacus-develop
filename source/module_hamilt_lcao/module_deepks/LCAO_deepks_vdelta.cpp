//This file contains subroutines related to V_delta, which is the deepks contribution to Hamiltonian
//defined as |alpha>V(D)<alpha|
#include "module_parameter/parameter.h"
//as well as subroutines for printing them for checking
//It also contains subroutine related to calculating e_delta_bands, which is basically
//tr (rho * V_delta)

//Four subroutines are contained in the file:
//5. cal_e_delta_band : calculates e_delta_bands

#ifdef __DEEPKS

#include "LCAO_deepks.h"
#include "module_base/vector3.h"
#include "module_base/parallel_reduce.h"
#include "module_base/timer.h"

//calculating sum of correction band energies
template <typename TK>
void LCAO_Deepks::cal_e_delta_band(const std::vector<std::vector<TK>>& dm, const int nks)
{
    ModuleBase::TITLE("LCAO_Deepks", "cal_e_delta_band");
    ModuleBase::timer::tick("LCAO_Deepks","cal_e_delta_band");
    TK e_delta_band_tmp = TK(0);
    for (int i = 0; i < PARAM.globalv.nlocal; ++i)
    {
        for (int j = 0; j < PARAM.globalv.nlocal; ++j)
        {
            const int mu = pv->global2local_row(j);
            const int nu = pv->global2local_col(i);
            
            if (mu >= 0 && nu >= 0)
            {                
                int iic;
                if (ModuleBase::GlobalFunc::IS_COLUMN_MAJOR_KS_SOLVER(PARAM.inp.ks_solver))
                {
                    iic = mu + nu * pv->nrow;
                }
                else
                {
                    iic = mu * pv->ncol + nu;
                }
                if constexpr (std::is_same<TK, double>::value)
                {
                    for (int is = 0; is < dm.size(); ++is)  //dm.size() == PARAM.inp.nspin
                    {
                        e_delta_band_tmp += dm[is][nu * this->pv->nrow + mu] * this->H_V_delta[0][iic];
                    }
                }
                else
                {
                    for(int ik = 0; ik < nks; ik++)
                    {
                        e_delta_band_tmp += dm[ik][nu * this->pv->nrow + mu] * this->H_V_delta_k[ik][iic];
                    }
                }
                
            }
        }
    }
    if constexpr (std::is_same<TK, double>::value)
    {
        this->e_delta_band = e_delta_band_tmp;
    }
    else
    {
        this->e_delta_band = e_delta_band_tmp.real();
    }
#ifdef __MPI
    Parallel_Reduce::reduce_all(this->e_delta_band);
#endif
    ModuleBase::timer::tick("LCAO_Deepks","cal_e_delta_band");
    return;
}

template void LCAO_Deepks::cal_e_delta_band<double>(const std::vector<std::vector<double>>& dm, const int nks);
template void LCAO_Deepks::cal_e_delta_band<std::complex<double>>(const std::vector<std::vector<std::complex<double>>>& dm, const int nks);

#endif
