//QO 2022-1-7
//This file contains subroutines for calculating O_delta, i.e., corrections of the bandgap,
#include "module_parameter/parameter.h"
//which is defind as sum_mu,nu rho^{hl}_mu,nu <chi_mu|alpha>V(D)<alpha|chi_nu>
//where rho^{hl}_mu,nu = C_{L\mu}C_{L\nu} - C_{H\mu}C_{H\nu}, L for LUMO, H for HOMO

#ifdef __DEEPKS

#include "LCAO_deepks.h"
#include "module_base/parallel_reduce.h"

template <typename TK, typename TH>
void LCAO_Deepks::cal_o_delta(const std::vector<std::vector<TH>>& dm_hl, const int nks)
{
    ModuleBase::TITLE("LCAO_Deepks", "cal_o_delta");

    this->o_delta.zero_out();
    for (int ik = 0; ik < nks; ik++)
    {
        for (int hl = 0; hl < 1; ++hl)
        {
            TK o_delta_tmp = TK(0.0);
            for (int i = 0; i < PARAM.globalv.nlocal; ++i)
            {
                for (int j = 0; j < PARAM.globalv.nlocal; ++j)
                {
                    const int mu = pv->global2local_row(j);
                    const int nu = pv->global2local_col(i);
                
                    if (mu >= 0 && nu >= 0)
                    {                
                        int iic;
                        if(PARAM.inp.ks_solver=="genelpa" || PARAM.inp.ks_solver=="scalapack_gvx" || PARAM.inp.ks_solver=="pexsi")  // save the matrix as column major format
                        {
                            iic = mu + nu * pv->nrow;
                        }
                        else
                        {
                            iic = mu * pv->ncol + nu;
                        }
                        if constexpr (std::is_same<TK, double>::value)
                        {
                            for (int is = 0; is < PARAM.inp.nspin; ++is)
                            {
                                o_delta_tmp += dm_hl[hl][is](nu, mu) * this->H_V_delta[0][iic];
                            }
                        }
                        else
                        {
                            o_delta_tmp += dm_hl[hl][ik](nu, mu) * this->H_V_delta_k[ik][iic];
                        }
                    }
                }
            }
            Parallel_Reduce::reduce_all(o_delta_tmp);
            if constexpr (std::is_same<TK, double>::value)
            {
                this->o_delta(ik,hl) = o_delta_tmp;
            }
            else
            {
                this->o_delta(ik,hl) = o_delta_tmp.real();
            }
        }
    }
    return;
}

template void LCAO_Deepks::cal_o_delta<double, ModuleBase::matrix>(
    const std::vector<std::vector<ModuleBase::matrix>>& dm_hl, 
    const int nks);

template void LCAO_Deepks::cal_o_delta<std::complex<double>, ModuleBase::ComplexMatrix>(
    const std::vector<std::vector<ModuleBase::ComplexMatrix>>& dm_hl, 
    const int nks);

#endif
