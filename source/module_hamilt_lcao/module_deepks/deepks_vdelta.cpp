// This file contains subroutines related to V_delta, which is the deepks contribution to Hamiltonian
// defined as |alpha>V(D)<alpha|
#include "module_parameter/parameter.h"
// as well as subroutines for printing them for checking
// It also contains subroutine related to calculating e_delta_bands, which is basically
// tr (rho * V_delta)

// One subroutine is contained in the file:
// 1. cal_e_delta_band : calculates e_delta_bands

#ifdef __DEEPKS

#include "deepks_vdelta.h"
#include "module_base/global_function.h"
#include "module_base/parallel_reduce.h"
#include "module_base/timer.h"
#include "module_base/tool_title.h"

// calculating sum of correction band energies
template <typename TK>
void DeePKS_domain::cal_e_delta_band(const std::vector<std::vector<TK>>& dm,
                                     const std::vector<std::vector<TK>>& H_V_delta,
                                     const int nks,
                                     const Parallel_Orbitals* pv,
                                     double& e_delta_band)
{
    ModuleBase::TITLE("DeePKS_domain", "cal_e_delta_band");
    ModuleBase::timer::tick("DeePKS_domain", "cal_e_delta_band");
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
                    for (int is = 0; is < dm.size(); ++is) // dm.size() == PARAM.inp.nspin
                    {
                        e_delta_band_tmp += dm[is][nu * pv->nrow + mu] * H_V_delta[0][iic];
                    }
                }
                else
                {
                    for (int ik = 0; ik < nks; ik++)
                    {
                        e_delta_band_tmp += dm[ik][nu * pv->nrow + mu] * H_V_delta[ik][iic];
                    }
                }
            }
        }
    }
    if constexpr (std::is_same<TK, double>::value)
    {
        e_delta_band = e_delta_band_tmp;
    }
    else
    {
        e_delta_band = e_delta_band_tmp.real();
    }
#ifdef __MPI
    Parallel_Reduce::reduce_all(e_delta_band);
#endif
    ModuleBase::timer::tick("DeePKS_domain", "cal_e_delta_band");
    return;
}

template void DeePKS_domain::cal_e_delta_band<double>(const std::vector<std::vector<double>>& dm,
                                                      const std::vector<std::vector<double>>& H_V_delta,
                                                      const int nks,
                                                      const Parallel_Orbitals* pv,
                                                      double& e_delta_band);
template void DeePKS_domain::cal_e_delta_band<std::complex<double>>(
    const std::vector<std::vector<std::complex<double>>>& dm,
    const std::vector<std::vector<std::complex<double>>>& H_V_delta,
    const int nks,
    const Parallel_Orbitals* pv,
    double& e_delta_band);

#endif
