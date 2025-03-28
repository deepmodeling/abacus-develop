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
                                     const std::vector<std::vector<TK>>& V_delta,
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
                        e_delta_band_tmp += dm[is][nu * pv->nrow + mu] * V_delta[0][iic];
                    }
                }
                else
                {
                    for (int ik = 0; ik < nks; ik++)
                    {
                        e_delta_band_tmp += dm[ik][nu * pv->nrow + mu] * V_delta[ik][iic];
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

template <typename TK, typename TH>
void DeePKS_domain::collect_h_mat(const Parallel_Orbitals& pv,
                                  const std::vector<std::vector<TK>>& h_in,
                                  std::vector<TH>& h_out,
                                  const int nlocal,
                                  const int nks)
{
    ModuleBase::TITLE("DeePKS_domain", "collect_h_tot");

    // construct the total H matrix
    for (int k = 0; k < nks; k++)
    {
#ifdef __MPI
        int ir = 0;
        int ic = 0;
        for (int i = 0; i < nlocal; i++)
        {
            std::vector<TK> lineH(nlocal - i, TK(0.0));

            ir = pv.global2local_row(i);
            if (ir >= 0)
            {
                // data collection
                for (int j = i; j < nlocal; j++)
                {
                    ic = pv.global2local_col(j);
                    if (ic >= 0)
                    {
                        int iic = 0;
                        if (ModuleBase::GlobalFunc::IS_COLUMN_MAJOR_KS_SOLVER(PARAM.inp.ks_solver))
                        {
                            iic = ir + ic * pv.nrow;
                        }
                        else
                        {
                            iic = ir * pv.ncol + ic;
                        }
                        lineH[j - i] = h_in[k][iic];
                    }
                }
            }
            else
            {
                // do nothing
            }

            Parallel_Reduce::reduce_all(lineH.data(), nlocal - i);

            for (int j = i; j < nlocal; j++)
            {
                h_out[k](i, j) = lineH[j - i];
                h_out[k](j, i) = h_out[k](i, j); // H is a symmetric matrix
            }
        }
#else
        for (int i = 0; i < nlocal; i++)
        {
            for (int j = i; j < nlocal; j++)
            {
                h_out[k](i, j) = h_in[k][i * nlocal + j];
                h_out[k](j, i) = h_out[k](i, j); // H is a symmetric matrix
            }
        }
#endif
    }
}

template void DeePKS_domain::cal_e_delta_band<double>(const std::vector<std::vector<double>>& dm,
                                                      const std::vector<std::vector<double>>& V_delta,
                                                      const int nks,
                                                      const Parallel_Orbitals* pv,
                                                      double& e_delta_band);
template void DeePKS_domain::cal_e_delta_band<std::complex<double>>(
    const std::vector<std::vector<std::complex<double>>>& dm,
    const std::vector<std::vector<std::complex<double>>>& V_delta,
    const int nks,
    const Parallel_Orbitals* pv,
    double& e_delta_band);

template void DeePKS_domain::collect_h_mat<double, ModuleBase::matrix>(
    const Parallel_Orbitals& pv,
    const std::vector<std::vector<double>>& h_in,
    std::vector<ModuleBase::matrix>& h_out,
    const int nlocal,
    const int nks);

template void DeePKS_domain::collect_h_mat<std::complex<double>, ModuleBase::ComplexMatrix>(
    const Parallel_Orbitals& pv,
    const std::vector<std::vector<std::complex<double>>>& h_in,
    std::vector<ModuleBase::ComplexMatrix>& h_out,
    const int nlocal,
    const int nks);

#endif
