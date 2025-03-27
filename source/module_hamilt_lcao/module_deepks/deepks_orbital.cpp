#include "module_parameter/parameter.h"

#ifdef __DEEPKS

#include "deepks_orbital.h"
#include "module_base/parallel_reduce.h"
#include "module_base/timer.h"

template <typename TK, typename TH>
void DeePKS_domain::cal_o_delta(const std::vector<TH>& dm_hl,
                                const std::vector<std::vector<TK>>& h_delta,
                                // std::vector<double>& o_delta,
                                ModuleBase::matrix& o_delta,
                                const Parallel_Orbitals& pv,
                                const int nks)
{
    ModuleBase::TITLE("DeePKS_domain", "cal_o_delta");
    ModuleBase::timer::tick("DeePKS_domain", "cal_o_delta");

    for (int ik = 0; ik < nks; ik++)
    {
        TK o_delta_tmp = TK(0.0);
        for (int i = 0; i < PARAM.globalv.nlocal; ++i)
        {
            for (int j = 0; j < PARAM.globalv.nlocal; ++j)
            {
                const int mu = pv.global2local_row(j);
                const int nu = pv.global2local_col(i);

                if (mu >= 0 && nu >= 0)
                {
                    int iic;
                    if (PARAM.inp.ks_solver == "genelpa" || PARAM.inp.ks_solver == "scalapack_gvx"
                        || PARAM.inp.ks_solver == "pexsi") // save the matrix as column major format
                    {
                        iic = mu + nu * pv.nrow;
                    }
                    else
                    {
                        iic = mu * pv.ncol + nu;
                    }
                    if constexpr (std::is_same<TK, double>::value)
                    {
                        for (int is = 0; is < PARAM.inp.nspin; ++is)
                        {
                            o_delta_tmp += dm_hl[is](nu, mu) * h_delta[ik][iic];
                        }
                    }
                    else
                    {
                        o_delta_tmp += dm_hl[ik](nu, mu) * h_delta[ik][iic];
                    }
                }
            }
        }
        Parallel_Reduce::reduce_all(o_delta_tmp);
        if constexpr (std::is_same<TK, double>::value)
        {
            // o_delta[ik] = o_delta_tmp;
            o_delta(ik, 0) = o_delta_tmp;
        }
        else
        {
            // o_delta[ik] = o_delta_tmp.real();
            o_delta(ik, 0) = o_delta_tmp.real();
        }
    }
    ModuleBase::timer::tick("DeePKS_domain", "cal_o_delta");
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

template void DeePKS_domain::cal_o_delta<double, ModuleBase::matrix>(const std::vector<ModuleBase::matrix>& dm_hl,
                                                                     const std::vector<std::vector<double>>& h_delta,
                                                                     //  std::vector<double>& o_delta,
                                                                     ModuleBase::matrix& o_delta,
                                                                     const Parallel_Orbitals& pv,
                                                                     const int nks);

template void DeePKS_domain::cal_o_delta<std::complex<double>, ModuleBase::ComplexMatrix>(
    const std::vector<ModuleBase::ComplexMatrix>& dm_hl,
    const std::vector<std::vector<std::complex<double>>>& h_delta,
    // std::vector<double>& o_delta,
    ModuleBase::matrix& o_delta,
    const Parallel_Orbitals& pv,
    const int nks);

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
