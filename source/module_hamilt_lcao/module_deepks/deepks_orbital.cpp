#include "module_parameter/parameter.h"

#ifdef __MLALGO

#include "deepks_orbital.h"
#include "module_base/parallel_reduce.h"
#include "module_base/timer.h"

template <typename TK, typename TH>
void DeePKS_domain::cal_o_delta(const std::vector<TH>& dm_hl,
                                const std::vector<std::vector<TK>>& h_delta,
                                // std::vector<double>& o_delta,
                                ModuleBase::matrix& o_delta,
                                const Parallel_Orbitals& pv,
                                const int nks,
                                const int nspin)
{
    ModuleBase::TITLE("DeePKS_domain", "cal_o_delta");
    ModuleBase::timer::tick("DeePKS_domain", "cal_o_delta");

    for (int ik = 0; ik < nks / nspin; ik++)
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
                    for (int is = 0; is < nspin; is++)
                    {
                        o_delta_tmp += dm_hl[ik + is * nks / nspin](nu, mu) * h_delta[ik][iic];
                    }
                }
            }
        }
        Parallel_Reduce::reduce_all(o_delta_tmp);
        
        const double* o_delta_ptr = reinterpret_cast<const double*>(&o_delta_tmp);
        o_delta(ik, 0) = o_delta_ptr[0]; // real part in complex case
    }
    ModuleBase::timer::tick("DeePKS_domain", "cal_o_delta");
    return;
}

void DeePKS_domain::check_o_delta(ModuleBase::matrix& o_delta)
{
    std::ofstream ofs("o_delta.dat");
    ofs << std::setprecision(10);
    for (int ik = 0; ik < o_delta.nr; ik++)
    {
        ofs << o_delta(ik, 0) << std::endl;
    }
}

template void DeePKS_domain::cal_o_delta<double, ModuleBase::matrix>(const std::vector<ModuleBase::matrix>& dm_hl,
                                                                     const std::vector<std::vector<double>>& h_delta,
                                                                     //  std::vector<double>& o_delta,
                                                                     ModuleBase::matrix& o_delta,
                                                                     const Parallel_Orbitals& pv,
                                                                     const int nks,
                                                                     const int nspin);

template void DeePKS_domain::cal_o_delta<std::complex<double>, ModuleBase::ComplexMatrix>(
    const std::vector<ModuleBase::ComplexMatrix>& dm_hl,
    const std::vector<std::vector<std::complex<double>>>& h_delta,
    // std::vector<double>& o_delta,
    ModuleBase::matrix& o_delta,
    const Parallel_Orbitals& pv,
    const int nks,
    const int nspin);

#endif
