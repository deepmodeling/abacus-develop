//QO 2022-1-7
//This file contains subroutines for calculating O_delta, i.e., corrections of the bandgap,
//which is defind as sum_mu,nu rho^{hl}_mu,nu <chi_mu|alpha>V(D)<alpha|chi_nu>
//where rho^{hl}_mu,nu = C_{L\mu}C_{L\nu} - C_{H\mu}C_{H\nu}, L for LUMO, H for HOMO

//There are two subroutines in this file:
//1. cal_o_delta, which is used for gamma point calculation
//2. cal_o_delta_k, which is used for multi-k calculation

#ifdef __DEEPKS

#include "LCAO_deepks.h"
#include "../src_parallel/parallel_reduce.h"

void LCAO_Deepks::cal_o_delta(const std::vector<ModuleBase::matrix> &dm_hl, const Parallel_Orbitals &ParaO)
{
    ModuleBase::TITLE("LCAO_Deepks", "cal_o_delta");
    this->o_delta = 0;
    for (int i = 0; i < GlobalV::NLOCAL; ++i)
    {
        for (int j = 0; j < GlobalV::NLOCAL; ++j)
        {
            const int mu = ParaO.trace_loc_row[j];
            const int nu = ParaO.trace_loc_col[i];
            
            if (mu >= 0 && nu >= 0)
            {                
                const int index = nu*ParaO.nrow + mu;
                for (int is = 0; is < GlobalV::NSPIN; ++is)
                {
                    this->o_delta += dm_hl[is](nu, mu) * this->H_V_delta[index];
                }
            }
        }
    }
    Parallel_Reduce::reduce_double_all(this->o_delta);
    return;
}


//calculating the correction of (LUMO-HOMO) energies, i.e., band gap corrections
//for multi_k calculations
void LCAO_Deepks::cal_o_delta_k(const std::vector<ModuleBase::ComplexMatrix> &dm_hl,
    const Parallel_Orbitals &ParaO,
    const int nks)
{
    ModuleBase::TITLE("LCAO_Deepks", "cal_o_delta_k");
    std::complex<double> o_delta_k=std::complex<double>(0.0,0.0);
    for (int i = 0; i < GlobalV::NLOCAL; ++i)
    {
        for (int j = 0; j < GlobalV::NLOCAL; ++j)
        {
            const int mu = ParaO.trace_loc_row[j];
            const int nu = ParaO.trace_loc_col[i];
            
            if (mu >= 0 && nu >= 0)
            {                
                int iic;
                if(GlobalV::KS_SOLVER=="genelpa" || GlobalV::KS_SOLVER=="scalapack_gvx")  // save the matrix as column major format
                {
                    iic = mu + nu*ParaO.nrow;
                }
                else
                {
                    iic = mu*ParaO.ncol + nu;
                }
                for(int ik=0; ik<nks; ik++)
                {
                    o_delta_k += dm_hl[ik](nu, mu) * this->H_V_delta_k[ik][iic];
                }
            }
        }
    }
    Parallel_Reduce::reduce_complex_double_all(o_delta_k);
    if(o_delta_k.imag()>1e-12)
    {
        GlobalV::ofs_running << "o_delta_k : " << o_delta_k << std::endl;
        //ModuleBase::WARNING_QUIT("o_delta_k","energy should be real!");
    }
    this->o_delta = o_delta_k.real();
    return;
}

#endif
