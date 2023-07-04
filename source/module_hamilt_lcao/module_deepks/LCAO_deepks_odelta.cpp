//QO 2022-1-7
//This file contains subroutines for calculating O_delta, i.e., corrections of the bandgap,
//which is defind as sum_mu,nu rho^{hl}_mu,nu <chi_mu|alpha>V(D)<alpha|chi_nu>
//where rho^{hl}_mu,nu = C_{L\mu}C_{L\nu} - C_{H\mu}C_{H\nu}, L for LUMO, H for HOMO

//There are two subroutines in this file:
//1. cal_o_delta, which is used for gamma point calculation
//2. cal_o_delta_k, which is used for multi-k calculation

#ifdef __DEEPKS

#include "LCAO_deepks.h"
#include "module_base/parallel_reduce.h"

void LCAO_Deepks::cal_o_delta(const std::vector<std::vector<ModuleBase::matrix>> &dm_hl, const Parallel_Orbitals &ParaO)
{
    ModuleBase::TITLE("LCAO_Deepks", "cal_o_delta");
    this->o_delta.zero_out();
    for (int hl = 0; hl < 1; ++hl)
    {
        for (int i = 0; i < GlobalV::NLOCAL; ++i)
        {
            for (int j = 0; j < GlobalV::NLOCAL; ++j)
            {
                const int mu = ParaO.global2local_row[j];
                const int nu = ParaO.global2local_col[i];
            
                if (mu >= 0 && nu >= 0)
                {                
                    const int index = nu*ParaO.nrow + mu;   
                    for (int is = 0; is < GlobalV::NSPIN; ++is)
                    {
                        this->o_delta(0,hl) += dm_hl[hl][is](nu, mu) * this->H_V_delta[index];
                    }
                }
            }
        }
        Parallel_Reduce::reduce_double_all(this->o_delta(0,hl));
    }
    return;
}


//calculating the correction of (LUMO-HOMO) energies, i.e., band gap corrections
//for multi_k calculations
void LCAO_Deepks::cal_o_delta_k(const std::vector<std::vector<ModuleBase::ComplexMatrix>> &dm_hl,
    const Parallel_Orbitals &ParaO,
    const int nks)
{
    ModuleBase::TITLE("LCAO_Deepks", "cal_o_delta_k");
    
    for(int ik=0; ik<nks; ik++)
    {
        for (int hl=0; hl<1; hl++)
        {
            std::complex<double> o_delta_k=std::complex<double>(0.0,0.0);
            for (int i = 0; i < GlobalV::NLOCAL; ++i)
            {
                for (int j = 0; j < GlobalV::NLOCAL; ++j)
                {
                    const int mu = ParaO.global2local_row[j];
                    const int nu = ParaO.global2local_col[i];
            
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
                        o_delta_k += dm_hl[hl][ik](nu, mu) * this->H_V_delta_k[ik][iic];
                    }
                } //end j
            } //end i
            Parallel_Reduce::reduce_complex_double_all(o_delta_k);
            this->o_delta(ik,hl) = o_delta_k.real();
        }// end hl
    }// end nks
    
    return;
}

#endif
