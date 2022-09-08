#include "cal_r_overlap_R.h"
#include "../src_pw/global.h"
#include "write_HS.h"
#include "../module_base/timer.h"
#include "module_esolver/esolver_ks_lcao.h"

namespace ModuleESolver
{
// if 'binary=true', output binary file.
// The 'sparse_threshold' is the accuracy of the sparse matrix. 
// If the absolute value of the matrix element is less than or equal to the 'sparse_threshold', it will be ignored.
void ESolver_KS_LCAO::output_HS_R(
    const std::string &SR_filename,
    const std::string &HR_filename_up,
    const std::string HR_filename_down,
    const bool &binary, 
    const double &sparse_threshold)
{
    ModuleBase::TITLE("ESolver_KS_LCAO","output_HS_R"); 
    ModuleBase::timer::tick("ESolver_KS_LCAO","output_HS_R"); 
    
    // add by jingan for out r_R matrix 2019.8.14
    if(INPUT.out_mat_r)
    {
        cal_r_overlap_R r_matrix;
        r_matrix.init(*this->LOWF.ParaV);
        r_matrix.out_r_overlap_R(GlobalV::NSPIN);
    }

    if(GlobalV::NSPIN==1||GlobalV::NSPIN==4)
    {
        // jingan add 2021-6-4, modify 2021-12-2
        this->UHM.calculate_HSR_sparse(0, sparse_threshold);
    }
    ///*
    else if(GlobalV::NSPIN==2)
    {
        // jingan add 2021-6-4
        for(int ik = 0; ik < GlobalC::kv.nks; ik++)
        {
            if(ik == 0 || ik == GlobalC::kv.nks/2)
            {
                if(GlobalV::NSPIN == 2)
                {
                    GlobalV::CURRENT_SPIN = GlobalC::kv.isk[ik];
                }

                for(int ir = 0; ir < GlobalC::rhopw->nrxx; ir++)
                {
                    GlobalC::pot.vr_eff1[ir] = GlobalC::pot.vr_eff( GlobalV::CURRENT_SPIN, ir);
                }
                    
                if(!GlobalV::GAMMA_ONLY_LOCAL)
                {
                    if(GlobalV::VL_IN_H)
                    {
                        Gint_inout inout(GlobalC::pot.vr_eff1, GlobalV::CURRENT_SPIN, Gint_Tools::job_type::vlocal);
                        this->UHM.GK.cal_gint(&inout);
                    }
                }

                this->UHM.calculate_HSR_sparse(GlobalV::CURRENT_SPIN, sparse_threshold);
            }
        }
    }

    HS_Matrix::save_HSR_sparse(*this->UHM.LM, sparse_threshold, binary, SR_filename, HR_filename_up, HR_filename_down);
    this->UHM.destroy_all_HSR_sparse();

    if(!GlobalV::GAMMA_ONLY_LOCAL) //LiuXh 20181011
    {
        this->UHM.GK.destroy_pvpR();
    } //LiuXh 20181011

    ModuleBase::timer::tick("ESolver_KS_LCAO","output_HS_R"); 
    return;
}


void ESolver_KS_LCAO::output_SR(const std::string &SR_filename, const bool &binary, const double &sparse_threshold)
{
    ModuleBase::TITLE("ESolver_KS_LCAO","output_SR");
    ModuleBase::timer::tick("ESolver_KS_LCAO","output_SR"); 

    this->UHM.calculate_SR_sparse(sparse_threshold);
    HS_Matrix::save_SR_sparse(*this->UHM.LM, sparse_threshold, binary, SR_filename);
    this->UHM.destroy_all_HSR_sparse();

    ModuleBase::timer::tick("ESolver_KS_LCAO","output_SR");
    return;
}
}