#include "write_HS_R.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "write_HS.h"
#include "module_base/timer.h"

// if 'binary=true', output binary file.
// The 'sparse_threshold' is the accuracy of the sparse matrix. 
// If the absolute value of the matrix element is less than or equal to the 'sparse_threshold', it will be ignored.
void ModuleIO::output_HS_R(
    const int &istep,
    const ModuleBase::matrix& v_eff,
    LCAO_Hamilt &UHM,
    const std::string &SR_filename,
    const std::string &HR_filename_up,
    const std::string HR_filename_down,
    const bool &binary, 
    const double &sparse_threshold)
{
    ModuleBase::TITLE("ModuleIO","output_HS_R"); 
    ModuleBase::timer::tick("ModuleIO","output_HS_R"); 

    if(GlobalV::NSPIN==1||GlobalV::NSPIN==4)
    {
        // jingan add 2021-6-4, modify 2021-12-2
        UHM.calculate_HSR_sparse(0, sparse_threshold);
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

                const double* vr_eff1 = &(v_eff(GlobalV::CURRENT_SPIN, 0));
                    
                if(!GlobalV::GAMMA_ONLY_LOCAL)
                {
                    if(GlobalV::VL_IN_H)
                    {
                        Gint_inout inout(vr_eff1, GlobalV::CURRENT_SPIN, Gint_Tools::job_type::vlocal);
                        UHM.GK.cal_gint(&inout);
                    }
                }

                UHM.calculate_HSR_sparse(GlobalV::CURRENT_SPIN, sparse_threshold);
            }
        }
    }

    ModuleIO::save_HSR_sparse(istep, *UHM.LM, sparse_threshold, binary, SR_filename, HR_filename_up, HR_filename_down);
    UHM.destroy_all_HSR_sparse();

    if(!GlobalV::GAMMA_ONLY_LOCAL) //LiuXh 20181011
    {
        UHM.GK.destroy_pvpR();
    } //LiuXh 20181011

    ModuleBase::timer::tick("ModuleIO","output_HS_R"); 
    return;
}


void ModuleIO::output_SR(const std::string &SR_filename,
    LCAO_Hamilt &UHM,
    const bool &binary,
    const double &sparse_threshold)
{
    ModuleBase::TITLE("ModuleIO","output_SR");
    ModuleBase::timer::tick("ModuleIO","output_SR"); 

    UHM.calculate_SR_sparse(sparse_threshold);
    ModuleIO::save_SR_sparse(*UHM.LM, sparse_threshold, binary, SR_filename);
    UHM.destroy_all_HSR_sparse();

    ModuleBase::timer::tick("ModuleIO","output_SR");
    return;
}
