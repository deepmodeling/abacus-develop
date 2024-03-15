#include "write_HS_R.h"

#include "module_base/timer.h"
#include "write_HS_sparse.h"

// if 'binary=true', output binary file.
// The 'sparse_threshold' is the accuracy of the sparse matrix. 
// If the absolute value of the matrix element is less than or equal to the 'sparse_threshold', it will be ignored.
void ModuleIO::output_HS_R(const int& istep,
                           const ModuleBase::matrix& v_eff,
                           LCAO_Hamilt& UHM,
                           const K_Vectors& kv,
                           hamilt::Hamilt<std::complex<double>>* p_ham,
                           const std::string& SR_filename,
                           const std::string& HR_filename_up,
                           const std::string HR_filename_down,
                           const bool& binary,
                           const double& sparse_threshold)
{
    ModuleBase::TITLE("ModuleIO","output_HS_R"); 
    ModuleBase::timer::tick("ModuleIO","output_HS_R"); 

    if(GlobalV::NSPIN==1||GlobalV::NSPIN==4)
    {
        // jingan add 2021-6-4, modify 2021-12-2
        UHM.calculate_HSR_sparse(0, sparse_threshold, kv.nmp, p_ham);
    }
    else if(GlobalV::NSPIN==2)
    {
        // save HR of current_spin first
        UHM.calculate_HSR_sparse(GlobalV::CURRENT_SPIN, sparse_threshold, kv.nmp, p_ham);
        // calculate HR of the other spin
        if(GlobalV::VL_IN_H)
        {
            int ik = 0;
            if(GlobalV::CURRENT_SPIN == 1)
            {
                ik = 0;
                GlobalV::CURRENT_SPIN = 0;
            } 
            else
            {
                ik = kv.nks / 2;
                GlobalV::CURRENT_SPIN = 1;
            }
            p_ham->refresh();
            p_ham->updateHk(ik);
        }
        UHM.calculate_HSR_sparse(GlobalV::CURRENT_SPIN, sparse_threshold, kv.nmp, p_ham);
    }

    ModuleIO::save_HSR_sparse(istep, *UHM.LM, sparse_threshold, binary, SR_filename, HR_filename_up, HR_filename_down);
    UHM.destroy_all_HSR_sparse();

    ModuleBase::timer::tick("ModuleIO","output_HS_R"); 
    return;
}

void ModuleIO::output_dH_R(const int& istep,
                           const ModuleBase::matrix& v_eff,
                           LCAO_Hamilt& UHM,
                           const K_Vectors& kv,
                           const bool& binary,
                           const double& sparse_threshold)
{
    ModuleBase::TITLE("ModuleIO","output_dH_R"); 
    ModuleBase::timer::tick("ModuleIO","output_dH_R"); 

    UHM.LM->Hloc_fixedR.resize(UHM.LM->ParaV->nnr);
    UHM.GK.allocate_pvdpR();
    if(GlobalV::NSPIN==1||GlobalV::NSPIN==4)
    {
        UHM.calculate_dH_sparse(0, sparse_threshold);
    }
    else if(GlobalV::NSPIN==2)
    {
        for (int ik = 0; ik < kv.nks; ik++)
        {
            if (ik == 0 || ik == kv.nks / 2)
            {
                if(GlobalV::NSPIN == 2)
                {
                    GlobalV::CURRENT_SPIN = kv.isk[ik];
                }

                //note: some MPI process will not have grids when MPI cores is too many, v_eff in these processes are empty
                const double* vr_eff1 = v_eff.nc * v_eff.nr > 0? &(v_eff(GlobalV::CURRENT_SPIN, 0)):nullptr;
                    
                if(!GlobalV::GAMMA_ONLY_LOCAL)
                {
                    if(GlobalV::VL_IN_H)
                    {
                        Gint_inout inout(vr_eff1, GlobalV::CURRENT_SPIN, Gint_Tools::job_type::dvlocal);
                        UHM.GK.cal_gint(&inout);
                    }
                }

                UHM.calculate_dH_sparse(GlobalV::CURRENT_SPIN, sparse_threshold);
            }
        }
    }

    ModuleIO::save_dH_sparse(istep, *UHM.LM, sparse_threshold, binary);
    UHM.destroy_dH_R_sparse();

    UHM.GK.destroy_pvdpR();

    ModuleBase::timer::tick("ModuleIO","output_HS_R"); 
    return;
}

void ModuleIO::output_S_R(
    LCAO_Hamilt &UHM,
    hamilt::Hamilt<std::complex<double>>* p_ham,
    const std::string &SR_filename,
    const bool &binary,
    const double &sparse_threshold)
{
    ModuleBase::TITLE("ModuleIO","output_S_R");
    ModuleBase::timer::tick("ModuleIO","output_S_R"); 

    UHM.calculate_SR_sparse(sparse_threshold, p_ham);
    ModuleIO::save_sparse(UHM.LM->SR_sparse, UHM.LM->all_R_coor, sparse_threshold, binary, SR_filename, *UHM.LM->ParaV, "S", 0);
    UHM.destroy_all_HSR_sparse();

    ModuleBase::timer::tick("ModuleIO","output_S_R");
    return;
}

void ModuleIO::output_T_R(
    const int istep,
    LCAO_Hamilt &UHM,
    const std::string &TR_filename,
    const bool &binary,
    const double &sparse_threshold
)
{
    ModuleBase::TITLE("ModuleIO","output_T_R");
    ModuleBase::timer::tick("ModuleIO","output_T_R"); 

    std::stringstream sst;
    if(GlobalV::CALCULATION == "md" && !GlobalV::out_app_flag)
    {
        sst << GlobalV::global_matrix_dir << istep << "_" << TR_filename;
    }
    else
    {
        sst << GlobalV::global_out_dir << TR_filename;
    }

    UHM.calculate_TR_sparse(sparse_threshold);
    ModuleIO::save_sparse(UHM.LM->TR_sparse, UHM.LM->all_R_coor, sparse_threshold, binary, sst.str().c_str(), *UHM.LM->ParaV, "T", istep);
    UHM.destroy_TR_sparse();

    ModuleBase::timer::tick("ModuleIO","output_T_R");
    return;
}