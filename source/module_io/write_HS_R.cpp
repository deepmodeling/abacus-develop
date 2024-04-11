#include "module_base/timer.h"
#include "write_HS_R.h"
#include "write_HS_sparse.h"
#include "module_hamilt_lcao/hamilt_lcaodft/spar_hsr.h"
#include "module_hamilt_lcao/hamilt_lcaodft/spar_dh.h"
#include "module_hamilt_lcao/hamilt_lcaodft/spar_st.h"


// if 'binary=true', output binary file.
// The 'sparse_thr' is the accuracy of the sparse matrix. 
// If the absolute value of the matrix element is less than or equal to the 'sparse_thr', it will be ignored.
void ModuleIO::output_HSR(const int& istep,
                           const ModuleBase::matrix& v_eff,
                           const Parallel_Orbitals &pv,
                           LCAO_Matrix& lm,
                           Grid_Driver &grid, // mohan add 2024-04-06
                           const K_Vectors& kv,
                           hamilt::Hamilt<std::complex<double>>* p_ham,
                           const std::string& SR_filename,
                           const std::string& HR_filename_up,
                           const std::string HR_filename_down,
                           const bool& binary,
                           const double& sparse_thr)
{
    ModuleBase::TITLE("ModuleIO","output_HSR"); 
    ModuleBase::timer::tick("ModuleIO","output_HSR"); 

    const int nspin = GlobalV::NSPIN;

    if(nspin==1||nspin==4)
    {
        const int spin_now = 0;
        // jingan add 2021-6-4, modify 2021-12-2
		sparse_format::cal_HSR(
				pv,
				lm,
                grid,
				spin_now, 
				sparse_thr, 
				kv.nmp, 
				p_ham);
	}
    else if(nspin==2)
    {
        const int spin_now = GlobalV::CURRENT_SPIN;

        // save HR of current_spin first
		sparse_format::cal_HSR(
				pv,
				lm,
				grid,
				spin_now, 
				sparse_thr, 
				kv.nmp, 
				p_ham);

        // cal HR of the other spin
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

		sparse_format::cal_HSR(
				pv,
				lm,
				grid,
				GlobalV::CURRENT_SPIN, 
				sparse_thr, 
				kv.nmp, 
				p_ham);
    }

	ModuleIO::save_HSR_sparse(
			istep, 
			lm, 
			sparse_thr, 
			binary, 
			SR_filename, 
			HR_filename_up, 
			HR_filename_down);

	lm.destroy_HS_R_sparse();

    ModuleBase::timer::tick("ModuleIO","output_HSR"); 
    return;
}


void ModuleIO::output_dHR(const int &istep,
                           const ModuleBase::matrix &v_eff,
                           LCAO_gen_fixedH &gen_h, // mohan add 2024-04-02
                           Gint_k &gint_k,  // mohan add 2024-04-01
                           LCAO_Matrix &lm,  // mohan add 2024-04-01
                           Grid_Driver &grid, // mohan add 2024-04-06
                           const K_Vectors& kv,
                           const bool &binary,
                           const double &sparse_thr)
{
    ModuleBase::TITLE("ModuleIO","output_dHR"); 
    ModuleBase::timer::tick("ModuleIO","output_dHR"); 

    lm.Hloc_fixedR.resize(lm.ParaV->nnr);

    gint_k.allocate_pvdpR();

    const int nspin = GlobalV::NSPIN;

    if(nspin==1||nspin==4)
    {
        // mohan add 2024-04-01
        const int cspin = GlobalV::CURRENT_SPIN; 

		sparse_format::cal_dH(
                lm,
                grid,
				gen_h,
				cspin, 
				sparse_thr, 
				gint_k);
	}
    else if(nspin==2)
    {
        for (int ik = 0; ik < kv.nks; ik++)
        {
            if (ik == 0 || ik == kv.nks / 2)
            {
                if(nspin == 2)
                {
                    GlobalV::CURRENT_SPIN = kv.isk[ik];
                }

                // note: some MPI process will not have grids when MPI cores are too many, 
                // v_eff in these processes are empty
                const double* vr_eff1 = v_eff.nc * v_eff.nr > 0? &(v_eff(GlobalV::CURRENT_SPIN, 0)):nullptr;
                    
                if(!GlobalV::GAMMA_ONLY_LOCAL)
                {
                    if(GlobalV::VL_IN_H)
                    {
                        Gint_inout inout(vr_eff1, GlobalV::CURRENT_SPIN, Gint_Tools::job_type::dvlocal);
                        gint_k.cal_gint(&inout);
                    }
                }

                const int cspin = GlobalV::CURRENT_SPIN;

				sparse_format::cal_dH(
                        lm,
                        grid,
						gen_h,
						cspin, 
						sparse_thr, 
						gint_k);
			}
        }
    }

    // mohan update 2024-04-01
    ModuleIO::save_dH_sparse(istep, lm, sparse_thr, binary);

    lm.destroy_dH_R_sparse();

    gint_k.destroy_pvdpR();

    ModuleBase::timer::tick("ModuleIO","output_dHR"); 
    return;
}

void ModuleIO::output_SR(
    Parallel_Orbitals &pv, 
    LCAO_Matrix &lm,
    Grid_Driver &grid,
    hamilt::Hamilt<std::complex<double>>* p_ham,
    const std::string &SR_filename,
    const bool &binary,
    const double &sparse_thr)
{
    ModuleBase::TITLE("ModuleIO","output_SR");
    ModuleBase::timer::tick("ModuleIO","output_SR"); 

	sparse_format::cal_SR(
            pv,
			lm.all_R_coor,
			lm.SR_sparse,
			lm.SR_soc_sparse,
			grid,
			sparse_thr, 
			p_ham);

    const int istep=0;

	ModuleIO::save_sparse(
			lm.SR_sparse, 
			lm.all_R_coor,
			sparse_thr, 
			binary, 
			SR_filename,
			*lm.ParaV, 
			"S", 
			istep
			);

    lm.destroy_HS_R_sparse();

    ModuleBase::timer::tick("ModuleIO","output_SR");
    return;
}

void ModuleIO::output_TR(
    const int istep,
    const UnitCell &ucell,
    const Parallel_Orbitals &pv,
    LCAO_Matrix &lm,
    Grid_Driver &grid,
    LCAO_gen_fixedH &gen_h, // mohan add 2024-04-02
    const std::string &TR_filename,
    const bool &binary,
    const double &sparse_thr
)
{
    ModuleBase::TITLE("ModuleIO","output_TR");
    ModuleBase::timer::tick("ModuleIO","output_TR"); 

    std::stringstream sst;
    if(GlobalV::CALCULATION == "md" && !GlobalV::out_app_flag)
    {
        sst << GlobalV::global_matrix_dir << istep << "_" << TR_filename;
    }
    else
    {
        sst << GlobalV::global_out_dir << TR_filename;
    }

	sparse_format::cal_TR(
			ucell,
			pv,
			lm,
			grid,
			gen_h, 
			sparse_thr);

	ModuleIO::save_sparse(
			lm.TR_sparse, 
			lm.all_R_coor,
			sparse_thr, 
			binary, 
		    sst.str().c_str(),
			*(lm.ParaV), 
			"T", 
			istep
			);

    lm.destroy_T_R_sparse();

    ModuleBase::timer::tick("ModuleIO","output_TR");
    return;
}
