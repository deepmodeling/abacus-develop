#include "write_HS_R.h"

#include "module_base/timer.h"
#include "module_hamilt_lcao/hamilt_lcaodft/LCAO_HS_arrays.hpp"
#include "module_hamilt_lcao/hamilt_lcaodft/spar_dh.h"
#include "module_hamilt_lcao/hamilt_lcaodft/spar_hsr.h"
#include "module_hamilt_lcao/hamilt_lcaodft/spar_st.h"
#include "write_HS_sparse.h"

// if 'binary=true', output binary file.
// The 'sparse_thr' is the accuracy of the sparse matrix.
// If the absolute value of the matrix element is less than or equal to the 'sparse_thr', it will be ignored.
void ModuleIO::output_HSR(const int& istep,
                          const ModuleBase::matrix& v_eff,
                          const Parallel_Orbitals& pv,
                          LCAO_Matrix& lm,
                          Grid_Driver& grid, // mohan add 2024-04-06
                          const K_Vectors& kv,
                          hamilt::Hamilt<std::complex<double>>* p_ham,
                          const std::string& SR_filename,
                          const std::string& HR_filename_up,
                          const std::string HR_filename_down,
                          const bool& binary,
                          const double& sparse_thr)
{
    ModuleBase::TITLE("ModuleIO", "output_HSR");
    ModuleBase::timer::tick("ModuleIO", "output_HSR");

    const int nspin = GlobalV::NSPIN;

    if (nspin == 1 || nspin == 4)
    {
        const int spin_now = 0;
        // jingan add 2021-6-4, modify 2021-12-2
        sparse_format::cal_HSR(pv, lm, grid, spin_now, sparse_thr, kv.nmp, p_ham);
    }
    else if (nspin == 2)
    {
        int spin_now = 1;

        // save HR of spin down first (the current spin always be down)
        sparse_format::cal_HSR(pv, lm, grid, spin_now, sparse_thr, kv.nmp, p_ham);

        // cal HR of the spin up
        if (GlobalV::VL_IN_H)
        {
            const int ik = 0;
            p_ham->refresh();
            p_ham->updateHk(ik);
            spin_now = 0;
        }

        sparse_format::cal_HSR(pv, lm, grid, spin_now, sparse_thr, kv.nmp, p_ham);
    }

    ModuleIO::save_HSR_sparse(istep, lm, sparse_thr, binary, SR_filename, HR_filename_up, HR_filename_down);

    lm.destroy_HS_R_sparse();

    ModuleBase::timer::tick("ModuleIO", "output_HSR");
    return;
}

void ModuleIO::output_dHR(const int& istep,
                          const ModuleBase::matrix& v_eff,
                          Gint_k& gint_k,    // mohan add 2024-04-01
                          LCAO_Matrix& lm,   // mohan add 2024-04-01
                          Grid_Driver& grid, // mohan add 2024-04-06
                          const ORB_gen_tables* uot,
                          const K_Vectors& kv,
                          const bool& binary,
                          const double& sparse_thr)
{
    ModuleBase::TITLE("ModuleIO", "output_dHR");
    ModuleBase::timer::tick("ModuleIO", "output_dHR");

    gint_k.allocate_pvdpR();

    const int nspin = GlobalV::NSPIN;

    if (nspin == 1 || nspin == 4)
    {
        // mohan add 2024-04-01
        const int cspin = 0;

        sparse_format::cal_dH(lm, grid, uot, cspin, sparse_thr, gint_k);
    }
    else if (nspin == 2)
    {
        for (int cspin = 0; cspin < 2; cspin++)
        {
            // note: some MPI process will not have grids when MPI cores are too many,
            // v_eff in these processes are empty
            const double* vr_eff1 = v_eff.nc * v_eff.nr > 0 ? &(v_eff(cspin, 0)) : nullptr;

            if (!GlobalV::GAMMA_ONLY_LOCAL)
            {
                if (GlobalV::VL_IN_H)
                {
                    Gint_inout inout(vr_eff1, cspin, Gint_Tools::job_type::dvlocal);
                    gint_k.cal_gint(&inout);
                }
            }

            sparse_format::cal_dH(lm, grid, uot, cspin, sparse_thr, gint_k);
        }
    }
    // mohan update 2024-04-01
    ModuleIO::save_dH_sparse(istep, lm, sparse_thr, binary);

    lm.destroy_dH_R_sparse();

    gint_k.destroy_pvdpR();

    ModuleBase::timer::tick("ModuleIO", "output_dHR");
    return;
}

void ModuleIO::output_SR(Parallel_Orbitals& pv,
                         LCAO_Matrix& lm,
                         Grid_Driver& grid,
                         hamilt::Hamilt<std::complex<double>>* p_ham,
                         const std::string& SR_filename,
                         const bool& binary,
                         const double& sparse_thr)
{
    ModuleBase::TITLE("ModuleIO", "output_SR");
    ModuleBase::timer::tick("ModuleIO", "output_SR");

    sparse_format::cal_SR(pv, lm.all_R_coor, lm.SR_sparse, lm.SR_soc_sparse, grid, sparse_thr, p_ham);

    const int istep = 0;

    ModuleIO::save_sparse(lm.SR_sparse, lm.all_R_coor, sparse_thr, binary, SR_filename, *lm.ParaV, "S", istep);

    lm.destroy_HS_R_sparse();

    ModuleBase::timer::tick("ModuleIO", "output_SR");
    return;
}

void ModuleIO::output_TR(const int istep,
                         const UnitCell& ucell,
                         const Parallel_Orbitals& pv,
                         LCAO_Matrix& lm,
                         Grid_Driver& grid,
                         const ORB_gen_tables* uot,
                         const std::string& TR_filename,
                         const bool& binary,
                         const double& sparse_thr)
{
    ModuleBase::TITLE("ModuleIO", "output_TR");
    ModuleBase::timer::tick("ModuleIO", "output_TR");

    std::stringstream sst;
    if (GlobalV::CALCULATION == "md" && !GlobalV::out_app_flag)
    {
        sst << GlobalV::global_matrix_dir << istep << "_" << TR_filename;
    }
    else
    {
        sst << GlobalV::global_out_dir << TR_filename;
    }

    // need Hloc_fixedR
    LCAO_HS_Arrays HS_arrays;

    sparse_format::cal_TR(ucell, pv, lm, HS_arrays, grid, uot, sparse_thr);

    ModuleIO::save_sparse(lm.TR_sparse, lm.all_R_coor, sparse_thr, binary, sst.str().c_str(), *(lm.ParaV), "T", istep);

    lm.destroy_T_R_sparse();

    ModuleBase::timer::tick("ModuleIO", "output_TR");
    return;
}
