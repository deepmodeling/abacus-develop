#include "hamilt_lcao.h"

#include "module_base/global_variable.h"
#include "module_base/timer.h"
#include "src_lcao/global_fp.h"
#include "src_pw/global.h"
#include "src_lcao/dftu.h"
#include "module_hsolver/diago_elpa.h"
#ifdef __DEEPKS
#include "module_deepks/LCAO_deepks.h"
#endif
#include "module_hsolver/hsolver_lcao.h"
#include "module_xc/xc_functional.h"

namespace hamilt
{
// case for nspin<4, gamma-only k-point
template class HamiltLCAO<double>;
// case for nspin<4, multi-k-points
// case for nspin == 4, non-collinear spin case
template class HamiltLCAO<std::complex<double>>;

template <typename T> void HamiltLCAO<T>::getMatrix(MatrixBlock<T> &hk_in, MatrixBlock<T> &sk_in)
{
    hk_in = MatrixBlock<T>{hmatrix_k,
                           (size_t)this->LM->ParaV->nrow,
                           (size_t)this->LM->ParaV->ncol,
                           this->LM->ParaV->desc};
    sk_in = MatrixBlock<T>{smatrix_k,
                           (size_t)this->LM->ParaV->nrow,
                           (size_t)this->LM->ParaV->ncol,
                           this->LM->ParaV->desc};
}

// case for nspin==4
/*template <>
void HamiltLCAO<std::complex<double>, std::complex<double>>::matrix(MatrixBlock<std::complex<double>> &hk_in,
                                                      MatrixBlock<std::complex<double>> &sk_in)
{
    this->getMatrix(hk_in, sk_in);
}*/
// case for nspin<4, multi-k-points
template <>
void HamiltLCAO<std::complex<double>>::matrix(MatrixBlock<std::complex<double>> &hk_in,
                                                      MatrixBlock<std::complex<double>> &sk_in)
{
    this->getMatrix(hk_in, sk_in);
}

// case for nspin<4, gamma_only
template <> void HamiltLCAO<double>::matrix(MatrixBlock<double> &hk_in, MatrixBlock<double> &sk_in)
{
    this->getMatrix(hk_in, sk_in);
}

template <> void HamiltLCAO<double>::updateHk(const int ik)
{
    ModuleBase::TITLE("HamiltLCAO", "updateHk");
    ModuleBase::timer::tick("HamiltLCAO", "updateHk");
    //this->hk_fixed_mock(ik);
    //this->hk_update_mock(ik);

    //-----------------------------------------
    //(1) prepare data for this k point.
    // copy the local potential from array.
    //-----------------------------------------
    if (GlobalV::NSPIN == 2)
    {
        GlobalV::CURRENT_SPIN = GlobalC::kv.isk[ik];
    }

    for (int ir = 0; ir < GlobalC::rhopw->nrxx; ir++)
    {
        GlobalC::pot.vr_eff1[ir] = GlobalC::pot.vr_eff(GlobalV::CURRENT_SPIN, ir);
        if(XC_Functional::get_func_type()==3)
        {
            GlobalC::pot.vofk_eff1[ir] = GlobalC::pot.vofk(GlobalV::CURRENT_SPIN, ir);
        }
    }

    if (!this->uhm->init_s)
    {
        ModuleBase::WARNING_QUIT("Hamilt_Linear::solve_using_cg", "Need init S matrix firstly");
    }

    //--------------------------------------------
    // (3) folding matrix, 
    // and diagonalize the H matrix (T+Vl+Vnl).
    //--------------------------------------------

    // Peize Lin add ik 2016-12-03
    this->uhm->calculate_Hgamma(ik, this->loc->dm_gamma);

    // Effective potential of DFT+U is added to total Hamiltonian here; Quxin adds on 20201029
    if (INPUT.dft_plus_u)
    {
        std::vector<double> eff_pot(this->lowf->ParaV->nloc);
        GlobalC::dftu.cal_eff_pot_mat_real(ik, this->non_first_scf, &eff_pot[0]);

        const int spin = GlobalC::kv.isk[ik];
        for (int irc = 0; irc < this->lowf->ParaV->nloc; irc++)
            this->LM->Hloc[irc] += eff_pot[irc];

    }

    // Peize Lin add at 2020.04.04
    if (GlobalC::restart.info_load.load_H && !GlobalC::restart.info_load.load_H_finish)
    {
        GlobalC::restart.load_disk(this->LM[0], "H", ik);
        GlobalC::restart.info_load.load_H_finish = true;
    }
    if (GlobalC::restart.info_save.save_H)
    {
        GlobalC::restart.save_disk(this->LM[0], "H", ik);
    }

    this->hmatrix_k = this->LM->Hloc.data();
    if( (GlobalC::CHR.get_new_e_iteration() && ik==0) || hsolver::HSolverLCAO::out_mat_hs)
    {
        if(this->smatrix_k==nullptr)
        {
            this->smatrix_k = new double[this->LM->Sloc.size()];
            this->allocated_smatrix = true;
        }
        const int inc = 1;
        BlasConnector::copy(this->LM->Sloc.size(), this->LM->Sloc.data(), inc, this->smatrix_k, inc);
        hsolver::DiagoElpa::is_already_decomposed = false;
    }
    ModuleBase::timer::tick("HamiltLCAO", "updateHk");
    return;
}

// gamma_only case, do nothing before loop-k
template<>
void HamiltLCAO<double>::constructHamilt()
{
    assert(GlobalV::NSPIN == GlobalC::kv.nks);
    return;
}

// multi-k case, calculate hamiltonian matrix k
template <> void HamiltLCAO<std::complex<double>>::updateHk(const int ik)
{
    ModuleBase::TITLE("HamiltLCAO", "updateHk");
    ModuleBase::timer::tick("HamiltLCAO", "each_k");
    //-----------------------------------------
    //(1) prepare data for this k point.
    // copy the local potential from array.
    //-----------------------------------------
    if (GlobalV::NSPIN == 2)
    {
        GlobalV::CURRENT_SPIN = GlobalC::kv.isk[ik];
    }
    for (int ir = 0; ir < GlobalC::rhopw->nrxx; ir++)
    {
        GlobalC::pot.vr_eff1[ir] = GlobalC::pot.vr_eff(GlobalV::CURRENT_SPIN, ir);
        if(XC_Functional::get_func_type()==3)
        {
            GlobalC::pot.vofk_eff1[ir] = GlobalC::pot.vofk(GlobalV::CURRENT_SPIN, ir);
        }
    }

    //--------------------------------------------
    //(2) check if we need to calculate 
    // pvpR = < phi0 | v(spin) | phiR> for a new spin.
    //--------------------------------------------
    if (GlobalV::CURRENT_SPIN == this->GK->get_spin())
    {
        //GlobalV::ofs_running << " Same spin, same vlocal integration." << std::endl;
    }
    else
    {
        //GlobalV::ofs_running << " (spin change)" << std::endl;
        this->GK->reset_spin(GlobalV::CURRENT_SPIN);

        // if you change the place of the following code,
        // rememeber to delete the #include	
        if (GlobalV::VL_IN_H)
        {
            if(XC_Functional::get_func_type()==3)
            {
                Gint_inout inout(GlobalC::pot.vr_eff1, GlobalC::pot.vofk_eff1, 0, Gint_Tools::job_type::vlocal_meta);
                this->uhm->GK.cal_gint(&inout);
            }
            else
            {
                // vlocal = Vh[rho] + Vxc[rho] + Vl(pseudo)
                Gint_inout inout(GlobalC::pot.vr_eff1, 0, Gint_Tools::job_type::vlocal);
                this->uhm->GK.cal_gint(&inout);
            }

            // added by zhengdy-soc, for non-collinear case
            // integral 4 times, is there any method to simplify?
            if (GlobalV::NSPIN == 4)
            {
                for (int is = 1;is < 4;is++)
                {
                    for (int ir = 0; ir < GlobalC::rhopw->nrxx; ir++)
                    {
                        GlobalC::pot.vr_eff1[ir] = GlobalC::pot.vr_eff(is, ir);
                        if(XC_Functional::get_func_type()==3)
                        {
                            GlobalC::pot.vofk_eff1[ir] = GlobalC::pot.vofk(is, ir);
                        }
                    }
                    
                    if(XC_Functional::get_func_type()==3)
                    {
                        Gint_inout inout(GlobalC::pot.vr_eff1, GlobalC::pot.vofk_eff1, is, Gint_Tools::job_type::vlocal_meta);
                        this->uhm->GK.cal_gint(&inout);
                    }
                    else
                    {
                        Gint_inout inout(GlobalC::pot.vr_eff1, is, Gint_Tools::job_type::vlocal);
                        this->uhm->GK.cal_gint(&inout);
                    }
                }
            }
        }
    }


    if (!this->uhm->init_s)
    {
        ModuleBase::WARNING_QUIT("Hamilt_Linear::solve_using_cg", "Need init S matrix firstly");
    }

    //--------------------------------------------
    // (3) folding matrix, 
    // and diagonalize the H matrix (T+Vl+Vnl).
    //--------------------------------------------

    // with k points
    ModuleBase::timer::tick("HamiltLCAO", "each_k");
    ModuleBase::timer::tick("HamiltLCAO", "H_k");
    this->uhm->calculate_Hk(ik);

    // Effective potential of DFT+U is added to total Hamiltonian here; Quxin adds on 20201029
    if (INPUT.dft_plus_u)
    {
        std::vector<std::complex<double>> eff_pot(this->lowf->ParaV->nloc);
        GlobalC::dftu.cal_eff_pot_mat_complex(ik, this->non_first_scf, &eff_pot[0]);

        for (int irc = 0; irc < this->lowf->ParaV->nloc; irc++)
            this->LM->Hloc2[irc] += eff_pot[irc];
    }

    ModuleBase::timer::tick("Efficience", "H_k");

    // Peize Lin add at 2020.04.04
    if (GlobalC::restart.info_load.load_H && !GlobalC::restart.info_load.load_H_finish)
    {
        GlobalC::restart.load_disk(this->LM[0], "H", ik);
        GlobalC::restart.info_load.load_H_finish = true;
    }
    if (GlobalC::restart.info_save.save_H)
    {
        GlobalC::restart.save_disk(this->LM[0], "H", ik);
    }

    this->hmatrix_k = this->LM->Hloc2.data();
    this->smatrix_k = this->LM->Sloc2.data();
    return;
}

//multi-k case, do something before loop-k
template<>
void HamiltLCAO<std::complex<double>>::constructHamilt()
{
    const Parallel_Orbitals* pv = this->lowf->ParaV;

    int start_spin = -1;
    GK->reset_spin(start_spin);
    GK->destroy_pvpR();
    GK->allocate_pvpR();

#ifdef __DEEPKS
    if (GlobalV::deepks_scf)
    {
        GlobalC::ld.cal_projected_DM_k(this->loc->dm_k,
            GlobalC::ucell,
            GlobalC::ORB,
            GlobalC::GridD,
            pv->trace_loc_row,
            pv->trace_loc_col,
            GlobalC::kv.nks,
            GlobalC::kv.kvec_d);
        GlobalC::ld.cal_descriptor();
        //calculate dE/dD
        GlobalC::ld.cal_gedm(GlobalC::ucell.nat);

        //calculate H_V_deltaR from saved <alpha(0)|psi(R)>
        GlobalC::ld.add_v_delta_k(GlobalC::ucell,
            GlobalC::ORB,
            GlobalC::GridD,
            pv->trace_loc_row,
            pv->trace_loc_col,
            pv->nnr);
    }
#endif
}

} // namespace hamilt