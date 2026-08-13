#include "esolver_ks_lcaopw.h"
#include "source_pw/module_pwdft/elecond.h"
#include "source_io/module_parameter/input_conv.h"
#include <iostream>

//--------------temporary----------------------------
#include "source_estate/module_charge/symm_rho.h"
#include "source_estate/occupy.h"
#include "source_hamilt/module_ewald/h_ewald_pw.h"
//-----force-------------------
#include "source_pw/module_pwdft/force_pw.h"
//-----stress------------------
#include "source_pw/module_pwdft/stress_pw.h"
//---------------------------------------------------
#include "source_estate/elecstate_pw.h"
#include "source_pw/module_pwdft/hamilt_lcaopw.h"
#include "source_pw/module_pwdft/hamilt_pw.h"
#include "source_hsolver/diago_iter_assist.h"
#include "source_hsolver/hsolver_lcaopw.h"
#include "source_hsolver/kernels/hegvd_op.h"
#include "source_base/kernels/math_kernel_op.h"
#include "source_io/module_parameter/parameter.h"
#include "source_hamilt/module_xc/xc_functional.h"

#include <ATen/kernels/blas.h>
#include <ATen/kernels/lapack.h>
#include <sys/time.h>
#ifdef __LCAO
#include "source_io/module_hs/write_vxc_lip.hpp"
#endif

namespace ModuleESolver
{

    template <typename T>
    ESolver_KS_LIP<T>::ESolver_KS_LIP()
    {
        this->classname = "ESolver_KS_LIP";
        this->basisname = "LIP";
    }
    template <typename T>
    ESolver_KS_LIP<T>::~ESolver_KS_LIP()
    {
        //****************************************************
        // do not add any codes in this deconstructor funcion
        //****************************************************
        delete this->psi_local;
        // delete Hamilt
        if (this->p_hamilt != nullptr)
        {
            delete this->p_hamilt;
            this->p_hamilt = nullptr;
        }
    }

    template <typename T>
    void ESolver_KS_LIP<T>::allocate_hamilt(const UnitCell& ucell)
    {
        this->p_hamilt = new hamilt::HamiltLIP<T>(this->pelec->pot, this->pw_wfc, &this->kv, &this->ppcell, &ucell
#ifdef __EXX
            , *this->exx_lip
#endif
        );
    }

    template <typename T>
    void ESolver_KS_LIP<T>::before_scf(UnitCell& ucell, const int istep)
    {
        ESolver_KS_PW<T>::before_scf(ucell, istep);
        auto* p_psi_init = static_cast<psi::PSIPrepare<T>*>(this->stp.p_psi_init);
        p_psi_init->initialize_lcao_in_pw(this->psi_local, GlobalV::ofs_running);
    }

    template <typename T>
    void ESolver_KS_LIP<T>::before_all_runners(BaseCell& basecell, const Input_para& inp)
    {
        basecell.require_kind(BaseCell::Kind::unit_cell, __FUNCTION__);
        UnitCell& ucell = static_cast<UnitCell&>(basecell);
        ESolver_KS_PW<T>::before_all_runners(basecell, inp);

        // Initialize LIP-specific info_lip_ from general_exx_info_ and input
        this->info_lip_.ccp_type = this->general_exx_info_.ccp_type;
        this->info_lip_.hse_omega = this->general_exx_info_.hse_omega;
        if (!inp.exx_fock_lambda.empty())
        {
            this->info_lip_.lambda = std::stod(inp.exx_fock_lambda[0]);
        }

        auto* p_psi_init = static_cast<psi::PSIPrepare<T>*>(this->stp.p_psi_init);
        delete this->psi_local;
        this->psi_local = new psi::Psi<T>(this->stp.psi_cpu->get_nk(),
                                          p_psi_init->psi_initer->nbands_start(),
                                          this->stp.psi_cpu->get_nbasis(),
                                          this->kv.ngk,
                                          true);
#ifdef __EXX
        if (inp.calculation == "scf" || inp.calculation == "relax"
            || inp.calculation == "cell-relax"
            || inp.calculation == "md") {
            if (this->general_exx_info_.cal_exx)
            {
                XC_Functional::set_xc_first_loop(ucell);
                this->exx_lip = std::unique_ptr<Exx_Lip<T>>(new Exx_Lip<T>(this->info_lip_,
                                                                           &this->kv,
                                                                           this->psi_local,
                                                                           this->stp.template get_psi_t<T, base_device::DEVICE_CPU>(),
                                                                           this->pw_wfc,
                                                                           this->pw_rho,
                                                                           &ucell,
                                                                           this->pelec));
            }
}
#endif
    }

    template <typename T>
    void ESolver_KS_LIP<T>::iter_init(UnitCell& ucell, const int istep, const int iter)
    {
        ESolver_KS_PW<T>::iter_init(ucell, istep, iter);
#ifdef __EXX
        if (this->general_exx_info_.cal_exx && !this->general_exx_info_.separate_loop && this->two_level_step) {
            this->exx_lip->cal_exx();
}
#endif
    }

    template <typename T>
    void ESolver_KS_LIP<T>::hamilt2rho_single(UnitCell& ucell, const int istep, const int iter, const double ethr)
    {
        ModuleBase::TITLE("ESolver_KS_LIP", "hamilt2rho_single");
        ModuleBase::timer::start("ESolver_KS_LIP", "hamilt2rho_single");

        // reset energy
        this->pelec->f_en.eband = 0.0;
        this->pelec->f_en.demet = 0.0;
        // choose if psi should be diag in subspace
        // be careful that istep start from 0 and iter start from 1
        // if (iter == 1)
        hsolver::DiagoIterAssist<T>::need_subspace = ((istep == 0 || istep == 1) && iter == 1) ? false : true;
        hsolver::DiagoIterAssist<T>::SCF_ITER = iter;
        hsolver::DiagoIterAssist<T>::PW_DIAG_THR = ethr;
        hsolver::DiagoIterAssist<T>::PW_DIAG_NMAX = this->inp_->pw_diag_nmax;
        bool skip_charge = this->inp_->calculation == "nscf" ? true : false;

        hsolver::HSolverLIP<T> hsolver_lip_obj(this->pw_wfc,
                                               PARAM.globalv.use_uspp,
                                               this->inp_->basis_type,
                                               this->inp_->calculation);
        hsolver_lip_obj.solve(static_cast<hamilt::Hamilt<T>*>(this->p_hamilt), *this->stp.template get_psi_t<T, base_device::DEVICE_CPU>(), this->pelec,
          *this->psi_local, skip_charge,ucell.tpiba,ucell.nat, this->general_exx_info_);

        // add exx
#ifdef __EXX
        bool cal_exx = this->general_exx_info_.cal_exx;
        double hybrid_alpha = this->general_exx_info_.hybrid_alpha;
        if (cal_exx)
        {
            this->pelec->set_exx(this->exx_lip->get_exx_energy(), cal_exx, hybrid_alpha); // Peize Lin add 2019-03-09
        }
#endif

        Symmetry_rho::symmetrize_rho(this->inp_->nspin, this->chr, this->pw_rhod, ucell.symm);

        // deband is calculated from "output" charge density calculated
        // in sum_band
        // need 'rho(out)' and 'vr (v_h(in) and v_xc(in))'
        this->pelec->f_en.deband = this->pelec->cal_delta_eband(ucell);

        ModuleBase::timer::end("ESolver_KS_LIP", "hamilt2rho_single");
    }

    template <typename T>
    void ESolver_KS_LIP<T>::iter_finish(UnitCell& ucell, const int istep, int& iter, bool& conv_esolver)
    {
        ESolver_KS_PW<T>::iter_finish(ucell, istep, iter, conv_esolver);

#ifdef __EXX
        if (this->general_exx_info_.cal_exx && conv_esolver)
        {
            // no separate_loop case
            if (!this->general_exx_info_.separate_loop)
            {
                this->general_exx_info_.hybrid_step = 1;

                // in no_separate_loop case, scf loop only did twice
                // in first scf loop, exx updated once in beginning,
                // in second scf loop, exx updated every iter

                if (!this->two_level_step)
                {
                    // update exx and redo scf
                    XC_Functional::set_xc_type(ucell.atoms[0].ncpp.xc_func);
                    iter = 0;
                    std::cout << " Entering 2nd SCF, where EXX is updated" << std::endl;
                    this->two_level_step++;
                    conv_esolver = false;
                }
            }
            // has separate_loop case
            // exx converged or get max exx steps
            else if (this->two_level_step == this->general_exx_info_.hybrid_step
                     || (iter == 1 && this->two_level_step != 0))
            {
                conv_esolver = true;
            }
            else
            {
                // update exx and redo scf
                if (this->two_level_step == 0)
                {
                    XC_Functional::set_xc_type(ucell.atoms[0].ncpp.xc_func);
                }

                std::cout << " Updating EXX " << std::flush;
                timeval t_start;
                gettimeofday(&t_start, nullptr);

                this->exx_lip->cal_exx();
                iter = 0;
                this->two_level_step++;

                timeval t_end;
                gettimeofday(&t_end, nullptr);
                std::cout << "and rerun SCF\t" << std::setprecision(3) << std::setiosflags(std::ios::scientific)
                          << (double)(t_end.tv_sec - t_start.tv_sec)
                                 + (double)(t_end.tv_usec - t_start.tv_usec) / 1000000.0
                          << std::defaultfloat << " (s)" << std::endl;
                conv_esolver = false;
            }
        }
#endif
    }

    template <typename T>
    void ESolver_KS_LIP<T>::after_all_runners(BaseCell& basecell)
    {
        basecell.require_kind(BaseCell::Kind::unit_cell, __FUNCTION__);
        UnitCell& ucell = static_cast<UnitCell&>(basecell);
        ESolver_KS_PW<T>::after_all_runners(basecell);

#ifdef __LCAO
        if (this->inp_->out_mat_xc)
        {
#ifdef __EXX
            bool cal_exx = this->general_exx_info_.cal_exx;
            double hybrid_alpha = this->general_exx_info_.hybrid_alpha;
#else
            bool cal_exx = false;
            double hybrid_alpha = 0.0;
#endif
            ModuleIO::write_Vxc(this->inp_->nspin,
                                PARAM.globalv.nlocal,
                                GlobalV::DRANK,
                                *this->stp.template get_psi_t<T, base_device::DEVICE_CPU>(),
                                ucell,
                                this->sf,
                                this->solvent,
                                *this->pw_wfc,
                                *this->pw_rho,
                                *this->pw_rhod,
                                this->locpp.vloc,
                                this->chr,
                                this->kv,
                                this->pelec->wg,
                                cal_exx,
                                hybrid_alpha
#ifdef __EXX
                                ,
                                *this->exx_lip
#endif
            );
        }
#endif
    }
    template class ESolver_KS_LIP<std::complex<float>>;
    template class ESolver_KS_LIP<std::complex<double>>;
    // LIP is not supported on GPU yet.
} // namespace ModuleESolver
