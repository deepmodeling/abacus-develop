#include "esolver_ks_lcao_tddft.h"

#include "module_io/cal_r_overlap_R.h"
#include "module_io/dipole_io.h"
#include "module_io/rho_io.h"
#include "module_io/td_current_io.h"
#include "module_io/write_HS.h"
#include "module_io/write_HS_R.h"
#include "module_io/write_wfc_nao.h"

//--------------temporary----------------------------
#include "module_base/blas_connector.h"
#include "module_base/global_function.h"
#include "module_base/lapack_connector.h"
#include "module_base/scalapack_connector.h"
#include "module_elecstate/module_charge/symmetry_rho.h"
#include "module_elecstate/occupy.h"
#include "module_hamilt_lcao/hamilt_lcaodft/LCAO_domain.h" // need divide_HS_in_frag
#include "module_hamilt_lcao/module_tddft/evolve_elec.h"
#include "module_hamilt_lcao/module_tddft/td_velocity.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_io/print_info.h"

//-----HSolver ElecState Hamilt--------
#include "module_elecstate/elecstate_lcao.h"
#include "module_elecstate/elecstate_lcao_tddft.h"
#include "module_hamilt_lcao/hamilt_lcaodft/hamilt_lcao.h"
#include "module_hsolver/hsolver_lcao.h"
#include "module_parameter/parameter.h"
#include "module_psi/psi.h"

//-----force& stress-------------------
#include "module_hamilt_lcao/hamilt_lcaodft/FORCE_STRESS.h"

//---------------------------------------------------

namespace ModuleESolver
{

ESolver_KS_LCAO_TDDFT::ESolver_KS_LCAO_TDDFT()
{
    classname = "ESolver_KS_LCAO_TDDFT";
    basisname = "LCAO";
}

ESolver_KS_LCAO_TDDFT::~ESolver_KS_LCAO_TDDFT()
{
    delete psi_laststep;
    if (Hk_laststep != nullptr)
    {
        for (int ik = 0; ik < kv.get_nks(); ++ik)
        {
            delete[] Hk_laststep[ik];
        }
        delete[] Hk_laststep;
    }
    if (Sk_laststep != nullptr)
    {
        for (int ik = 0; ik < kv.get_nks(); ++ik)
        {
            delete[] Sk_laststep[ik];
        }
        delete[] Sk_laststep;
    }
}

void ESolver_KS_LCAO_TDDFT::before_all_runners(const Input_para& inp, UnitCell& ucell)
{
    // 1) run "before_all_runners" in ESolver_KS
    ESolver_KS::before_all_runners(inp, ucell);

    // 2) initialize the local pseudopotential with plane wave basis set
    GlobalC::ppcell.init_vloc(GlobalC::ppcell.vloc, pw_rho);

    // 3) initialize the electronic states for TDDFT
    if (this->pelec == nullptr)
    {
        this->pelec = new elecstate::ElecStateLCAO_TDDFT(&this->chr,
                                                         &kv,
                                                         kv.get_nks(),
                                                         &this->GK, // mohan add 2024-04-01
                                                         this->pw_rho,
                                                         pw_big);
    }

    // 4) read the local orbitals and construct the interpolation tables.
    // initialize the pv
    LCAO_domain::init_basis_lcao(this->pv,
                                 inp.onsite_radius,
                                 inp.lcao_ecut,
                                 inp.lcao_dk,
                                 inp.lcao_dr,
                                 inp.lcao_rmax,
                                 ucell,
                                 two_center_bundle_,
                                 orb_);

    // 5) allocate H and S matrices according to computational resources
    LCAO_domain::divide_HS_in_frag(PARAM.globalv.gamma_only_local, this->pv, kv.get_nks(), orb_);

    // 6) initialize Density Matrix
    dynamic_cast<elecstate::ElecStateLCAO<std::complex<double>>*>(this->pelec)
        ->init_DM(&kv, &this->pv, PARAM.inp.nspin);

    // 8) initialize the charge density
    this->pelec->charge->allocate(PARAM.inp.nspin);
    this->pelec->omega = GlobalC::ucell.omega; // this line is very odd.

    // 9) initializee the potential
    this->pelec->pot = new elecstate::Potential(pw_rhod,
                                                pw_rho,
                                                &GlobalC::ucell,
                                                &(GlobalC::ppcell.vloc),
                                                &(sf),
                                                &(pelec->f_en.etxc),
                                                &(pelec->f_en.vtxc));

    // this line should be optimized
    this->pelec_td = dynamic_cast<elecstate::ElecStateLCAO_TDDFT*>(this->pelec);
}

void ESolver_KS_LCAO_TDDFT::hamilt2density(const int istep, const int iter, const double ethr)
{
    pelec->charge->save_rho_before_sum_band();

    if (wf.init_wfc == "file")
    {
        if (istep >= 1)
        {
            module_tddft::Evolve_elec::solve_psi(istep,
                                                 PARAM.inp.nbands,
                                                 PARAM.globalv.nlocal,
                                                 this->p_hamilt,
                                                 this->pv,
                                                 this->psi,
                                                 this->psi_laststep,
                                                 this->Hk_laststep,
                                                 this->Sk_laststep,
                                                 this->pelec_td->ekb,
                                                 td_htype,
                                                 PARAM.inp.propagator,
                                                 kv.get_nks());
            this->pelec_td->psiToRho_td(this->psi[0]);
        }
        this->pelec_td->psiToRho_td(this->psi[0]);
    }
    else if (istep >= 2)
    {
        module_tddft::Evolve_elec::solve_psi(istep,
                                             PARAM.inp.nbands,
                                             PARAM.globalv.nlocal,
                                             this->p_hamilt,
                                             this->pv,
                                             this->psi,
                                             this->psi_laststep,
                                             this->Hk_laststep,
                                             this->Sk_laststep,
                                             this->pelec_td->ekb,
                                             td_htype,
                                             PARAM.inp.propagator,
                                             kv.get_nks());
        this->pelec_td->psiToRho_td(this->psi[0]);
    }
    else
    {
        // reset energy
        this->pelec->f_en.eband = 0.0;
        this->pelec->f_en.demet = 0.0;
        if (this->psi != nullptr)
        {
            hsolver::HSolverLCAO<std::complex<double>> hsolver_lcao_obj(&this->pv, PARAM.inp.ks_solver);
            hsolver_lcao_obj.solve(this->p_hamilt, this->psi[0], this->pelec_td, false);
        }
    }
    // else
    // {
    //     ModuleBase::WARNING_QUIT("ESolver_KS_LCAO", "HSolver has not been initialed!");
    // }

    // print occupation of each band
    if (iter == 1 && istep <= 2)
    {
        GlobalV::ofs_running << "---------------------------------------------------------------"
                                "---------------------------------"
                             << std::endl;
        GlobalV::ofs_running << "occupation : " << std::endl;
        GlobalV::ofs_running << "ik  iband     occ " << std::endl;
        GlobalV::ofs_running << std::setprecision(6);
        GlobalV::ofs_running << std::setiosflags(std::ios::showpoint);
        for (int ik = 0; ik < kv.get_nks(); ik++)
        {
            for (int ib = 0; ib < PARAM.inp.nbands; ib++)
            {
                std::setprecision(6);
                GlobalV::ofs_running << ik + 1 << "     " << ib + 1 << "      " << this->pelec_td->wg(ik, ib)
                                     << std::endl;
            }
        }
        GlobalV::ofs_running << std::endl;
        GlobalV::ofs_running << "---------------------------------------------------------------"
                                "---------------------------------"
                             << std::endl;
    }

    for (int ik = 0; ik < kv.get_nks(); ++ik)
    {
        this->pelec_td->print_band(ik, PARAM.inp.printe, iter);
    }

    // using new charge density.
    this->pelec->cal_energies(1);

    // symmetrize the charge density only for ground state
    if (istep <= 1)
    {
        Symmetry_rho srho;
        for (int is = 0; is < PARAM.inp.nspin; is++)
        {
            srho.begin(is, *(pelec->charge), pw_rho, GlobalC::ucell.symm);
        }
    }

    // (6) compute magnetization, only for spin==2
    GlobalC::ucell.magnet.compute_magnetization(this->pelec->charge->nrxx,
                                                this->pelec->charge->nxyz,
                                                this->pelec->charge->rho,
                                                pelec->nelec_spin.data());

    // (7) calculate delta energy
    this->pelec->f_en.deband = this->pelec->cal_delta_eband();
}

void ESolver_KS_LCAO_TDDFT::update_pot(const int istep, const int iter)
{
    // print Hamiltonian and Overlap matrix
    if (this->conv_esolver)
    {
        if (!PARAM.globalv.gamma_only_local)
        {
            this->GK.renew(true);
        }
        for (int ik = 0; ik < kv.get_nks(); ++ik)
        {
            if (PARAM.inp.out_mat_hs[0])
            {
                this->p_hamilt->updateHk(ik);
            }
            bool bit = false; // LiuXh, 2017-03-21
            // if set bit = true, there would be error in soc-multi-core
            // calculation, noted by zhengdy-soc
            if (this->psi != nullptr && (istep % PARAM.inp.out_interval == 0))
            {
                hamilt::MatrixBlock<complex<double>> h_mat, s_mat;
                this->p_hamilt->matrix(h_mat, s_mat);
                if (PARAM.inp.out_mat_hs[0])
                {
                    ModuleIO::save_mat(istep,
                                       h_mat.p,
                                       PARAM.globalv.nlocal,
                                       bit,
                                       PARAM.inp.out_mat_hs[1],
                                       1,
                                       PARAM.inp.out_app_flag,
                                       "H",
                                       "data-" + std::to_string(ik),
                                       this->pv,
                                       GlobalV::DRANK);

                    ModuleIO::save_mat(istep,
                                       s_mat.p,
                                       PARAM.globalv.nlocal,
                                       bit,
                                       PARAM.inp.out_mat_hs[1],
                                       1,
                                       PARAM.inp.out_app_flag,
                                       "S",
                                       "data-" + std::to_string(ik),
                                       this->pv,
                                       GlobalV::DRANK);
                }
            }
        }
    }

    if (elecstate::ElecStateLCAO<std::complex<double>>::out_wfc_lcao
        && (this->conv_esolver || iter == PARAM.inp.scf_nmax) && (istep % PARAM.inp.out_interval == 0))
    {
        ModuleIO::write_wfc_nao(elecstate::ElecStateLCAO<std::complex<double>>::out_wfc_lcao,
                                this->psi[0],
                                this->pelec->ekb,
                                this->pelec->wg,
                                this->pelec->klist->kvec_c,
                                this->pv,
                                istep);
    }

    // Calculate new potential according to new Charge Density
    if (!this->conv_esolver)
    {
        if (PARAM.inp.nspin == 4)
        {
            GlobalC::ucell.cal_ux();
        }
        this->pelec->pot->update_from_charge(this->pelec->charge, &GlobalC::ucell);
        this->pelec->f_en.descf = this->pelec->cal_delta_escf();
    }
    else
    {
        this->pelec->cal_converged();
    }

    const int nloc = this->pv.nloc;
    const int ncol_nbands = this->pv.ncol_bands;
    const int nrow = this->pv.nrow;
    const int nbands = PARAM.inp.nbands;
    const int nlocal = PARAM.globalv.nlocal;

    // store wfc and Hk laststep
    if (istep >= (wf.init_wfc == "file" ? 0 : 1) && this->conv_esolver)
    {
        if (this->psi_laststep == nullptr)
        {
#ifdef __MPI
            this->psi_laststep = new psi::Psi<std::complex<double>>(kv.get_nks(), ncol_nbands, nrow, nullptr);
#else
            this->psi_laststep = new psi::Psi<std::complex<double>>(kv.get_nks(), nbands, nlocal, nullptr);
#endif
        }

        if (td_htype == 1)
        {
            if (this->Hk_laststep == nullptr)
            {
                this->Hk_laststep = new std::complex<double>*[kv.get_nks()];
                for (int ik = 0; ik < kv.get_nks(); ++ik)
                {
                    this->Hk_laststep[ik] = new std::complex<double>[nloc];
                    ModuleBase::GlobalFunc::ZEROS(Hk_laststep[ik], nloc);
                }
            }
            if (this->Sk_laststep == nullptr)
            {
                this->Sk_laststep = new std::complex<double>*[kv.get_nks()];
                for (int ik = 0; ik < kv.get_nks(); ++ik)
                {
                    this->Sk_laststep[ik] = new std::complex<double>[nloc];
                    ModuleBase::GlobalFunc::ZEROS(Sk_laststep[ik], nloc);
                }
            }
        }

        for (int ik = 0; ik < kv.get_nks(); ++ik)
        {
            this->psi->fix_k(ik);
            this->psi_laststep->fix_k(ik);
            int size0 = psi->get_nbands() * psi->get_nbasis();
            for (int index = 0; index < size0; ++index)
            {
                psi_laststep[0].get_pointer()[index] = psi[0].get_pointer()[index];
            }

            // store Hamiltonian
            if (td_htype == 1)
            {
                this->p_hamilt->updateHk(ik);
                hamilt::MatrixBlock<complex<double>> h_mat, s_mat;
                this->p_hamilt->matrix(h_mat, s_mat);
                BlasConnector::copy(nloc, h_mat.p, 1, Hk_laststep[ik], 1);
                BlasConnector::copy(nloc, s_mat.p, 1, Sk_laststep[ik], 1);
            }
        }

        // calculate energy density matrix for tddft
        if (istep >= (wf.init_wfc == "file" ? 0 : 2) && module_tddft::Evolve_elec::td_edm == 0)
        {
            this->cal_edm_tddft();
        }
    }

    // print "eigen value" for tddft
    if (this->conv_esolver)
    {
        GlobalV::ofs_running << "---------------------------------------------------------------"
                                "---------------------------------"
                             << std::endl;
        GlobalV::ofs_running << "Eii : " << std::endl;
        GlobalV::ofs_running << "ik  iband    Eii (eV)" << std::endl;
        GlobalV::ofs_running << std::setprecision(6);
        GlobalV::ofs_running << std::setiosflags(std::ios::showpoint);

        for (int ik = 0; ik < kv.get_nks(); ik++)
        {
            for (int ib = 0; ib < PARAM.inp.nbands; ib++)
            {
                GlobalV::ofs_running << ik + 1 << "     " << ib + 1 << "      "
                                     << this->pelec_td->ekb(ik, ib) * ModuleBase::Ry_to_eV << std::endl;
            }
        }
        GlobalV::ofs_running << std::endl;
        GlobalV::ofs_running << "---------------------------------------------------------------"
                                "---------------------------------"
                             << std::endl;
    }
}

void ESolver_KS_LCAO_TDDFT::after_scf(const int istep)
{
    for (int is = 0; is < PARAM.inp.nspin; is++)
    {
        if (module_tddft::Evolve_elec::out_dipole == 1)
        {
            std::stringstream ss_dipole;
            ss_dipole << PARAM.globalv.global_out_dir << "SPIN" << is + 1 << "_DIPOLE";
            ModuleIO::write_dipole(pelec->charge->rho_save[is], pelec->charge->rhopw, is, istep, ss_dipole.str());
        }
    }
    if (TD_Velocity::out_current == true)
    {
        elecstate::DensityMatrix<std::complex<double>, double>* tmp_DM
            = dynamic_cast<elecstate::ElecStateLCAO<std::complex<double>>*>(this->pelec)->get_DM();

        ModuleIO::write_current(istep,
                                this->psi,
                                pelec,
                                kv,
                                two_center_bundle_.overlap_orb.get(),
                                tmp_DM->get_paraV_pointer(),
                                orb_,
                                this->RA);
    }
    ESolver_KS_LCAO<std::complex<double>, double>::after_scf(istep);
}

} // namespace ModuleESolver
