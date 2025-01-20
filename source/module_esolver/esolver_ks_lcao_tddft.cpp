#include "esolver_ks_lcao_tddft.h"

#include "module_io/cal_r_overlap_R.h"
#include "module_io/dipole_io.h"
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
#include "module_elecstate/module_dm/cal_dm_psi.h"
#include "module_elecstate/module_dm/cal_edm_tddft.h"
#include "module_elecstate/module_dm/density_matrix.h"
#include "module_elecstate/occupy.h"
#include "module_hamilt_lcao/module_tddft/evolve_elec.h"
#include "module_hamilt_lcao/module_tddft/td_velocity.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_io/print_info.h"

//-----HSolver ElecState Hamilt--------
#include "module_elecstate/cal_ux.h"
#include "module_elecstate/elecstate_lcao.h"
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

void ESolver_KS_LCAO_TDDFT::before_all_runners(UnitCell& ucell, const Input_para& inp)
{
    // 1) run before_all_runners in ESolver_KS_LCAO
    ESolver_KS_LCAO<std::complex<double>, double>::before_all_runners(ucell, inp);

    // this line should be optimized
    // this->pelec = dynamic_cast<elecstate::ElecStateLCAO_TDDFT*>(this->pelec);
}

void ESolver_KS_LCAO_TDDFT::hamilt2density_single(UnitCell& ucell, const int istep, const int iter, const double ethr)
{
    if (PARAM.inp.init_wfc == "file")
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
                                                 this->pelec->ekb,
                                                 td_htype,
                                                 PARAM.inp.propagator,
                                                 kv.get_nks());
            this->weight_dm_rho();
        }
        this->weight_dm_rho();
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
                                             this->pelec->ekb,
                                             td_htype,
                                             PARAM.inp.propagator,
                                             kv.get_nks());
        this->weight_dm_rho();
    }
    else
    {
        // reset energy
        this->pelec->f_en.eband = 0.0;
        this->pelec->f_en.demet = 0.0;
        if (this->psi != nullptr)
        {
            bool skip_charge = PARAM.inp.calculation == "nscf" ? true : false;
            hsolver::HSolverLCAO<std::complex<double>> hsolver_lcao_obj(&this->pv, PARAM.inp.ks_solver);
            hsolver_lcao_obj.solve(this->p_hamilt, this->psi[0], this->pelec, skip_charge);
        }
    }

    // symmetrize the charge density only for ground state
    if (istep <= 1)
    {
        Symmetry_rho srho;
        for (int is = 0; is < PARAM.inp.nspin; is++)
        {
            srho.begin(is, *(pelec->charge), pw_rho, ucell.symm);
        }
    }

    // (7) calculate delta energy
    this->pelec->f_en.deband = this->pelec->cal_delta_eband(ucell);
}

void ESolver_KS_LCAO_TDDFT::iter_finish(UnitCell& ucell, const int istep, int& iter)
{
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
                GlobalV::ofs_running << ik + 1 << "     " << ib + 1 << "      " << this->pelec->wg(ik, ib) << std::endl;
            }
        }
        GlobalV::ofs_running << std::endl;
        GlobalV::ofs_running << "---------------------------------------------------------------"
                                "---------------------------------"
                             << std::endl;
    }

    ESolver_KS_LCAO<std::complex<double>, double>::iter_finish(ucell, istep, iter);
}

void ESolver_KS_LCAO_TDDFT::update_pot(UnitCell& ucell, const int istep, const int iter)
{
    // Calculate new potential according to new Charge Density
    if (!this->conv_esolver)
    {
        elecstate::cal_ux(ucell);
        this->pelec->pot->update_from_charge(this->pelec->charge, &ucell);
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
    if (istep >= (PARAM.inp.init_wfc == "file" ? 0 : 1) && this->conv_esolver)
    {
        if (this->psi_laststep == nullptr)
        {
#ifdef __MPI
            this->psi_laststep = new psi::Psi<std::complex<double>>(kv.get_nks(), ncol_nbands, nrow, kv.ngk, true);
#else
            this->psi_laststep = new psi::Psi<std::complex<double>>(kv.get_nks(), nbands, nlocal, kv.ngk, true);
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
        if (istep >= (PARAM.inp.init_wfc == "file" ? 0 : 2) && module_tddft::Evolve_elec::td_edm == 0)
        {
            elecstate::cal_edm_tddft(this->pv, this->pelec, this->kv, this->p_hamilt);
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
                                     << this->pelec->ekb(ik, ib) * ModuleBase::Ry_to_eV << std::endl;
            }
        }
        GlobalV::ofs_running << std::endl;
        GlobalV::ofs_running << "---------------------------------------------------------------"
                                "---------------------------------"
                             << std::endl;
    }
}

void ESolver_KS_LCAO_TDDFT::after_scf(UnitCell& ucell, const int istep)
{
    ModuleBase::TITLE("ESolver_KS_LCAO_TDDFT", "after_scf");
    ModuleBase::timer::tick("ESolver_KS_LCAO_TDDFT", "after_scf");

    for (int is = 0; is < PARAM.inp.nspin; is++)
    {
        if (module_tddft::Evolve_elec::out_dipole == 1)
        {
            std::stringstream ss_dipole;
            ss_dipole << PARAM.globalv.global_out_dir << "SPIN" << is + 1 << "_DIPOLE";
            ModuleIO::write_dipole(ucell,
                                   pelec->charge->rho_save[is],
                                   pelec->charge->rhopw,
                                   is,
                                   istep,
                                   ss_dipole.str());
        }
    }
    if (TD_Velocity::out_current == true)
    {
        elecstate::DensityMatrix<std::complex<double>, double>* tmp_DM
            = dynamic_cast<elecstate::ElecStateLCAO<std::complex<double>>*>(this->pelec)->get_DM();

        ModuleIO::write_current(ucell,
                                this->gd,
                                istep,
                                this->psi,
                                pelec,
                                kv,
                                two_center_bundle_.overlap_orb.get(),
                                tmp_DM->get_paraV_pointer(),
                                orb_,
                                this->RA);
    }
    ESolver_KS_LCAO<std::complex<double>, double>::after_scf(ucell, istep);

    ModuleBase::timer::tick("ESolver_KS_LCAO_TDDFT", "after_scf");
}

void ESolver_KS_LCAO_TDDFT::weight_dm_rho()
{
    if (PARAM.inp.ocp == 1)
    {
        this->pelec->fixed_weights(PARAM.inp.ocp_kb, PARAM.inp.nbands, PARAM.inp.nelec);
    }
    this->pelec->calEBand();

    ModuleBase::GlobalFunc::NOTE("Calculate the density matrix.");

    auto _pes = dynamic_cast<elecstate::ElecStateLCAO<std::complex<double>>*>(this->pelec);
    elecstate::cal_dm_psi(_pes->DM->get_paraV_pointer(), _pes->wg, this->psi[0], *(_pes->DM));
    _pes->DM->cal_DMR();

    this->pelec->psiToRho(this->psi[0]);
}

} // namespace ModuleESolver
