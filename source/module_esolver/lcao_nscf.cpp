#include "module_elecstate/module_charge/symmetry_rho.h"
#include "module_esolver/esolver_ks_lcao.h"
#include "module_hamilt_lcao/hamilt_lcaodft/hamilt_lcao.h"
#include "module_hamilt_lcao/module_dftu/dftu.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
//
#include "module_base/timer.h"
#include "module_cell/module_neighbor/sltk_atom_arrange.h"
#include "module_cell/module_neighbor/sltk_grid_driver.h"
#include "module_io/berryphase.h"
#include "module_io/cube_io.h"
#include "module_io/get_pchg_lcao.h"
#include "module_io/get_wf_lcao.h"
#include "module_io/to_wannier90_lcao.h"
#include "module_io/to_wannier90_lcao_in_pw.h"
#include "module_io/write_HS_R.h"
#include "module_io/write_elecstat_pot.h"
#include "module_parameter/parameter.h"
#ifdef __DEEPKS
#include "module_hamilt_lcao/module_deepks/LCAO_deepks.h"
#endif
#include "module_base/formatter.h"
#include "module_elecstate/elecstate_lcao.h"
#include "module_elecstate/module_dm/cal_dm_psi.h"
#include "module_hamilt_general/module_ewald/H_Ewald_pw.h"
#include "module_hamilt_general/module_vdw/vdw.h"
#include "module_hamilt_lcao/hamilt_lcaodft/LCAO_domain.h"
#include "module_hamilt_lcao/hamilt_lcaodft/operator_lcao/op_exx_lcao.h"
#include "module_hamilt_lcao/hamilt_lcaodft/operator_lcao/operator_lcao.h"
#include "module_hamilt_lcao/module_deltaspin/spin_constrain.h"
#include "module_io/read_wfc_nao.h"
#include "module_io/write_elecstat_pot.h"
#include "module_io/write_wfc_nao.h"
#ifdef __EXX
#include "module_io/restart_exx_csr.h"
#endif

namespace ModuleESolver
{

template <typename TK, typename TR>
void ESolver_KS_LCAO<TK, TR>::nscf() {
    ModuleBase::TITLE("ESolver_KS_LCAO", "nscf");

    std::cout << " NON-SELF CONSISTENT CALCULATIONS" << std::endl;

    time_t time_start = std::time(nullptr);

    // mohan add 2021-02-09
    // in ions, istep starts from 1,
    // then when the istep is a variable of scf or nscf,
    // istep becomes istep-1, this should be fixed in future
    int istep = 0;
    hsolver::HSolverLCAO<TK> hsolver_lcao_obj(&(this->pv), PARAM.inp.ks_solver);
    hsolver_lcao_obj.solve(this->p_hamilt, this->psi[0], this->pelec, true);

    time_t time_finish = std::time(nullptr);
    ModuleBase::GlobalFunc::OUT_TIME("cal_bands", time_start, time_finish);

    GlobalV::ofs_running << " end of band structure calculation " << std::endl;
    GlobalV::ofs_running << " band eigenvalue in this processor (eV) :" << std::endl;

    const int nspin = PARAM.inp.nspin;
    const int nbands = PARAM.inp.nbands;

    for (int ik = 0; ik < this->kv.get_nks(); ++ik)
    {
        if (nspin == 2)
        {
            if (ik == 0)
            {
                GlobalV::ofs_running << " spin up :" << std::endl;
            }
            if (ik == (this->kv.get_nks() / 2))
            {
                GlobalV::ofs_running << " spin down :" << std::endl;
            }
        }

        GlobalV::ofs_running << " k-points" << ik + 1 << "(" << this->kv.get_nkstot() << "): " << this->kv.kvec_c[ik].x
                             << " " << this->kv.kvec_c[ik].y << " " << this->kv.kvec_c[ik].z << std::endl;

        for (int ib = 0; ib < nbands; ++ib)
        {
            GlobalV::ofs_running << " spin" << this->kv.isk[ik] + 1 << "final_state " << ib + 1 << " "
                                 << this->pelec->ekb(ik, ib) * ModuleBase::Ry_to_eV << " "
                                 << this->pelec->wg(ik, ib) * this->kv.get_nks() << std::endl;
        }
        GlobalV::ofs_running << std::endl;
    }
    if (PARAM.inp.out_bandgap) {
        std::cout << FmtCore::format("\n * * * * * *\n << Start %s.\n", "band gap calculation");
        if (!PARAM.globalv.two_fermi) {
            this->pelec->cal_bandgap();
            GlobalV::ofs_running << " E_bandgap " << this->pelec->bandgap * ModuleBase::Ry_to_eV << " eV" << std::endl;
        }
        else
        {
            this->pelec->cal_bandgap_updw();
            GlobalV::ofs_running << " E_bandgap_up " << this->pelec->bandgap_up * ModuleBase::Ry_to_eV << " eV"
                                 << std::endl;
            GlobalV::ofs_running << " E_bandgap_dw " << this->pelec->bandgap_dw * ModuleBase::Ry_to_eV << " eV"
                                 << std::endl;
        }
        std::cout << FmtCore::format(" >> Finish %s.\n * * * * * *\n", "band gap calculation");
    }

    // add by jingan in 2018.11.7
    if (PARAM.inp.calculation == "nscf" && PARAM.inp.towannier90)
    {
#ifdef __LCAO
        std::cout << FmtCore::format("\n * * * * * *\n << Start %s.\n", "Wave function to Wannier90");
        if (PARAM.inp.wannier_method == 1) {
            toWannier90_LCAO_IN_PW myWannier(PARAM.inp.out_wannier_mmn,
                                             PARAM.inp.out_wannier_amn,
                                             PARAM.inp.out_wannier_unk,
                                             PARAM.inp.out_wannier_eig,
                                             PARAM.inp.out_wannier_wvfn_formatted,
                                             PARAM.inp.nnkpfile,
                                             PARAM.inp.wannier_spin);

            myWannier.calculate(this->pelec->ekb,
                                this->pw_wfc,
                                this->pw_big,
                                this->sf,
                                this->kv,
                                this->psi,
                                &(this->pv));
        }
        else if (PARAM.inp.wannier_method == 2)
        {
            toWannier90_LCAO myWannier(PARAM.inp.out_wannier_mmn,
                                       PARAM.inp.out_wannier_amn,
                                       PARAM.inp.out_wannier_unk,
                                       PARAM.inp.out_wannier_eig,
                                       PARAM.inp.out_wannier_wvfn_formatted,
                                       PARAM.inp.nnkpfile,
                                       PARAM.inp.wannier_spin,
                                       orb_);

            myWannier.calculate(this->pelec->ekb, this->kv, *(this->psi), &(this->pv));
        }
        std::cout << FmtCore::format(" >> Finish %s.\n * * * * * *\n", "Wave function to Wannier90");
#endif
    }

    // add by jingan
    if (berryphase::berry_phase_flag
        && ModuleSymmetry::Symmetry::symm_flag != 1) {
        std::cout << FmtCore::format("\n * * * * * *\n << Start %s.\n", "Berry phase calculation");
        berryphase bp(&(this->pv));
        bp.lcao_init(this->kv,
                     this->GridT,
                     orb_); // additional step before calling
                                   // macroscopic_polarization (why capitalize
                                   // the function name?)
        bp.Macroscopic_polarization(this->pw_wfc->npwk_max,
                                    this->psi,
                                    this->pw_rho,
                                    this->pw_wfc,
                                    this->kv);
        std::cout << FmtCore::format(" >> Finish %s.\n * * * * * *\n", "Berry phase calculation");
    }

    // below is for DeePKS NSCF calculation
#ifdef __DEEPKS
    if (PARAM.inp.deepks_out_labels || PARAM.inp.deepks_scf) {
        std::cout << FmtCore::format("\n * * * * * *\n << Start %s.\n", "DeepKS output");
        const elecstate::DensityMatrix<TK, double>* dm
            = dynamic_cast<const elecstate::ElecStateLCAO<TK>*>(this->pelec)->get_DM();
        this->dpks_cal_projected_DM(dm);
        GlobalC::ld.cal_descriptor(GlobalC::ucell.nat); // final descriptor
        GlobalC::ld.cal_gedm(GlobalC::ucell.nat);
        std::cout << FmtCore::format(" >> Finish %s.\n * * * * * *\n", "DeepKS output");
    }
#endif

    this->create_Output_Mat_Sparse(0).write();

    // mulliken charge analysis
    if (PARAM.inp.out_mul) {
        std::cout << FmtCore::format("\n * * * * * *\n << Start %s.\n", "Mulliken charge analysis");
        elecstate::ElecStateLCAO<TK>* pelec_lcao
            = dynamic_cast<elecstate::ElecStateLCAO<TK>*>(this->pelec);
        this->pelec->calculate_weights();
        this->pelec->calEBand();
        elecstate::cal_dm_psi(&(this->pv), pelec_lcao->wg, *(this->psi), *(pelec_lcao->get_DM()));
        this->cal_mag(istep, true);
        std::cout << FmtCore::format(" >> Finish %s.\n * * * * * *\n", "Mulliken charge analysis");
    }

    /// write potential
    if (PARAM.inp.out_pot == 1 || PARAM.inp.out_pot == 3)
    {
        for (int is = 0; is < PARAM.inp.nspin; is++)
        {
            std::string fn = PARAM.globalv.global_out_dir + "/SPIN" + std::to_string(is + 1) + "_POT.cube";

            ModuleIO::write_vdata_palgrid(GlobalC::Pgrid,
                this->pelec->pot->get_effective_v(is),
                is,
                PARAM.inp.nspin,
                0,
                fn,
                0.0, // efermi
                &(GlobalC::ucell),
                3,  // precision
                0); // out_fermi
        }
    }
    else if (PARAM.inp.out_pot == 2)
    {
        std::string fn = PARAM.globalv.global_out_dir + "/ElecStaticPot.cube";
        ModuleIO::write_elecstat_pot(
#ifdef __MPI
            this->pw_big->bz,
            this->pw_big->nbz,
#endif
            fn,
            0, // istep
            this->pw_rhod,
            this->pelec->charge,
            &(GlobalC::ucell),
            this->pelec->pot->get_fixed_v());
    }

    // write wfc
    if (PARAM.inp.out_wfc_lcao)
    {
        std::cout << FmtCore::format("\n * * * * * *\n << Start %s.\n", "writing wave function");
        ModuleIO::write_wfc_nao(PARAM.inp.out_wfc_lcao,
                                *this->psi,
                                this->pelec->ekb,
                                this->pelec->wg,
                                this->pelec->klist->kvec_c,
                                this->pv,
                                istep);
        std::cout << FmtCore::format(" >> Finish %s.\n * * * * * *\n", "writing wave function");
    }
}

template class ESolver_KS_LCAO<double, double>;
template class ESolver_KS_LCAO<std::complex<double>, double>;
template class ESolver_KS_LCAO<std::complex<double>, std::complex<double>>;
} // namespace ModuleESolver
