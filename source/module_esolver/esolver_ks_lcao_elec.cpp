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
#include "module_io/istate_charge.h"
#include "module_io/istate_envelope.h"
#include "module_io/to_wannier90_lcao.h"
#include "module_io/to_wannier90_lcao_in_pw.h"
#include "module_io/write_HS_R.h"
#ifdef __DEEPKS
#include "module_hamilt_lcao/module_deepks/LCAO_deepks.h"
#endif
#include "module_elecstate/elecstate_lcao.h"
#include "module_hamilt_general/module_ewald/H_Ewald_pw.h"
#include "module_hamilt_general/module_vdw/vdw.h"
#include "module_hamilt_lcao/hamilt_lcaodft/operator_lcao/op_exx_lcao.h"
#include "module_io/dm_io.h"

#include "module_hamilt_lcao/module_deltaspin/spin_constrain.h"

namespace ModuleESolver
{

template <typename TK, typename TR>
void ESolver_KS_LCAO<TK, TR>::set_matrix_grid(Record_adj& ra)
{
    ModuleBase::TITLE("ESolver_KS_LCAO", "set_matrix_grid");
    ModuleBase::timer::tick("ESolver_KS_LCAO", "set_matrix_grid");

    // (1) Find adjacent atoms for each atom.
    GlobalV::SEARCH_RADIUS = atom_arrange::set_sr_NL(GlobalV::ofs_running,
                                                     GlobalV::OUT_LEVEL,
                                                     GlobalC::ORB.get_rcutmax_Phi(),
                                                     GlobalC::ucell.infoNL.get_rcutmax_Beta(),
                                                     GlobalV::GAMMA_ONLY_LOCAL);

    atom_arrange::search(GlobalV::SEARCH_PBC,
                         GlobalV::ofs_running,
                         GlobalC::GridD,
                         GlobalC::ucell,
                         GlobalV::SEARCH_RADIUS,
                         GlobalV::test_atom_input);

    // ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running,"SEARCH ADJACENT ATOMS");

    // (3) Periodic condition search for each grid.
    this->GridT.set_pbc_grid(this->pw_rho->nx,
                             this->pw_rho->ny,
                             this->pw_rho->nz,
                             this->pw_big->bx,
                             this->pw_big->by,
                             this->pw_big->bz,
                             this->pw_big->nbx,
                             this->pw_big->nby,
                             this->pw_big->nbz,
                             this->pw_big->nbxx,
                             this->pw_big->nbzp_start,
                             this->pw_big->nbzp,
                             this->pw_rho->ny,
                             this->pw_rho->nplane,
                             this->pw_rho->startz_current);

    // (2)For each atom, calculate the adjacent atoms in different cells
    // and allocate the space for H(R) and S(R).
    // If k point is used here, allocate HlocR after atom_arrange.
    Parallel_Orbitals* pv = this->UHM.LM->ParaV;
    ra.for_2d(*pv, GlobalV::GAMMA_ONLY_LOCAL);
    if (!GlobalV::GAMMA_ONLY_LOCAL)
    {
        // need to first calculae lgd.
        // using GridT.init.
        this->GridT.cal_nnrg(pv);
    }

    ModuleBase::timer::tick("ESolver_KS_LCAO", "set_matrix_grid");
    return;
}

template <typename TK, typename TR>
void ESolver_KS_LCAO<TK, TR>::beforesolver(const int istep)
{
    ModuleBase::TITLE("ESolver_KS_LCAO", "beforesolver");
    ModuleBase::timer::tick("ESolver_KS_LCAO", "beforesolver");

    // 1. prepare HS matrices, prepare grid integral
    this->set_matrix_grid(this->RA);

    // 2. density matrix extrapolation

    // set the augmented orbitals index.
    // after ParaO and GridT,
    // this information is used to calculate
    // the force.

    // init psi
    if (this->psi == nullptr)
    {
        int nsk;
        int ncol;
        if (GlobalV::GAMMA_ONLY_LOCAL)
        {
            nsk = GlobalV::NSPIN;
            ncol = this->LOWF.ParaV->ncol_bands;
            if (GlobalV::KS_SOLVER == "genelpa" || GlobalV::KS_SOLVER == "lapack_gvx"
                || GlobalV::KS_SOLVER == "cusolver")
            {
                ncol = this->LOWF.ParaV->ncol;
            }
        }
        else
        {
            nsk = this->kv.nks;
#ifdef __MPI
            ncol = this->LOWF.ParaV->ncol_bands;
#else
            ncol = GlobalV::NBANDS;
#endif
        }
        this->psi = new psi::Psi<TK>(nsk, ncol, this->LOWF.ParaV->nrow, nullptr);
    }

    // prepare grid in Gint
    this->UHM.grid_prepare(this->GridT, *this->pw_rho, *this->pw_big);

    // init Hamiltonian
    if (this->p_hamilt != nullptr)
    {
        delete this->p_hamilt;
        this->p_hamilt = nullptr;
    }
    if (this->p_hamilt == nullptr)
    {
        elecstate::DensityMatrix<TK, double>* DM = dynamic_cast<elecstate::ElecStateLCAO<TK>*>(this->pelec)->get_DM();
        this->p_hamilt = new hamilt::HamiltLCAO<TK, TR>(GlobalV::GAMMA_ONLY_LOCAL ? &(this->UHM.GG) : nullptr,
                                                        GlobalV::GAMMA_ONLY_LOCAL ? nullptr : &(this->UHM.GK),
                                                        &(this->UHM.genH),
                                                        &(this->LM),
                                                        &(this->LOC),
                                                        this->pelec->pot,
                                                        this->kv,
                                                        DM);
    }
    // init density kernel and wave functions.
    this->LOC.allocate_dm_wfc(this->GridT, this->pelec, this->LOWF, this->psi, this->kv, istep);

    //======================================
    // do the charge extrapolation before the density matrix is regenerated.
    // mohan add 2011-04-08
    // because once atoms are moving out of this processor,
    // the density matrix will not std::map the new atomic configuration,
    //======================================
    // THIS IS A BUG, BECAUSE THE INDEX GlobalC::GridT.trace_lo
    // HAS BEEN REGENERATED, SO WE NEED TO
    // REALLOCATE DENSITY MATRIX FIRST, THEN READ IN DENSITY MATRIX,
    // AND USE DENSITY MATRIX TO DO RHO GlobalV::CALCULATION.-- mohan 2013-03-31
    //======================================
    if (GlobalV::chg_extrap == "dm" && istep > 1) // xiaohui modify 2015-02-01
    {
        for (int is = 0; is < GlobalV::NSPIN; is++)
        {
            ModuleBase::GlobalFunc::ZEROS(this->pelec->charge->rho[is], this->pw_rho->nrxx);
            std::stringstream ssd;
            ssd << GlobalV::global_out_dir << "SPIN" << is + 1 << "_DM";
            // reading density matrix,
            double& ef_tmp = this->pelec->eferm.get_ef(is);
            ModuleIO::read_dm(
#ifdef __MPI
                this->GridT.nnrg,
                this->GridT.trace_lo,
#endif
                is,
                ssd.str(),
                this->LOC.DM,
                this->LOC.DM_R,
                ef_tmp,
                &(GlobalC::ucell));
        }

        // calculate the charge density
        if (GlobalV::GAMMA_ONLY_LOCAL)
        {
            Gint_inout inout(this->LOC.DM, this->pelec->charge->rho, Gint_Tools::job_type::rho);
            this->UHM.GG.cal_gint(&inout);
            if (XC_Functional::get_func_type() == 3 || XC_Functional::get_func_type() == 5)
            {
                for (int is = 0; is < GlobalV::NSPIN; is++)
                {
                    ModuleBase::GlobalFunc::ZEROS(this->pelec->charge->kin_r[0], this->pw_rho->nrxx);
                }
                Gint_inout inout1(this->LOC.DM, this->pelec->charge->kin_r, Gint_Tools::job_type::tau);
                this->UHM.GG.cal_gint(&inout1);
            }
        }
        else
        {
            Gint_inout inout(this->LOC.DM_R, this->pelec->charge->rho, Gint_Tools::job_type::rho);
            this->UHM.GK.cal_gint(&inout);
            if (XC_Functional::get_func_type() == 3 || XC_Functional::get_func_type() == 5)
            {
                for (int is = 0; is < GlobalV::NSPIN; is++)
                {
                    ModuleBase::GlobalFunc::ZEROS(this->pelec->charge->kin_r[0], this->pw_rho->nrxx);
                }
                Gint_inout inout1(this->LOC.DM_R, this->pelec->charge->kin_r, Gint_Tools::job_type::tau);
                this->UHM.GK.cal_gint(&inout1);
            }
        }

        // renormalize the charge density
        this->pelec->charge->renormalize_rho();
    }

#ifdef __DEEPKS
    // for each ionic step, the overlap <psi|alpha> must be rebuilt
    // since it depends on ionic positions
    if (GlobalV::deepks_setorb)
    {
        const Parallel_Orbitals* pv = this->UHM.LM->ParaV;
        // build and save <psi(0)|alpha(R)> at beginning
        GlobalC::ld.build_psialpha(GlobalV::CAL_FORCE, GlobalC::ucell, GlobalC::ORB, GlobalC::GridD, GlobalC::UOT);

        if (GlobalV::deepks_out_unittest)
        {
            GlobalC::ld.check_psialpha(GlobalV::CAL_FORCE, GlobalC::ucell, GlobalC::ORB, GlobalC::GridD, GlobalC::UOT);
        }
    }
#endif
    if (GlobalV::sc_mag_switch == 1)
    {
        SpinConstrain<TK, psi::DEVICE_CPU>& sc = SpinConstrain<TK, psi::DEVICE_CPU>::getScInstance();
        sc.init_sc(GlobalV::sc_thr,
                   GlobalV::nsc,
                   GlobalV::nsc_min,
                   GlobalV::alpha_trial,
                   GlobalV::sccut,
                   GlobalV::decay_grad_switch,
                   GlobalC::ucell,
                   GlobalV::sc_file,
                   GlobalV::NPOL,
                   &(this->orb_con.ParaV),
                   GlobalV::NSPIN,
                   this->kv,
                   GlobalV::KS_SOLVER,
                   &(this->LM),
                   this->phsol,
                   this->p_hamilt,
                   this->psi,
                   this->pelec);
    }
    //=========================================================
    // cal_ux should be called before init_scf because
    // the direction of ux is used in noncoline_rho
    //=========================================================
    if(GlobalV::NSPIN == 4 && GlobalV::DOMAG) GlobalC::ucell.cal_ux();
    ModuleBase::timer::tick("ESolver_KS_LCAO", "beforesolver");
}

template <typename TK, typename TR>
void ESolver_KS_LCAO<TK, TR>::beforescf(int istep)
{
    ModuleBase::TITLE("ESolver_KS_LCAO", "beforescf");
    ModuleBase::timer::tick("ESolver_KS_LCAO", "beforescf");

    if (GlobalC::ucell.cell_parameter_updated)
    {
        this->init_after_vc(INPUT, GlobalC::ucell);
    }
    if (GlobalC::ucell.ionic_position_updated)
    {
        this->CE.update_all_dis(GlobalC::ucell);
        this->CE.extrapolate_charge(
#ifdef __MPI
            &(GlobalC::Pgrid),
#endif
            GlobalC::ucell,
            this->pelec->charge,
            &(this->sf));
    }

    //----------------------------------------------------------
    // about vdw, jiyy add vdwd3 and linpz add vdwd2
    //----------------------------------------------------------
    auto vdw_solver = vdw::make_vdw(GlobalC::ucell, INPUT);
    if (vdw_solver != nullptr)
    {
        this->pelec->f_en.evdw = vdw_solver->get_energy();
    }

    this->beforesolver(istep);
    this->pelec->init_scf(istep, this->sf.strucFac);
    // initalize DMR
    // DMR should be same size with Hamiltonian(R)
    dynamic_cast<elecstate::ElecStateLCAO<TK>*>(this->pelec)
        ->get_DM()
        ->init_DMR(*(dynamic_cast<hamilt::HamiltLCAO<TK, TR>*>(this->p_hamilt)->getHR()));

    // the electron charge density should be symmetrized,
    // here is the initialization
    Symmetry_rho srho;
    for (int is = 0; is < GlobalV::NSPIN; is++)
    {
        srho.begin(is, *(this->pelec->charge), this->pw_rho, GlobalC::Pgrid, GlobalC::ucell.symm);
    }
// Peize Lin add 2016-12-03
#ifdef __EXX
    if (GlobalC::exx_info.info_ri.real_number)
        this->exd->exx_beforescf(this->kv, *this->p_chgmix);
    else
        this->exc->exx_beforescf(this->kv, *this->p_chgmix);
#endif // __EXX
       // 1. calculate ewald energy.
       // mohan update 2021-02-25
    if (!GlobalV::test_skip_ewald)
    {
        this->pelec->f_en.ewald_energy = H_Ewald_pw::compute_ewald(GlobalC::ucell, this->pw_rho, this->sf.strucFac);
    }

    this->p_hamilt->non_first_scf = istep;

    ModuleBase::timer::tick("ESolver_KS_LCAO", "beforescf");
    return;
}

template <typename TK, typename TR>
void ESolver_KS_LCAO<TK, TR>::othercalculation(const int istep)
{
    ModuleBase::TITLE("ESolver_KS_LCAO", "othercalculation");
    ModuleBase::timer::tick("ESolver_KS_LCAO", "othercalculation");
    if (GlobalV::CALCULATION == "get_S")
    {
        this->get_S();
        ModuleBase::QUIT();
    }

    if (GlobalV::CALCULATION == "test_memory")
    {
        Cal_Test::test_memory(this->pw_rho,
                              this->pw_wfc,
                              this->p_chgmix->get_mixing_mode(),
                              this->p_chgmix->get_mixing_ndim());
        return;
    }

    if (GlobalV::CALCULATION == "test_neighbour")
    {
        // test_search_neighbor();
        if (GlobalV::SEARCH_RADIUS < 0)
        {
            std::cout << " SEARCH_RADIUS : " << GlobalV::SEARCH_RADIUS << std::endl;
            std::cout << " please make sure search_radius > 0" << std::endl;
        }

        atom_arrange::search(GlobalV::SEARCH_PBC,
                             GlobalV::ofs_running,
                             GlobalC::GridD,
                             GlobalC::ucell,
                             GlobalV::SEARCH_RADIUS,
                             GlobalV::test_atom_input,
                             1);
        return;
    }

    this->beforesolver(istep);
    // pelec should be initialized before these calculations
    this->pelec->init_scf(istep, this->sf.strucFac);
    // self consistent calculations for electronic ground state
    if (GlobalV::CALCULATION == "nscf")
    {
        this->nscf();
    }
    else if (GlobalV::CALCULATION == "get_pchg")
    {
        IState_Charge ISC(this->psi, this->LOC);
        ISC.begin(this->UHM.GG,
                  this->pelec,
                  this->pw_rho,
                  this->pw_big,
                  GlobalV::GAMMA_ONLY_LOCAL,
                  GlobalV::NBANDS_ISTATE,
                  GlobalV::NBANDS,
                  GlobalV::nelec,
                  GlobalV::NSPIN,
                  GlobalV::global_out_dir,
                  GlobalV::MY_RANK,
                  GlobalV::ofs_warning);
    }
    else if (GlobalV::CALCULATION == "get_wf")
    {
        IState_Envelope IEP(this->pelec);
        if (GlobalV::GAMMA_ONLY_LOCAL)
            IEP.begin(this->psi,
                      this->pw_rho,
                      this->pw_wfc,
                      this->pw_big,
                      this->LOWF,
                      this->UHM.GG,
                      INPUT.out_wfc_pw,
                      this->wf.out_wfc_r,
                      this->kv);
        else
            IEP.begin(this->psi,
                      this->pw_rho,
                      this->pw_wfc,
                      this->pw_big,
                      this->LOWF,
                      this->UHM.GK,
                      INPUT.out_wfc_pw,
                      this->wf.out_wfc_r,
                      this->kv);
    }
    else
    {
        ModuleBase::WARNING_QUIT("ESolver_KS_LCAO<TK, TR>::othercalculation", "CALCULATION type not supported");
    }

    ModuleBase::timer::tick("ESolver_KS_LCAO", "othercalculation");
    return;
}
template <>
void ESolver_KS_LCAO<double, double>::get_S()
{
    ModuleBase::TITLE("ESolver_KS_LCAO", "get_S");
    ModuleBase::WARNING_QUIT("ESolver_KS_LCAO<TK, TR>::get_S", "not implemented for");
}
template <>
void ESolver_KS_LCAO<std::complex<double>, double>::get_S()
{
    ModuleBase::TITLE("ESolver_KS_LCAO", "get_S");
    // (1) Find adjacent atoms for each atom.
    GlobalV::SEARCH_RADIUS = atom_arrange::set_sr_NL(GlobalV::ofs_running,
                                                     GlobalV::OUT_LEVEL,
                                                     GlobalC::ORB.get_rcutmax_Phi(),
                                                     GlobalC::ucell.infoNL.get_rcutmax_Beta(),
                                                     GlobalV::GAMMA_ONLY_LOCAL);

    atom_arrange::search(GlobalV::SEARCH_PBC,
                         GlobalV::ofs_running,
                         GlobalC::GridD,
                         GlobalC::ucell,
                         GlobalV::SEARCH_RADIUS,
                         GlobalV::test_atom_input);

    this->RA.for_2d(this->orb_con.ParaV, GlobalV::GAMMA_ONLY_LOCAL);
    this->UHM.genH.LM->ParaV = &this->orb_con.ParaV;
    if (this->p_hamilt == nullptr)
    {
        this->p_hamilt = new hamilt::HamiltLCAO<std::complex<double>, double>(this->UHM.genH.LM, this->kv);
        dynamic_cast<hamilt::OperatorLCAO<std::complex<double>, double>*>(this->p_hamilt->ops)->contributeHR();
    }
    ModuleIO::output_S_R(this->UHM, this->p_hamilt, "SR.csr");
}
template <>
void ESolver_KS_LCAO<std::complex<double>, std::complex<double>>::get_S()
{
    ModuleBase::TITLE("ESolver_KS_LCAO", "get_S");
    // (1) Find adjacent atoms for each atom.
    GlobalV::SEARCH_RADIUS = atom_arrange::set_sr_NL(GlobalV::ofs_running,
                                                     GlobalV::OUT_LEVEL,
                                                     GlobalC::ORB.get_rcutmax_Phi(),
                                                     GlobalC::ucell.infoNL.get_rcutmax_Beta(),
                                                     GlobalV::GAMMA_ONLY_LOCAL);

    atom_arrange::search(GlobalV::SEARCH_PBC,
                         GlobalV::ofs_running,
                         GlobalC::GridD,
                         GlobalC::ucell,
                         GlobalV::SEARCH_RADIUS,
                         GlobalV::test_atom_input);

    this->RA.for_2d(this->orb_con.ParaV, GlobalV::GAMMA_ONLY_LOCAL);
    this->UHM.genH.LM->ParaV = &this->orb_con.ParaV;
    if (this->p_hamilt == nullptr)
    {
        this->p_hamilt
            = new hamilt::HamiltLCAO<std::complex<double>, std::complex<double>>(this->UHM.genH.LM, this->kv);
        dynamic_cast<hamilt::OperatorLCAO<std::complex<double>, std::complex<double>>*>(this->p_hamilt->ops)
            ->contributeHR();
    }
    ModuleIO::output_S_R(this->UHM, this->p_hamilt, "SR.csr");
}

template <typename TK, typename TR>
void ESolver_KS_LCAO<TK, TR>::nscf()
{
    ModuleBase::TITLE("ESolver_KS_LCAO", "nscf");

    std::cout << " NON-SELF CONSISTENT CALCULATIONS" << std::endl;

    time_t time_start = std::time(NULL);

#ifdef __EXX
#ifdef __MPI
    // Peize Lin add 2018-08-14
    if (GlobalC::exx_info.info_global.cal_exx)
    {
        // GlobalC::exx_lcao.cal_exx_elec_nscf(this->LOWF.ParaV[0]);
        const std::string file_name_exx = GlobalV::global_out_dir + "HexxR_" + std::to_string(GlobalV::MY_RANK);
        if (GlobalC::exx_info.info_ri.real_number)
            this->exd->read_Hexxs(file_name_exx);
        else
            this->exc->read_Hexxs(file_name_exx);

        hamilt::HamiltLCAO<TK, TR>* hamilt_lcao = dynamic_cast<hamilt::HamiltLCAO<TK, TR>*>(this->p_hamilt);
        auto exx = new hamilt::OperatorEXX<hamilt::OperatorLCAO<TK, TR>>(&this->LM,
                                                                         hamilt_lcao->getHR(),
                                                                         &(hamilt_lcao->getHk(&this->LM)),
                                                                         this->kv);
        hamilt_lcao->getOperator()->add(exx);
    }
#endif // __MPI
#endif // __EXX

    // mohan add 2021-02-09
    // in ions, istep starts from 1,
    // then when the istep is a variable of scf or nscf,
    // istep becomes istep-1, this should be fixed in future
    int istep = 0;
    if (this->phsol != nullptr)
    {
        this->phsol->solve(this->p_hamilt, this->psi[0], this->pelec, GlobalV::KS_SOLVER, true);
    }
    else
    {
        ModuleBase::WARNING_QUIT("ESolver_KS_PW", "HSolver has not been initialed!");
    }

    time_t time_finish = std::time(NULL);
    ModuleBase::GlobalFunc::OUT_TIME("cal_bands", time_start, time_finish);

    GlobalV::ofs_running << " end of band structure calculation " << std::endl;
    GlobalV::ofs_running << " band eigenvalue in this processor (eV) :" << std::endl;

    for (int ik = 0; ik < this->kv.nks; ik++)
    {
        if (GlobalV::NSPIN == 2)
        {
            if (ik == 0)
            {
                GlobalV::ofs_running << " spin up :" << std::endl;
            }
            if (ik == (this->kv.nks / 2))
            {
                GlobalV::ofs_running << " spin down :" << std::endl;
            }
        }

        GlobalV::ofs_running << " k-points" << ik + 1 << "(" << this->kv.nkstot << "): " << this->kv.kvec_c[ik].x << " "
                             << this->kv.kvec_c[ik].y << " " << this->kv.kvec_c[ik].z << std::endl;

        for (int ib = 0; ib < GlobalV::NBANDS; ib++)
        {
            GlobalV::ofs_running << " spin" << this->kv.isk[ik] + 1 << "final_state " << ib + 1 << " "
                                 << this->pelec->ekb(ik, ib) * ModuleBase::Ry_to_eV << " "
                                 << this->pelec->wg(ik, ib) * this->kv.nks << std::endl;
        }
        GlobalV::ofs_running << std::endl;
    }
    if (GlobalV::out_bandgap)
    {
        if (!GlobalV::TWO_EFERMI)
        {
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
    }

    // add by jingan in 2018.11.7
    if (GlobalV::CALCULATION == "nscf" && INPUT.towannier90)
    {
#ifdef __LCAO
        if (INPUT.wannier_method == 1)
        {
            toWannier90_LCAO_IN_PW myWannier(
                INPUT.out_wannier_mmn,
                INPUT.out_wannier_amn,
                INPUT.out_wannier_unk, 
                INPUT.out_wannier_eig,
                INPUT.out_wannier_wvfn_formatted,
                INPUT.nnkpfile,
                INPUT.wannier_spin
            );

            myWannier.calculate(this->pelec->ekb, this->pw_wfc, this->pw_big, this->sf, this->kv, this->psi, this->LOWF.ParaV);
        }
        else if (INPUT.wannier_method == 2)
        {
            toWannier90_LCAO myWannier(
                INPUT.out_wannier_mmn,
                INPUT.out_wannier_amn,
                INPUT.out_wannier_unk, 
                INPUT.out_wannier_eig,
                INPUT.out_wannier_wvfn_formatted,
                INPUT.nnkpfile,
                INPUT.wannier_spin
            );

            myWannier.calculate(this->pelec->ekb, this->kv, *(this->psi), this->LOWF.ParaV);
        }
#endif
    }

    // add by jingan
    if (berryphase::berry_phase_flag && ModuleSymmetry::Symmetry::symm_flag != 1)
    {
        berryphase bp(this->LOWF);
        bp.Macroscopic_polarization(this->pw_wfc->npwk_max, this->psi, this->pw_rho, this->pw_wfc, this->kv);
    }

    // below is for DeePKS NSCF calculation
#ifdef __DEEPKS
    const Parallel_Orbitals* pv = this->LOWF.ParaV;
    if (GlobalV::deepks_out_labels || GlobalV::deepks_scf)
    {
        const elecstate::DensityMatrix<TK, double>* dm
            = dynamic_cast<const elecstate::ElecStateLCAO<TK>*>(this->pelec)->get_DM();
        this->dpks_cal_projected_DM(dm);
        GlobalC::ld.cal_descriptor(); // final descriptor
        GlobalC::ld.cal_gedm(GlobalC::ucell.nat);
    }
#endif
    return;
}

template class ESolver_KS_LCAO<double, double>;
template class ESolver_KS_LCAO<std::complex<double>, double>;
template class ESolver_KS_LCAO<std::complex<double>, std::complex<double>>;
} // namespace ModuleESolver
