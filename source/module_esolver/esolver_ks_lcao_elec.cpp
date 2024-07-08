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
#include "module_elecstate/module_dm/cal_dm_psi.h"
#include "module_hamilt_general/module_ewald/H_Ewald_pw.h"
#include "module_hamilt_general/module_vdw/vdw.h"
#include "module_hamilt_lcao/hamilt_lcaodft/LCAO_domain.h"
#include "module_hamilt_lcao/hamilt_lcaodft/operator_lcao/op_exx_lcao.h"
#include "module_hamilt_lcao/hamilt_lcaodft/operator_lcao/operator_lcao.h"
#include "module_hamilt_lcao/module_deltaspin/spin_constrain.h"
#include "module_io/rho_io.h"
#include "module_io/write_pot.h"
#include "module_io/write_wfc_nao.h"
#include "module_base/formatter.h"
#ifdef __EXX
#include "module_io/restart_exx_csr.h"
#endif

namespace ModuleESolver {

template <typename TK, typename TR>
void ESolver_KS_LCAO<TK, TR>::set_matrix_grid(Record_adj& ra) {
    ModuleBase::TITLE("ESolver_KS_LCAO", "set_matrix_grid");
    ModuleBase::timer::tick("ESolver_KS_LCAO", "set_matrix_grid");

    // (1) Find adjacent atoms for each atom.
    GlobalV::SEARCH_RADIUS
        = atom_arrange::set_sr_NL(GlobalV::ofs_running,
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

    // ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running,"SEARCH ADJACENT
    // ATOMS");

    // (3) Periodic condition search for each grid.
    double dr_uniform = 0.001;
    std::vector<double> rcuts;
    std::vector<std::vector<double>> psi_u;
    std::vector<std::vector<double>> dpsi_u;
    std::vector<std::vector<double>> d2psi_u;

    Gint_Tools::init_orb(dr_uniform,
                         rcuts,
                         GlobalC::ucell,
                         psi_u,
                         dpsi_u,
                         d2psi_u);

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
                             this->pw_rho->startz_current,
                             GlobalC::ucell,
                             dr_uniform,
                             rcuts,
                             psi_u,
                             dpsi_u,
                             d2psi_u,
                             GlobalV::NUM_STREAM);
    psi_u.clear();
    psi_u.shrink_to_fit();
    dpsi_u.clear();
    dpsi_u.shrink_to_fit();
    d2psi_u.clear();
    d2psi_u.shrink_to_fit();
    // (2)For each atom, calculate the adjacent atoms in different cells
    // and allocate the space for H(R) and S(R).
    // If k point is used here, allocate HlocR after atom_arrange.
    Parallel_Orbitals* pv = this->LM.ParaV;
    ra.for_2d(*pv, GlobalV::GAMMA_ONLY_LOCAL);
    if (!GlobalV::GAMMA_ONLY_LOCAL) {
        // need to first calculae lgd.
        // using GridT.init.
        this->GridT.cal_nnrg(pv);
    }

    ModuleBase::timer::tick("ESolver_KS_LCAO", "set_matrix_grid");
    return;
}

template <typename TK, typename TR>
void ESolver_KS_LCAO<TK, TR>::beforesolver(const int istep) {
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
    if (this->psi == nullptr) {
        int nsk = 0;
        int ncol = 0;
        if (GlobalV::GAMMA_ONLY_LOCAL) {
            nsk = GlobalV::NSPIN;
            ncol = this->orb_con.ParaV.ncol_bands;
            if (GlobalV::KS_SOLVER == "genelpa"
                || GlobalV::KS_SOLVER == "lapack"
                || GlobalV::KS_SOLVER == "pexsi"
                || GlobalV::KS_SOLVER == "cusolver"
                || GlobalV::KS_SOLVER == "cusolvermp") {
                ncol = this->orb_con.ParaV.ncol;
            }
        } else {
            nsk = this->kv.get_nks();
#ifdef __MPI
            ncol = this->orb_con.ParaV.ncol_bands;
#else
            ncol = GlobalV::NBANDS;
#endif
        }
        this->psi
            = new psi::Psi<TK>(nsk, ncol, this->orb_con.ParaV.nrow, nullptr);
    }

    // prepare grid in Gint
    LCAO_domain::grid_prepare(this->GridT,
                              this->GG,
                              this->GK,
                              *this->pw_rho,
                              *this->pw_big);

    // init Hamiltonian
    if (this->p_hamilt != nullptr) {
        delete this->p_hamilt;
        this->p_hamilt = nullptr;
    }
    if (this->p_hamilt == nullptr) {
        elecstate::DensityMatrix<TK, double>* DM
            = dynamic_cast<elecstate::ElecStateLCAO<TK>*>(this->pelec)
                  ->get_DM();
        this->p_hamilt = new hamilt::HamiltLCAO<TK, TR>(
            GlobalV::GAMMA_ONLY_LOCAL ? &(this->GG) : nullptr,
            GlobalV::GAMMA_ONLY_LOCAL ? nullptr : &(this->GK),
            &(this->LM),
            &this->orb_con.ParaV,
            this->pelec->pot,
            this->kv,
            two_center_bundle_,
#ifdef __EXX
            DM,
            GlobalC::exx_info.info_ri.real_number ? &this->exd->two_level_step
                                                  : &this->exc->two_level_step);
#else
            DM);
#endif
    }
    // init density kernel and wave functions.
    this->LOC.allocate_dm_wfc(this->GridT,
                              this->pelec,
                              this->psi,
                              this->kv,
                              istep);

#ifdef __DEEPKS
    // for each ionic step, the overlap <psi|alpha> must be rebuilt
    // since it depends on ionic positions
    if (GlobalV::deepks_setorb) {
        const Parallel_Orbitals* pv = this->LM.ParaV;
        // build and save <psi(0)|alpha(R)> at beginning
        GlobalC::ld.build_psialpha(GlobalV::CAL_FORCE,
                                   GlobalC::ucell,
                                   GlobalC::ORB,
                                   GlobalC::GridD,
                                   *(two_center_bundle_.overlap_orb_alpha));

        if (GlobalV::deepks_out_unittest) {
            GlobalC::ld.check_psialpha(GlobalV::CAL_FORCE,
                                       GlobalC::ucell,
                                       GlobalC::ORB,
                                       GlobalC::GridD);
        }
    }
#endif
    if (GlobalV::sc_mag_switch) {
        SpinConstrain<TK, base_device::DEVICE_CPU>& sc
            = SpinConstrain<TK, base_device::DEVICE_CPU>::getScInstance();
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
    if (GlobalV::NSPIN == 4 && GlobalV::DOMAG) {
        GlobalC::ucell.cal_ux();
    }
    ModuleBase::timer::tick("ESolver_KS_LCAO", "beforesolver");
}

template <typename TK, typename TR>
void ESolver_KS_LCAO<TK, TR>::before_scf(const int istep) {
    ModuleBase::TITLE("ESolver_KS_LCAO", "before_scf");
    ModuleBase::timer::tick("ESolver_KS_LCAO", "before_scf");

    if (GlobalC::ucell.cell_parameter_updated) {
        this->init_after_vc(INPUT, GlobalC::ucell);
    }
    if (GlobalC::ucell.ionic_position_updated) {
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
    if (vdw_solver != nullptr) {
        this->pelec->f_en.evdw = vdw_solver->get_energy();
    }

    this->beforesolver(istep);

    // Peize Lin add 2016-12-03
#ifdef __EXX // set xc type before the first cal of xc in pelec->init_scf
    if (GlobalC::exx_info.info_ri.real_number) {
        this->exd->exx_beforescf(this->kv, *this->p_chgmix);
    } else {
        this->exc->exx_beforescf(this->kv, *this->p_chgmix);
    }
#endif // __EXX

    this->pelec->init_scf(istep, this->sf.strucFac);
    if (GlobalV::out_chg == 2) {
        for (int is = 0; is < GlobalV::NSPIN; is++) {
            std::stringstream ss;
            ss << GlobalV::global_out_dir << "SPIN" << is + 1
               << "_CHG_INI.cube";
            ModuleIO::write_rho(
#ifdef __MPI
                this->pw_big->bz, // bz first, then nbz
                this->pw_big->nbz,
                this->pw_rho->nplane,
                this->pw_rho->startz_current,
#endif
                this->pelec->charge->rho[is],
                is,
                GlobalV::NSPIN,
                0,
                ss.str(),
                this->pw_rho->nx,
                this->pw_rho->ny,
                this->pw_rho->nz,
                this->pelec->eferm.ef,
                &(GlobalC::ucell),
                11);
        }
    }

    ModuleIO::write_pot(GlobalV::out_pot,
                        GlobalV::NSPIN,
                        GlobalV::global_out_dir,
#ifdef __MPI
                        this->pw_big->bz,
                        this->pw_big->nbz,
                        this->pw_rho->nplane,
                        this->pw_rho->startz_current,
#endif
                        this->pw_rho->nx,
                        this->pw_rho->ny,
                        this->pw_rho->nz,
                        this->pelec->pot->get_effective_v());

    // initalize DMR
    // DMR should be same size with Hamiltonian(R)
    dynamic_cast<elecstate::ElecStateLCAO<TK>*>(this->pelec)
        ->get_DM()
        ->init_DMR(*(dynamic_cast<hamilt::HamiltLCAO<TK, TR>*>(this->p_hamilt)
                         ->getHR()));

    if (GlobalV::dm_to_rho) {
        std::string zipname = "output_DM0.npz";
        elecstate::DensityMatrix<TK, double>* dm
            = dynamic_cast<const elecstate::ElecStateLCAO<TK>*>(this->pelec)
                  ->get_DM();
        this->read_mat_npz(zipname, *(dm->get_DMR_pointer(1)));
        if (GlobalV::NSPIN == 2) {
            zipname = "output_DM1.npz";
            this->read_mat_npz(zipname, *(dm->get_DMR_pointer(2)));
        }

        this->pelec->psiToRho(*this->psi);

        this->create_Output_Rho(0, istep).write();
        if (GlobalV::NSPIN == 2) {
            this->create_Output_Rho(1, istep).write();
        }

        return;
    }

    // the electron charge density should be symmetrized,
    // here is the initialization
    Symmetry_rho srho;
    for (int is = 0; is < GlobalV::NSPIN; is++) {
        srho.begin(is,
                   *(this->pelec->charge),
                   this->pw_rho,
                   GlobalC::Pgrid,
                   GlobalC::ucell.symm);
    }

    // 1. calculate ewald energy.
    // mohan update 2021-02-25
    if (!GlobalV::test_skip_ewald) {
        this->pelec->f_en.ewald_energy
            = H_Ewald_pw::compute_ewald(GlobalC::ucell,
                                        this->pw_rho,
                                        this->sf.strucFac);
    }

    this->p_hamilt->non_first_scf = istep;

    ModuleBase::timer::tick("ESolver_KS_LCAO", "before_scf");
    return;
}

template <typename TK, typename TR>
void ESolver_KS_LCAO<TK, TR>::others(const int istep) {
    ModuleBase::TITLE("ESolver_KS_LCAO", "others");
    ModuleBase::timer::tick("ESolver_KS_LCAO", "others");

    const std::string cal_type = GlobalV::CALCULATION;

    if (cal_type == "get_S") {
        std::cout << FmtCore::format("\n * * * * * *\n << Start %s.\n", "writing the overlap matrix");
        this->get_S();
        std::cout << FmtCore::format(" >> Finish %s.\n * * * * * *\n", "writing the overlap matrix");

        ModuleBase::QUIT();

        // return; // use 'return' will cause segmentation fault. by mohan
        // 2024-06-09
    } else if (cal_type == "test_memory") {
        std::cout << FmtCore::format("\n * * * * * *\n << Start %s.\n", "testing memory");
        Cal_Test::test_memory(this->pw_rho,
                              this->pw_wfc,
                              this->p_chgmix->get_mixing_mode(),
                              this->p_chgmix->get_mixing_ndim());
        std::cout << FmtCore::format(" >> Finish %s.\n * * * * * *\n", "testing memory");
        return;
    } else if (cal_type == "test_neighbour") {
        // test_search_neighbor();
        std::cout << FmtCore::format("\n * * * * * *\n << Start %s.\n", "testing neighbour");
        if (GlobalV::SEARCH_RADIUS < 0) {
            std::cout << " SEARCH_RADIUS : " << GlobalV::SEARCH_RADIUS
                      << std::endl;
            std::cout << " please make sure search_radius > 0" << std::endl;
        }

        atom_arrange::search(GlobalV::SEARCH_PBC,
                             GlobalV::ofs_running,
                             GlobalC::GridD,
                             GlobalC::ucell,
                             GlobalV::SEARCH_RADIUS,
                             GlobalV::test_atom_input,
                             true);
        std::cout << FmtCore::format(" >> Finish %s.\n * * * * * *\n", "testing neighbour");
        return;
    }

    this->beforesolver(istep);
    // pelec should be initialized before these calculations
    this->pelec->init_scf(istep, this->sf.strucFac);
    // self consistent calculations for electronic ground state
    if (GlobalV::CALCULATION == "nscf") {
        this->nscf();
    } else if (cal_type == "get_pchg") {
        std::cout << FmtCore::format("\n * * * * * *\n << Start %s.\n", "getting partial charge");
        IState_Charge ISC(this->psi, &(this->orb_con.ParaV));
        ISC.begin(this->GG,
                  this->pelec->charge->rho,
                  this->pelec->wg,
                  this->pelec->eferm.get_all_ef(),
                  this->pw_rho->nrxx,
                  this->pw_rho->nplane,
                  this->pw_rho->startz_current,
                  this->pw_rho->nx,
                  this->pw_rho->ny,
                  this->pw_rho->nz,
                  this->pw_big->bz,
                  this->pw_big->nbz,
                  GlobalV::GAMMA_ONLY_LOCAL,
                  GlobalV::NBANDS_ISTATE,
                  INPUT.get_out_band_kb(),
                  GlobalV::NBANDS,
                  GlobalV::nelec,
                  GlobalV::NSPIN,
                  GlobalV::NLOCAL,
                  GlobalV::global_out_dir,
                  GlobalV::MY_RANK,
                  GlobalV::ofs_warning,
                  &GlobalC::ucell,
                  &GlobalC::GridD);
        std::cout << FmtCore::format(" >> Finish %s.\n * * * * * *\n", "getting partial charge");
    } else if (cal_type == "get_wf") {
        std::cout << FmtCore::format("\n * * * * * *\n << Start %s.\n", "getting wave function");
        IState_Envelope IEP(this->pelec);
        if (GlobalV::GAMMA_ONLY_LOCAL) {
            IEP.begin(this->psi,
                      this->pw_rho,
                      this->pw_wfc,
                      this->pw_big,
                      this->orb_con.ParaV,
                      this->GG,
                      INPUT.out_wfc_pw,
                      this->wf.out_wfc_r,
                      this->kv,
                      GlobalV::nelec,
                      GlobalV::NBANDS_ISTATE,
                      INPUT.get_out_band_kb(),
                      GlobalV::NBANDS,
                      GlobalV::NSPIN,
                      GlobalV::NLOCAL,
                      GlobalV::global_out_dir);
        } else {
            IEP.begin(this->psi,
                      this->pw_rho,
                      this->pw_wfc,
                      this->pw_big,
                      this->orb_con.ParaV,
                      this->GK,
                      INPUT.out_wfc_pw,
                      this->wf.out_wfc_r,
                      this->kv,
                      GlobalV::nelec,
                      GlobalV::NBANDS_ISTATE,
                      INPUT.get_out_band_kb(),
                      GlobalV::NBANDS,
                      GlobalV::NSPIN,
                      GlobalV::NLOCAL,
                      GlobalV::global_out_dir);
        }
        std::cout << FmtCore::format(" >> Finish %s.\n * * * * * *\n", "getting wave function");
    } else {
        ModuleBase::WARNING_QUIT("ESolver_KS_LCAO<TK, TR>::others",
                                 "CALCULATION type not supported");
    }

    ModuleBase::timer::tick("ESolver_KS_LCAO", "others");
    return;
}

template <>
void ESolver_KS_LCAO<double, double>::get_S(void) {
    ModuleBase::TITLE("ESolver_KS_LCAO", "get_S");
    ModuleBase::WARNING_QUIT("ESolver_KS_LCAO<double,double>::get_S",
                             "not implemented for");
}

template <>
void ESolver_KS_LCAO<std::complex<double>, double>::get_S(void) {
    ModuleBase::TITLE("ESolver_KS_LCAO", "get_S");
    // (1) Find adjacent atoms for each atom.
    GlobalV::SEARCH_RADIUS
        = atom_arrange::set_sr_NL(GlobalV::ofs_running,
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

    if (this->p_hamilt == nullptr) {
        this->p_hamilt = new hamilt::HamiltLCAO<std::complex<double>, double>(
            &this->orb_con.ParaV,
            this->kv,
            *(two_center_bundle_.overlap_orb));
        dynamic_cast<hamilt::OperatorLCAO<std::complex<double>, double>*>(
            this->p_hamilt->ops)
            ->contributeHR();
    }

    // mohan add 2024-06-09
    const std::string fn = GlobalV::global_out_dir + "SR.csr";

    std::cout << " The file is saved in " << fn << std::endl;

    ModuleIO::output_SR(orb_con.ParaV,
                        this->LM,
                        GlobalC::GridD,
                        this->p_hamilt,
                        fn);

    return;
}

template <>
void ESolver_KS_LCAO<std::complex<double>, std::complex<double>>::get_S(void) {
    ModuleBase::TITLE("ESolver_KS_LCAO", "get_S");
    // (1) Find adjacent atoms for each atom.
    GlobalV::SEARCH_RADIUS
        = atom_arrange::set_sr_NL(GlobalV::ofs_running,
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
    this->LM.ParaV = &this->orb_con.ParaV;
    if (this->p_hamilt == nullptr) {
        this->p_hamilt = new hamilt::HamiltLCAO<std::complex<double>,
                                                std::complex<double>>(
            &this->orb_con.ParaV,
            this->kv,
            *(two_center_bundle_.overlap_orb));
        dynamic_cast<
            hamilt::OperatorLCAO<std::complex<double>, std::complex<double>>*>(
            this->p_hamilt->ops)
            ->contributeHR();
    }

    // mohan add 2024-06-09
    const std::string fn = GlobalV::global_out_dir + "SR.csr";

    std::cout << " The file is saved in " << fn << std::endl;

    ModuleIO::output_SR(orb_con.ParaV,
                        this->LM,
                        GlobalC::GridD,
                        this->p_hamilt,
                        fn);

    return;
}

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
    if (this->phsol != nullptr) {
        this->phsol->solve(this->p_hamilt,
                           this->psi[0],
                           this->pelec,
                           GlobalV::KS_SOLVER,
                           true);
    } else {
        ModuleBase::WARNING_QUIT("ESolver_KS_PW",
                                 "HSolver has not been initialed!");
    }

    time_t time_finish = std::time(nullptr);
    ModuleBase::GlobalFunc::OUT_TIME("cal_bands", time_start, time_finish);

    GlobalV::ofs_running << " end of band structure calculation " << std::endl;
    GlobalV::ofs_running << " band eigenvalue in this processor (eV) :"
                         << std::endl;

    const int nspin = GlobalV::NSPIN;
    const int nbands = GlobalV::NBANDS;

    for (int ik = 0; ik < this->kv.get_nks(); ++ik) {
        if (nspin == 2) {
            if (ik == 0) {
                GlobalV::ofs_running << " spin up :" << std::endl;
            }
            if (ik == (this->kv.get_nks() / 2)) {
                GlobalV::ofs_running << " spin down :" << std::endl;
            }
        }

        GlobalV::ofs_running
            << " k-points" << ik + 1 << "(" << this->kv.get_nkstot()
            << "): " << this->kv.kvec_c[ik].x << " " << this->kv.kvec_c[ik].y
            << " " << this->kv.kvec_c[ik].z << std::endl;

        for (int ib = 0; ib < nbands; ++ib) {
            GlobalV::ofs_running
                << " spin" << this->kv.isk[ik] + 1 << "final_state " << ib + 1
                << " " << this->pelec->ekb(ik, ib) * ModuleBase::Ry_to_eV << " "
                << this->pelec->wg(ik, ib) * this->kv.get_nks() << std::endl;
        }
        GlobalV::ofs_running << std::endl;
    }
    if (GlobalV::out_bandgap) {
        std::cout << FmtCore::format("\n * * * * * *\n << Start %s.\n", "band gap calculation");
        if (!GlobalV::TWO_EFERMI) {
            this->pelec->cal_bandgap();
            GlobalV::ofs_running << " E_bandgap "
                                 << this->pelec->bandgap * ModuleBase::Ry_to_eV
                                 << " eV" << std::endl;
        } else {
            this->pelec->cal_bandgap_updw();
            GlobalV::ofs_running
                << " E_bandgap_up "
                << this->pelec->bandgap_up * ModuleBase::Ry_to_eV << " eV"
                << std::endl;
            GlobalV::ofs_running
                << " E_bandgap_dw "
                << this->pelec->bandgap_dw * ModuleBase::Ry_to_eV << " eV"
                << std::endl;
        }
        std::cout << FmtCore::format(" >> Finish %s.\n * * * * * *\n", "band gap calculation");
    }

    // add by jingan in 2018.11.7
    if (GlobalV::CALCULATION == "nscf" && INPUT.towannier90) {
#ifdef __LCAO
        std::cout << FmtCore::format("\n * * * * * *\n << Start %s.\n", "Wave function to Wannier90");
        if (INPUT.wannier_method == 1) {
            toWannier90_LCAO_IN_PW myWannier(INPUT.out_wannier_mmn,
                                             INPUT.out_wannier_amn,
                                             INPUT.out_wannier_unk,
                                             INPUT.out_wannier_eig,
                                             INPUT.out_wannier_wvfn_formatted,
                                             INPUT.nnkpfile,
                                             INPUT.wannier_spin);

            myWannier.calculate(this->pelec->ekb,
                                this->pw_wfc,
                                this->pw_big,
                                this->sf,
                                this->kv,
                                this->psi,
                                &(this->orb_con.ParaV));
        } else if (INPUT.wannier_method == 2) {
            toWannier90_LCAO myWannier(INPUT.out_wannier_mmn,
                                       INPUT.out_wannier_amn,
                                       INPUT.out_wannier_unk,
                                       INPUT.out_wannier_eig,
                                       INPUT.out_wannier_wvfn_formatted,
                                       INPUT.nnkpfile,
                                       INPUT.wannier_spin);

            myWannier.calculate(this->pelec->ekb,
                                this->kv,
                                *(this->psi),
                                &(this->orb_con.ParaV));
        }
        std::cout << FmtCore::format(" >> Finish %s.\n * * * * * *\n", "Wave function to Wannier90");
#endif
    }

    // add by jingan
    if (berryphase::berry_phase_flag
        && ModuleSymmetry::Symmetry::symm_flag != 1) {
        std::cout << FmtCore::format("\n * * * * * *\n << Start %s.\n", "Berry phase calculation");
        berryphase bp(this->orb_con.ParaV);
        bp.lcao_init(this->kv,
                     this->GridT); // additional step before calling
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
    if (GlobalV::deepks_out_labels || GlobalV::deepks_scf) {
        std::cout << FmtCore::format("\n * * * * * *\n << Start %s.\n", "DeepKS output");
        const elecstate::DensityMatrix<TK, double>* dm
            = dynamic_cast<const elecstate::ElecStateLCAO<TK>*>(this->pelec)
                  ->get_DM();
        this->dpks_cal_projected_DM(dm);
        GlobalC::ld.cal_descriptor(GlobalC::ucell.nat); // final descriptor
        GlobalC::ld.cal_gedm(GlobalC::ucell.nat);
        std::cout << FmtCore::format(" >> Finish %s.\n * * * * * *\n", "DeepKS output");
    }
#endif

    this->create_Output_Mat_Sparse(0).write();

    // mulliken charge analysis
    if (GlobalV::out_mul) {
        std::cout << FmtCore::format("\n * * * * * *\n << Start %s.\n", "Mulliken charge analysis");
        elecstate::ElecStateLCAO<TK>* pelec_lcao
            = dynamic_cast<elecstate::ElecStateLCAO<TK>*>(this->pelec);
        this->pelec->calculate_weights();
        this->pelec->calEBand();
        elecstate::cal_dm_psi(&(this->orb_con.ParaV),
                              pelec_lcao->wg,
                              *(this->psi),
                              *(pelec_lcao->get_DM()));
        this->cal_mag(istep, true);
        std::cout << FmtCore::format(" >> Finish %s.\n * * * * * *\n", "Mulliken charge analysis");
    }

    /// write potential
    this->create_Output_Potential(0).write();

    // write wfc
    if (INPUT.out_wfc_lcao)
    {
        std::cout << FmtCore::format("\n * * * * * *\n << Start %s.\n", "writing wave function");
        ModuleIO::write_wfc_nao(INPUT.out_wfc_lcao,
                                *this->psi,
                                this->pelec->ekb,
                                this->pelec->wg,
                                this->pelec->klist->kvec_c,
                                this->orb_con.ParaV,
                                istep);
        std::cout << FmtCore::format(" >> Finish %s.\n * * * * * *\n", "writing wave function");
    }
}

template class ESolver_KS_LCAO<double, double>;
template class ESolver_KS_LCAO<std::complex<double>, double>;
template class ESolver_KS_LCAO<std::complex<double>, std::complex<double>>;
} // namespace ModuleESolver
