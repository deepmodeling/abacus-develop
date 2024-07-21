#include "esolver_ks_lcao.h"

#include "module_base/global_variable.h"
#include "module_base/tool_title.h"
#include "module_io/dos_nao.h"
#include "module_io/nscf_band.h"
#include "module_io/output_dmk.h"
#include "module_io/output_log.h"
#include "module_io/output_mulliken.h"
#include "module_io/output_sk.h"
#include "module_io/to_qo.h"
#include "module_io/write_HS.h"
#include "module_io/write_Vxc.hpp"
#include "module_io/write_istate_info.h"
#include "module_io/write_proj_band_lcao.h"
#include "module_parameter/parameter.h"

//--------------temporary----------------------------
#include <memory>

#include "module_base/global_function.h"
#include "module_cell/module_neighbor/sltk_grid_driver.h"
#include "module_elecstate/module_charge/symmetry_rho.h"
#include "module_elecstate/occupy.h"
#include "module_hamilt_lcao/hamilt_lcaodft/LCAO_domain.h" // need divide_HS_in_frag
#include "module_hamilt_lcao/module_dftu/dftu.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_io/print_info.h"
#ifdef __EXX
#include "module_ri/RPA_LRI.h"
#include "module_io/restart_exx_csr.h"
#endif

#ifdef __DEEPKS
#include "module_hamilt_lcao/module_deepks/LCAO_deepks.h"
#include "module_hamilt_lcao/module_deepks/LCAO_deepks_interface.h"
#endif
//-----force& stress-------------------
#include "module_hamilt_lcao/hamilt_lcaodft/FORCE_STRESS.h"

//-----HSolver ElecState Hamilt--------
#include "module_elecstate/elecstate_lcao.h"
#include "module_hamilt_lcao/hamilt_lcaodft/hamilt_lcao.h"
#include "module_hsolver/hsolver_lcao.h"
// function used by deepks
#include "module_elecstate/cal_dm.h"
//---------------------------------------------------

#include "module_hamilt_lcao/module_deltaspin/spin_constrain.h"
#include "module_io/io_dmk.h"
#include "module_io/write_dmr.h"
#include "module_io/write_wfc_nao.h"

namespace ModuleESolver
{

//------------------------------------------------------------------------------
//! the 1st function of ESolver_KS_LCAO: constructor
//! mohan add 2024-05-11
//------------------------------------------------------------------------------
template <typename TK, typename TR>
ESolver_KS_LCAO<TK, TR>::ESolver_KS_LCAO()
{
    this->classname = "ESolver_KS_LCAO";
    this->basisname = "LCAO";
#ifdef __EXX
    // 1. currently this initialization must be put in constructor rather than `before_all_runners()`
    //  because the latter is not reused by ESolver_LCAO_TDDFT, 
    //  which cause the failure of the subsequent procedure reused by ESolver_LCAO_TDDFT
    // 2. always construct but only initialize when if(cal_exx) is true
    //  because some members like two_level_step are used outside if(cal_exx)
    if (GlobalC::exx_info.info_ri.real_number)
    {
        this->exx_lri_double = std::make_shared<Exx_LRI<double>>(GlobalC::exx_info.info_ri);
        this->exd = std::make_shared<Exx_LRI_Interface<TK, double>>(exx_lri_double);
    }
    else
    {
        this->exx_lri_complex = std::make_shared<Exx_LRI<std::complex<double>>>(GlobalC::exx_info.info_ri);
        this->exc = std::make_shared<Exx_LRI_Interface<TK, std::complex<double>>>(exx_lri_complex);
    }
#endif
}

//------------------------------------------------------------------------------
//! the 2nd function of ESolver_KS_LCAO: deconstructor
//! mohan add 2024-05-11
//------------------------------------------------------------------------------
template <typename TK, typename TR>
ESolver_KS_LCAO<TK, TR>::~ESolver_KS_LCAO()
{
}

//------------------------------------------------------------------------------
//! the 3rd function of ESolver_KS_LCAO: init
//! mohan add 2024-05-11
//! 1) calculate overlap matrix S or initialize
//! 2) init ElecState
//! 3) init LCAO basis
//! 4) redundant ParaV and LM pointers
//! 5) initialize Density Matrix
//! 6) initialize Hamilt in LCAO
//! 7) initialize exx
//! 8) Quxin added for DFT+U
//! 9) ppcell
//! 10) init HSolver
//! 11) inititlize the charge density.
//! 12) initialize the potential.
//! 13) initialize deepks
//! 14) set occupations?
//------------------------------------------------------------------------------
template <typename TK, typename TR>
void ESolver_KS_LCAO<TK, TR>::before_all_runners(const Input_para& inp, UnitCell& ucell)
{
    ModuleBase::TITLE("ESolver_KS_LCAO", "before_all_runners");
    ModuleBase::timer::tick("ESolver_KS_LCAO", "before_all_runners");

    // 1) calculate overlap matrix S
    if (GlobalV::CALCULATION == "get_S")
    {
        // 1.1) read pseudopotentials
        ucell.read_pseudo(GlobalV::ofs_running);

        // 1.2) symmetrize things
        if (ModuleSymmetry::Symmetry::symm_flag == 1)
        {
            ucell.symm.analy_sys(ucell.lat, ucell.st, ucell.atoms, GlobalV::ofs_running);
            ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "SYMMETRY");
        }

        // 1.3) Setup k-points according to symmetry.
        this->kv
            .set(ucell.symm, GlobalV::global_kpoint_card, GlobalV::NSPIN, ucell.G, ucell.latvec, GlobalV::ofs_running);
        ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "INIT K-POINTS");

        Print_Info::setup_parameters(ucell, this->kv);
    }
    else
    {
        // 1) else, call before_all_runners() in ESolver_KS
        ESolver_KS<TK>::before_all_runners(inp, ucell);
    } // end ifnot get_S

    // 2) init ElecState
    // autoset nbands in ElecState, it should before basis_init (for Psi 2d
    // divid)
    if (this->pelec == nullptr)
    {
        // TK stands for double and complex<double>?
        this->pelec = new elecstate::ElecStateLCAO<TK>(
            &(this->chr), // use which parameter?
            &(this->kv),
            this->kv.get_nks(),
            &(this->GG),  // mohan add 2024-04-01
            &(this->GK),  // mohan add 2024-04-01
            this->pw_rho,
            this->pw_big);
    }

    // 3) init LCAO basis
    // reading the localized orbitals/projectors
    // construct the interpolation tables.
    this->init_basis_lcao(inp, ucell);
    //------------------init Basis_lcao----------------------

    // 5) initialize density matrix
    // DensityMatrix is allocated here, DMK is also initialized here
    // DMR is not initialized here, it will be constructed in each before_scf
    dynamic_cast<elecstate::ElecStateLCAO<TK>*>(this->pelec)
        ->init_DM(&this->kv, &(this->ParaV), GlobalV::NSPIN);

    // this function should be removed outside of the function
    if (GlobalV::CALCULATION == "get_S")
    {
        ModuleBase::timer::tick("ESolver_KS_LCAO", "init");
        return;
    }

    // 6) initialize Hamilt in LCAO
    // * allocate H and S matrices according to computational resources
    // * set the 'trace' between local H/S and global H/S
    LCAO_domain::divide_HS_in_frag(GlobalV::GAMMA_ONLY_LOCAL, ParaV, this->kv.get_nks());

#ifdef __EXX
    // 7) initialize exx
    // PLEASE simplify the Exx_Global interface
    if (GlobalV::CALCULATION == "scf" || GlobalV::CALCULATION == "relax"
        || GlobalV::CALCULATION == "cell-relax"
        || GlobalV::CALCULATION == "md")
    {
        if (GlobalC::exx_info.info_global.cal_exx)
        {
            XC_Functional::set_xc_first_loop(ucell);
            // initialize 2-center radial tables for EXX-LRI
            if (GlobalC::exx_info.info_ri.real_number) { this->exx_lri_double->init(MPI_COMM_WORLD, this->kv); }
            else { this->exx_lri_complex->init(MPI_COMM_WORLD, this->kv); }
        }
    }
#endif

    // 8) initialize DFT+U
    if (GlobalV::dft_plus_u) {
        GlobalC::dftu.init(ucell, &this->ParaV, this->kv.get_nks());
    }

    // 9) initialize ppcell
    GlobalC::ppcell.init_vloc(GlobalC::ppcell.vloc, this->pw_rho);

    // 10) initialize the HSolver
    if (this->phsol == nullptr)
    {
        this->phsol = new hsolver::HSolverLCAO<TK>(&(this->ParaV));
        this->phsol->method = GlobalV::KS_SOLVER;
    }

    // 11) inititlize the charge density
    this->pelec->charge->allocate(GlobalV::NSPIN);
    this->pelec->omega = GlobalC::ucell.omega;

    // 12) initialize the potential
    if (this->pelec->pot == nullptr)
    {
        this->pelec->pot = new elecstate::Potential(this->pw_rhod,
                                                    this->pw_rho,
                                                    &GlobalC::ucell,
                                                    &(GlobalC::ppcell.vloc),
                                                    &(this->sf),
                                                    &(this->pelec->f_en.etxc),
                                                    &(this->pelec->f_en.vtxc));
    }

#ifdef __DEEPKS
    // 13) initialize deepks
    if (GlobalV::deepks_scf)
    {
        // load the DeePKS model from deep neural network
        GlobalC::ld.load_model(PARAM.inp.deepks_model);
        // read pdm from file for NSCF or SCF-restart, do it only once in whole calculation
        GlobalC::ld.read_projected_DM((GlobalV::init_chg == "file"), GlobalV::deepks_equiv, *GlobalC::ORB.Alpha);
    }
#endif

    // 14) set occupations
    if (PARAM.inp.ocp)
    {
        this->pelec->fixed_weights(PARAM.inp.ocp_kb, GlobalV::NBANDS, GlobalV::nelec);
    }

    // 15) if kpar is not divisible by nks, print a warning
    if (GlobalV::KPAR_LCAO > 1)
    {
        if (this->kv.get_nks() % GlobalV::KPAR_LCAO != 0)
        {
            ModuleBase::WARNING("ESolver_KS_LCAO::before_all_runners",
                                 "nks is not divisible by kpar.");
            std::cout << "\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
                        "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
                        "%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl;
            std::cout << " Warning: nks (" << this->kv.get_nks() << ") is not divisible by kpar ("
                      << GlobalV::KPAR_LCAO << ")." << std::endl;
            std::cout << " This may lead to poor load balance. It is strongly suggested to" << std::endl;
            std::cout << " set nks to be divisible by kpar, but if this is really what" << std::endl;
            std::cout << " you want, please ignore this warning." << std::endl;
            std::cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
                             "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
                             "%%%%%%%%%%%%\n";
        }
    }

    ModuleBase::timer::tick("ESolver_KS_LCAO", "before_all_runners");
    return;
}

//------------------------------------------------------------------------------
//! the 4th function of ESolver_KS_LCAO: init_after_vc
//! mohan add 2024-05-11
//------------------------------------------------------------------------------
template <typename TK, typename TR>
void ESolver_KS_LCAO<TK, TR>::init_after_vc(const Input_para& inp, UnitCell& ucell)
{
    ModuleBase::TITLE("ESolver_KS_LCAO", "init_after_vc");
    ModuleBase::timer::tick("ESolver_KS_LCAO", "init_after_vc");

    ESolver_KS<TK>::init_after_vc(inp, ucell);
    if (inp.mdp.md_prec_level == 2)
    {
        delete this->pelec;
        this->pelec = new elecstate::ElecStateLCAO<TK>(
            &(this->chr),
            &(this->kv),
            this->kv.get_nks(),
            &(this->GG), // mohan add 2024-04-01
            &(this->GK), // mohan add 2024-04-01
            this->pw_rho,
            this->pw_big);

        dynamic_cast<elecstate::ElecStateLCAO<TK>*>(this->pelec)
            ->init_DM(&this->kv, &this->ParaV, GlobalV::NSPIN);

        GlobalC::ppcell.init_vloc(GlobalC::ppcell.vloc, this->pw_rho);

        this->pelec->charge->allocate(GlobalV::NSPIN);
        this->pelec->omega = GlobalC::ucell.omega;

        // Initialize the potential.
        if (this->pelec->pot == nullptr)
        {
            this->pelec->pot = new elecstate::Potential(this->pw_rhod,
                                                        this->pw_rho,
                                                        &GlobalC::ucell,
                                                        &(GlobalC::ppcell.vloc),
                                                        &(this->sf),
                                                        &(this->pelec->f_en.etxc),
                                                        &(this->pelec->f_en.vtxc));
        }
    }

    ModuleBase::timer::tick("ESolver_KS_LCAO", "init_after_vc");
    return;
}

//------------------------------------------------------------------------------
//! the 5th function of ESolver_KS_LCAO: cal_energy
//! mohan add 2024-05-11
//------------------------------------------------------------------------------
template <typename TK, typename TR>
double ESolver_KS_LCAO<TK, TR>::cal_energy()
{
    ModuleBase::TITLE("ESolver_KS_LCAO", "cal_energy");

    return this->pelec->f_en.etot;
}

//------------------------------------------------------------------------------
//! the 6th function of ESolver_KS_LCAO: cal_force
//! mohan add 2024-05-11
//------------------------------------------------------------------------------
template <typename TK, typename TR>
void ESolver_KS_LCAO<TK, TR>::cal_force(ModuleBase::matrix& force)
{
    ModuleBase::TITLE("ESolver_KS_LCAO", "cal_force");
    ModuleBase::timer::tick("ESolver_KS_LCAO", "cal_force");

    Force_Stress_LCAO<TK> fsl(this->RA, GlobalC::ucell.nat);

    fsl.getForceStress(GlobalV::CAL_FORCE,
                       GlobalV::CAL_STRESS,
                       GlobalV::TEST_FORCE,
                       GlobalV::TEST_STRESS,
                       this->ParaV,
                       this->pelec,
                       this->psi,
                       this->GG, // mohan add 2024-04-01
                       this->GK, // mohan add 2024-04-01
                       two_center_bundle_,
                       force,
                       this->scs,
                       this->sf,
                       this->kv,
                       this->pw_rho,
#ifdef __EXX
                       *this->exx_lri_double,
                       *this->exx_lri_complex,
#endif
                       &GlobalC::ucell.symm);

    // delete RA after cal_force

    this->RA.delete_grid();

    this->have_force = true;

    ModuleBase::timer::tick("ESolver_KS_LCAO", "cal_force");
}

//------------------------------------------------------------------------------
//! the 7th function of ESolver_KS_LCAO: cal_stress
//! mohan add 2024-05-11
//------------------------------------------------------------------------------
template <typename TK, typename TR>
void ESolver_KS_LCAO<TK, TR>::cal_stress(ModuleBase::matrix& stress)
{
    ModuleBase::TITLE("ESolver_KS_LCAO", "cal_stress");
    ModuleBase::timer::tick("ESolver_KS_LCAO", "cal_stress");

    if (!this->have_force)
    {
        ModuleBase::matrix fcs;
        this->cal_force(fcs);
    }
    stress = this->scs; // copy the stress
    this->have_force = false;

    ModuleBase::timer::tick("ESolver_KS_LCAO", "cal_stress");
}

//------------------------------------------------------------------------------
//! the 8th function of ESolver_KS_LCAO: after_all_runners
//! mohan add 2024-05-11
//------------------------------------------------------------------------------
template <typename TK, typename TR>
void ESolver_KS_LCAO<TK, TR>::after_all_runners()
{
    ModuleBase::TITLE("ESolver_KS_LCAO", "after_all_runners");
    ModuleBase::timer::tick("ESolver_KS_LCAO", "after_all_runners");

    GlobalV::ofs_running << "\n\n --------------------------------------------" << std::endl;
    GlobalV::ofs_running << std::setprecision(16);
    GlobalV::ofs_running << " !FINAL_ETOT_IS " << this->pelec->f_en.etot * ModuleBase::Ry_to_eV << " eV" << std::endl;
    GlobalV::ofs_running << " --------------------------------------------\n\n" << std::endl;

    if (PARAM.inp.out_dos != 0 || PARAM.inp.out_band[0] != 0 || PARAM.inp.out_proj_band != 0)
    {
        GlobalV::ofs_running << "\n\n\n\n";
        GlobalV::ofs_running << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
                                ">>>>>>>>>>>>>>>>>>>>>>>>>"
                             << std::endl;
        GlobalV::ofs_running << " |                                            "
                                "                        |"
                             << std::endl;
        GlobalV::ofs_running << " | Post-processing of data:                   "
                                "                        |"
                             << std::endl;
        GlobalV::ofs_running << " | DOS (density of states) and bands will be "
                                "output here.             |"
                             << std::endl;
        GlobalV::ofs_running << " | If atomic orbitals are used, Mulliken "
                                "charge analysis can be done. |"
                             << std::endl;
        GlobalV::ofs_running << " | Also the .bxsf file containing fermi "
                                "surface information can be    |"
                             << std::endl;
        GlobalV::ofs_running << " | done here.                                 "
                                "                        |"
                             << std::endl;
        GlobalV::ofs_running << " |                                            "
                                "                        |"
                             << std::endl;
        GlobalV::ofs_running << " <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"
                                "<<<<<<<<<<<<<<<<<<<<<<<<<"
                             << std::endl;
        GlobalV::ofs_running << "\n\n\n\n";
    }
    // qianrui modify 2020-10-18
    if (GlobalV::CALCULATION == "scf" || GlobalV::CALCULATION == "md" || GlobalV::CALCULATION == "relax")
    {
        ModuleIO::write_istate_info(this->pelec->ekb, this->pelec->wg, this->kv, &(GlobalC::Pkpoints));
    }

    const int nspin0 = (GlobalV::NSPIN == 2) ? 2 : 1;

    if (PARAM.inp.out_band[0])
    {
        for (int is = 0; is < nspin0; is++)
        {
            std::stringstream ss2;
            ss2 << GlobalV::global_out_dir << "BANDS_" << is + 1 << ".dat";
            GlobalV::ofs_running << "\n Output bands in file: " << ss2.str() << std::endl;
            ModuleIO::nscf_band(is,
                                ss2.str(),
                                GlobalV::NBANDS,
                                0.0,
                                PARAM.inp.out_band[1],
                                this->pelec->ekb,
                                this->kv,
                                &(GlobalC::Pkpoints));
        }
    } // out_band

    if (PARAM.inp.out_proj_band) // Projeced band structure added by jiyy-2022-4-20
    {
        ModuleIO::write_proj_band_lcao(this->psi, this->ParaV, this->pelec, this->kv, GlobalC::ucell, this->p_hamilt);
    }

    if (PARAM.inp.out_dos)
    {
        ModuleIO::out_dos_nao(this->psi,
                              this->ParaV,
                              this->pelec->ekb,
                              this->pelec->wg,
                              PARAM.inp.dos_edelta_ev,
                              PARAM.inp.dos_scale,
                              PARAM.inp.dos_sigma,
                              *(this->pelec->klist),
                              GlobalC::Pkpoints,
                              GlobalC::ucell,
                              this->pelec->eferm,
                              GlobalV::NBANDS,
                              this->p_hamilt);
    }

    if (PARAM.inp.out_mat_xc)
    {
        ModuleIO::write_Vxc<TK, TR>(GlobalV::NSPIN,
            GlobalV::NLOCAL,
            GlobalV::DRANK,
            &this->ParaV,
            *this->psi,
            GlobalC::ucell,
            this->sf,
            *this->pw_rho,
            *this->pw_rhod,
            GlobalC::ppcell.vloc,
            *this->pelec->charge,
            this->GG,
            this->GK,
            this->kv,
            this->pelec->wg,
            GlobalC::GridD
#ifdef __EXX
            , this->exx_lri_double ? &this->exx_lri_double->Hexxs : nullptr
            , this->exx_lri_complex ? &this->exx_lri_complex->Hexxs : nullptr
#endif
        );
    }

    ModuleBase::timer::tick("ESolver_KS_LCAO", "after_all_runners");
}

//------------------------------------------------------------------------------
//! the 9th function of ESolver_KS_LCAO: init_basis_lcao
//! mohan add 2024-05-11
//------------------------------------------------------------------------------
template <typename TK, typename TR>
void ESolver_KS_LCAO<TK, TR>::init_basis_lcao(const Input_para& inp, UnitCell& ucell)
{
    ModuleBase::TITLE("ESolver_KS_LCAO", "init_basis_lcao");

    // autoset NB2D first
    if (GlobalV::NB2D == 0)
    {
        if (GlobalV::NLOCAL > 0)
        {
            GlobalV::NB2D = (GlobalV::NSPIN == 4) ? 2 : 1;
        }
        if (GlobalV::NLOCAL > 500)
        {
            GlobalV::NB2D = 32;
        }
        if (GlobalV::NLOCAL > 1000)
        {
            GlobalV::NB2D = 64;
        }
    }

    // * reading the localized orbitals/projectors
    // * construct the interpolation tables.

    two_center_bundle_.build_orb(ucell.ntype, ucell.orbital_fn);
    two_center_bundle_.build_alpha(GlobalV::deepks_setorb, &ucell.descriptor_file);
    two_center_bundle_.build_orb_onsite(PARAM.inp.onsite_radius);
    // currently deepks only use one descriptor file, so cast bool to int is
    // fine

    // TODO Due to the omnipresence of GlobalC::ORB, we still have to rely
    // on the old interface for now.
    two_center_bundle_.to_LCAO_Orbitals(GlobalC::ORB, inp.lcao_ecut, inp.lcao_dk, inp.lcao_dr, inp.lcao_rmax);

    ucell.infoNL.setupNonlocal(ucell.ntype, ucell.atoms, GlobalV::ofs_running, GlobalC::ORB);

    two_center_bundle_.build_beta(ucell.ntype, ucell.infoNL.Beta);

    int Lmax = 0;
#ifdef __EXX
    Lmax = GlobalC::exx_info.info_ri.abfs_Lmax;
#endif

#ifdef USE_NEW_TWO_CENTER
    two_center_bundle_.tabulate();
#else
    two_center_bundle_.tabulate(inp.lcao_ecut, inp.lcao_dk, inp.lcao_dr, inp.lcao_rmax);
#endif

    // setup_2d_division
#ifdef __MPI
    // storage form of H and S matrices on each processor
    // is determined in 'divide_HS_2d' subroutine

    int try_nb = ParaV.init(GlobalV::NLOCAL, GlobalV::NLOCAL, GlobalV::NB2D, DIAG_WORLD);
    try_nb += ParaV.set_nloc_wfc_Eij(GlobalV::NBANDS, GlobalV::ofs_running, GlobalV::ofs_warning);
    if (try_nb != 0)
    {
        ParaV.set(GlobalV::NLOCAL, GlobalV::NLOCAL, 1, ParaV.blacs_ctxt);
        try_nb = ParaV.set_nloc_wfc_Eij(GlobalV::NBANDS, GlobalV::ofs_running, GlobalV::ofs_warning);
    }

    // init blacs context for genelpa
    ParaV.set_desc_wfc_Eij(GlobalV::NLOCAL, GlobalV::NBANDS, ParaV.nrow);

#else
    ParaV.set_serial(GlobalV::NLOCAL, GlobalV::NLOCAL);
    ParaV.nrow_bands = GlobalV::NLOCAL;
    ParaV.ncol_bands = GlobalV::NBANDS;
    // Zhang Xiaoyang enable the serial version of LCAO and recovered this function usage. 2024-07-06
#endif

    ParaV.set_atomic_trace(GlobalC::ucell.get_iat2iwt(), GlobalC::ucell.nat, GlobalV::NLOCAL);

    return;
}

//------------------------------------------------------------------------------
//! the 10th function of ESolver_KS_LCAO: iter_init
//! mohan add 2024-05-11
//------------------------------------------------------------------------------
template <typename TK, typename TR>
void ESolver_KS_LCAO<TK, TR>::iter_init(const int istep, const int iter)
{
    ModuleBase::TITLE("ESolver_KS_LCAO", "iter_init");

    if (iter == 1)
    {
        this->p_chgmix->init_mixing(); // init mixing
        this->p_chgmix->mixing_restart_step = GlobalV::SCF_NMAX + 1;
        this->p_chgmix->mixing_restart_count = 0;
        // this output will be removed once the feeature is stable
        if (GlobalC::dftu.uramping > 0.01)
        {
            std::cout << " U-Ramping! Current U = ";
            for (int i = 0; i < GlobalC::dftu.U0.size(); i++)
            {
                std::cout << GlobalC::dftu.U[i] * ModuleBase::Ry_to_eV << " ";
            }
            std::cout << " eV " << std::endl;
        }
    }
    // for mixing restart
    if (iter == this->p_chgmix->mixing_restart_step && GlobalV::MIXING_RESTART > 0.0)
    {
        this->p_chgmix->init_mixing();
        this->p_chgmix->mixing_restart_count++;
        if (GlobalV::dft_plus_u)
        {
            GlobalC::dftu.uramping_update(); // update U by uramping if uramping > 0.01
            if (GlobalC::dftu.uramping > 0.01)
            {
                std::cout << " U-Ramping! Current U = ";
                for (int i = 0; i < GlobalC::dftu.U0.size(); i++)
                {
                    std::cout << GlobalC::dftu.U[i] * ModuleBase::Ry_to_eV << " ";
                }
                std::cout << " eV " << std::endl;
            }
            if (GlobalC::dftu.uramping > 0.01 && !GlobalC::dftu.u_converged())
            {
                this->p_chgmix->mixing_restart_step = GlobalV::SCF_NMAX + 1;
            }
        }
        if (GlobalV::MIXING_DMR) // for mixing_dmr
        {
            // allocate memory for dmr_mdata
            const elecstate::DensityMatrix<TK, double>* dm
                = dynamic_cast<const elecstate::ElecStateLCAO<TK>*>(this->pelec)->get_DM();
            int nnr_tmp = dm->get_DMR_pointer(1)->get_nnr();
            this->p_chgmix->allocate_mixing_dmr(nnr_tmp);
        }
    }

    // mohan update 2012-06-05
    this->pelec->f_en.deband_harris = this->pelec->cal_delta_eband();

    // mohan move it outside 2011-01-13
    // first need to calculate the weight according to
    // electrons number.
    if (istep == 0 && this->wf.init_wfc == "file")
    {
        if (iter == 1)
        {
            std::cout << " WAVEFUN -> CHARGE " << std::endl;

            // calculate the density matrix using read in wave functions
            // and the ncalculate the charge density on grid.

            this->pelec->skip_weights = true;
            this->pelec->psiToRho(*this->psi);
            this->pelec->skip_weights = false;

            // calculate the local potential(rho) again.
            // the grid integration will do in later grid integration.

            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            // a puzzle remains here.
            // if I don't renew potential,
            // The scf_thr is very small.
            // OneElectron, Hartree and
            // Exc energy are all correct
            // except the band energy.
            //
            // solved by mohan 2010-09-10
            // there are there rho here:
            // rho1: formed by read in orbitals.
            // rho2: atomic rho, used to construct H
            // rho3: generated by after diagonalize
            // here converged because rho3 and rho1
            // are very close.
            // so be careful here, make sure
            // rho1 and rho2 are the same rho.
            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            if (GlobalV::NSPIN == 4)
            {
                GlobalC::ucell.cal_ux();
            }

            //! update the potentials by using new electron charge density
            this->pelec->pot->update_from_charge(this->pelec->charge, &GlobalC::ucell);

            //! compute the correction energy for metals
            this->pelec->f_en.descf = this->pelec->cal_delta_escf();
        }
    }

#ifdef __EXX
    // calculate exact-exchange
    if (GlobalC::exx_info.info_ri.real_number)
    {
        this->exd->exx_eachiterinit(*dynamic_cast<const elecstate::ElecStateLCAO<TK>*>(this->pelec)->get_DM(), iter);
    }
    else
    {
        this->exc->exx_eachiterinit(*dynamic_cast<const elecstate::ElecStateLCAO<TK>*>(this->pelec)->get_DM(), iter);
    }
#endif

    if (GlobalV::dft_plus_u)
    {
        if (istep != 0 || iter != 1)
        {
            GlobalC::dftu.set_dmr(dynamic_cast<elecstate::ElecStateLCAO<TK>*>(this->pelec)->get_DM());
        }
        // Calculate U and J if Yukawa potential is used
        GlobalC::dftu.cal_slater_UJ(this->pelec->charge->rho, this->pw_rho->nrxx);
    }

#ifdef __DEEPKS
    // the density matrixes of DeePKS have been updated in each iter
    GlobalC::ld.set_hr_cal(true);

    // HR in HamiltLCAO should be recalculate
    if (GlobalV::deepks_scf)
    {
        this->p_hamilt->refresh();
    }
#endif

    if (PARAM.inp.vl_in_h)
    {
        // update Gint_K
        if (!GlobalV::GAMMA_ONLY_LOCAL)
        {
            this->GK.renew();
        }
        // update real space Hamiltonian
        this->p_hamilt->refresh();
    }

    // run the inner lambda loop to contrain atomic moments with the DeltaSpin
    // method
    if (PARAM.inp.sc_mag_switch && iter > PARAM.inp.sc_scf_nmin)
    {
        SpinConstrain<TK, base_device::DEVICE_CPU>& sc = SpinConstrain<TK, base_device::DEVICE_CPU>::getScInstance();
        sc.run_lambda_loop(iter - 1);
    }
}

//------------------------------------------------------------------------------
//! the 11th function of ESolver_KS_LCAO: hamilt2density
//! mohan add 2024-05-11
//! 1) save input rho
//! 2) save density matrix DMR for mixing
//! 3) solve the Hamiltonian and output band gap
//! 4) print bands for each k-point and each band
//! 5) EXX:
//! 6) DFT+U: compute local occupation number matrix and energy correction
//! 7) DeePKS: compute delta_e
//! 8) DeltaSpin:
//! 9) use new charge density to calculate energy
//! 10) symmetrize the charge density
//! 11) compute magnetization, only for spin==2
//! 12) calculate delta energy
//------------------------------------------------------------------------------
template <typename TK, typename TR>
void ESolver_KS_LCAO<TK, TR>::hamilt2density(int istep, int iter, double ethr)
{
    ModuleBase::TITLE("ESolver_KS_LCAO", "hamilt2density");

    // 1) save input rho
    this->pelec->charge->save_rho_before_sum_band();

    // 2) save density matrix DMR for mixing
    if (GlobalV::MIXING_RESTART > 0 && GlobalV::MIXING_DMR && this->p_chgmix->mixing_restart_count > 0)
    {
        elecstate::DensityMatrix<TK, double>* dm = dynamic_cast<elecstate::ElecStateLCAO<TK>*>(this->pelec)->get_DM();
        dm->save_DMR();
    }

    // 3) solve the Hamiltonian and output band gap
    if (this->phsol != nullptr)
    {
        // reset energy
        this->pelec->f_en.eband = 0.0;
        this->pelec->f_en.demet = 0.0;

        this->phsol->solve(this->p_hamilt, this->psi[0], this->pelec, GlobalV::KS_SOLVER);

        if (PARAM.inp.out_bandgap)
        {
            if (!GlobalV::TWO_EFERMI)
            {
                this->pelec->cal_bandgap();
            }
            else
            {
                this->pelec->cal_bandgap_updw();
            }
        }
    }
    else
    {
        ModuleBase::WARNING_QUIT("ESolver_KS_PW", "HSolver has not been initialed!");
    }

    // 4) print bands for each k-point and each band
    for (int ik = 0; ik < this->kv.get_nks(); ++ik)
    {
        this->pelec->print_band(ik, PARAM.inp.printe, iter);
    }

    // 5) what's the exd used for?
#ifdef __EXX
    if (GlobalC::exx_info.info_ri.real_number)
    {
        this->exd->exx_hamilt2density(*this->pelec, this->ParaV, iter);
    }
    else
    {
        this->exc->exx_hamilt2density(*this->pelec, this->ParaV, iter);
    }
#endif

    // 6) calculate the local occupation number matrix and energy correction in
    // DFT+U
    if (GlobalV::dft_plus_u)
    {
        // only old DFT+U method should calculated energy correction in esolver,
        // new DFT+U method will calculate energy in calculating Hamiltonian
        if (GlobalV::dft_plus_u == 2)
        {
            if (GlobalC::dftu.omc != 2)
            {
                const std::vector<std::vector<TK>>& tmp_dm
                    = dynamic_cast<elecstate::ElecStateLCAO<TK>*>(this->pelec)->get_DM()->get_DMK_vector();
                this->dftu_cal_occup_m(iter, tmp_dm);
            }
            GlobalC::dftu.cal_energy_correction(istep);
        }
        GlobalC::dftu.output();
    }

    // (7) for deepks, calculate delta_e
#ifdef __DEEPKS
    if (GlobalV::deepks_scf)
    {
        const std::vector<std::vector<TK>>& dm
            = dynamic_cast<const elecstate::ElecStateLCAO<TK>*>(this->pelec)->get_DM()->get_DMK_vector();

        this->dpks_cal_e_delta_band(dm);
    }
#endif

    // 8) for delta spin
    if (PARAM.inp.sc_mag_switch)
    {
        SpinConstrain<TK, base_device::DEVICE_CPU>& sc = SpinConstrain<TK, base_device::DEVICE_CPU>::getScInstance();
        sc.cal_MW(iter, this->p_hamilt);
    }

    // 9) use new charge density to calculate energy
    this->pelec->cal_energies(1);

    // 10) symmetrize the charge density
    Symmetry_rho srho;
    for (int is = 0; is < GlobalV::NSPIN; is++)
    {
        srho.begin(is, *(this->pelec->charge), this->pw_rho, GlobalC::Pgrid, GlobalC::ucell.symm);
    }

    // 11) compute magnetization, only for spin==2
    GlobalC::ucell.magnet.compute_magnetization(this->pelec->charge->nrxx,
                                                this->pelec->charge->nxyz,
                                                this->pelec->charge->rho,
                                                this->pelec->nelec_spin.data());

    // 12) calculate delta energy
    this->pelec->f_en.deband = this->pelec->cal_delta_eband();
}

//------------------------------------------------------------------------------
//! the 12th function of ESolver_KS_LCAO: update_pot
//! mohan add 2024-05-11
//! 1) print Hamiltonian and Overlap matrix (why related to update_pot()?)
//! 2) print wavefunctions (why related to update_pot()?)
//! 3) print potential
//------------------------------------------------------------------------------
template <typename TK, typename TR>
void ESolver_KS_LCAO<TK, TR>::update_pot(const int istep, const int iter)
{
    ModuleBase::TITLE("ESolver_KS_LCAO", "update_pot");

    // 1) print Hamiltonian and Overlap matrix
    if (this->conv_elec || iter == GlobalV::SCF_NMAX)
    {
        if (!GlobalV::GAMMA_ONLY_LOCAL && (hsolver::HSolverLCAO<TK>::out_mat_hs[0] || GlobalV::deepks_v_delta))
        {
            this->GK.renew(true);
        }
        for (int ik = 0; ik < this->kv.get_nks(); ++ik)
        {
            if (hsolver::HSolverLCAO<TK>::out_mat_hs[0]|| GlobalV::deepks_v_delta)
            {
                this->p_hamilt->updateHk(ik);
            }
            bool bit = false; // LiuXh, 2017-03-21
            // if set bit = true, there would be error in soc-multi-core
            // calculation, noted by zhengdy-soc
            if (this->psi != nullptr && (istep % PARAM.inp.out_interval == 0))
            {
                hamilt::MatrixBlock<TK> h_mat, s_mat;
                this->p_hamilt->matrix(h_mat, s_mat);
                if (hsolver::HSolverLCAO<TK>::out_mat_hs[0])
                {
                    ModuleIO::save_mat(istep,
                                       h_mat.p,
                                       GlobalV::NLOCAL,
                                       bit,
                                       hsolver::HSolverLCAO<TK>::out_mat_hs[1],
                                       1,
                                       GlobalV::out_app_flag,
                                       "H",
                                       "data-" + std::to_string(ik),
                                       this->ParaV,
                                       GlobalV::DRANK);
                    ModuleIO::save_mat(istep,
                                       s_mat.p,
                                       GlobalV::NLOCAL,
                                       bit,
                                       hsolver::HSolverLCAO<TK>::out_mat_hs[1],
                                       1,
                                       GlobalV::out_app_flag,
                                       "S",
                                       "data-" + std::to_string(ik),
                                       this->ParaV,
                                       GlobalV::DRANK);
                }
#ifdef __DEEPKS
                if(GlobalV::deepks_v_delta)
                {
                    GlobalC::ld.save_h_mat(h_mat.p,this->ParaV.nloc);
                }
#endif
            }
        }
    }

    // 2) print wavefunctions
    if (elecstate::ElecStateLCAO<TK>::out_wfc_lcao && (this->conv_elec || iter == GlobalV::SCF_NMAX)
        && (istep % PARAM.inp.out_interval == 0))
    {
        ModuleIO::write_wfc_nao(elecstate::ElecStateLCAO<TK>::out_wfc_lcao,
                                this->psi[0],
                                this->pelec->ekb,
                                this->pelec->wg,
                                this->pelec->klist->kvec_c,
                                this->ParaV,
                                istep);
    }

    // 3) print potential
    if (this->conv_elec || iter == GlobalV::SCF_NMAX)
    {
        if (GlobalV::out_pot < 0) // mohan add 2011-10-10
        {
            GlobalV::out_pot = -2;
        }
    }

    if (!this->conv_elec)
    {
        if (GlobalV::NSPIN == 4)
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
}

//------------------------------------------------------------------------------
//! the 13th function of ESolver_KS_LCAO: iter_finish
//! mohan add 2024-05-11
//! 1) mix density matrix
//! 2) output charge density
//! 3) output exx matrix
//! 4) output charge density and density matrix
//! 5) cal_MW? (why put it here?)
//! 6) calculate the total energy?
//------------------------------------------------------------------------------
template <typename TK, typename TR>
void ESolver_KS_LCAO<TK, TR>::iter_finish(int iter)
{
    ModuleBase::TITLE("ESolver_KS_LCAO", "iter_finish");

    // 1) mix density matrix if mixing_restart + mixing_dmr + not first
    // mixing_restart at every iter
    if (GlobalV::MIXING_RESTART > 0 && this->p_chgmix->mixing_restart_count > 0 && GlobalV::MIXING_DMR)
    {
        elecstate::DensityMatrix<TK, double>* dm = dynamic_cast<elecstate::ElecStateLCAO<TK>*>(this->pelec)->get_DM();
        this->p_chgmix->mix_dmr(dm);
    }

    // 2) save charge density
    // Peize Lin add 2020.04.04
    if (GlobalC::restart.info_save.save_charge)
    {
        for (int is = 0; is < GlobalV::NSPIN; ++is)
        {
            GlobalC::restart.save_disk("charge", is, this->pelec->charge->nrxx, this->pelec->charge->rho[is]);
        }
    }

#ifdef __EXX
    // 3) save exx matrix
    int two_level_step = GlobalC::exx_info.info_ri.real_number ? this->exd->two_level_step : this->exc->two_level_step;

    if (GlobalC::restart.info_save.save_H && two_level_step > 0
        && (!GlobalC::exx_info.info_global.separate_loop || iter == 1)) // to avoid saving the same value repeatedly
    {
        ////////// for Add_Hexx_Type::k
        /*
        hamilt::HS_Matrix_K<TK> Hexxk_save(&this->ParaV, 1);
        for (int ik = 0; ik < this->kv.get_nks(); ++ik) {
            Hexxk_save.set_zero_hk();

            hamilt::OperatorEXX<hamilt::OperatorLCAO<TK, TR>> opexx_save(&Hexxk_save,
                                                                         nullptr, 
                                                                         this->kv);

            opexx_save.contributeHk(ik);

            GlobalC::restart.save_disk("Hexx",
                                       ik,
                                       this->ParaV.get_local_size(),
                                       Hexxk_save.get_hk());
        }*/
        ////////// for Add_Hexx_Type:R
        const std::string& restart_HR_path = GlobalC::restart.folder + "HexxR" + std::to_string(GlobalV::MY_RANK);
        if (GlobalC::exx_info.info_ri.real_number)
        {
            ModuleIO::write_Hexxs_csr(restart_HR_path, GlobalC::ucell, this->exd->get_Hexxs());
        }
        else
        {
            ModuleIO::write_Hexxs_csr(restart_HR_path, GlobalC::ucell, this->exc->get_Hexxs());
        }
        if (GlobalV::MY_RANK == 0)
        {
            GlobalC::restart.save_disk("Eexx", 0, 1, &this->pelec->f_en.exx);
        }
    }
#endif

    // 4) output charge density and density matrix
    bool print = false;
    if (this->out_freq_elec && iter % this->out_freq_elec == 0)
    {
        print = true;
    }

    if (print)
    {
        for (int is = 0; is < GlobalV::NSPIN; is++)
        {
            this->create_Output_Rho(is, iter, "tmp_").write();
            if (XC_Functional::get_func_type() == 3 || XC_Functional::get_func_type() == 5)
            {
                this->create_Output_Kin(is, iter, "tmp_").write();
            }
        }
    }

    // 5) cal_MW?
    // escon: energy of spin constraint depends on Mi, so cal_energies should be
    // called after cal_MW
    if (PARAM.inp.sc_mag_switch)
    {
        SpinConstrain<TK, base_device::DEVICE_CPU>& sc = SpinConstrain<TK, base_device::DEVICE_CPU>::getScInstance();
        sc.cal_MW(iter, this->p_hamilt);
    }

    // 6) calculate the total energy.
    this->pelec->cal_energies(2);
}

//------------------------------------------------------------------------------
//! the 14th function of ESolver_KS_LCAO: after_scf
//! mohan add 2024-05-11
//! 1) write charge difference into files for charge extrapolation
//! 2) write density matrix for sparse matrix
//! 3) write charge density
//! 4) write density matrix
//! 5) write Vxc
//! 6) write Exx matrix
//! 7) write potential
//! 8) write convergence
//! 9) write fermi energy
//! 10) write eigenvalues
//! 11) write deepks information
//! 12) write rpa information
//! 13) write HR in npz format
//! 14) write dm in npz format
//! 15) write md related
//! 16) write spin constrian MW?
//! 17) delete grid
//! 18) write quasi-orbitals
//------------------------------------------------------------------------------
template <typename TK, typename TR>
void ESolver_KS_LCAO<TK, TR>::after_scf(const int istep)
{
    ModuleBase::TITLE("ESolver_KS_LCAO", "after_scf");

    // 1) write charge difference into files for charge extrapolation
    if (GlobalV::CALCULATION != "scf")
    {
        this->CE.save_files(istep,
                            GlobalC::ucell,
#ifdef __MPI
                            this->pw_big,
#endif
                            this->pelec->charge,
                            &this->sf);
    }

    // 2) write density matrix for sparse matrix
    ModuleIO::write_dmr(dynamic_cast<const elecstate::ElecStateLCAO<TK>*>(this->pelec)->get_DM()->get_DMR_vector(),
                        this->ParaV,
                        PARAM.inp.out_dm1,
                        false,
                        GlobalV::out_app_flag,
                        istep);

    // 3) write charge density
    if (PARAM.inp.out_chg)
    {
        for (int is = 0; is < GlobalV::NSPIN; is++)
        {
            this->create_Output_Rho(is, istep).write();
            if (XC_Functional::get_func_type() == 3 || XC_Functional::get_func_type() == 5)
            {
                this->create_Output_Kin(is, istep).write();
            }
        }
    }

    // 4) write density matrix
    if (PARAM.inp.out_dm)
    {
        std::vector<double> efermis(GlobalV::NSPIN == 2 ? 2 : 1);
        for (int ispin = 0; ispin < efermis.size(); ispin++)
        {
            efermis[ispin] = this->pelec->eferm.get_efval(ispin);
        }
        const int precision = 3;
        ModuleIO::write_dmk(dynamic_cast<const elecstate::ElecStateLCAO<TK>*>(this->pelec)->get_DM()->get_DMK_vector(),
                            precision,
                            efermis,
                            &(GlobalC::ucell),
                            this->ParaV);
    }

#ifdef __EXX
    // 5) write Exx matrix
    if (GlobalC::exx_info.info_global.cal_exx) // Peize Lin add if 2022.11.14
    {
        const std::string file_name_exx = GlobalV::global_out_dir + "HexxR"
                                          + std::to_string(GlobalV::MY_RANK);
        if (GlobalC::exx_info.info_ri.real_number) {
            ModuleIO::write_Hexxs_csr(file_name_exx, GlobalC::ucell, this->exd->get_Hexxs());
        } else {
            ModuleIO::write_Hexxs_csr(file_name_exx, GlobalC::ucell, this->exc->get_Hexxs());
        }
    }
#endif

    // 6) write potential
    this->create_Output_Potential(istep).write();

    // 7) write convergence
    ModuleIO::output_convergence_after_scf(this->conv_elec, this->pelec->f_en.etot);

    // 8) write fermi energy
    ModuleIO::output_efermi(this->conv_elec, this->pelec->eferm.ef);

    // 9) write eigenvalues
    if (GlobalV::OUT_LEVEL != "m")
    {
        this->pelec->print_eigenvalue(GlobalV::ofs_running);
    }

    // 10) write deepks information
#ifdef __DEEPKS
    std::shared_ptr<LCAO_Deepks> ld_shared_ptr(&GlobalC::ld, [](LCAO_Deepks*) {});
    LCAO_Deepks_Interface LDI = LCAO_Deepks_Interface(ld_shared_ptr);
    ModuleBase::timer::tick("ESolver_KS_LCAO", "out_deepks_labels");
    LDI.out_deepks_labels(this->pelec->f_en.etot,
                          this->pelec->klist->get_nks(),
                          GlobalC::ucell.nat,
                          GlobalV::NLOCAL,
                          this->pelec->ekb,
                          this->pelec->klist->kvec_d,
                          GlobalC::ucell,
                          GlobalC::ORB,
                          GlobalC::GridD,
                          &(this->ParaV),
                          *(this->psi),
                          dynamic_cast<const elecstate::ElecStateLCAO<TK>*>(this->pelec)->get_DM(),
                          GlobalV::deepks_v_delta);

    ModuleBase::timer::tick("ESolver_KS_LCAO", "out_deepks_labels");
#endif

#ifdef __EXX
    // 11) write rpa information
    if (PARAM.inp.rpa)
    {
        // ModuleRPA::DFT_RPA_interface
        // rpa_interface(GlobalC::exx_info.info_global);
        // rpa_interface.rpa_exx_lcao().info.files_abfs = GlobalV::rpa_orbitals;
        RPA_LRI<TK, double> rpa_lri_double(GlobalC::exx_info.info_ri);
        rpa_lri_double.cal_postSCF_exx(*dynamic_cast<const elecstate::ElecStateLCAO<TK>*>(this->pelec)->get_DM(),
                                       MPI_COMM_WORLD,
                                       this->kv);
        rpa_lri_double.init(MPI_COMM_WORLD, this->kv);
        rpa_lri_double.out_for_RPA(this->ParaV, *(this->psi), this->pelec);
    }
#endif

    // 12) write HR in npz format
    if (PARAM.inp.out_hr_npz)
    {
        this->p_hamilt->updateHk(0); // first k point, up spin
        hamilt::HamiltLCAO<std::complex<double>, double>* p_ham_lcao
            = dynamic_cast<hamilt::HamiltLCAO<std::complex<double>, double>*>(this->p_hamilt);
        std::string zipname = "output_HR0.npz";
        this->output_mat_npz(zipname, *(p_ham_lcao->getHR()));

        if (GlobalV::NSPIN == 2)
        {
            this->p_hamilt->updateHk(this->kv.get_nks() / 2); // the other half of k points, down spin
            hamilt::HamiltLCAO<std::complex<double>, double>* p_ham_lcao
                = dynamic_cast<hamilt::HamiltLCAO<std::complex<double>, double>*>(this->p_hamilt);
            zipname = "output_HR1.npz";
            this->output_mat_npz(zipname, *(p_ham_lcao->getHR()));
        }
    }

    // 13) write dm in npz format
    if (PARAM.inp.out_dm_npz)
    {
        const elecstate::DensityMatrix<TK, double>* dm
            = dynamic_cast<const elecstate::ElecStateLCAO<TK>*>(this->pelec)->get_DM();
        std::string zipname = "output_DM0.npz";
        this->output_mat_npz(zipname, *(dm->get_DMR_pointer(1)));

        if (GlobalV::NSPIN == 2)
        {
            zipname = "output_DM1.npz";
            this->output_mat_npz(zipname, *(dm->get_DMR_pointer(2)));
        }
    }

    // 14) write md related
    if (!md_skip_out(GlobalV::CALCULATION, istep, PARAM.inp.out_interval))
    {
        this->create_Output_Mat_Sparse(istep).write();
        // mulliken charge analysis
        if (PARAM.inp.out_mul)
        {
            this->cal_mag(istep, true);
        }
    }

    // 15) write spin constrian MW?
    // spin constrain calculations, added by Tianqi Zhao.
    if (PARAM.inp.sc_mag_switch) {
        SpinConstrain<TK, base_device::DEVICE_CPU>& sc
            = SpinConstrain<TK, base_device::DEVICE_CPU>::getScInstance();
        sc.cal_MW(istep, true);
        sc.print_Mag_Force();
    }

    // 16) delete grid
    if (!GlobalV::CAL_FORCE && !GlobalV::CAL_STRESS)
    {
        RA.delete_grid();
    }

    // 17) write quasi-orbitals, added by Yike Huang.
    if (PARAM.inp.qo_switch)
    {
        toQO tqo(PARAM.inp.qo_basis, PARAM.inp.qo_strategy, GlobalV::qo_thr, GlobalV::qo_screening_coeff);
        tqo.initialize(GlobalV::global_out_dir,
                       PARAM.inp.pseudo_dir,
                       PARAM.inp.orbital_dir,
                       &GlobalC::ucell,
                       this->kv.kvec_d,
                       GlobalV::ofs_running,
                       GlobalV::MY_RANK,
                       GlobalV::NPROC);
        tqo.calculate();
    }
}

//------------------------------------------------------------------------------
//! the 15th function of ESolver_KS_LCAO: do_after_converge
//! mohan add 2024-05-11
//------------------------------------------------------------------------------
template <typename TK, typename TR>
bool ESolver_KS_LCAO<TK, TR>::do_after_converge(int& iter)
{
    ModuleBase::TITLE("ESolver_KS_LCAO", "do_after_converge");

    if (GlobalV::dft_plus_u)
    {
        // use the converged occupation matrix for next MD/Relax SCF calculation
        GlobalC::dftu.initialed_locale = true;
    }
    // FIXME: for developer who want to test restarting DeePKS with same Descriptor/PDM in last MD step
    // RUN: " GlobalC::ld.set_init_pdm(true); " can skip the calculation of PDM in the next iter_init

#ifdef __EXX
    if (GlobalC::exx_info.info_global.cal_exx)
    {
        if (GlobalC::exx_info.info_ri.real_number) {
            return this->exd->exx_after_converge(
                *this->p_hamilt,
                *dynamic_cast<const elecstate::ElecStateLCAO<TK>*>(this->pelec)
                ->get_DM(),
                this->kv,
                iter);
        }
        else {
            return this->exc->exx_after_converge(
                *this->p_hamilt,
                *dynamic_cast<const elecstate::ElecStateLCAO<TK>*>(this->pelec)
                ->get_DM(),
                this->kv,
                iter);
        }
    }
#endif // __EXX

    return true;
}

//------------------------------------------------------------------------------
//! the 16th function of ESolver_KS_LCAO: create_Output_DM
//! mohan add 2024-05-11
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//! the 17th function of ESolver_KS_LCAO: create_Output_DM1
//! mohan add 2024-05-11
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//! the 18th function of ESolver_KS_LCAO: create_Output_Mat_Sparse
//! mohan add 2024-05-11
//------------------------------------------------------------------------------
template <typename TK, typename TR>
ModuleIO::Output_Mat_Sparse<TK> ESolver_KS_LCAO<TK, TR>::create_Output_Mat_Sparse(int istep)
{
    return ModuleIO::Output_Mat_Sparse<TK>(hsolver::HSolverLCAO<TK>::out_mat_hsR,
        hsolver::HSolverLCAO<TK>::out_mat_dh,
        hsolver::HSolverLCAO<TK>::out_mat_t,
        PARAM.inp.out_mat_r,
        istep,
        this->pelec->pot->get_effective_v(),
        this->ParaV,
        this->GK, // mohan add 2024-04-01
        two_center_bundle_,
        GlobalC::GridD, // mohan add 2024-04-06
        this->kv,
        this->p_hamilt);
}

//------------------------------------------------------------------------------
//! the 19th function of ESolver_KS_LCAO: md_skip_out
//! mohan add 2024-05-11
//------------------------------------------------------------------------------
template <typename TK, typename TR>
bool ESolver_KS_LCAO<TK, TR>::md_skip_out(std::string calculation, int istep, int interval)
{
    if (calculation == "md")
    {
        if (istep % interval != 0)
        {
            return true;
        }
    }
    return false;
}

template <typename TK, typename TR>
void ESolver_KS_LCAO<TK, TR>::cal_mag(const int istep, const bool print)
{
    auto cell_index = CellIndex(GlobalC::ucell.get_atomLabels(),
                                GlobalC::ucell.get_atomCounts(),
                                GlobalC::ucell.get_lnchiCounts(),
                                GlobalV::NSPIN);
    auto out_sk = ModuleIO::Output_Sk<TK>(this->p_hamilt,
                                          &(this->ParaV),
                                          GlobalV::NSPIN,
                                          this->kv.get_nks());
    auto out_dmk = ModuleIO::Output_DMK<TK>(dynamic_cast<const elecstate::ElecStateLCAO<TK>*>(this->pelec)->get_DM(),
                                            &(this->ParaV),
                                            GlobalV::NSPIN,
                                            this->kv.get_nks());
    auto mulp = ModuleIO::Output_Mulliken<TK>(&(out_sk),
                                              &(out_dmk),
                                              &(this->ParaV),
                                              &cell_index,
                                              this->kv.isk,
                                              GlobalV::NSPIN);
    auto atom_chg = mulp.get_atom_chg();
    /// used in updating mag info in STRU file
    GlobalC::ucell.atom_mulliken = mulp.get_atom_mulliken(atom_chg);
    if (print && GlobalV::MY_RANK == 0)
    {
        /// write the Orbital file
        cell_index.write_orb_info(GlobalV::global_out_dir);
        /// write mulliken.txt
        mulp.write(istep, GlobalV::global_out_dir);
        /// write atomic mag info in running log file
        mulp.print_atom_mag(atom_chg, GlobalV::ofs_running);
    }
}

//------------------------------------------------------------------------------
//! the 20th,21th,22th functions of ESolver_KS_LCAO
//! mohan add 2024-05-11
//------------------------------------------------------------------------------
template class ESolver_KS_LCAO<double, double>;
template class ESolver_KS_LCAO<std::complex<double>, double>;
template class ESolver_KS_LCAO<std::complex<double>, std::complex<double>>;
} // namespace ModuleESolver
