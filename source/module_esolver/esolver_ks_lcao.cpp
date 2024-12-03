#include "esolver_ks_lcao.h"

#include "module_base/formatter.h"
#include "module_base/global_variable.h"
#include "module_base/tool_title.h"
#include "module_elecstate/module_dm/cal_dm_psi.h"
#include "module_hamilt_lcao/module_deltaspin/spin_constrain.h"
#include "module_hamilt_lcao/module_dftu/dftu.h"
#include "module_io/berryphase.h"
#include "module_io/cube_io.h"
#include "module_io/dos_nao.h"
#include "module_io/io_dmk.h"
#include "module_io/io_npz.h"
#include "module_io/nscf_band.h"
#include "module_io/output_dmk.h"
#include "module_io/output_log.h"
#include "module_io/output_mat_sparse.h"
#include "module_io/output_mulliken.h"
#include "module_io/output_sk.h"
#include "module_io/to_qo.h"
#include "module_io/to_wannier90_lcao.h"
#include "module_io/to_wannier90_lcao_in_pw.h"
#include "module_io/write_HS.h"
#include "module_io/write_dmr.h"
#include "module_io/write_eband_terms.hpp"
#include "module_io/write_elecstat_pot.h"
#include "module_io/write_istate_info.h"
#include "module_io/write_proj_band_lcao.h"
#include "module_io/write_vxc.hpp"
#include "module_io/write_wfc_nao.h"
#include "module_parameter/parameter.h"

//--------------temporary----------------------------
#include "module_base/global_function.h"
#include "module_cell/module_neighbor/sltk_grid_driver.h"
#include "module_elecstate/module_charge/symmetry_rho.h"
#include "module_elecstate/occupy.h"
#include "module_hamilt_lcao/hamilt_lcaodft/LCAO_domain.h" // need divide_HS_in_frag
#include "module_hamilt_lcao/hamilt_lcaodft/hs_matrix_k.hpp"
#include "module_hamilt_lcao/module_dftu/dftu.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_io/print_info.h"
#include "module_elecstate/cal_ux.h"

#include <memory>
#ifdef __EXX
#include "module_io/restart_exx_csr.h"
#include "module_ri/RPA_LRI.h"
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
// #include "module_elecstate/cal_dm.h"
//---------------------------------------------------

// test RDMFT
#include "module_rdmft/rdmft.h"
#include <iostream>

namespace ModuleESolver
{

//------------------------------------------------------------------------------
//! the 1st function of ESolver_KS_LCAO: constructor
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
//------------------------------------------------------------------------------
template <typename TK, typename TR>
ESolver_KS_LCAO<TK, TR>::~ESolver_KS_LCAO()
{
}

//------------------------------------------------------------------------------
//! the 3rd function of ESolver_KS_LCAO: init
//! 1) calculate overlap matrix S or initialize
//! 2) init ElecState
//! 3) init LCAO basis
//! 4) initialize the density matrix
//! 5) initialize Hamilt in LCAO
//! 6) initialize exx
//! 7) initialize DFT+U
//! 8) ppcell
//! 9) inititlize the charge density
//! 10) initialize the potential.
//! 11) initialize deepks
//! 12) set occupations
//! 13) print a warning if needed
//------------------------------------------------------------------------------
template <typename TK, typename TR>
void ESolver_KS_LCAO<TK, TR>::before_all_runners(UnitCell& ucell, const Input_para& inp)
{
    ModuleBase::TITLE("ESolver_KS_LCAO", "before_all_runners");
    ModuleBase::timer::tick("ESolver_KS_LCAO", "before_all_runners");

    ESolver_KS<TK>::before_all_runners(ucell, inp);

    // 2) init ElecState
    // autoset nbands in ElecState, it should before basis_init (for Psi 2d division)
    if (this->pelec == nullptr)
    {
        // TK stands for double and complex<double>?
        this->pelec = new elecstate::ElecStateLCAO<TK>(&(this->chr), // use which parameter?
                                                       &(this->kv),
                                                       this->kv.get_nks(),
                                                       &(this->GG), // mohan add 2024-04-01
                                                       &(this->GK), // mohan add 2024-04-01
                                                       this->pw_rho,
                                                       this->pw_big);
    }

    // 3) init LCAO basis
    // reading the localized orbitals/projectors
    // construct the interpolation tables.
    LCAO_domain::init_basis_lcao(this->pv,
                                 inp.onsite_radius,
                                 inp.lcao_ecut,
                                 inp.lcao_dk,
                                 inp.lcao_dr,
                                 inp.lcao_rmax,
                                 ucell,
                                 two_center_bundle_,
                                 orb_);

    // 4) initialize the density matrix
    // DensityMatrix is allocated here, DMK is also initialized here
    // DMR is not initialized here, it will be constructed in each before_scf
    dynamic_cast<elecstate::ElecStateLCAO<TK>*>(this->pelec)->init_DM(&this->kv, &(this->pv), PARAM.inp.nspin);

    // 5) initialize Hamilt in LCAO
    // * allocate H and S matrices according to computational resources
    // * set the 'trace' between local H/S and global H/S
    LCAO_domain::divide_HS_in_frag(PARAM.globalv.gamma_only_local, ucell , pv, this->kv.get_nks(), orb_);

#ifdef __EXX
    // 6) initialize exx
    // PLEASE simplify the Exx_Global interface
    if (PARAM.inp.calculation == "scf" || PARAM.inp.calculation == "relax" || PARAM.inp.calculation == "cell-relax"
        || PARAM.inp.calculation == "md")
    {
        if (GlobalC::exx_info.info_global.cal_exx)
        {
            XC_Functional::set_xc_first_loop(ucell);
            // initialize 2-center radial tables for EXX-LRI
            if (GlobalC::exx_info.info_ri.real_number)
            {
                this->exx_lri_double->init(MPI_COMM_WORLD, this->kv, orb_);
                this->exd->exx_before_all_runners(this->kv, ucell, this->pv);
            }
            else
            {
                this->exx_lri_complex->init(MPI_COMM_WORLD, this->kv, orb_);
                this->exc->exx_before_all_runners(this->kv, ucell, this->pv);
            }
        }
    }
#endif

    // 7) initialize DFT+U
    if (PARAM.inp.dft_plus_u)
    {
        GlobalC::dftu.init(ucell, &this->pv, this->kv.get_nks(), orb_);
    }

    // 8) initialize ppcell
    this->ppcell.init_vloc(ucell,this->pw_rho);
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "LOCAL POTENTIAL");

    // 9) inititlize the charge density
    this->pelec->charge->allocate(PARAM.inp.nspin);
    this->pelec->omega = ucell.omega;

    // 10) initialize the potential
    if (this->pelec->pot == nullptr)
    {
        this->pelec->pot = new elecstate::Potential(this->pw_rhod,
                                                    this->pw_rho,
                                                    &ucell,
                                                    &(this->ppcell.vloc),
                                                    &(this->sf),
                                                    &(this->pelec->f_en.etxc),
                                                    &(this->pelec->f_en.vtxc));
    }

#ifdef __DEEPKS
    // 11) initialize deepks
    if (PARAM.inp.deepks_scf)
    {
        // load the DeePKS model from deep neural network
        GlobalC::ld.load_model(PARAM.inp.deepks_model);
        // read pdm from file for NSCF or SCF-restart, do it only once in whole calculation
        GlobalC::ld.read_projected_DM((PARAM.inp.init_chg == "file"), PARAM.inp.deepks_equiv, *orb_.Alpha);
    }
#endif

    // 12) set occupations
    // tddft does not need to set occupations in the first scf
    if (PARAM.inp.ocp && inp.esolver_type != "tddft")
    {
        this->pelec->fixed_weights(PARAM.inp.ocp_kb, PARAM.inp.nbands, PARAM.inp.nelec);
    }

    // 13) if kpar is not divisible by nks, print a warning
    if (GlobalV::KPAR_LCAO > 1)
    {
        if (this->kv.get_nks() % GlobalV::KPAR_LCAO != 0)
        {
            ModuleBase::WARNING("ESolver_KS_LCAO::before_all_runners", "nks is not divisible by kpar.");
            std::cout << "\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
                         "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
                         "%%%%%%%%%%%%%%%%%%%%%%%%%%"
                      << std::endl;
            std::cout << " Warning: nks (" << this->kv.get_nks() << ") is not divisible by kpar (" << GlobalV::KPAR_LCAO
                      << ")." << std::endl;
            std::cout << " This may lead to poor load balance. It is strongly suggested to" << std::endl;
            std::cout << " set nks to be divisible by kpar, but if this is really what" << std::endl;
            std::cout << " you want, please ignore this warning." << std::endl;
            std::cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
                         "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
                         "%%%%%%%%%%%%\n";
        }
    }

    // 14) initialize rdmft, added by jghan
    if( PARAM.inp.rdmft == true )
    {
        rdmft_solver.init( this->GG, this->GK, this->pv, ucell, this->kv, *(this->pelec),
                                this->orb_, two_center_bundle_, PARAM.inp.dft_functional, PARAM.inp.rdmft_power_alpha);
    }

    ModuleBase::timer::tick("ESolver_KS_LCAO", "before_all_runners");
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
void ESolver_KS_LCAO<TK, TR>::cal_force(UnitCell& ucell, ModuleBase::matrix& force)
{
    ModuleBase::TITLE("ESolver_KS_LCAO", "cal_force");
    ModuleBase::timer::tick("ESolver_KS_LCAO", "cal_force");

    Force_Stress_LCAO<TK> fsl(this->RA, ucell.nat);

    fsl.getForceStress(PARAM.inp.cal_force,
                       PARAM.inp.cal_stress,
                       PARAM.inp.test_force,
                       PARAM.inp.test_stress,
                       ucell,
                       this->pv,
                       this->pelec,
                       this->psi,
                       this->GG, // mohan add 2024-04-01
                       this->GK, // mohan add 2024-04-01
                       two_center_bundle_,
                       orb_,
                       force,
                       this->scs,
                       this->ppcell,
                       this->sf,
                       this->kv,
                       this->pw_rho,
#ifdef __EXX
                       *this->exx_lri_double,
                       *this->exx_lri_complex,
#endif
                       &ucell.symm);

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
void ESolver_KS_LCAO<TK, TR>::cal_stress(UnitCell& ucell, ModuleBase::matrix& stress)
{
    ModuleBase::TITLE("ESolver_KS_LCAO", "cal_stress");
    ModuleBase::timer::tick("ESolver_KS_LCAO", "cal_stress");

    if (!this->have_force)
    {
        ModuleBase::matrix fcs;
        this->cal_force(ucell, fcs);
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
void ESolver_KS_LCAO<TK, TR>::after_all_runners(UnitCell& ucell)
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
    if (PARAM.inp.calculation == "scf" || PARAM.inp.calculation == "md" || PARAM.inp.calculation == "relax")
    {
        ModuleIO::write_istate_info(this->pelec->ekb, this->pelec->wg, this->kv, &(GlobalC::Pkpoints));
    }

    const int nspin0 = (PARAM.inp.nspin == 2) ? 2 : 1;

    if (PARAM.inp.out_band[0])
    {
        for (int is = 0; is < nspin0; is++)
        {
            std::stringstream ss2;
            ss2 << PARAM.globalv.global_out_dir << "BANDS_" << is + 1 << ".dat";
            GlobalV::ofs_running << "\n Output bands in file: " << ss2.str() << std::endl;
            ModuleIO::nscf_band(is,
                                ss2.str(),
                                PARAM.inp.nbands,
                                0.0,
                                PARAM.inp.out_band[1],
                                this->pelec->ekb,
                                this->kv,
                                &(GlobalC::Pkpoints));
        }
    } // out_band

    if (PARAM.inp.out_proj_band) // Projeced band structure added by jiyy-2022-4-20
    {
        ModuleIO::write_proj_band_lcao(this->psi, this->pv, this->pelec, this->kv, ucell, this->p_hamilt);
    }

    if (PARAM.inp.out_dos)
    {
        ModuleIO::out_dos_nao(this->psi,
                              this->pv,
                              this->pelec->ekb,
                              this->pelec->wg,
                              PARAM.inp.dos_edelta_ev,
                              PARAM.inp.dos_scale,
                              PARAM.inp.dos_sigma,
                              *(this->pelec->klist),
                              GlobalC::Pkpoints,
                              ucell,
                              this->pelec->eferm,
                              PARAM.inp.nbands,
                              this->p_hamilt);
    }

    if (PARAM.inp.out_mat_xc)
    {
        ModuleIO::write_Vxc<TK, TR>(PARAM.inp.nspin,
                                    PARAM.globalv.nlocal,
                                    GlobalV::DRANK,
                                    &this->pv,
                                    *this->psi,
                                    ucell,
                                    this->sf,
                                    *this->pw_rho,
                                    *this->pw_rhod,
                                    this->ppcell.vloc,
                                    *this->pelec->charge,
                                    this->GG,
                                    this->GK,
                                    this->kv,
                                    orb_.cutoffs(),
                                    this->pelec->wg,
                                    GlobalC::GridD
#ifdef __EXX
                                    ,
                                    this->exx_lri_double ? &this->exx_lri_double->Hexxs : nullptr,
                                    this->exx_lri_complex ? &this->exx_lri_complex->Hexxs : nullptr
#endif
        );
    }

    if (PARAM.inp.out_eband_terms)
    {
        ModuleIO::write_eband_terms<TK, TR>(PARAM.inp.nspin,
                                            PARAM.globalv.nlocal,
                                            GlobalV::DRANK,
                                            &this->pv,
                                            *this->psi,
                                            ucell,
                                            this->sf,
                                            *this->pw_rho,
                                            *this->pw_rhod,
                                            this->ppcell.vloc,
                                            *this->pelec->charge,
                                            this->GG,
                                            this->GK,
                                            this->kv,
                                            this->pelec->wg,
                                            GlobalC::GridD,
                                            orb_.cutoffs(),
                                            this->two_center_bundle_
#ifdef __EXX
                                            ,
                                            this->exx_lri_double ? &this->exx_lri_double->Hexxs : nullptr,
                                            this->exx_lri_complex ? &this->exx_lri_complex->Hexxs : nullptr
#endif
        );
    }

    ModuleBase::timer::tick("ESolver_KS_LCAO", "after_all_runners");
}

//------------------------------------------------------------------------------
//! the 10th function of ESolver_KS_LCAO: iter_init
//! mohan add 2024-05-11
//------------------------------------------------------------------------------
template <typename TK, typename TR>
void ESolver_KS_LCAO<TK, TR>::iter_init(UnitCell& ucell, const int istep, const int iter)
{
    ModuleBase::TITLE("ESolver_KS_LCAO", "iter_init");

    // call iter_init() of ESolver_KS
    ESolver_KS<TK>::iter_init(ucell, istep, iter);

    if (iter == 1)
    {
        this->p_chgmix->init_mixing(); // init mixing
        this->p_chgmix->mixing_restart_step = PARAM.inp.scf_nmax + 1;
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
    if (iter == this->p_chgmix->mixing_restart_step && PARAM.inp.mixing_restart > 0.0)
    {
        this->p_chgmix->init_mixing();
        this->p_chgmix->mixing_restart_count++;
        if (PARAM.inp.dft_plus_u)
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
                this->p_chgmix->mixing_restart_step = PARAM.inp.scf_nmax + 1;
            }
        }
        if (PARAM.inp.mixing_dmr) // for mixing_dmr
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
    if (istep == 0 && PARAM.inp.init_wfc == "file")
    {
        if (iter == 1)
        {
            std::cout << " WAVEFUN -> CHARGE " << std::endl;

            // calculate the density matrix using read in wave functions
            // and the ncalculate the charge density on grid.

            this->pelec->skip_weights = true;
            this->pelec->calculate_weights();
            if (!PARAM.inp.dm_to_rho)
            {
                auto _pelec = dynamic_cast<elecstate::ElecStateLCAO<TK>*>(this->pelec);
                _pelec->calEBand();
                elecstate::cal_dm_psi(_pelec->DM->get_paraV_pointer(), _pelec->wg, *this->psi, *(_pelec->DM));
                _pelec->DM->cal_DMR();
            }
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

            elecstate::cal_ux(ucell);

            //! update the potentials by using new electron charge density
            this->pelec->pot->update_from_charge(this->pelec->charge, &ucell);

            //! compute the correction energy for metals
            this->pelec->f_en.descf = this->pelec->cal_delta_escf();
        }
    }

#ifdef __EXX
    // calculate exact-exchange
    if (PARAM.inp.calculation != "nscf")
    {
        if (GlobalC::exx_info.info_ri.real_number)
        {
            this->exd->exx_eachiterinit(istep,
                                        *dynamic_cast<const elecstate::ElecStateLCAO<TK>*>(this->pelec)->get_DM(),
                                        this->kv,
                                        iter);
        }
        else
        {
            this->exc->exx_eachiterinit(istep,
                                        *dynamic_cast<const elecstate::ElecStateLCAO<TK>*>(this->pelec)->get_DM(),
                                        this->kv,
                                        iter);
        }
    }
#endif

    if (PARAM.inp.dft_plus_u)
    {
        if (istep != 0 || iter != 1)
        {
            GlobalC::dftu.set_dmr(dynamic_cast<elecstate::ElecStateLCAO<TK>*>(this->pelec)->get_DM());
        }
        // Calculate U and J if Yukawa potential is used
        GlobalC::dftu.cal_slater_UJ(ucell,this->pelec->charge->rho, this->pw_rho->nrxx);
    }

#ifdef __DEEPKS
    // the density matrixes of DeePKS have been updated in each iter
    GlobalC::ld.set_hr_cal(true);

    // HR in HamiltLCAO should be recalculate
    if (PARAM.inp.deepks_scf)
    {
        this->p_hamilt->refresh();
    }
#endif

    if (PARAM.inp.vl_in_h)
    {
        // update Gint_K
        if (!PARAM.globalv.gamma_only_local)
        {
            this->GK.renew();
        }
        // update real space Hamiltonian
        this->p_hamilt->refresh();
    }

    // save density matrix DMR for mixing
    if (PARAM.inp.mixing_restart > 0 && PARAM.inp.mixing_dmr && this->p_chgmix->mixing_restart_count > 0)
    {
        elecstate::DensityMatrix<TK, double>* dm = dynamic_cast<elecstate::ElecStateLCAO<TK>*>(this->pelec)->get_DM();
        dm->save_DMR();
    }
}

//------------------------------------------------------------------------------
//! the 11th function of ESolver_KS_LCAO: hamilt2density_single
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
void ESolver_KS_LCAO<TK, TR>::hamilt2density_single(UnitCell& ucell, int istep, int iter, double ethr)
{
    ModuleBase::TITLE("ESolver_KS_LCAO", "hamilt2density_single");

    // reset energy
    this->pelec->f_en.eband = 0.0;
    this->pelec->f_en.demet = 0.0;
    bool skip_charge = PARAM.inp.calculation == "nscf" ? true : false;

    // run the inner lambda loop to contrain atomic moments with the DeltaSpin method
    bool skip_solve = false;
    if (PARAM.inp.sc_mag_switch)
    {
        spinconstrain::SpinConstrain<TK>& sc = spinconstrain::SpinConstrain<TK>::getScInstance();
        if(!sc.mag_converged() && this->drho>0 && this->drho < PARAM.inp.sc_scf_thr)
        {
            // optimize lambda to get target magnetic moments, but the lambda is not near target
            sc.run_lambda_loop(iter-1);
            sc.set_mag_converged(true);
            skip_solve = true;
        }
        else if(sc.mag_converged())
        {
            // optimize lambda to get target magnetic moments, but the lambda is not near target
            sc.run_lambda_loop(iter-1);
            skip_solve = true;
        }
    }
    if(!skip_solve)
    {
        hsolver::HSolverLCAO<TK> hsolver_lcao_obj(&(this->pv), PARAM.inp.ks_solver);
        hsolver_lcao_obj.solve(this->p_hamilt, this->psi[0], this->pelec, skip_charge);
    }

    // 5) what's the exd used for?
#ifdef __EXX
    if (PARAM.inp.calculation != "nscf")
    {
        if (GlobalC::exx_info.info_ri.real_number)
        {
            this->exd->exx_hamilt2density(*this->pelec, this->pv, iter);
        }
        else
        {
            this->exc->exx_hamilt2density(*this->pelec, this->pv, iter);
        }
    }
#endif

    // 10) symmetrize the charge density
    Symmetry_rho srho;
    for (int is = 0; is < PARAM.inp.nspin; is++)
    {
        srho.begin(is, *(this->pelec->charge), this->pw_rho, ucell.symm);
    }

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
void ESolver_KS_LCAO<TK, TR>::update_pot(UnitCell& ucell, const int istep, const int iter)
{
    ModuleBase::TITLE("ESolver_KS_LCAO", "update_pot");

    // 1) print Hamiltonian and Overlap matrix
    if (this->conv_esolver || iter == PARAM.inp.scf_nmax)
    {
        if (!PARAM.globalv.gamma_only_local && (PARAM.inp.out_mat_hs[0] || PARAM.inp.deepks_v_delta))
        {
            this->GK.renew(true);
        }
        for (int ik = 0; ik < this->kv.get_nks(); ++ik)
        {
            if (PARAM.inp.out_mat_hs[0] || PARAM.inp.deepks_v_delta)
            {
                this->p_hamilt->updateHk(ik);
            }
            bool bit = false; // LiuXh, 2017-03-21
            // if set bit = true, there would be error in soc-multi-core
            // calculation, noted by zhengdy-soc
            if (this->psi != nullptr && (istep % PARAM.inp.out_interval == 0))
            {
                hamilt::MatrixBlock<TK> h_mat;
                hamilt::MatrixBlock<TK> s_mat;

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
#ifdef __DEEPKS
                if (PARAM.inp.deepks_out_labels && PARAM.inp.deepks_v_delta)
                {
                    DeePKS_domain::save_h_mat(h_mat.p, this->pv.nloc);
                }
#endif
            }
        }
    }

    // 2) print wavefunctions
    if (elecstate::ElecStateLCAO<TK>::out_wfc_lcao && (this->conv_esolver || iter == PARAM.inp.scf_nmax)
        && (istep % PARAM.inp.out_interval == 0))
    {
        ModuleIO::write_wfc_nao(elecstate::ElecStateLCAO<TK>::out_wfc_lcao,
                                this->psi[0],
                                this->pelec->ekb,
                                this->pelec->wg,
                                this->pelec->klist->kvec_c,
                                this->pv,
                                istep);
    }

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

}

//------------------------------------------------------------------------------
//! the 13th function of ESolver_KS_LCAO: iter_finish
//! mohan add 2024-05-11
//! 1) mix density matrix
//! 2) output charge density
//! 3) output exx matrix
//! 4) output charge density and density matrix
//------------------------------------------------------------------------------
template <typename TK, typename TR>
void ESolver_KS_LCAO<TK, TR>::iter_finish(UnitCell& ucell, const int istep, int& iter)
{
    ModuleBase::TITLE("ESolver_KS_LCAO", "iter_finish");

    // 6) calculate the local occupation number matrix and energy correction in
    // DFT+U
    if (PARAM.inp.dft_plus_u)
    {
        // only old DFT+U method should calculated energy correction in esolver,
        // new DFT+U method will calculate energy in calculating Hamiltonian
        if (PARAM.inp.dft_plus_u == 2)
        {
            if (GlobalC::dftu.omc != 2)
            {
                const std::vector<std::vector<TK>>& tmp_dm
                    = dynamic_cast<elecstate::ElecStateLCAO<TK>*>(this->pelec)->get_DM()->get_DMK_vector();
                ModuleDFTU::dftu_cal_occup_m(iter, ucell,tmp_dm, this->kv, this->p_chgmix->get_mixing_beta(), this->p_hamilt);
            }
            GlobalC::dftu.cal_energy_correction(ucell,istep);
        }
        GlobalC::dftu.output(ucell);
    }

    // (7) for deepks, calculate delta_e
#ifdef __DEEPKS
    if (PARAM.inp.deepks_scf)
    {
        const std::vector<std::vector<TK>>& dm
            = dynamic_cast<const elecstate::ElecStateLCAO<TK>*>(this->pelec)->get_DM()->get_DMK_vector();

        GlobalC::ld.dpks_cal_e_delta_band(dm, this->kv.get_nks());
    }
#endif

    // 8) for delta spin
    if (PARAM.inp.sc_mag_switch)
    {
        spinconstrain::SpinConstrain<TK>& sc = spinconstrain::SpinConstrain<TK>::getScInstance();
        sc.cal_mi_lcao(iter);
    }

    // call iter_finish() of ESolver_KS
    ESolver_KS<TK>::iter_finish(ucell, istep, iter);

    // 1) mix density matrix if mixing_restart + mixing_dmr + not first
    // mixing_restart at every iter
    if (PARAM.inp.mixing_restart > 0 && this->p_chgmix->mixing_restart_count > 0 && PARAM.inp.mixing_dmr)
    {
        elecstate::DensityMatrix<TK, double>* dm = dynamic_cast<elecstate::ElecStateLCAO<TK>*>(this->pelec)->get_DM();
        this->p_chgmix->mix_dmr(dm);
    }

    // 2) save charge density
    // Peize Lin add 2020.04.04
    if (GlobalC::restart.info_save.save_charge)
    {
        for (int is = 0; is < PARAM.inp.nspin; ++is)
        {
            GlobalC::restart.save_disk("charge", is, this->pelec->charge->nrxx, this->pelec->charge->rho[is]);
        }
    }

#ifdef __EXX
    // 3) save exx matrix
    if (PARAM.inp.calculation != "nscf")
    {
        if (GlobalC::exx_info.info_global.cal_exx)
        {
            GlobalC::exx_info.info_ri.real_number ? this->exd->exx_iter_finish(this->kv,
                                                                               ucell,
                                                                               *this->p_hamilt,
                                                                               *this->pelec,
                                                                               *this->p_chgmix,
                                                                               this->scf_ene_thr,
                                                                               iter,
                                                                               istep,
                                                                               this->conv_esolver)
                                                  : this->exc->exx_iter_finish(this->kv,
                                                                               ucell,
                                                                               *this->p_hamilt,
                                                                               *this->pelec,
                                                                               *this->p_chgmix,
                                                                               this->scf_ene_thr,
                                                                               iter,
                                                                               istep,
                                                                               this->conv_esolver);
        }
    }
#endif

    // 6) use the converged occupation matrix for next MD/Relax SCF calculation
    if (PARAM.inp.dft_plus_u && this->conv_esolver)
    {
        GlobalC::dftu.initialed_locale = true;
    }
}

//------------------------------------------------------------------------------
//! the 14th function of ESolver_KS_LCAO: after_scf
//! mohan add 2024-05-11
//! 1) call after_scf() of ESolver_KS
//! 2) write density matrix for sparse matrix
//! 4) write density matrix
//! 6) write Exx matrix
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
void ESolver_KS_LCAO<TK, TR>::after_scf(UnitCell& ucell, const int istep)
{
    ModuleBase::TITLE("ESolver_KS_LCAO", "after_scf");
    // 1) calculate the kinetic energy density tau, sunliang 2024-09-18
    if (PARAM.inp.out_elf[0] > 0)
    {
        this->pelec->cal_tau(*(this->psi));
    }
    
    //! 5) call after_scf() of ESolver_KS
    ESolver_KS<TK>::after_scf(ucell, istep);

    //! 6) write density matrix for sparse matrix
    ModuleIO::write_dmr(dynamic_cast<const elecstate::ElecStateLCAO<TK>*>(this->pelec)->get_DM()->get_DMR_vector(),
                        this->pv,
                        PARAM.inp.out_dm1,
                        false,
                        PARAM.inp.out_app_flag,
                        istep);

    //! 7) write density matrix
    if (PARAM.inp.out_dm)
    {
        std::vector<double> efermis(PARAM.inp.nspin == 2 ? 2 : 1);
        for (int ispin = 0; ispin < efermis.size(); ispin++)
        {
            efermis[ispin] = this->pelec->eferm.get_efval(ispin);
        }
        const int precision = 3;
        ModuleIO::write_dmk(dynamic_cast<const elecstate::ElecStateLCAO<TK>*>(this->pelec)->get_DM()->get_DMK_vector(),
                            precision,
                            efermis,
                            &(ucell),
                            this->pv);
    }

#ifdef __EXX
    //! 8) write Hexx matrix for NSCF (see `out_chg` in docs/advanced/input_files/input-main.md)
    if (PARAM.inp.calculation != "nscf")
    {
        if (GlobalC::exx_info.info_global.cal_exx && PARAM.inp.out_chg[0]
            && istep % PARAM.inp.out_interval == 0) // Peize Lin add if 2022.11.14
        {
            const std::string file_name_exx = PARAM.globalv.global_out_dir + "HexxR" + std::to_string(GlobalV::MY_RANK);
            if (GlobalC::exx_info.info_ri.real_number)
            {
                ModuleIO::write_Hexxs_csr(file_name_exx, ucell, this->exd->get_Hexxs());
            }
            else
            {
                ModuleIO::write_Hexxs_csr(file_name_exx, ucell, this->exc->get_Hexxs());
            }
        }
    }
#endif

    //! 9) Write DeePKS information
#ifdef __DEEPKS
    std::shared_ptr<LCAO_Deepks> ld_shared_ptr(&GlobalC::ld, [](LCAO_Deepks*) {});
    LCAO_Deepks_Interface LDI = LCAO_Deepks_Interface(ld_shared_ptr);
    ModuleBase::timer::tick("ESolver_KS_LCAO", "out_deepks_labels");
    LDI.out_deepks_labels(this->pelec->f_en.etot,
                          this->pelec->klist->get_nks(),
                          ucell.nat,
                          PARAM.globalv.nlocal,
                          this->pelec->ekb,
                          this->pelec->klist->kvec_d,
                          ucell,
                          orb_,
                          GlobalC::GridD,
                          &(this->pv),
                          *(this->psi),
                          dynamic_cast<const elecstate::ElecStateLCAO<TK>*>(this->pelec)->get_DM(),
                          PARAM.inp.deepks_v_delta);

    ModuleBase::timer::tick("ESolver_KS_LCAO", "out_deepks_labels");
#endif

    //! 10) Perform RDMFT calculations
    /******** test RDMFT *********/
    if ( PARAM.inp.rdmft == true ) // rdmft, added by jghan, 2024-10-17
    {
        ModuleBase::matrix occ_number_ks(this->pelec->wg);
        for(int ik=0; ik < occ_number_ks.nr; ++ik) 
        { 
            for(int inb=0; inb < occ_number_ks.nc; ++inb)
            {
                occ_number_ks(ik, inb) /= this->kv.wk[ik];
            }
        }
        this->rdmft_solver.update_elec(occ_number_ks, *(this->psi));

        //! initialize the gradients of Etotal with respect to occupation numbers and wfc, 
        //! and set all elements to 0. 
        ModuleBase::matrix dE_dOccNum(this->pelec->wg.nr, this->pelec->wg.nc, true);
        psi::Psi<TK> dE_dWfc(this->psi->get_nk(), this->psi->get_nbands(), this->psi->get_nbasis()); 
        dE_dWfc.zero_out();

        double Etotal_RDMFT = this->rdmft_solver.run(dE_dOccNum, dE_dWfc);
    }
    /******** test RDMFT *********/


#ifdef __EXX
    // 11) Write RPA information.
    if (PARAM.inp.rpa)
    {
        // ModuleRPA::DFT_RPA_interface
        // rpa_interface(GlobalC::exx_info.info_global);
        RPA_LRI<TK, double> rpa_lri_double(GlobalC::exx_info.info_ri);
        rpa_lri_double.cal_postSCF_exx(*dynamic_cast<const elecstate::ElecStateLCAO<TK>*>(this->pelec)->get_DM(),
                                       MPI_COMM_WORLD,
                                       this->kv,
                                       orb_);
        rpa_lri_double.init(MPI_COMM_WORLD, this->kv, orb_.cutoffs());
        rpa_lri_double.out_for_RPA(this->pv, *(this->psi), this->pelec);
    }
#endif

    // 12) write HR in npz format.
    if (PARAM.inp.out_hr_npz)
    {
        this->p_hamilt->updateHk(0); // first k point, up spin
        hamilt::HamiltLCAO<std::complex<double>, double>* p_ham_lcao
            = dynamic_cast<hamilt::HamiltLCAO<std::complex<double>, double>*>(this->p_hamilt);
        std::string zipname = "output_HR0.npz";
        ModuleIO::output_mat_npz(ucell, zipname, *(p_ham_lcao->getHR()));

        if (PARAM.inp.nspin == 2)
        {
            this->p_hamilt->updateHk(this->kv.get_nks() / 2); // the other half of k points, down spin
            hamilt::HamiltLCAO<std::complex<double>, double>* p_ham_lcao
                = dynamic_cast<hamilt::HamiltLCAO<std::complex<double>, double>*>(this->p_hamilt);
            zipname = "output_HR1.npz";
            ModuleIO::output_mat_npz(ucell, zipname, *(p_ham_lcao->getHR()));
        }
    }

    // 13) write density matrix in the 'npz' format.
    if (PARAM.inp.out_dm_npz)
    {
        const elecstate::DensityMatrix<TK, double>* dm
            = dynamic_cast<const elecstate::ElecStateLCAO<TK>*>(this->pelec)->get_DM();
        std::string zipname = "output_DM0.npz";
        ModuleIO::output_mat_npz(ucell, zipname, *(dm->get_DMR_pointer(1)));

        if (PARAM.inp.nspin == 2)
        {
            zipname = "output_DM1.npz";
            ModuleIO::output_mat_npz(ucell, zipname, *(dm->get_DMR_pointer(2)));
        }
    }

    //! 14) Print out information every 'out_interval' steps.
    if (PARAM.inp.calculation != "md" || istep % PARAM.inp.out_interval == 0)
    {
        //! Print out sparse matrix
        ModuleIO::output_mat_sparse(PARAM.inp.out_mat_hs2,
                                    PARAM.inp.out_mat_dh,
                                    PARAM.inp.out_mat_t,
                                    PARAM.inp.out_mat_r,
                                    istep,
                                    this->pelec->pot->get_effective_v(),
                                    this->pv,
                                    this->GK,
                                    two_center_bundle_,
                                    orb_,
                                    ucell,
                                    GlobalC::GridD,
                                    this->kv,
                                    this->p_hamilt);
        
        //! Perform Mulliken charge analysis
        if (PARAM.inp.out_mul)
        {
            ModuleIO::cal_mag(&(this->pv), this->p_hamilt, this->kv, this->pelec, ucell, istep, true);
        }
    }

    //! 15) Print out atomic magnetization only when 'spin_constraint' is on.
    if (PARAM.inp.sc_mag_switch) 
    {
        spinconstrain::SpinConstrain<TK>& sc = spinconstrain::SpinConstrain<TK>::getScInstance();
        sc.cal_mi_lcao(istep);
        sc.print_Mi(GlobalV::ofs_running);
        sc.print_Mag_Force(GlobalV::ofs_running);
    }

    //! 16) Clean up RA. 
    //! this should be last function and put it in the end, mohan request 2024-11-28
    if (!PARAM.inp.cal_force && !PARAM.inp.cal_stress)
    {
        RA.delete_grid();
    }

    //! 17) Print out quasi-orbitals.
    if (PARAM.inp.qo_switch)
    {
        toQO tqo(PARAM.inp.qo_basis, PARAM.inp.qo_strategy, PARAM.inp.qo_thr, PARAM.inp.qo_screening_coeff);
        tqo.initialize(PARAM.globalv.global_out_dir,
                       PARAM.inp.pseudo_dir,
                       PARAM.inp.orbital_dir,
                       &ucell,
                       this->kv.kvec_d,
                       GlobalV::ofs_running,
                       GlobalV::MY_RANK,
                       GlobalV::NPROC);
        tqo.calculate();
    }

    //! 18) Print out kinetic matrix.
    if (PARAM.inp.out_mat_tk[0])
    {
        hamilt::HS_Matrix_K<TK> hsk(&pv, true);
        hamilt::HContainer<TR> hR(&pv);
        hamilt::Operator<TK>* ekinetic
            = new hamilt::EkineticNew<hamilt::OperatorLCAO<TK, TR>>(&hsk,
                                                                    this->kv.kvec_d,
                                                                    &hR,
                                                                    &ucell,
                                                                    orb_.cutoffs(),
                                                                    &GlobalC::GridD,
                                                                    two_center_bundle_.kinetic_orb.get());

        const int nspin_k = (PARAM.inp.nspin == 2 ? 2 : 1);
        for (int ik = 0; ik < this->kv.get_nks() / nspin_k; ++ik)
        {
            ekinetic->init(ik);
            ModuleIO::save_mat(0,
                               hsk.get_hk(),
                               PARAM.globalv.nlocal,
                               false,
                               PARAM.inp.out_mat_tk[1],
                               1,
                               PARAM.inp.out_app_flag,
                               "T",
                               "data-" + std::to_string(ik),
                               this->pv,
                               GlobalV::DRANK);
        }

        // where is new? mohan ask 2024-11-28
        delete ekinetic;
    }

    //! 19) Wannier 90 function, added by jingan in 2018.11.7
    if (PARAM.inp.calculation == "nscf" && PARAM.inp.towannier90)
    {
        std::cout << FmtCore::format("\n * * * * * *\n << Start %s.\n", "Wave function to Wannier90");
        if (PARAM.inp.wannier_method == 1)
        {
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
    }

    //! 20) berry phase calculations, added by jingan
    if (PARAM.inp.calculation == "nscf" && 
        berryphase::berry_phase_flag && 
        ModuleSymmetry::Symmetry::symm_flag != 1)
    {
        std::cout << FmtCore::format("\n * * * * * *\n << Start %s.\n", "Berry phase calculation");
        berryphase bp(&(this->pv));
        bp.lcao_init(this->kv,
                     this->GridT,
                     orb_); // additional step before calling
                            // macroscopic_polarization (why capitalize
                            // the function name?)
        bp.Macroscopic_polarization(this->pw_wfc->npwk_max, this->psi, this->pw_rho, this->pw_wfc, this->kv);
        std::cout << FmtCore::format(" >> Finish %s.\n * * * * * *\n", "Berry phase calculation");
    }
}

template class ESolver_KS_LCAO<double, double>;
template class ESolver_KS_LCAO<std::complex<double>, double>;
template class ESolver_KS_LCAO<std::complex<double>, std::complex<double>>;
} // namespace ModuleESolver
