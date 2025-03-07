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
#include "module_io/write_elecstat_pot.h"
#include "module_io/write_istate_info.h"
#include "module_io/write_proj_band_lcao.h"
#include "module_io/write_wfc_nao.h"
#include "module_parameter/parameter.h"


//be careful of hpp, there may be multiple definitions of functions, 20250302, mohan
#include "module_io/write_eband_terms.hpp"
#include "module_io/write_vxc.hpp"
#include "module_hamilt_lcao/hamilt_lcaodft/hs_matrix_k.hpp"

//--------------temporary----------------------------
#include "module_base/global_function.h"
#include "module_cell/module_neighbor/sltk_grid_driver.h"
#include "module_elecstate/cal_ux.h"
#include "module_elecstate/module_charge/symmetry_rho.h"
#include "module_elecstate/occupy.h"
#include "module_hamilt_lcao/hamilt_lcaodft/LCAO_domain.h" // need DeePKS_init
#include "module_hamilt_lcao/module_dftu/dftu.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_io/print_info.h"

#include <memory>

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
//! 5) initialize exx
//! 6) initialize DFT+U
//! 7) ppcell
//! 8) inititlize the charge density
//! 9) initialize the potential.
//! 10) initialize deepks
//! 11) set occupations
//! 12) print a warning if needed
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

#ifdef __EXX
    // 5) initialize exx
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
                this->exx_lri_double->init(MPI_COMM_WORLD, ucell, this->kv, orb_);
                this->exd->exx_before_all_runners(this->kv, ucell, this->pv);
            }
            else
            {
                this->exx_lri_complex->init(MPI_COMM_WORLD, ucell, this->kv, orb_);
                this->exc->exx_before_all_runners(this->kv, ucell, this->pv);
            }
        }
    }
#endif

    // 6) initialize DFT+U
    if (PARAM.inp.dft_plus_u)
    {
        auto* dftu = ModuleDFTU::DFTU::get_instance();
        dftu->init(ucell, &this->pv, this->kv.get_nks(), &orb_);
    }

    // 7) initialize ppcell
    this->locpp.init_vloc(ucell, this->pw_rho);
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "LOCAL POTENTIAL");

    // 8) inititlize the charge density
    this->chr.allocate(PARAM.inp.nspin);
    this->pelec->omega = ucell.omega;

    // 9) initialize the potential
    if (this->pelec->pot == nullptr)
    {
        this->pelec->pot = new elecstate::Potential(this->pw_rhod,
                                                    this->pw_rho,
                                                    &ucell,
                                                    &(this->locpp.vloc),
                                                    &(this->sf),
                                                    &(this->solvent),
                                                    &(this->pelec->f_en.etxc),
                                                    &(this->pelec->f_en.vtxc));
    }

#ifdef __DEEPKS
    // 10) initialize deepks
    LCAO_domain::DeePKS_init(ucell, pv, this->kv.get_nks(), orb_, this->ld, GlobalV::ofs_running);
    if (PARAM.inp.deepks_scf)
    {
        // load the DeePKS model from deep neural network
        DeePKS_domain::load_model(PARAM.inp.deepks_model, ld.model_deepks);
        // read pdm from file for NSCF or SCF-restart, do it only once in whole calculation
        DeePKS_domain::read_pdm((PARAM.inp.init_chg == "file"),
                                PARAM.inp.deepks_equiv,
                                ld.init_pdm,
                                ucell.nat,
                                orb_.Alpha[0].getTotal_nchi() * ucell.nat,
                                ld.lmaxd,
                                ld.inl_l,
                                *orb_.Alpha,
                                ld.pdm);
    }
#endif

    // 11) set occupations
    // tddft does not need to set occupations in the first scf
    if (PARAM.inp.ocp && inp.esolver_type != "tddft")
    {
        this->pelec->fixed_weights(PARAM.inp.ocp_kb, PARAM.inp.nbands, PARAM.inp.nelec);
    }

    // 12) if kpar is not divisible by nks, print a warning
    if (PARAM.globalv.kpar_lcao > 1)
    {
        if (this->kv.get_nks() % PARAM.globalv.kpar_lcao != 0)
        {
            ModuleBase::WARNING("ESolver_KS_LCAO::before_all_runners", "nks is not divisible by kpar.");
            std::cout << "\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
                         "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
                         "%%%%%%%%%%%%%%%%%%%%%%%%%%"
                      << std::endl;
            std::cout << " Warning: nks (" << this->kv.get_nks() << ") is not divisible by kpar ("
                      << PARAM.globalv.kpar_lcao << ")." << std::endl;
            std::cout << " This may lead to poor load balance. It is strongly suggested to" << std::endl;
            std::cout << " set nks to be divisible by kpar, but if this is really what" << std::endl;
            std::cout << " you want, please ignore this warning." << std::endl;
            std::cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
                         "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
                         "%%%%%%%%%%%%\n";
        }
    }

    // 13) initialize rdmft, added by jghan
    if (PARAM.inp.rdmft == true)
    {
        rdmft_solver.init(this->GG,
                          this->GK,
                          this->pv,
                          ucell,
                          this->gd,
                          this->kv,
                          *(this->pelec),
                          this->orb_,
                          two_center_bundle_,
                          PARAM.inp.dft_functional,
                          PARAM.inp.rdmft_power_alpha);
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

    fsl.getForceStress(ucell,
                       PARAM.inp.cal_force,
                       PARAM.inp.cal_stress,
                       PARAM.inp.test_force,
                       PARAM.inp.test_stress,
                       this->gd,
                       this->pv,
                       this->pelec,
                       this->psi,
                       this->GG, // mohan add 2024-04-01
                       this->GK, // mohan add 2024-04-01
                       two_center_bundle_,
                       orb_,
                       force,
                       this->scs,
                       this->locpp,
                       this->sf,
                       this->kv,
                       this->pw_rho,
                       this->solvent,
#ifdef __DEEPKS
                       this->ld,
#endif
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

    // if the users do not want to calculate forces but want stress,
    // we call cal_force
    if (!this->have_force)
    {
        ModuleBase::matrix fcs;
        this->cal_force(ucell, fcs);
    }

    // the 'scs' stress has already been calculated in 'cal_force'
    stress = this->scs;
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
        ModuleIO::write_istate_info(this->pelec->ekb, this->pelec->wg, this->kv);
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
                                this->kv);
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
                                    this->solvent,
                                    *this->pw_rho,
                                    *this->pw_rhod,
                                    this->locpp.vloc,
                                    this->chr,
                                    this->GG,
                                    this->GK,
                                    this->kv,
                                    orb_.cutoffs(),
                                    this->pelec->wg,
                                    this->gd
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
                                            this->solvent,
                                            *this->pw_rho,
                                            *this->pw_rhod,
                                            this->locpp.vloc,
                                            this->chr,
                                            this->GG,
                                            this->GK,
                                            this->kv,
                                            this->pelec->wg,
                                            this->gd,
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
        this->p_chgmix->mix_reset(); // init mixing
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
    this->pelec->f_en.deband_harris = this->pelec->cal_delta_eband(ucell);

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
            this->pelec->pot->update_from_charge(&this->chr, &ucell);

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
                                        ucell,
                                        *dynamic_cast<const elecstate::ElecStateLCAO<TK>*>(this->pelec)->get_DM(),
                                        this->kv,
                                        iter);
        }
        else
        {
            this->exc->exx_eachiterinit(istep,
                                        ucell,
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
        GlobalC::dftu.cal_slater_UJ(ucell, this->chr.rho, this->pw_rho->nrxx);
    }

#ifdef __DEEPKS
    // the density matrixes of DeePKS have been updated in each iter
    ld.set_hr_cal(true);

    // HR in HamiltLCAO should be recalculate
    if (PARAM.inp.deepks_scf)
    {
        this->p_hamilt->refresh();
    }
#endif

    if (PARAM.inp.vl_in_h)
    {
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
        if (!sc.mag_converged() && this->drho > 0 && this->drho < PARAM.inp.sc_scf_thr)
        {
            // optimize lambda to get target magnetic moments, but the lambda is not near target
            sc.run_lambda_loop(iter - 1);
            sc.set_mag_converged(true);
            skip_solve = true;
        }
        else if (sc.mag_converged())
        {
            // optimize lambda to get target magnetic moments, but the lambda is not near target
            sc.run_lambda_loop(iter - 1);
            skip_solve = true;
        }
    }
    if (!skip_solve)
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
        srho.begin(is, this->chr, this->pw_rho, ucell.symm);
    }

    // 12) calculate delta energy
    this->pelec->f_en.deband = this->pelec->cal_delta_eband(ucell);
}

//------------------------------------------------------------------------------
//! the 12th function of ESolver_KS_LCAO: update_pot
//! mohan add 2024-05-11
//! 1) print potential
//------------------------------------------------------------------------------
template <typename TK, typename TR>
void ESolver_KS_LCAO<TK, TR>::update_pot(UnitCell& ucell, const int istep, const int iter, const bool conv_esolver)
{
    ModuleBase::TITLE("ESolver_KS_LCAO", "update_pot");

    if (!conv_esolver)
    {
        elecstate::cal_ux(ucell);
        this->pelec->pot->update_from_charge(&this->chr, &ucell);
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
void ESolver_KS_LCAO<TK, TR>::iter_finish(UnitCell& ucell, const int istep, int& iter, bool& conv_esolver)
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
                ModuleDFTU::dftu_cal_occup_m(iter,
                                             ucell,
                                             tmp_dm,
                                             this->kv,
                                             this->p_chgmix->get_mixing_beta(),
                                             this->p_hamilt);
            }
            GlobalC::dftu.cal_energy_correction(ucell, istep);
        }
        GlobalC::dftu.output(ucell);
    }

    // (7) for deepks, calculate delta_e
#ifdef __DEEPKS
    if (PARAM.inp.deepks_scf)
    {
        const std::vector<std::vector<TK>>& dm
            = dynamic_cast<const elecstate::ElecStateLCAO<TK>*>(this->pelec)->get_DM()->get_DMK_vector();

        ld.dpks_cal_e_delta_band(dm, this->kv.get_nks());
        this->pelec->f_en.edeepks_scf = ld.E_delta - ld.e_delta_band;
        this->pelec->f_en.edeepks_delta = ld.E_delta;
    }
#endif

    // 8) for delta spin
    if (PARAM.inp.sc_mag_switch)
    {
        spinconstrain::SpinConstrain<TK>& sc = spinconstrain::SpinConstrain<TK>::getScInstance();
        sc.cal_mi_lcao(iter);
    }

    // call iter_finish() of ESolver_KS
    ESolver_KS<TK>::iter_finish(ucell, istep, iter, conv_esolver);

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
            GlobalC::restart.save_disk("charge", is, this->chr.nrxx, this->chr.rho[is]);
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
                                                                               conv_esolver)
                                                  : this->exc->exx_iter_finish(this->kv,
                                                                               ucell,
                                                                               *this->p_hamilt,
                                                                               *this->pelec,
                                                                               *this->p_chgmix,
                                                                               this->scf_ene_thr,
                                                                               iter,
                                                                               istep,
                                                                               conv_esolver);
        }
    }
#endif

    // 6) use the converged occupation matrix for next MD/Relax SCF calculation
    if (PARAM.inp.dft_plus_u && conv_esolver)
    {
        GlobalC::dftu.initialed_locale = true;
    }
}

template class ESolver_KS_LCAO<double, double>;
template class ESolver_KS_LCAO<std::complex<double>, double>;
template class ESolver_KS_LCAO<std::complex<double>, std::complex<double>>;
} // namespace ModuleESolver
