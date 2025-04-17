#include "esolver_ks_lcao.h"

#include "module_io/write_dos_lcao.h"       // write DOS and PDOS
#include "module_io/write_proj_band_lcao.h" // projcted band structure

#include "module_base/formatter.h"
#include "module_base/global_variable.h"
#include "module_base/tool_title.h"
#include "module_elecstate/elecstate_tools.h"

#include "module_elecstate/module_dm/cal_dm_psi.h"
#include "module_hamilt_lcao/module_deltaspin/spin_constrain.h"
#include "module_hamilt_lcao/module_dftu/dftu.h"
#include "module_io/berryphase.h"
#include "module_io/cal_ldos.h"
#include "module_io/cube_io.h"
#include "module_io/io_dmk.h"
#include "module_io/io_npz.h"
#include "module_io/output_dmk.h"
#include "module_io/output_log.h"
#include "module_io/output_mat_sparse.h"
#include "module_io/output_mulliken.h"
#include "module_io/output_sk.h"
#include "module_io/read_wfc_nao.h"
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

// be careful of hpp, there may be multiple definitions of functions, 20250302, mohan
#include "module_hamilt_lcao/hamilt_lcaodft/hs_matrix_k.hpp"
#include "module_io/write_eband_terms.hpp"
#include "module_io/write_vxc.hpp"
#include "module_io/write_vxc_r.hpp"

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

template <typename TK, typename TR>
ESolver_KS_LCAO<TK, TR>::~ESolver_KS_LCAO()
{
}

template <typename TK, typename TR>
void ESolver_KS_LCAO<TK, TR>::before_all_runners(UnitCell& ucell, const Input_para& inp)
{
    ModuleBase::TITLE("ESolver_KS_LCAO", "before_all_runners");
    ModuleBase::timer::tick("ESolver_KS_LCAO", "before_all_runners");

    // 1) before_all_runners in ESolver_KS
    ESolver_KS<TK>::before_all_runners(ucell, inp);

    // 2) init ElecState
    // autoset nbands in ElecState before basis_init (for Psi 2d division)
    if (this->pelec == nullptr)
    {
        // TK stands for double and complex<double>?
        this->pelec = new elecstate::ElecStateLCAO<TK>(&(this->chr), // use which parameter?
                                                       &(this->kv),
                                                       this->kv.get_nks(),
                                                       &(this->GG),
                                                       &(this->GK),
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

    // 4) initialize electronic wave function psi
    if (this->psi == nullptr)
    {
        int nsk = 0;
        int ncol = 0;
        if (PARAM.globalv.gamma_only_local)
        {
            nsk = PARAM.inp.nspin;
            ncol = this->pv.ncol_bands;
            if (PARAM.inp.ks_solver == "genelpa" || PARAM.inp.ks_solver == "elpa" || PARAM.inp.ks_solver == "lapack"
                || PARAM.inp.ks_solver == "pexsi" || PARAM.inp.ks_solver == "cusolver"
                || PARAM.inp.ks_solver == "cusolvermp")
            {
                ncol = this->pv.ncol;
            }
        }
        else
        {
            nsk = this->kv.get_nks();
#ifdef __MPI
            ncol = this->pv.ncol_bands;
#else
            ncol = PARAM.inp.nbands;
#endif
        }
        this->psi = new psi::Psi<TK>(nsk, ncol, this->pv.nrow, this->kv.ngk, true);
    }

    // 5) read psi from file
    if (PARAM.inp.init_wfc == "file")
    {
        if (!ModuleIO::read_wfc_nao(PARAM.globalv.global_readin_dir, this->pv, *(this->psi), this->pelec))
        {
            ModuleBase::WARNING_QUIT("ESolver_KS_LCAO", "read electronic wave functions failed");
        }
    }

    // 6) initialize the density matrix
    // DensityMatrix is allocated here, DMK is also initialized here
    // DMR is not initialized here, it will be constructed in each before_scf
    dynamic_cast<elecstate::ElecStateLCAO<TK>*>(this->pelec)->init_DM(&this->kv, &(this->pv), PARAM.inp.nspin);

    // 7) initialize exact exchange calculations
#ifdef __EXX
    if (PARAM.inp.calculation == "scf" || PARAM.inp.calculation == "relax" || PARAM.inp.calculation == "cell-relax"
        || PARAM.inp.calculation == "md")
    {
        if (GlobalC::exx_info.info_global.cal_exx)
        {
            if (PARAM.inp.init_wfc != "file")
            { // if init_wfc==file, directly enter the EXX loop
                XC_Functional::set_xc_first_loop(ucell);
            }

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

    // 8) initialize DFT+U
    if (PARAM.inp.dft_plus_u)
    {
        auto* dftu = ModuleDFTU::DFTU::get_instance();
        dftu->init(ucell, &this->pv, this->kv.get_nks(), &orb_);
    }

    // 9) initialize local pseudopotentials
    this->locpp.init_vloc(ucell, this->pw_rho);
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "LOCAL POTENTIAL");

    // 10) inititlize the charge density
    this->chr.allocate(PARAM.inp.nspin);
    this->pelec->omega = ucell.omega;

    // 11) initialize the potential
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

    // 12) initialize deepks
#ifdef __DEEPKS
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
                                ld.inl2l,
                                *orb_.Alpha,
                                ld.pdm);
    }
#endif

    // 13) set occupations
    // tddft does not need to set occupations in the first scf
    if (PARAM.inp.ocp && inp.esolver_type != "tddft")
    {
        elecstate::fixed_weights(PARAM.inp.ocp_kb,
                                 PARAM.inp.nbands,
                                 PARAM.inp.nelec,
                                 this->pelec->klist,
                                 this->pelec->wg,
                                 this->pelec->skip_weights);
    }

    // 14) if kpar is not divisible by nks, print a warning
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

    // 15) initialize rdmft, added by jghan
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

template <typename TK, typename TR>
double ESolver_KS_LCAO<TK, TR>::cal_energy()
{
    ModuleBase::TITLE("ESolver_KS_LCAO", "cal_energy");

    return this->pelec->f_en.etot;
}

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

    ESolver_KS<TK>::after_all_runners(ucell);

    const int nspin0 = (PARAM.inp.nspin == 2) ? 2 : 1;

    // 4) write projected band structure
    if (PARAM.inp.out_proj_band)
    {
        ModuleIO::write_proj_band_lcao(this->psi, this->pv, this->pelec, this->kv, ucell, this->p_hamilt);
    }

    // 5) print out density of states (DOS)
    if (PARAM.inp.out_dos)
	{
		ModuleIO::write_dos_lcao(this->psi,
				this->p_hamilt,
				this->pv,
				ucell,
				*(this->pelec->klist),
				PARAM.inp.nbands,
				this->pelec->eferm,
				this->pelec->ekb,
				this->pelec->wg,
				PARAM.inp.dos_edelta_ev,
				PARAM.inp.dos_scale,
				PARAM.inp.dos_sigma,
                GlobalV::ofs_running);
    }

    // out ldos
    if (PARAM.inp.out_ldos[0])
    {
        ModuleIO::Cal_ldos<TK>::cal_ldos_lcao(reinterpret_cast<elecstate::ElecStateLCAO<TK>*>(this->pelec),
                                              this->psi[0],
                                              this->Pgrid,
                                              ucell);
    }

    // 6) print out exchange-correlation potential
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
    if (PARAM.inp.out_mat_xc2)
    {
        ModuleIO::write_Vxc_R<TK, TR>(PARAM.inp.nspin,
                                      &this->pv,
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
                                      this->gd
#ifdef __EXX
                                      ,
                                      this->exx_lri_double ? &this->exx_lri_double->Hexxs : nullptr,
                                      this->exx_lri_complex ? &this->exx_lri_complex->Hexxs : nullptr
#endif
        );
    }

    // 7) write eband terms
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

    // first need to calculate the weight according to
    // electrons number.
    if (istep == 0 && PARAM.inp.init_wfc == "file")
    {
        int exx_two_level_step = 0;
#ifdef __EXX
        if (GlobalC::exx_info.info_global.cal_exx)
        {
            // the following steps are only needed in the first outer exx loop
            exx_two_level_step
                = GlobalC::exx_info.info_ri.real_number ? this->exd->two_level_step : this->exc->two_level_step;
        }
#endif
        if (iter == 1 && exx_two_level_step == 0)
        {
            std::cout << " WAVEFUN -> CHARGE " << std::endl;

            // calculate the density matrix using read in wave functions
            // and then calculate the charge density on grid.

            this->pelec->skip_weights = true;
            elecstate::calculate_weights(this->pelec->ekb,
                                         this->pelec->wg,
                                         this->pelec->klist,
                                         this->pelec->eferm,
                                         this->pelec->f_en,
                                         this->pelec->nelec_spin,
                                         this->pelec->skip_weights);

            auto _pelec = dynamic_cast<elecstate::ElecStateLCAO<TK>*>(this->pelec);
            elecstate::calEBand(_pelec->ekb, _pelec->wg, _pelec->f_en);
            elecstate::cal_dm_psi(_pelec->DM->get_paraV_pointer(), _pelec->wg, *this->psi, *(_pelec->DM));
            _pelec->DM->cal_DMR();

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
    if (iter == 1 && istep == 0)
    {
        // initialize DMR
        this->ld.init_DMR(ucell, orb_, this->pv, this->gd);
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

template <typename TK, typename TR>
void ESolver_KS_LCAO<TK, TR>::hamilt2rho_single(UnitCell& ucell, int istep, int iter, double ethr)
{
    ModuleBase::TITLE("ESolver_KS_LCAO", "hamilt2rho_single");

    // i1) reset energy
    this->pelec->f_en.eband = 0.0;
    this->pelec->f_en.demet = 0.0;
    bool skip_charge = PARAM.inp.calculation == "nscf" ? true : false;

    // 2) run the inner lambda loop to contrain atomic moments with the DeltaSpin method
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

    // 3) run Hsolver
    if (!skip_solve)
    {
        hsolver::HSolverLCAO<TK> hsolver_lcao_obj(&(this->pv), PARAM.inp.ks_solver);
        hsolver_lcao_obj.solve(this->p_hamilt, this->psi[0], this->pelec, skip_charge);
    }

    // 4) EXX
#ifdef __EXX
    if (PARAM.inp.calculation != "nscf")
    {
        if (GlobalC::exx_info.info_ri.real_number)
        {
            this->exd->exx_hamilt2rho(*this->pelec, this->pv, iter);
        }
        else
        {
            this->exc->exx_hamilt2rho(*this->pelec, this->pv, iter);
        }
    }
#endif

    // 5) symmetrize the charge density
    Symmetry_rho srho;
    for (int is = 0; is < PARAM.inp.nspin; is++)
    {
        srho.begin(is, this->chr, this->pw_rho, ucell.symm);
    }

    // 6) calculate delta energy
    this->pelec->f_en.deband = this->pelec->cal_delta_eband(ucell);
}

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

template <typename TK, typename TR>
void ESolver_KS_LCAO<TK, TR>::iter_finish(UnitCell& ucell, const int istep, int& iter, bool& conv_esolver)
{
    ModuleBase::TITLE("ESolver_KS_LCAO", "iter_finish");

    // 1) calculate the local occupation number matrix and energy correction
    // in DFT+U
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

    // 2) for deepks, calculate delta_e
#ifdef __DEEPKS
    if (PARAM.inp.deepks_scf)
    {
        const std::vector<std::vector<TK>>& dm
            = dynamic_cast<const elecstate::ElecStateLCAO<TK>*>(this->pelec)->get_DM()->get_DMK_vector();

        ld.dpks_cal_e_delta_band(dm, this->kv.get_nks());
        DeePKS_domain::update_dmr(this->kv.kvec_d, dm, ucell, orb_, this->pv, this->gd, ld.dm_r);
        this->pelec->f_en.edeepks_scf = ld.E_delta - ld.e_delta_band;
        this->pelec->f_en.edeepks_delta = ld.E_delta;
    }
#endif

    // 3) for delta spin
    if (PARAM.inp.sc_mag_switch)
    {
        spinconstrain::SpinConstrain<TK>& sc = spinconstrain::SpinConstrain<TK>::getScInstance();
        sc.cal_mi_lcao(iter);
    }

    // 4) call iter_finish() of ESolver_KS
    ESolver_KS<TK>::iter_finish(ucell, istep, iter, conv_esolver);

    // 5) mix density matrix if mixing_restart + mixing_dmr + not first
    // mixing_restart at every iter
    if (PARAM.inp.mixing_restart > 0 && this->p_chgmix->mixing_restart_count > 0 && PARAM.inp.mixing_dmr)
    {
        elecstate::DensityMatrix<TK, double>* dm = dynamic_cast<elecstate::ElecStateLCAO<TK>*>(this->pelec)->get_DM();
        this->p_chgmix->mix_dmr(dm);
    }

    // 6) save charge density
    // Peize Lin add 2020.04.04
    if (GlobalC::restart.info_save.save_charge)
    {
        for (int is = 0; is < PARAM.inp.nspin; ++is)
        {
            GlobalC::restart.save_disk("charge", is, this->chr.nrxx, this->chr.rho[is]);
        }
    }

#ifdef __EXX
    // 7) save exx matrix
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

    // 8) use the converged occupation matrix for next MD/Relax SCF calculation
    if (PARAM.inp.dft_plus_u && conv_esolver)
    {
        GlobalC::dftu.initialed_locale = true;
    }
}

template class ESolver_KS_LCAO<double, double>;
template class ESolver_KS_LCAO<std::complex<double>, double>;
template class ESolver_KS_LCAO<std::complex<double>, std::complex<double>>;
} // namespace ModuleESolver
