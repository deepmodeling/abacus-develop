#include "esolver_ks_lcao.h"
#include "source_estate/elecstate_tools.h"
#include "source_lcao/module_deltaspin/spin_constrain.h"
#include "source_io/read_wfc_nao.h"
#include "source_lcao/hs_matrix_k.hpp" // there may be multiple definitions if using hpp
#include "source_estate/cal_ux.h"
#include "source_estate/module_charge/symmetry_rho.h"
#include "source_lcao/LCAO_domain.h" // need DeePKS_init
#include "source_lcao/module_dftu/dftu.h"
#include "source_lcao/FORCE_STRESS.h"
#include "source_estate/elecstate_lcao.h"
#include "source_lcao/hamilt_lcao.h"
#include "source_hsolver/hsolver_lcao.h"
#ifdef __EXX
#include "../source_lcao/module_ri/exx_opt_orb.h"
#endif
#include "source_lcao/module_rdmft/rdmft.h"
#include "source_estate/module_charge/chgmixing.h" // use charge mixing, mohan add 20251006 
#include "source_estate/module_dm/setup_dm.h" // setup dm from electronic wave functions
#include "source_io/ctrl_runner_lcao.h" // use ctrl_runner_lcao() 
#include "source_io/ctrl_iter_lcao.h" // use ctrl_iter_lcao() 
#include "source_io/ctrl_scf_lcao.h" // use ctrl_scf_lcao()
#include "source_psi/setup_psi.h" // mohan add 20251019
#include "source_io/read_wfc_nao.h" 
#include "source_io/print_info.h"

namespace ModuleESolver
{

template <typename TK, typename TR>
ESolver_KS_LCAO<TK, TR>::ESolver_KS_LCAO()
{
    this->classname = "ESolver_KS_LCAO";
    this->basisname = "LCAO";
    this->exx_nao.init(); // mohan add 20251008
}

template <typename TK, typename TR>
ESolver_KS_LCAO<TK, TR>::~ESolver_KS_LCAO()
{
	//****************************************************
	// do not add any codes in this deconstructor funcion
	//****************************************************
}

template <typename TK, typename TR>
void ESolver_KS_LCAO<TK, TR>::before_all_runners(UnitCell& ucell, const Input_para& inp)
{
    ModuleBase::TITLE("ESolver_KS_LCAO", "before_all_runners");
    ModuleBase::timer::tick("ESolver_KS_LCAO", "before_all_runners");

    // 1) before_all_runners in ESolver_KS
    ESolver_KS<TK>::before_all_runners(ucell, inp);

    // 2) autoset nbands in ElecState before init_basis (for Psi 2d division)
    if (this->pelec == nullptr)
    {
        // TK stands for double and std::complex<double>?
        this->pelec = new elecstate::ElecStateLCAO<TK>(&(this->chr), &(this->kv),
          this->kv.get_nks(), this->pw_rho, this->pw_big);
    }

    // 3) read LCAO orbitals/projectors and construct the interpolation tables.
    LCAO_domain::init_basis_lcao(this->pv, inp.onsite_radius, inp.lcao_ecut,
      inp.lcao_dk, inp.lcao_dr, inp.lcao_rmax, ucell, two_center_bundle_, orb_);

    // 4) setup EXX calculations
    if (inp.calculation == "gen_opt_abfs")
    {
#ifdef __EXX
        Exx_Opt_Orb exx_opt_orb;
        exx_opt_orb.generate_matrix(GlobalC::exx_info.info_opt_abfs, this->kv, ucell, this->orb_);
#else
        ModuleBase::WARNING_QUIT("ESolver_KS_LCAO::before_all_runners", "calculation=gen_opt_abfs must compile __EXX");
#endif
        return;
    }

    // 5) init electronic wave function psi
    Setup_Psi<TK>::allocate_psi(this->psi, this->kv, this->pv, inp);

    //! read psi from file
    if (inp.init_wfc == "file" && inp.esolver_type != "tddft")
    {
        if (!ModuleIO::read_wfc_nao(PARAM.globalv.global_readin_dir,
             this->pv, *this->psi, this->pelec->ekb, this->pelec->wg, this->kv.ik2iktot,
             this->kv.get_nkstot(), inp.nspin))
        {
            ModuleBase::WARNING_QUIT("ESolver_KS_LCAO", "read electronic wave functions failed");
        }
    }


    // 7) init DMK, but DMR is constructed in before_scf()
    dynamic_cast<elecstate::ElecStateLCAO<TK>*>(this->pelec)->init_DM(&this->kv, &(this->pv), inp.nspin);

    // 8) init exact exchange calculations
    this->exx_nao.before_runner(ucell, this->kv, this->orb_, this->pv, inp);

    // 9) initialize DFT+U
    if (inp.dft_plus_u)
    {
        auto* dftu = ModuleDFTU::DFTU::get_instance();
        dftu->init(ucell, &this->pv, this->kv.get_nks(), &orb_);
    }

    // 10) init local pseudopotentials
    this->locpp.init_vloc(ucell, this->pw_rho);
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "LOCAL POTENTIAL");

    // 11) init charge density
    this->chr.allocate(inp.nspin);
    this->pelec->omega = ucell.omega;

    // 12) init potentials
    if (this->pelec->pot == nullptr)
    {
        this->pelec->pot = new elecstate::Potential(this->pw_rhod, this->pw_rho,
          &ucell, &(this->locpp.vloc), &(this->sf), &(this->solvent),
          &(this->pelec->f_en.etxc), &(this->pelec->f_en.vtxc));
    }

    // 13) init deepks
    this->deepks.before_runner(ucell, this->kv.get_nks(), this->orb_, this->pv, inp);

    // 14) set occupations, tddft does not need to set occupations in the first scf
    if (inp.ocp && inp.esolver_type != "tddft")
    {
        elecstate::fixed_weights(inp.ocp_kb, inp.nbands, inp.nelec,
          this->pelec->klist, this->pelec->wg, this->pelec->skip_weights);
    }

    // 15) if kpar is not divisible by nks, print a warning
    ModuleIO::print_kpar(this->kv.get_nks(), PARAM.globalv.kpar_lcao);

    // 16) init rdmft, added by jghan
    if (inp.rdmft == true)
    {
        rdmft_solver.init(this->pv, ucell,
          this->gd, this->kv, *(this->pelec), this->orb_,
          two_center_bundle_, inp.dft_functional, inp.rdmft_power_alpha);
    }

    ModuleBase::timer::tick("ESolver_KS_LCAO", "before_all_runners");
    return;
}


template <typename TK, typename TR>
void ESolver_KS_LCAO<TK, TR>::before_scf(UnitCell& ucell, const int istep)
{
    ModuleBase::TITLE("ESolver_KS_LCAO", "before_scf");
    ModuleBase::timer::tick("ESolver_KS_LCAO", "before_scf");

    //! 1) call before_scf() of ESolver_KS.
    ESolver_KS<TK>::before_scf(ucell, istep);

    auto* estate = dynamic_cast<elecstate::ElecStateLCAO<TK>*>(this->pelec);
    if(!estate)
    {
        ModuleBase::WARNING_QUIT("ESolver_KS_LCAO::before_scf","pelec does not exist");
    }

    //! 2) find search radius
    double search_radius = atom_arrange::set_sr_NL(GlobalV::ofs_running,
      PARAM.inp.out_level, orb_.get_rcutmax_Phi(), ucell.infoNL.get_rcutmax_Beta(),
      PARAM.globalv.gamma_only_local);

    //! 3) use search_radius to search adj atoms
    atom_arrange::search(PARAM.globalv.search_pbc, GlobalV::ofs_running,
      this->gd, ucell, search_radius, PARAM.inp.test_atom_input);

    //! 4) initialize NAO basis set
    // here new is a unique pointer, which will be deleted automatically
    gint_info_.reset(
        new ModuleGint::GintInfo(
        this->pw_big->nbx, this->pw_big->nby, this->pw_big->nbz,
        this->pw_rho->nx, this->pw_rho->ny, this->pw_rho->nz,
        0, 0, this->pw_big->nbzp_start,
        this->pw_big->nbx, this->pw_big->nby, this->pw_big->nbzp,
        orb_.Phi, ucell, this->gd));
    ModuleGint::Gint::set_gint_info(gint_info_.get());

    // 7) For each atom, calculate the adjacent atoms in different cells
    // and allocate the space for H(R) and S(R).
    // If k point is used here, allocate HlocR after atom_arrange.
    this->RA.for_2d(ucell, this->gd, this->pv, PARAM.globalv.gamma_only_local, orb_.cutoffs());

    // 8) initialize the Hamiltonian operators
    // if atom moves, then delete old pointer and add a new one
    if (this->p_hamilt != nullptr)
    {
        delete this->p_hamilt;
        this->p_hamilt = nullptr;
    }
    if (this->p_hamilt == nullptr)
    {
        elecstate::DensityMatrix<TK, double>* DM = estate->get_DM();

        this->p_hamilt = new hamilt::HamiltLCAO<TK, TR>(
            ucell, this->gd, &this->pv, this->pelec->pot, this->kv,
            two_center_bundle_, orb_, DM, this->deepks
#ifdef __EXX
            ,
            istep,
            GlobalC::exx_info.info_ri.real_number ? &this->exx_nao.exd->two_level_step : &this->exx_nao.exc->two_level_step,
            GlobalC::exx_info.info_ri.real_number ? &this->exx_nao.exd->get_Hexxs() : nullptr,
            GlobalC::exx_info.info_ri.real_number ? nullptr : &this->exx_nao.exc->get_Hexxs()
#endif
        );
    }

    // 9) for each ionic step, the overlap <phi|alpha> must be rebuilt
    // since it depends on ionic positions
    this->deepks.build_overlap(ucell, orb_, pv, gd, *(two_center_bundle_.overlap_orb_alpha), PARAM.inp);

    // 10) prepare sc calculation
    if (PARAM.inp.sc_mag_switch)
    {
        spinconstrain::SpinConstrain<TK>& sc = spinconstrain::SpinConstrain<TK>::getScInstance();
        sc.init_sc(PARAM.inp.sc_thr, PARAM.inp.nsc, PARAM.inp.nsc_min, PARAM.inp.alpha_trial,
                   PARAM.inp.sccut, PARAM.inp.sc_drop_thr, ucell, &(this->pv),
                   PARAM.inp.nspin, this->kv, this->p_hamilt, this->psi, this->pelec);
    }

    // 11) set xc type before the first cal of xc in pelec->init_scf, Peize Lin add 2016-12-03
#ifdef __EXX
    if (PARAM.inp.calculation != "nscf")
    {
        if (GlobalC::exx_info.info_ri.real_number)
        {
            this->exx_nao.exd->exx_beforescf(istep, this->kv, *this->p_chgmix, ucell, orb_);
        }
        else
        {
            this->exx_nao.exc->exx_beforescf(istep, this->kv, *this->p_chgmix, ucell, orb_);
        }
    }
#endif

    // 12) init_scf, should be before_scf? mohan add 2025-03-10
    this->pelec->init_scf(istep, ucell, this->Pgrid, this->sf.strucFac, this->locpp.numeric, ucell.symm);

    // 13) initalize DMR
    // DMR should be same size with Hamiltonian(R)
    auto* hamilt_lcao = dynamic_cast<hamilt::HamiltLCAO<TK, TR>*>(this->p_hamilt);
    if(!hamilt_lcao)
    {
        ModuleBase::WARNING_QUIT("ESolver_KS_LCAO::before_scf","p_hamilt does not exist");
    }
    estate->get_DM()->init_DMR(*hamilt_lcao->getHR());

#ifdef __MLALGO
    // 14) initialize DMR of DeePKS
    this->deepks.ld.init_DMR(ucell, orb_, this->pv, this->gd);
#endif

    // 15) two cases are considered:
    // 1. DMK in DensityMatrix is not empty (istep > 0), then DMR is initialized by DMK
    // 2. DMK in DensityMatrix is empty (istep == 0), then DMR is initialized by zeros
    if (istep > 0)
    {
        estate->get_DM()->cal_DMR();
    }

    // 16) the electron charge density should be symmetrized,
    Symmetry_rho srho;
    for (int is = 0; is < PARAM.inp.nspin; is++)
    {
        srho.begin(is, this->chr, this->pw_rho, ucell.symm);
    }

    // 17) why we need to set this sentence? mohan add 2025-03-10
    this->p_hamilt->non_first_scf = istep;

    // 18) update of RDMFT, added by jghan
    if (PARAM.inp.rdmft == true)
    {
        rdmft_solver.update_ion(ucell, *(this->pw_rho), this->locpp.vloc, this->sf.strucFac);
    }

    ModuleBase::timer::tick("ESolver_KS_LCAO", "before_scf");
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

    deepks.dpks_out_type = "tot";  // for deepks method

    fsl.getForceStress(ucell, PARAM.inp.cal_force, PARAM.inp.cal_stress, 
                       PARAM.inp.test_force, PARAM.inp.test_stress,
                       this->gd, this->pv, this->pelec, this->psi,
                       two_center_bundle_, orb_, force, this->scs,
                       this->locpp, this->sf, this->kv,
                       this->pw_rho, this->solvent, this->deepks,
                       this->exx_nao, &ucell.symm);

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

template <typename TK, typename TR>
void ESolver_KS_LCAO<TK, TR>::after_all_runners(UnitCell& ucell)
{
    ModuleBase::TITLE("ESolver_KS_LCAO", "after_all_runners");
    ModuleBase::timer::tick("ESolver_KS_LCAO", "after_all_runners");

    ESolver_KS<TK>::after_all_runners(ucell);

    const int nspin0 = (PARAM.inp.nspin == 2) ? 2 : 1;

	auto* estate = dynamic_cast<elecstate::ElecStateLCAO<TK>*>(this->pelec);
    auto* hamilt_lcao = dynamic_cast<hamilt::HamiltLCAO<TK, TR>*>(this->p_hamilt);

	if(!estate)
	{
		ModuleBase::WARNING_QUIT("ESolver_KS_LCAO::after_all_runners","pelec does not exist");
	}

	if(!hamilt_lcao)
	{
		ModuleBase::WARNING_QUIT("ESolver_KS_LCAO::after_all_runners","p_hamilt does not exist");
	}

	ModuleIO::ctrl_runner_lcao<TK, TR>(ucell,
		  PARAM.inp, this->kv, estate, this->pv, this->Pgrid, 
		  this->gd, this->psi, this->chr, hamilt_lcao,
          this->two_center_bundle_,
          this->orb_, this->pw_rho, this->pw_rhod,
          this->sf, this->locpp.vloc, this->exx_nao, this->solvent);

    ModuleBase::timer::tick("ESolver_KS_LCAO", "after_all_runners");
}

template <typename TK, typename TR>
void ESolver_KS_LCAO<TK, TR>::iter_init(UnitCell& ucell, const int istep, const int iter)
{
    ModuleBase::TITLE("ESolver_KS_LCAO", "iter_init");

    // call iter_init() of ESolver_KS
    ESolver_KS<TK>::iter_init(ucell, istep, iter);

    // cast pointers
	auto* estate = dynamic_cast<elecstate::ElecStateLCAO<TK>*>(this->pelec);
	if(!estate)
	{
		ModuleBase::WARNING_QUIT("ESolver_KS_LCAO::iter_init","pelec does not exist");
	}

	elecstate::DensityMatrix<TK, double>* dm = estate->get_DM();

    module_charge::chgmixing_ks_lcao(iter, this->p_chgmix, dm->get_DMR_pointer(1)->get_nnr(), PARAM.inp); 

    // mohan update 2012-06-05
    estate->f_en.deband_harris = estate->cal_delta_eband(ucell);

    if (istep == 0 && PARAM.inp.init_wfc == "file")
	{
		int exx_two_level_step = 0;
#ifdef __EXX
		if (GlobalC::exx_info.info_global.cal_exx)
		{
			// the following steps are only needed in the first outer exx loop
			exx_two_level_step
				= GlobalC::exx_info.info_ri.real_number ? this->exx_nao.exd->two_level_step : this->exx_nao.exc->two_level_step;
		}
#endif
		elecstate::setup_dm<TK>(ucell, estate, this->psi, this->chr, iter, exx_two_level_step);
	}

#ifdef __EXX
    // calculate exact-exchange
    if (PARAM.inp.calculation != "nscf")
    {
        if (GlobalC::exx_info.info_ri.real_number)
        {
            this->exx_nao.exd->exx_eachiterinit(istep, ucell, *dm, this->kv, iter);
        }
        else
        {
            this->exx_nao.exc->exx_eachiterinit(istep, ucell, *dm, this->kv, iter);
        }
    }
#endif

    if (PARAM.inp.dft_plus_u)
    {
        if (istep != 0 || iter != 1)
        {
            GlobalC::dftu.set_dmr(dm);
        }
        // Calculate U and J if Yukawa potential is used
        GlobalC::dftu.cal_slater_UJ(ucell, this->chr.rho, this->pw_rho->nrxx);
    }

#ifdef __MLALGO
    // the density matrixes of DeePKS have been updated in each iter
    this->deepks.ld.set_hr_cal(true);

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
        dm->save_DMR();
    }
}

template <typename TK, typename TR>
void ESolver_KS_LCAO<TK, TR>::hamilt2rho_single(UnitCell& ucell, int istep, int iter, double ethr)
{
    ModuleBase::TITLE("ESolver_KS_LCAO", "hamilt2rho_single");

    // 1) reset energy
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
            this->exx_nao.exd->exx_hamilt2rho(*this->pelec, this->pv, iter);
        }
        else
        {
            this->exx_nao.exc->exx_hamilt2rho(*this->pelec, this->pv, iter);
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
void ESolver_KS_LCAO<TK, TR>::iter_finish(UnitCell& ucell, const int istep, int& iter, bool& conv_esolver)
{
    ModuleBase::TITLE("ESolver_KS_LCAO", "iter_finish");

    auto* estate = dynamic_cast<elecstate::ElecStateLCAO<TK>*>(this->pelec);
    auto* hamilt_lcao = dynamic_cast<hamilt::HamiltLCAO<TK, TR>*>(this->p_hamilt);

    if(!estate)
    {
        ModuleBase::WARNING_QUIT("ESolver_KS_LCAO::iter_finish","pelec does not exist");
    }

    if(!hamilt_lcao)
    {
        ModuleBase::WARNING_QUIT("ESolver_KS_LCAO::iter_finish","p_hamilt does not exist");
    }

	const std::vector<std::vector<TK>>& dm_vec = estate->get_DM()->get_DMK_vector();

    // 1) calculate the local occupation number matrix and energy correction in DFT+U
    if (PARAM.inp.dft_plus_u)
    {
        // only old DFT+U method should calculated energy correction in esolver,
        // new DFT+U method will calculate energy in calculating Hamiltonian
        if (PARAM.inp.dft_plus_u == 2)
        {
            if (GlobalC::dftu.omc != 2)
            {
                ModuleDFTU::dftu_cal_occup_m(iter, ucell, dm_vec, this->kv,
                  this->p_chgmix->get_mixing_beta(), hamilt_lcao);
            }
            GlobalC::dftu.cal_energy_correction(ucell, istep);
        }
        GlobalC::dftu.output(ucell);
    }

    // 2) for deepks, calculate delta_e, output labels during electronic steps
#ifdef __MLALGO
    if (PARAM.inp.deepks_scf)
    {
        this->deepks.ld.dpks_cal_e_delta_band(dm_vec, this->kv.get_nks());
        DeePKS_domain::update_dmr(this->kv.kvec_d, dm_vec, ucell, orb_, this->pv, this->gd, this->deepks.ld.dm_r);
        estate->f_en.edeepks_scf = this->deepks.ld.E_delta - this->deepks.ld.e_delta_band;
        estate->f_en.edeepks_delta = this->deepks.ld.E_delta;
    }
#endif

    // 3) for delta spin
    if (PARAM.inp.sc_mag_switch)
    {
        spinconstrain::SpinConstrain<TK>& sc = spinconstrain::SpinConstrain<TK>::getScInstance();
        sc.cal_mi_lcao(iter);
    }

    // call iter_finish() of ESolver_KS, where band gap is printed,
    // eig and occ are printed, magnetization is calculated,
    // charge mixing is performed, potential is updated, 
    // HF and kS energies are computed, meta-GGA, Jason and restart
    ESolver_KS<TK>::iter_finish(ucell, istep, iter, conv_esolver);

    // mix density matrix if mixing_restart + mixing_dmr + not first
    // mixing_restart at every iter except the last iter
    if(iter != PARAM.inp.scf_nmax && !conv_esolver)
    {
        if (PARAM.inp.mixing_restart > 0 && this->p_chgmix->mixing_restart_count > 0 && PARAM.inp.mixing_dmr)
        {
            elecstate::DensityMatrix<TK, double>* dm = estate->get_DM();
            this->p_chgmix->mix_dmr(dm);
        }
    }

    // use the converged occupation matrix for next MD/Relax SCF calculation
    if (PARAM.inp.dft_plus_u && conv_esolver)
    {
        GlobalC::dftu.initialed_locale = true;
    }

    // control the output related to the finished iteration
    ModuleIO::ctrl_iter_lcao<TK, TR>(ucell, PARAM.inp, this->kv, estate,
      this->pv, this->gd, this->psi, this->chr, this->p_chgmix, 
      hamilt_lcao, this->orb_, this->deepks, 
      this->exx_nao, iter, istep, conv_esolver, this->scf_ene_thr);

}

template <typename TK, typename TR>
void ESolver_KS_LCAO<TK, TR>::after_scf(UnitCell& ucell, const int istep, const bool conv_esolver)
{
    ModuleBase::TITLE("ESolver_KS_LCAO", "after_scf");
    ModuleBase::timer::tick("ESolver_KS_LCAO", "after_scf");

    //! 1) call after_scf() of ESolver_KS
    ESolver_KS<TK>::after_scf(ucell, istep, conv_esolver);

    auto* estate = dynamic_cast<elecstate::ElecStateLCAO<TK>*>(this->pelec);
    auto* hamilt_lcao = dynamic_cast<hamilt::HamiltLCAO<TK, TR>*>(this->p_hamilt);

    if(!estate)
    {
        ModuleBase::WARNING_QUIT("ESolver_KS_LCAO::after_scf","pelec does not exist");
    }

    if(!hamilt_lcao)
    {
        ModuleBase::WARNING_QUIT("ESolver_KS_LCAO::after_scf","p_hamilt does not exist");
    }

    //! 2) output of lcao every few ionic steps
    ModuleIO::ctrl_scf_lcao<TK, TR>(ucell,
            PARAM.inp, this->kv, estate, this->pv,
            this->gd, this->psi, hamilt_lcao,
            this->two_center_bundle_,
            this->orb_, this->pw_wfc, this->pw_rho,
            this->pw_big, this->sf,
            this->rdmft_solver, this->deepks, this->exx_nao, 
            this->conv_esolver, this->scf_nmax_flag,
            istep);


    //! 3) Clean up RA, which is used to serach for adjacent atoms
    if (!PARAM.inp.cal_force && !PARAM.inp.cal_stress)
    {
        this->RA.delete_grid();
    }

    ModuleBase::timer::tick("ESolver_KS_LCAO", "after_scf");
}


template class ESolver_KS_LCAO<double, double>;
template class ESolver_KS_LCAO<std::complex<double>, double>;
template class ESolver_KS_LCAO<std::complex<double>, std::complex<double>>;
} // namespace ModuleESolver
