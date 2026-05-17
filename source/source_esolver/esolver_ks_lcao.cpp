#include "esolver_ks_lcao.h"
#include "source_base/module_external/blacs_connector.h"
#include "source_cell/module_neighbor/sltk_atom_arrange.h"
#include "source_estate/elecstate_tools.h"
#include "source_lcao/module_deltaspin/spin_constrain.h"
#include "source_lcao/module_deltaspin/deltaspin_lcao.h"
#include "source_lcao/dftu_lcao.h"
#include "source_lcao/hs_matrix_k.hpp" // there may be multiple definitions if using hpp
#include "source_estate/module_charge/symmetry_rho.h"
#include "source_lcao/LCAO_domain.h" // need DeePKS_init
#include "source_lcao/FORCE_STRESS.h"
#include "source_estate/elecstate_lcao.h"
#include "source_lcao/hamilt_lcao.h"
#include "source_hsolver/hsolver_lcao.h"
#ifdef __EXX
#include "../source_lcao/module_ri/exx_opt_orb.h"
#endif
#include "source_lcao/module_rdmft/rdmft.h"
#include "source_estate/module_charge/chgmixing.h" // use charge mixing, mohan add 20251006
#include "source_estate/module_dm/init_dm.h" // init dm from electronic wave functions
#include "source_io/module_ctrl/ctrl_runner_lcao.h" // use ctrl_runner_lcao() 
#include "source_io/module_ctrl/ctrl_iter_lcao.h" // use ctrl_iter_lcao() 
#include "source_io/module_ctrl/ctrl_scf_lcao.h" // use ctrl_scf_lcao()
#include "source_io/module_output/print_info.h"
#include "source_lcao/rho_tau_lcao.h" // mohan add 20251024
#include "source_lcao/LCAO_set.h" // mohan add 20251111
#include "source_psi/setup_psi.h" // use Setup_Psi for deallocate_psi

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
    Setup_Psi<TK>::deallocate_psi(this->psi);
}

template <typename TK, typename TR>
void ESolver_KS_LCAO<TK, TR>::before_all_runners(UnitCell& ucell, const Input_para& inp)
{
    ModuleBase::TITLE("ESolver_KS_LCAO", "before_all_runners");
    ModuleBase::timer::start("ESolver_KS_LCAO", "before_all_runners");

    // 1) before_all_runners in ESolver_KS
    ESolver_KS::before_all_runners(ucell, inp);

    // 2) autoset nbands in ElecState before init_basis (for Psi 2d division)
    if (this->pelec == nullptr)
    {
        // TK stands for double and std::complex<double>?
        this->pelec = new elecstate::ElecStateLCAO<TK>(&(this->chr), &(this->kv),
          this->kv.get_nks(), this->pw_big);
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

    LCAO_domain::set_psi_occ_dm_chg<TK>(this->kv, this->psi, this->pv, this->pelec,
      this->dmat, this->chr, inp);

    LCAO_domain::set_pot<TK>(ucell, this->kv, this->sf, *this->pw_rho, *this->pw_rhod,
      this->pelec, this->orb_, this->pv, this->locpp, this->dftu,
      this->solvent, this->exx_nao, this->deepks, inp);

    //! if kpar is not divisible by nks, print a warning
    ModuleIO::print_kpar(this->kv.get_nks(), PARAM.globalv.kpar_lcao);

    //! init rdmft, added by jghan
    if (inp.rdmft == true)
    {
        rdmft_solver.init(this->pv, ucell,
          this->gd, this->kv, *(this->pelec), this->orb_,
          two_center_bundle_, inp.dft_functional, inp.rdmft_power_alpha);
    }

    ModuleBase::timer::end("ESolver_KS_LCAO", "before_all_runners");
    return;
}


template <typename TK, typename TR>
void ESolver_KS_LCAO<TK, TR>::before_scf(UnitCell& ucell, const int istep)
{
    ModuleBase::TITLE("ESolver_KS_LCAO", "before_scf");
    ModuleBase::timer::start("ESolver_KS_LCAO", "before_scf");

    //! 1) call before_scf() of ESolver_KS.
    ESolver_KS::before_scf(ucell, istep);

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
        this->p_hamilt = new hamilt::HamiltLCAO<TK, TR>(
            ucell, this->gd, &this->pv, this->pelec->pot, this->kv,
            two_center_bundle_, orb_, this->dmat.dm, &this->dftu, this->deepks, istep, exx_nao);
    }

    // 9) for each ionic step, the overlap <phi|alpha> must be rebuilt
    // since it depends on ionic positions
    this->deepks.build_overlap(ucell, orb_, pv, gd, *(two_center_bundle_.overlap_orb_alpha), PARAM.inp);

    // 10) prepare sc calculation
    init_deltaspin_lcao<TK>(ucell, PARAM.inp, &(this->pv), this->kv, this->p_hamilt, this->psi, this->dmat.dm, this->pelec);

    // 11) set xc type before the first cal of xc in pelec->init_scf, Peize Lin add 2016-12-03
    this->exx_nao.before_scf(ucell, this->kv, orb_, this->p_chgmix, istep, PARAM.inp);

    // 12) initalize DM(R), which has the same size with Hamiltonian(R)
    auto* hamilt_lcao = dynamic_cast<hamilt::HamiltLCAO<TK, TR>*>(this->p_hamilt);

    if(!hamilt_lcao)
    {
        ModuleBase::WARNING_QUIT("ESolver_KS_LCAO::before_scf","p_hamilt does not exist");
    }
    this->dmat.dm->init_DMR(*hamilt_lcao->getHR());

    // 13.1) decide the strategy for initializing DMR and HR
    if(istep == 0)//if the first scf step, readin DMR from file,
    {
        //calculate or readin the density matrix DMR
        if(PARAM.inp.init_chg == "dm")
        {
            //! 13.1.1) init charge density from density matrix file
            LCAO_domain::init_chg_dm<TK>(PARAM.globalv.global_readin_dir, PARAM.inp.nspin,
                this->dmat, ucell, &(this->pv), this->pelec->charge);
        }
        if(PARAM.inp.init_chg == "hr")
        {
            //! 13.1.2) init charge density from Hamiltonian matrix file
            LCAO_domain::init_chg_hr<TK, TR>(PARAM.globalv.global_readin_dir, PARAM.inp.nspin,
                static_cast<hamilt::Hamilt<TK>*>(this->p_hamilt), ucell, &(this->pv), this->psi[0], this->pelec, *this->dmat.dm,
                this->chr, PARAM.inp.ks_solver);
        }
    }
    else if(PARAM.inp.esolver_type!="tddft")//if not, use the DMR calculated from last step
    {
        // 13.1.2) two cases are considered:
        // 1. DMK in DensityMatrix is not empty (istep > 0), then DMR is initialized by DMK
        // 2. DMK in DensityMatrix is empty (istep == 0), then DMR is initialized by zeros
        this->dmat.dm->cal_DMR();
    }
    // 13.2) init_scf, should be before_scf? mohan add 2025-03-10
    this->pelec->init_scf(ucell, this->Pgrid, this->sf.strucFac, this->locpp.numeric, ucell.symm);

#ifdef __MLALGO
    // 14) initialize DM2(R) of DeePKS, the DM2(R) is different from DM(R)
    this->deepks.ld.init_DMR(ucell, orb_, this->pv, this->gd);
#endif

    // 16) the electron charge density should be symmetrized,
    Symmetry_rho::symmetrize_rho(PARAM.inp.nspin, this->chr, this->pw_rho, ucell.symm);

    // 17) update of RDMFT, added by jghan
    if (PARAM.inp.rdmft == true)
    {
        rdmft_solver.update_ion(ucell, *(this->pw_rho), this->locpp.vloc, this->sf.strucFac);
    }

    ModuleBase::timer::end("ESolver_KS_LCAO", "before_scf");
    return;
}


template <typename TK, typename TR>
double ESolver_KS_LCAO<TK, TR>::cal_energy()
{
    return this->pelec->f_en.etot;
}

template <typename TK, typename TR>
void ESolver_KS_LCAO<TK, TR>::cal_force(UnitCell& ucell, ModuleBase::matrix& force)
{
    ModuleBase::TITLE("ESolver_KS_LCAO", "cal_force");
    ModuleBase::timer::start("ESolver_KS_LCAO", "cal_force");

    Force_Stress_LCAO<TK> fsl(this->RA, ucell.nat);

    deepks.dpks_out_type = "tot";  // for deepks method

    fsl.getForceStress(ucell, PARAM.inp.cal_force, PARAM.inp.cal_stress, 
                       PARAM.inp.test_force, PARAM.inp.test_stress,
                       this->gd, this->pv, this->pelec, this->dmat, this->psi,
                       two_center_bundle_, orb_, force, this->scs,
                       this->locpp, this->sf, this->kv,
                       this->pw_rho, this->solvent, this->dftu, this->deepks,
                       this->exx_nao, &ucell.symm);

    // delete RA after cal_force
    this->RA.delete_grid();

    this->have_force = true;

    ModuleBase::timer::end("ESolver_KS_LCAO", "cal_force");
}

template <typename TK, typename TR>
void ESolver_KS_LCAO<TK, TR>::cal_stress(UnitCell& ucell, ModuleBase::matrix& stress)
{
    ModuleBase::TITLE("ESolver_KS_LCAO", "cal_stress");
    ModuleBase::timer::start("ESolver_KS_LCAO", "cal_stress");

    if (!this->have_force)
    {
        ModuleBase::matrix fcs;
        this->cal_force(ucell, fcs);
    }

    // the stress has been calculated in 'cal_force'
    stress = this->scs;
    this->have_force = false;

    ModuleBase::timer::end("ESolver_KS_LCAO", "cal_stress");
}

template <typename TK, typename TR>
void ESolver_KS_LCAO<TK, TR>::after_all_runners(UnitCell& ucell)
{
    ModuleBase::TITLE("ESolver_KS_LCAO", "after_all_runners");
    ModuleBase::timer::start("ESolver_KS_LCAO", "after_all_runners");

    ESolver_KS::after_all_runners(ucell);

    auto* hamilt_lcao = dynamic_cast<hamilt::HamiltLCAO<TK, TR>*>(this->p_hamilt);
    if(!hamilt_lcao)
    {
	    ModuleBase::WARNING_QUIT("ESolver_KS_LCAO::after_all_runners","p_hamilt does not exist");
    }

    ModuleIO::ctrl_runner_lcao<TK, TR>(ucell,
		    PARAM.inp, this->kv, this->pelec, this->dmat, this->pv, this->Pgrid, 
		    this->gd, this->psi, this->chr, hamilt_lcao,
		    this->two_center_bundle_,
		    this->orb_, this->pw_rho, this->pw_rhod,
		    this->sf, this->locpp.vloc, this->exx_nao, this->solvent);


#ifdef __MPI
#ifdef __LCAO
    // Exit BLACS environment for LCAO calculations
    Cblacs_exit(1);
#endif
#endif

    ModuleBase::timer::end("ESolver_KS_LCAO", "after_all_runners");
}

template <typename TK, typename TR>
void ESolver_KS_LCAO<TK, TR>::iter_init(UnitCell& ucell, const int istep, const int iter)
{
    ModuleBase::TITLE("ESolver_KS_LCAO", "iter_init");

    // call iter_init() of ESolver_KS
    ESolver_KS::iter_init(ucell, istep, iter);

    module_charge::chgmixing_ks_lcao(iter, this->p_chgmix, this->dftu, 
      this->dmat.dm->get_DMR_pointer(1)->get_nnr(), PARAM.inp); 

    if (iter == 1)
    {
        this->gint_precision_controller_.set_mode(PARAM.inp.gint_precision);
        this->gint_precision_controller_.reset_for_new_scf();
        this->gint_info_->set_exec_precision(this->gint_precision_controller_.current_precision());
        if (PARAM.inp.gint_precision == "mix")
        {
            GlobalV::ofs_running << "\n >> Gint mixed-precision mode: starting SCF with fp32"
                                 << " (will switch to fp64 when drho is small enough)" << std::endl;
            std::cout << " >> NOTICE: Gint grid-integration starts with fp32 (mixed-precision mode)" << std::endl;
        }
        else if (PARAM.inp.gint_precision == "single")
        {
            GlobalV::ofs_running << "\n >> Gint single-precision mode: using fp32 throughout SCF" << std::endl;
            std::cout << " >> NOTICE: Gint grid-integration uses fp32 throughout SCF (single-precision mode)" << std::endl;
        }
    }

    // mohan update 2012-06-05
    this->pelec->f_en.deband_harris = this->pelec->cal_delta_eband(ucell);

    if (istep == 0 && PARAM.inp.init_wfc == "file")
	{
		int exx_two_level_step = 0;
#ifdef __EXX
		if (GlobalC::exx_info.info_global.cal_exx)
		{
			// the following steps are only needed in the first outer exx loop
			exx_two_level_step
				= GlobalC::exx_info.info_ri.real_number ? 
                  this->exx_nao.exd->two_level_step : this->exx_nao.exc->two_level_step;
		}
#endif
		elecstate::init_dm<TK>(ucell, this->pelec, this->dmat, this->psi, this->chr, iter, exx_two_level_step);
	}

#ifdef __EXX
    // calculate exact-exchange
    if (PARAM.inp.calculation != "nscf")
    {
        if (GlobalC::exx_info.info_ri.real_number)
        {
            this->exx_nao.exd->exx_eachiterinit(istep, ucell, *this->dmat.dm, this->kv, iter);
        }
        else
        {
            this->exx_nao.exc->exx_eachiterinit(istep, ucell, *this->dmat.dm, this->kv, iter);
        }
    }
#endif

    init_dftu_lcao<TK>(istep, iter, PARAM.inp, &(this->dftu), this->dmat.dm, ucell, this->chr.rho, this->pw_rho->nrxx);

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
        this->dmat.dm->save_DMR();
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

    // =====================================================================
    // 2) DeltaSpin: inner lambda loop to constrain atomic magnetic moments
    // =====================================================================
    // The DeltaSpin method implements constrained LSDA via Lagrange multipliers:
    //   E'[rho] = E[rho] - sum_i lambda_i . (M_i - M_target_i)
    //
    // The constrained energy functional adds a penalty term that drives each
    // atom's magnetic moment M_i toward its target value M_target_i.
    // The Lagrange multiplier lambda_i acts as a "magnetic force" (eV/uB).
    //
    // Code paths by spin type:
    // ---------------------------------------------------------------
    // nspin=2 (collinear):
    //   - Only z-component of magnetization is constrained (M_z per atom)
    //   - H_delta = lambda_z * sigma_z (diagonal, opposite sign per spin channel)
    //   - DMR: uses switch_dmr(2) -> spin-difference density (rho_up - rho_dn)
    //   - cal_coeff_lambda: coefficients[spin] = +/- lambda_z
    //
    // nspin=4 (non-collinear):
    //   - Full 3D magnetization vector constrained (Mx, My, Mz per atom)
    //   - H_delta = lambda . sigma (full 2x2 Pauli matrix with spin-flip terms)
    //   - DMR: spinor density matrix (2x2 blocks interleaved)
    //   - cal_coeff_lambda: 4 coeffs for 2x2 spinor block
    //
    // direction_only mode:
    //   - Designed for non-collinear: removes parallel lambda component so
    //     only transverse (directional) constraint remains
    //   - CRITICAL: for nspin=2, direction_only projects lambda to ZERO because
    //     the only constrained direction (z) IS the parallel direction.
    //     Therefore direction_only MUST be disabled during Phase 1 BFGS.
    //
    // sc_scf_thr parameter:
    //   - Threshold at which DeltaSpin lambda loop begins (default 1e-3)
    //   - Should be 10-100x larger than scf_thr so lambda loop starts when
    //     charge density is "reasonably stable" but not fully converged
    //   - Recommended: sc_scf_thr = 1e-4 to 1e-3, scf_thr = 1e-8
    //   - mixing_restart is auto-set to sc_scf_thr for DeltaSpin calculations
    // =====================================================================
    bool skip_solve = false;
    if (PARAM.inp.sc_mag_switch)
    {
        spinconstrain::SpinConstrain<TK>& sc = spinconstrain::SpinConstrain<TK>::getScInstance();

        if (PARAM.inp.sc_lambda_strategy == "linear_scan")
        {
            // Linear scan mode: sweep lambda values for E(lambda) mapping (debug/diagnostic)
            sc.run_lambda_linear_scan(iter - 1);
            skip_solve = true;
        }
        else if (PARAM.inp.sc_direction_only)
        {
            // ================================================================
            // direction_only mode: two-phase constraint strategy
            // ================================================================
            //
            // PROBLEM with original approach (fixed lambda = 500 eV/uB):
            //   - Over-polarizes s/p orbitals, giving unphysical Mi ~ 13 uB
            //   - Prevents natural relaxation to magnetic ground state
            //   - Cannot converge (lambda never decays)
            //
            // SOLUTION: two-phase BFGS + gradual lambda decay
            //
            // Phase 1 (iter 1-5): BFGS constrains moment MAGNITUDE to target
            //   - direction_only projection is DISABLED so BFGS can optimize
            //     the z-component of lambda (the only constrained direction
            //     in collinear mode)
            //   - BFGS finds physically reasonable lambda values (~0.3 eV/uB)
            //     that give Mi ≈ target (±2 uB for Fe)
            //   - skip_solve = true: BFGS inner loop does its own subspace
            //     diagonalization; outer HSolver is skipped
            //
            // Phase 2 (iter > 5): normal SCF with gradual lambda decay
            //   - Lambda decays by factor 0.5^(1/3) per step (~halves every 3 steps)
            //   - Broyden mixing history is reset at transition to avoid charge
            //     density discontinuity between BFGS inner loop state and outer SCF
            //   - skip_solve = false: outer HSolver runs, charge mixing converges
            //   - System naturally relaxes to correct magnetic ground state
            //     (Mi ≈ ±2.33 uB for this Fe AFM system)
            //
            // Phase transition threshold (default: 5 iterations):
            //   - 3-5 iterations is enough for BFGS to establish AFM order
            //   - Too few (< 3): lambda not optimized, poor initial state
            //   - Too many (> 10): waste of iterations, system already stable
            // ================================================================
            if (iter <= 5)
            {
                // Phase 1: BFGS constraint loop
                // CRITICAL: disable direction_only projection for collinear mode.
                // Without this, the projection removes the z-component of lambda
                // (since target direction = z), making lambda = 0 and disabling
                // the constraint entirely.
                sc.set_direction_only(false);
                sc.run_lambda_loop(iter - 1);
                sc.set_direction_only(true);
                skip_solve = true; // BFGS inner loop handles diagonalization
            }
            else
            {
                // Phase 2: gradual lambda decay, normal SCF
                if (iter == 6)
                {
                    // CRITICAL: reset Broyden mixing history at Phase 1->2 transition.
                    //
                    // Why this is needed:
                    //   - Phase 1 BFGS uses subspace diagonalization to update the
                    //     density matrix (DM) directly, without going through the
                    //     outer SCF charge mixing loop
                    //   - The Broyden mixer stores a history of previous density
                    //     differences and Jacobian approximations
                    //   - These history entries are based on BFGS-updated DM states,
                    //     which are incompatible with the outer SCF charge density
                    //     produced by HSolver
                    //   - If we don't reset, the mixer will combine irrelevant history
                    //     with the new SCF state, causing large oscillations
                    //
                    // mix_reset() clears:
                    //   - Broyden/Pulay mixing history (stored density differences)
                    //   - rho_mdata mixing data structure
                    //   - tau_mdata (if meta-GGA)
                    //   Effectively "forgets" all previous SCF history.
                    this->p_chgmix->mix_reset();
                }

                // Gradual lambda decay: factor = 0.5^(1/3) per step
                // This means lambda halves every 3 iterations.
                //
                // Why gradual (not instant zero)?
                //   - Instant removal of lambda causes discontinuous change in
                //     the Hamiltonian, leading to large charge density oscillations
                //   - Gradual decay lets the system adiabatically relax from the
                //     constrained state to the unconstrained ground state
                //
                // Decay rate tuning:
                //   - Halves every 2 steps: faster decay, may oscillate
                //   - Halves every 3 steps: recommended default (current)
                //   - Halves every 5 steps: slower, more stable but longer SCF
                int nat = sc.get_nat();
                auto lambda = sc.get_sc_lambda();
                const double DECAY = std::pow(0.5, 1.0 / 3.0); // ~0.7937 per step
                for (int ia = 0; ia < nat; ++ia)
                    for (int ic = 0; ic < 3; ++ic)
                        lambda[ia][ic] *= DECAY;
                sc.set_lambda(lambda);

                // Normal SCF: skip_solve = false
                // HSolver will run with the decayed lambda, and charge mixing
                // will converge the density naturally.
            }
        }
        else if (!sc.mag_converged() && this->drho > 0 && this->drho < PARAM.inp.sc_scf_thr)
        {
            // Standard DeltaSpin: first lambda loop invocation.
            // Triggered when charge density has converged to within sc_scf_thr.
            //
            // sc_scf_thr recommended settings:
            //   - Default: 1e-3 (100,000x larger than default scf_thr = 1e-8)
            //   - For difficult systems: 1e-4 to 1e-5
            //   - Should be large enough that charge density is "reasonably stable"
            //     but small enough that lambda optimization is meaningful
            //
            // mixing_restart is auto-set to sc_scf_thr for DeltaSpin calculations
            // (see read_input_item_elec_stru.cpp), ensuring clean mixing history
            // before the lambda loop starts.
            sc.run_lambda_loop(iter - 1);
            sc.set_mag_converged(true);
            skip_solve = true;
        }
        else if (sc.mag_converged())
        {
            // Standard DeltaSpin: subsequent lambda loop invocations.
            // Lambda is refined for the current charge density at each SCF step.
            sc.run_lambda_loop(iter - 1);
            skip_solve = true;
        }
    }

    // 3) run Hsolver
    if (!skip_solve)
    {
        hsolver::HSolverLCAO<TK> hsolver_lcao_obj(&(this->pv), PARAM.inp.ks_solver);
        hsolver_lcao_obj.solve(static_cast<hamilt::Hamilt<TK>*>(this->p_hamilt), this->psi[0], this->pelec, *this->dmat.dm, 
          this->chr, PARAM.inp.nspin, skip_charge);
    }
    else
    {
        // Lambda loop updated the density matrix (DM) but not the real-space charge density.
        // HSolver was skipped, so we need to sync rho from DM manually.
        LCAO_domain::dm2rho(this->dmat.dm->get_DMR_vector(), PARAM.inp.nspin, &this->chr);
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
    Symmetry_rho::symmetrize_rho(PARAM.inp.nspin, this->chr, this->pw_rho, ucell.symm);

    // 6) calculate delta energy
    this->pelec->f_en.deband = this->pelec->cal_delta_eband(ucell);
}


template <typename TK, typename TR>
void ESolver_KS_LCAO<TK, TR>::iter_finish(UnitCell& ucell, const int istep, int& iter, bool& conv_esolver)
{
    ModuleBase::TITLE("ESolver_KS_LCAO", "iter_finish");

    auto* hamilt_lcao = dynamic_cast<hamilt::HamiltLCAO<TK, TR>*>(this->p_hamilt);

    if(!hamilt_lcao)
    {
        ModuleBase::WARNING_QUIT("ESolver_KS_LCAO::iter_finish","p_hamilt does not exist");
    }

	const std::vector<std::vector<TK>>& dm_vec = this->dmat.dm->get_DMK_vector();

    // 1) calculate the local occupation number matrix and energy correction in DFT+U
    finish_dftu_lcao<TK>(iter, conv_esolver, PARAM.inp, &(this->dftu), ucell, dm_vec, this->kv, this->p_chgmix->get_mixing_beta(), hamilt_lcao);

    // 2) for deepks, calculate delta_e, output labels during electronic steps
    this->deepks.delta_e(ucell, this->kv, this->orb_, this->pv, this->gd, dm_vec, this->pelec->f_en, PARAM.inp);

    // 3) for delta spin
    cal_mi_lcao_wrapper<TK>(iter, PARAM.inp);

    // 3b) direction_only: report magnetic moment status
    if (PARAM.inp.sc_direction_only && PARAM.inp.sc_mag_switch)
    {
        spinconstrain::SpinConstrain<TK>& sc = spinconstrain::SpinConstrain<TK>::getScInstance();
        const int nat = sc.get_nat();
        const auto& Mi = sc.get_Mi();
        const auto& target = sc.get_target_mag();
        const auto& constrain = sc.get_constrain();
        auto lambda = sc.get_sc_lambda();

        double lambda_abs = 0;
        for (int ia = 0; ia < nat; ++ia)
            for (int ic = 0; ic < 3; ++ic)
                if (constrain[ia][ic] != 0)
                    lambda_abs += std::abs(lambda[ia][ic] * ModuleBase::Ry_to_eV);

        GlobalV::ofs_running << " [DS-dir] iter " << iter << "  |lambda|=" << lambda_abs << " eV/uB" << std::endl;
        for (int ia = 0; ia < nat; ++ia)
            for (int ic = 0; ic < 3; ++ic)
                if (constrain[ia][ic] != 0)
                    GlobalV::ofs_running << "   Atom " << ia << " comp " << ic << ": Mi=" << Mi[ia][ic]
                                         << "  T=" << target[ia][ic] << std::endl;
    }

    // call iter_finish() of ESolver_KS, where band gap is printed,
    // eig and occ are printed, magnetization is calculated,
    // charge mixing is performed, potential is updated, 
    // HF and kS energies are computed, meta-GGA, Jason and restart
    ESolver_KS::iter_finish(ucell, istep, iter, conv_esolver);
    const bool precision_switched = this->gint_precision_controller_.update_after_iteration(this->drho, this->scf_thr);
    this->gint_info_->set_exec_precision(this->gint_precision_controller_.current_precision());
    if (precision_switched)
    {
        GlobalV::ofs_running << "\n >> Gint precision switched: fp32 -> fp64 (drho = "
                             << this->drho << ")" << std::endl;
        std::cout << " >> NOTICE: Gint grid-integration precision switched from fp32 to fp64" << std::endl;
    }

    // mix density matrix if mixing_restart + mixing_dmr + not first
    // mixing_restart at every iter except the last iter
    if(iter != PARAM.inp.scf_nmax && !conv_esolver)
    {
        if (PARAM.inp.mixing_restart > 0 && this->p_chgmix->mixing_restart_count > 0 && PARAM.inp.mixing_dmr)
        {
            this->p_chgmix->mix_dmr(this->dmat.dm);
        }
    }

    // control the output related to the finished iteration
    ModuleIO::ctrl_iter_lcao<TK, TR>(ucell, PARAM.inp, this->kv, this->pelec, *this->dmat.dm,
      this->pv, this->gd, this->psi, this->chr, this->p_chgmix, 
      hamilt_lcao, this->orb_, this->deepks, 
      this->exx_nao, iter, istep, conv_esolver, this->scf_ene_thr);
}

template <typename TK, typename TR>
void ESolver_KS_LCAO<TK, TR>::after_scf(UnitCell& ucell, const int istep, const bool conv_esolver)
{
    ModuleBase::TITLE("ESolver_KS_LCAO", "after_scf");
    ModuleBase::timer::start("ESolver_KS_LCAO", "after_scf");

    auto* hamilt_lcao = dynamic_cast<hamilt::HamiltLCAO<TK, TR>*>(this->p_hamilt);

    if(!hamilt_lcao)
    {
        ModuleBase::WARNING_QUIT("ESolver_KS_LCAO::after_scf","p_hamilt does not exist");
    }

    if (PARAM.inp.out_elf[0] > 0)
	{
		LCAO_domain::dm2tau(this->dmat.dm->get_DMR_vector(), PARAM.inp.nspin, this->pelec->charge);
	}

    //! 1) call after_scf() of ESolver_KS
    ESolver_KS::after_scf(ucell, istep, conv_esolver);

    //! 2) output of lcao every few ionic steps
    ModuleIO::ctrl_scf_lcao<TK, TR>(ucell,
            PARAM.inp, this->kv, this->pelec, this->dmat.dm, this->pv,
            this->gd, this->psi, hamilt_lcao, this->dftu, this->two_center_bundle_,
            this->orb_, this->pw_wfc, this->pw_rho, this->pw_big, this->sf,
            this->rdmft_solver, this->deepks, this->exx_nao,
            this->conv_esolver, this->scf_nmax_flag, istep);

    //! 3) Clean up RA, which is used to serach for adjacent atoms
    if (!PARAM.inp.cal_force && !PARAM.inp.cal_stress)
    {
        this->RA.delete_grid();
    }

    ModuleBase::timer::end("ESolver_KS_LCAO", "after_scf");
}

template class ESolver_KS_LCAO<double, double>;
template class ESolver_KS_LCAO<std::complex<double>, double>;
template class ESolver_KS_LCAO<std::complex<double>, std::complex<double>>;
} // namespace ModuleESolver
