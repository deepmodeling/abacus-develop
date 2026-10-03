#include "esolver_lr_lcao_tddft.h"
#include <cstdlib>

#include <algorithm>
#include <iomanip>
#include "source_lcao/module_lr/utils/lr_io.h"
#include "source_lcao/module_lr/utils/lr_util.h"
#include "source_lcao/module_lr/hamilt_casida.h"
#include "source_lcao/module_lr/hamilt_ulr.hpp"
#include "source_lcao/module_lr/potentials/pot_hxc_lrtd.h"
#include "source_lcao/lcao_nonlocal_info.h"
#include "source_hamilt/module_xc/xc_functional.h"
#include "source_lcao/module_lr/hsolver_lrtd.hpp"
#include "source_lcao/module_lr/lr_spectrum.h"
#include "source_lcao/module_lr/lr_density.hpp"
#include "source_hamilt/module_gint/gint.h"
#include <memory>
#include "source_lcao/hamilt_lcao.h"
#include "source_io/module_wf/read_wfc_nao.h"
#include "source_io/module_output/cube_io.h"
#include "source_io/module_output/print_info.h"
#include "source_cell/module_neighbor/sltk_atom_arrange.h"
#include "source_lcao/module_lr/utils/lr_util_print.h"
#include "source_base/module_external/scalapack_connector.h"
#include "source_io/module_parameter/parameter.h"
#ifdef __JSON
#include "source_io/module_json/output_info.h"
#endif
#include "source_lcao/module_lr/ri_benchmark/ri_benchmark.h"
#include "source_lcao/module_lr/operator_casida/operator_lr_diag.h" // for precondition
#ifdef __EXX
#include "source_lcao/module_ri/exx_lri_interface.h"
#include "source_hamilt/module_xc/exx_info.h" // for init_exx_info
#include "source_hamilt/module_xc/xc_functional.h"   // for set_xc_type
#endif

// gradient
#include "source_lcao/module_lr/Grad/multipliers/zeq_solver.h"

#ifdef __EXX
namespace
{

    /// One `Exx_LRI` carries ONE Coulomb operator, and it may be needed for two different
    /// reasons: the LR kernel (when `xc_kernel` is a hybrid) and the ground-state force (when
    /// `dft_functional` is a hybrid). That operator is NOT chosen here -- `Exx_LRI` reads
    /// `info_ri.coulomb_param`, which `input_conv` builds from `dft_functional` alone. So when
    /// the two disagree, the kernel silently gets the ground state's screening, not its own.
    ///
    /// This used to be decided by a local `exx_ccp_type()` returning Erfc for "hse" and bare
    /// Hf otherwise, written onto `info_global.ccp_type`. With general range-separated hybrids
    /// that two-way split is not even expressible ($\alpha/r+\beta\,\mathrm{erfc}(\mu r)/r$
    /// is both at once), and the write was in any case dead for the RI path: only `Exx_LRI`
    /// runs here, and it never looks at `ccp_type`.
    void warn_if_kernel_differs_from_gs(const std::string& xc_kernel, const std::string& dft_functional,
        std::ofstream& ofs_running)
    {
        const bool k = LR::exx_kernel_list().count(xc_kernel) > 0;
        const bool g = LR::exx_kernel_list().count(dft_functional) > 0;
        if (k && xc_kernel != dft_functional)
        {
            ofs_running << " WARNING: xc_kernel (" << xc_kernel << ") and dft_functional ("
                << dft_functional << ") are not the same functional. The Coulomb operator of"
                " Exx_LRI follows dft_functional" << (g ? "" : ", which is not a hybrid at all"
                " (no coulomb_param, so the LR exchange kernel vanishes)") << "; set both to the"
                " same hybrid to get the kernel you asked for." << std::endl;
        }
    }
}
#endif

#ifdef __EXX
template<>
void ModuleESolver::ESolver_LR<double>::share_exx_lri(std::shared_ptr<Exx_LRI<double>>& exx_ks)
{
    ModuleBase::TITLE("ESolver_LR<double>", "share_exx_lri");
    this->exx_lri = exx_ks;
}
template<>
void ModuleESolver::ESolver_LR<std::complex<double>>::share_exx_lri(std::shared_ptr<Exx_LRI<std::complex<double>>>& exx_ks)
{
    ModuleBase::TITLE("ESolver_LR<complex>", "share_exx_lri");
    this->exx_lri = exx_ks;
}
template<>
void ModuleESolver::ESolver_LR<std::complex<double>>::share_exx_lri(std::shared_ptr<Exx_LRI<double>>& exx_ks)
{
    throw std::runtime_error("ESolver_LR<std::complex<double>>::share_exx_lri: cannot share double to std::complex<double>");
}
template<>
void ModuleESolver::ESolver_LR<double>::share_exx_lri(std::shared_ptr<Exx_LRI<std::complex<double>>>& exx_ks)
{
    throw std::runtime_error("ESolver_LR<double>::share_exx_lri: cannot share std::complex<double> to double");
}
#endif

using namespace LR;

template<typename T, typename TR>
int ModuleESolver::ESolver_LR<T, TR>::cal_nupdown_form_occ(const ModuleBase::matrix& wg)
{   // only for nspin=2
    const int& nk = wg.nr / 2;
    auto occ_sum_k = [&](const int& is, const int& ib)->double { double o = 0.0; for (int ik = 0;ik < nk;++ik) { o += wg(is * nk + ik, ib); } return o;};
    // Sum the occupations of each channel FIRST and round once, instead of rounding band by band
    // and summing the differences. A half-occupied degenerate frontier pair (OH's 2-Pi doublet
    // smears its odd electron as 0.5/0.5 over the two pi_down orbitals) otherwise makes the answer
    // a coin flip: the stored values are 0.5000000052 and 0.4999999947, so one rounds up and one
    // down, and which way they land is pure noise.
    double up = 0.0, dn = 0.0;
    for (int ib = 0;ib < wg.nc;++ib) { up += occ_sum_k(0, ib); dn += occ_sum_k(1, ib); }
    return static_cast<int>(std::lround(up) - std::lround(dn));
}

template<typename T, typename TR>
void ModuleESolver::ESolver_LR<T, TR>::setup_2center_table(TwoCenterBundle& two_center_bundle, LCAO_Orbitals& orb, UnitCell& ucell)
{
    if (this->inp_->vnl_in_h)
    {
        auto* lcao_nl = new LCAONonlocalInfo();
        lcao_nl->setupNonlocal(ucell.ntype, ucell.atoms, this->ofs_running_, orb,
                               this->inp_->basis_type, this->inp_->out_element_info,
                               this->inp_->lspinorb, this->inp_->nspin, this->my_rank_);
        ucell.infoNL.reset(lcao_nl);
        two_center_bundle.build_beta(ucell.ntype, lcao_nl->get_nonlocal().get_Beta_data());
    }
    // NOTE: tabulate() must be called AFTER build_beta(), otherwise the
    // nonlocal (beta) two-center tables are left empty.
#ifdef __FFT_TWO_CENTER
    two_center_bundle.tabulate();
#else
    two_center_bundle.tabulate(this->inp_->lcao_ecut, this->inp_->lcao_dk, this->inp_->lcao_dr, this->inp_->lcao_rmax);
#endif
}

template<typename T, typename TR>
void ModuleESolver::ESolver_LR<T, TR>::parameter_check()const
{
    const std::set<std::string> lr_solvers = { "dav", "lapack" , "spectrum", "dav_subspace", "cg", "elpa", "plot" };
    // "rpa" and "bse" have no xc kernel at all; everything else is either a (semi)local
    // functional or a hybrid, both of which `LR_Util` enumerates. Listing the names a third
    // time here is what used to make a newly supported hybrid fail at input parsing.
    const std::set<std::string> kernel_less = { "rpa", "bse" };
    const std::set<std::string> abs_gauge = { "velocity", "length" };
    if (lr_solvers.find(this->inp_->lr_solver) == lr_solvers.end()) {
        throw std::invalid_argument("ESolver_LR: unknown type of lr_solver");
    }
    if (!kernel_less.count(this->xc_kernel) && !LR_Util::has_local_xc(this->xc_kernel)
        && !LR_Util::hybrid_xc_list().count(this->xc_kernel)) {
        throw std::invalid_argument("ESolver_LR: unknown type of xc_kernel");
    }
    if (this->nspin != 1 && this->nspin != 2) {
        throw std::invalid_argument("LR-TDDFT only supports nspin = 1 or 2 now");
    }
    if (abs_gauge.find(this->inp_->abs_gauge) == abs_gauge.end()) {
        throw std::invalid_argument("ESolver_LR: unknown type of abs_gauge");
    }
    if (this->inp_->cal_force && LR_Util::has_local_xc(this->xc_kernel))
    {
        std::cout << "To calculate LR-TDDFT gradients, Libxc should be compiled with kxc, i.e. `-DDISABLE_KXC=OFF` with cmake." << std::endl;
    }
}

template<typename T, typename TR>
void ModuleESolver::ESolver_LR<T, TR>::set_dimension()
{
    this->nspin = this->inp_->nspin;
    this->nstates = this->inp_->lr_nstates;
    this->nbasis = PARAM.globalv.nlocal;
    int ks_nbands = this->inp_->nbands;
    this->nocc_max = LR_Util::cal_nocc(LR_Util::cal_nelec(*this->ucell_));
    if (this->inp_->ri_hartree_benchmark == "aims" || this->inp_->ri_hartree_benchmark == "aims-librpa"
        && !this->inp_->aims_nbasis.empty())
    {
        // calculate total number of basis funcs, see https://en.cppreference.com/w/cpp/algorithm/inner_product
        this->nbasis = std::inner_product(this->inp_->aims_nbasis.begin(), /* iterator1.begin */
                                          this->inp_->aims_nbasis.end(),  /* iterator1.end */
                                          this->ucell_->atoms,  /* iterator2.begin */
                                          0,  /* init value */
                                          std::plus<int>(), /* iter op1 */
                                          [](const int& a, const Atom& b) { return a * b.na; }); /* iter op2 */
        std::cout << "nbasis from aims: " << this->nbasis << std::endl;
        for (int it = 0; it < this->ucell_->ntype; ++it)
        {
            this->ucell_->atoms[it].nw = this->inp_->aims_nbasis[it];
        }
        const_cast<UnitCell*>(this->ucell_)->set_iat2iwt(1);  // update iat2iwt for aims_nbasis 25-05-23

        int nbands_file = 0;
        int nk_file = 0;
        int nspin_file = 0;
        int nocc_file = 0;
        LR_IO::parse_band_out_file(this->inp_->rpa_outdir, nbands_file, nk_file, nspin_file, nocc_file);
        std::cout << "nocc from band_out: " << nocc_file << std::endl;
        ks_nbands = nbands_file;
        this->nocc_max = nocc_file;
    }
    // calculate the number of occupied and unoccupied states
    // which determines the basis size of the excited states    
    this->nocc_in = std::max(1, std::min(this->inp_->nocc, this->nocc_max));
    this->nvirt_in = ks_nbands - this->nocc_max;   //nbands-nocc
    if (this->inp_->nvirt > this->nvirt_in) { this->ofs_running_ << "ESolver_LR: input nvirt is too large to cover by nbands, set nvirt = nbands - nocc = " << this->nvirt_in << std::endl; }
    else if (this->inp_->nvirt > 0) { this->nvirt_in = this->inp_->nvirt; }
    this->nbands = this->nocc_in + this->nvirt_in;
    this->nk = this->inp_->nspin == 2 ? this->kv.get_nks() / 2 : this->kv.get_nks();
    this->nocc.resize(nspin, nocc_in);
    this->nvirt.resize(nspin, nvirt_in);
    if (this->nstates <= 0) { 
        this->nstates = nk * nocc_in * nvirt_in;
        this->ofs_running_ << "ESolver_LR: lr_nstates <= 0, set nstates = nk * nocc * nvirt = " << this->nstates << std::endl;
    }
    for (int is = 0;is < nspin;++is) { this->npairs.push_back(nocc[is] * nvirt[is]); }
    this->ofs_running_ << "Setting LR-TDDFT parameters: " << std::endl;
    this->ofs_running_ << "number of occupied bands: " << nocc_in << std::endl;
    this->ofs_running_ << "number of virtual bands: " << nvirt_in << std::endl;
    this->ofs_running_ << "number of Atom orbitals (LCAO-basis size): " << this->nbasis << std::endl;
    this->ofs_running_ << "number of KS bands: " << this->eig_ks.nc << std::endl;
    this->ofs_running_ << "number of excited states to be solved: " << this->nstates << std::endl;
}

template<typename T, typename TR>
void ModuleESolver::ESolver_LR<T, TR>::reset_dim_spin2()
{
	if (nspin != 2) 
	{ 
		return; 
	}
	if (nupdown != 0)
    {
        this->openshell = true;
        nupdown > 0 ? ((nocc[1] -= nupdown) && (nvirt[1] += nupdown)) : ((nocc[0] += nupdown) && (nvirt[0] -= nupdown));
        npairs = { nocc[0] * nvirt[0], nocc[1] * nvirt[1] };
        std::cout << "** Solve the spin-up and spin-down states separately for open-shell system. **" << std::endl;
    }
	for (int is : {0, 1}) 
	{ 
		if (npairs[is] <= 0) 
		{ 
			throw std::invalid_argument(std::string("ESolver_LR: npairs (nocc*nvirt) <= 0 for spin") + std::string(is == 0 ? "up" : "down")); 
		} 
	}

	if (nstates > (npairs[0] + npairs[1]) * nk) 
	{ 
		throw std::invalid_argument("ESolver_LR: nstates > nocc*nvirt*nk"); 
	}
	if (this->inp_->lr_unrestricted) 
	{ 
		this->openshell = true; 
	}
    if (!this->openshell) 
	{ 
		std::cout << " ** Assuming degenerate spin-up and spin-down states  **" << std::endl; 
	}
}

template <typename T, typename TR>
ModuleESolver::ESolver_LR<T, TR>::ESolver_LR(const Input_para& inp,
                                            const std::string& in_dir,
                                            const std::string& out_dir)
    : in_dir(in_dir), out_dir(out_dir)
{
#ifdef __EXX
    init_exx_info(this->exx_info, inp);
#endif
}

template <typename T, typename TR>
void ModuleESolver::ESolver_LR<T, TR>::before_all_runners(BaseCell& basecell, const Input_para& inp)
{
    basecell.require_kind(BaseCell::Kind::unitcell, __FUNCTION__);
    UnitCell& ucell = static_cast<UnitCell&>(basecell);
    this->ucell_ = &ucell;
    this->inp_ = &inp;
    if (inp.esolver_type == "ks-lr")
    {
#ifdef __JSON
        // The embedded KS run happens before Relax_Driver starts its first step.
        Json::init_output_array_obj();
#endif
        // the ground-state solver is a member, not a temporary: `runner` re-runs its SCF on every
        // ionic step, and the objects aliased from it have to stay alive as long as it does.
        // Its SCF is deliberately NOT run here -- `before_all_runners` is called once, outside the
        // relaxation loop, so the ground state has to be recomputed from `runner(istep)` instead.
        this->ks_ = LR_Util::make_unique<ModuleESolver::ESolver_KS_LCAO<T, TR>>();
        this->ks_->before_all_runners(basecell, inp);
    }
    else
    {
        this->initialize_from_unitcell_(ucell, inp);
    }
}

template <typename T, typename TR>
void ModuleESolver::ESolver_LR<T, TR>::bind_ground_state_aliases_()
{
    if (this->ks_)
    {
        // `ESolver_FP::before_all_runners` was never called on this object, so its own `pw_rhod`,
        // `Pgrid`, `sf` and `locpp` are empty -- the excited-state force needs all four. Point at
        // the ground-state solver's, which `before_scf` refreshes for the current geometry every
        // ionic step. `pw_rho_flag` stays false so the destructor does not free what it borrowed.
        this->pw_rho = this->ks_->pw_rho;
        this->pw_rhod = this->ks_->pw_rhod;
        this->pw_big = this->ks_->pw_big;
        this->pw_rho_flag = false;
        this->pgrid_ptr_ = &this->ks_->Pgrid;
        this->sf_ptr_ = &this->ks_->sf;
        this->locpp_ptr_ = &this->ks_->locpp;
        this->gd_ptr_ = &this->ks_->gd;
        // the two-center tables depend on the orbitals only, so the ground-state solver's single
        // build serves every geometry. This used to be moved over only for `abs_gauge velocity`,
        // leaving `LR_Force` and `cal_hs_grad` with an empty bundle in the length gauge.
        this->tcb_ptr_ = &this->ks_->two_center_bundle_;
    }
    else
    {
        this->pgrid_ptr_ = &this->Pgrid;
        this->sf_ptr_ = &this->sf;
        this->locpp_ptr_ = &this->locpp;
        this->gd_ptr_ = &this->gd_own_;
        this->tcb_ptr_ = &this->two_center_bundle_own_;
    }
}

template <typename T, typename TR>
void ModuleESolver::ESolver_LR<T, TR>::initialize_from_ks_(UnitCell& ucell, const Input_para& inp)
{
    ModuleBase::TITLE("ESolver_LR", "ESolver_LR(KS)");
    ModuleESolver::ESolver_KS_LCAO<T, TR>& ks_sol = *this->ks_;
    this->bind_ground_state_aliases_();

	if (this->inp_->lr_solver == "spectrum") 
	{
		throw std::invalid_argument("when lr_solver==spectrum, esolver_type must be `lr` to skip KS calculation.");
	}

    // xc kernel
    this->xc_kernel = LR_Util::tolower(inp.xc_kernel);
    //kv
    this->kv = ks_sol.kv;   // copy: cheap, and the KS solver keeps using its own

    this->parameter_check();

    this->set_dimension();

    // setup_2d_division is not need to be covered in #ifdef __MPI, see its implementation
    LR_Util::setup_2d_division(this->paraMat_, 1, this->nbasis, this->nbasis);
    this->set_parallel_orbitals_band(this->paraMat_, this->nbands);
    if (this->inp_->cal_force)
    {
        LR_Util::setup_2d_division(this->paraMat_all_, 1, this->nbasis, this->nbasis);
        this->set_parallel_orbitals_band(this->paraMat_all_, this->inp_->nbands);
    }

    this->paraMat_.atom_begin_row = ks_sol.pv.atom_begin_row;
    this->paraMat_.atom_begin_col = ks_sol.pv.atom_begin_col;
    this->paraMat_.iat2iwt_ = ucell.get_iat2iwt();

    LR_Util::setup_2d_division(this->paraC_, 1, this->nbasis, this->nbands
#ifdef __MPI
        , this->paraMat_.blacs_ctxt
#endif
    );

    // allocate psi_ks and eig_ks in the [nocc, nvirt] window
#ifdef __MPI
    this->psi_ks.reset(new psi::Psi<T>(this->kv.get_nks(),
        this->paraC_.get_col_size(),
        this->paraC_.get_row_size(),
        this->kv.ngk,
        true));
#else
    this->psi_ks.reset(new psi::Psi<T>(this->kv.get_nks(), this->nbands, this->nbasis, this->kv.ngk, true));
#endif
    this->eig_ks.create(this->kv.get_nks(), this->nbands);
    this->pelec = new elecstate::ElecStateLCAO<T>();
    orb_cutoff_ = ks_sol.orb_.cutoffs();

#ifdef __EXX
    // Two independent reasons to need an Exx_LRI: the LR exchange kernel, and the ground-state
    // EXX terms of the gradient (the H_gs[T+Z] and W multipliers, gated on gs_is_hybrid).
    // `initialize_from_unitcell_` has always covered both; this path used to test only the first,
    // so a local kernel on top of a hybrid ground state had no exx_lri at all when asked for forces.
    if (exx_kernel_list().count(xc_kernel) || (this->inp_->cal_force && gs_is_hybrid(this->inp_->dft_functional)))
    {
        std::string dft_functional = LR_Util::tolower(this->inp_->dft_functional);
        // Either object would be built from the same `info_ri.coulomb_param`, which `input_conv`
        // derives from dft_functional alone -- so whenever the ground-state solver has one of the
        // right type it is the same object we would construct, already up to date for this
        // geometry. Sharing it also skips a `cal_exx_ions` per ionic step.
        const bool share = (ks_sol.exx_nao.exd && std::is_same<T, double>::value)
                        || (ks_sol.exx_nao.exc && std::is_same<T, std::complex<double>>::value);
        warn_if_kernel_differs_from_gs(xc_kernel, dft_functional, this->ofs_running_);
        if (share) { this->exx_owned_ = false; }   // `refresh_from_ks_` re-binds it every step
        else    // construct C, V from scratch
        {
            // `input_conv` already filled `info_ri.coulomb_param` from INPUT.
            exx_info.sync_from_global();
            // populate ABFs/JLE file lists from UnitCell; keep in sync with Exx_NAO::init
            exx_info.info_ri.files_abfs = ucell.abfs_orbital_files;
            exx_info.info_opt_abfs.files_abfs = ucell.abfs_orbital_files;
            exx_info.info_opt_abfs.files_jles = ucell.jle_orbital_files;
            this->exx_lri = std::make_shared<Exx_LRI<T>>(exx_info.info_ri);
            this->exx_lri->init(MPI_COMM_WORLD, ucell, this->kv, ks_sol.orb_);
            this->exx_owned_ = true;   // the position-dependent `cal_exx_ions` is left to `refresh_from_ks_`
        }
    }
#endif

    refresh_from_ks_(ucell);
    this->ks_initialized_ = true;
}

template <typename T, typename TR>
void ModuleESolver::ESolver_LR<T, TR>::refresh_from_ks_(UnitCell& ucell)
{
    ModuleBase::TITLE("ESolver_LR", "refresh_from_ks_");
    ModuleESolver::ESolver_KS_LCAO<T, TR>& ks_sol = *this->ks_;
    // The ground-state solver owns `psi` and reuses it across ionic steps, so it cannot be stolen.
    // `eig_ks_all` / `wg_ks_all` are nspin x nbands matrices -- a few kB, copied rather than aliased
    // so that they survive the KS solver overwriting `pelec` on the next step.
    this->psi_ks_all_ = ks_sol.psi;
    this->eig_ks_all = ks_sol.pelec->ekb;
    this->wg_ks_all = ks_sol.pelec->wg;
    const int start_band = this->nocc_max - *std::max_element(nocc.begin(), nocc.end());

    for (int ik = 0;ik < this->kv.get_nks();++ik)
    {
        // copy the KS orbitals in the [nocc, nvirt] window
#ifdef __MPI
        Cpxgemr2d(this->nbasis, this->nbands, &(*this->psi_ks_all_)(ik, 0, 0), 1, start_band + 1, ks_sol.pv.desc_wfc,
            &(*this->psi_ks)(ik, 0, 0), 1, 1, this->paraC_.desc, this->paraC_.blacs_ctxt);
#else
        // serial: each band is `nbasis` contiguous coefficients, so the window is a plain
        // band-by-band copy (this loop used to compute the two pointers and copy nothing,
        // leaving `psi_ks` uninitialized in every non-MPI build)
        for (int ib = 0;ib < this->nbands;++ib)
        {
            const auto* start = &(*this->psi_ks_all_)(ik, start_band + ib, 0);
            auto* to = &(*this->psi_ks)(ik, ib, 0);
            std::copy(start, start + this->nbasis, to);
        }
#endif
        // copy the KS bands in the [nocc, nvirt] window
        for (int ib = 0;ib < this->nbands;++ib) { this->eig_ks(ik, ib) = this->eig_ks_all(ik, start_band + ib); }
    }

    if (nspin == 2)
    {
        const int nupdown_now = cal_nupdown_form_occ(ks_sol.pelec->wg);
        if (!this->ks_initialized_)
        {
            this->nupdown = nupdown_now;
            reset_dim_spin2();   // shifts nocc/nvirt between the spin channels: must run exactly once
        }
        else if (nupdown_now != this->nupdown)
        {
            ModuleBase::WARNING_QUIT("ESolver_LR::refresh_from_ks_",
                "the ground-state spin population changed between ionic steps, but nocc/nvirt and"
                " the distributions built from them were fixed at the first step.");
        }
    }
    // `gint_info_` stays owned by the KS solver: its `before_scf` rebuilds and re-publishes it
    //init potential and calculate kernels using ground state charge
    init_pot(*ks_sol.pelec->charge);

#ifdef __EXX
    if (exx_kernel_list().count(xc_kernel) || (this->inp_->cal_force && gs_is_hybrid(this->inp_->dft_functional)))
    {
        if (this->exx_owned_)
        {   // Cs/Vs follow the atoms, so they are rebuilt for every geometry
            this->exx_lri->cal_exx_ions(ucell, this->inp_->out_ri_cv);
        }
        else if (ks_sol.exx_nao.exd) { this->share_exx_lri(ks_sol.exx_nao.exd->exx_ptr); }
        else if (ks_sol.exx_nao.exc) { this->share_exx_lri(ks_sol.exx_nao.exc->exx_ptr); }
    }
#endif
    // the grid-integration tables hang off a static pointer that the ground-state solver
    // re-publishes in its `before_scf`; make sure it names the object we integrate on
    ModuleGint::Gint::set_gint_info(this->ks_->gint_info_.get());

    // the Z-vector window: after `reset_dim_spin2`, so nocc/nvirt/openshell are final
#ifdef __MPI
    this->fill_z_window_(ks_sol.pv.desc_wfc);
#else
    this->fill_z_window_(nullptr);
#endif
}


template <typename T, typename TR>
void ModuleESolver::ESolver_LR<T, TR>::initialize_from_unitcell_(UnitCell& ucell, const Input_para& inp)
{
    ModuleBase::TITLE("ESolver_LR", "ESolver_LR(from scratch)");
    this->bind_ground_state_aliases_();
    // xc kernel
    this->xc_kernel = LR_Util::tolower(inp.xc_kernel);

    // necessary steps in ESolver_FP
    ESolver_FP::before_all_runners(ucell, inp);
    this->pelec = new elecstate::ElecStateLCAO<T>();

    // necessary steps in ESolver_KS::before_all_runners : symmetry and k-points
    if (ModuleSymmetry::Symmetry::symm_flag == 1)
    {
        const int cal_symm_repr[2] = {this->inp_->cal_symm_repr[0], this->inp_->cal_symm_repr[1]};
        ucell.symm.analy_sys(ucell.lat, ucell.st, ucell.atoms, this->ofs_running_,
                             this->inp_->symmetry_prec, this->inp_->nspin, this->inp_->calculation, cal_symm_repr);
        ModuleBase::GlobalFunc::DONE(this->ofs_running_, "SYMMETRY");
    }
    const bool use_ibz = false;
    const bool gamma_only_local = PARAM.globalv.gamma_only_local;
    const double kspacing[3] = {this->inp_->kspacing[0], this->inp_->kspacing[1], this->inp_->kspacing[2]};
    const double koffset[3] = {this->inp_->koffset[0], this->inp_->koffset[1], this->inp_->koffset[2]};
    this->kv.set(ucell, ucell.symm, this->inp_->kpoint_file, this->inp_->nspin, ucell.G, ucell.latvec, this->ofs_running_, this->ofs_warning_, use_ibz, this->out_dir, gamma_only_local, kspacing, this->inp_->kmesh_type, koffset);
    ModuleBase::GlobalFunc::DONE(this->ofs_running_, "INIT K-POINTS");
    ModuleIO::print_parameters(ucell, this->kv, inp);

    this->parameter_check();

    /// read orbitals and build the interpolation table
    two_center_bundle_own_.build_orb(ucell.ntype, ucell.orbital_fn.data(), inp.orbital_dir);

    LCAO_Orbitals orb;
    two_center_bundle_own_.to_LCAO_Orbitals(orb, inp.lcao_ecut, inp.lcao_dk, inp.lcao_dr, inp.lcao_rmax,
                                        inp.out_element_info, inp.cal_force);
    orb_cutoff_ = orb.cutoffs();
    if (LR_Util::tolower(this->inp_->abs_gauge) == "velocity")
    {
        setup_2center_table(this->two_center_bundle_own_, orb, ucell);
    }

    this->set_dimension();
    //  setup 2d-block distribution for AO-matrix and KS wfc
    LR_Util::setup_2d_division(this->paraMat_, 1, this->nbasis, this->nbasis);
    this->set_parallel_orbitals_band(this->paraMat_, this->nbands);
    if (this->inp_->cal_force)
    {
        LR_Util::setup_2d_division(this->paraMat_all_, 1, this->nbasis, this->nbasis);
        this->set_parallel_orbitals_band(this->paraMat_all_, this->inp_->nbands);
    }

    // read the ground state info
    // now ModuleIO::read_wfc_nao needs `Parallel_Orbitals` and can only read all the bands
    // it need improvement to read only the bands needed
    this->psi_ks.reset(new psi::Psi<T>(this->kv.get_nks(),
                                       this->paraMat_.ncol_bands,
                                       this->paraMat_.get_row_size(),
                                       this->kv.ngk,
                                       true));
    this->read_ks_wfc();


    if (nspin == 2)
    {   // `read_ks_wfc` fills `wg_ks`, not `pelec->wg` -- reading the latter here meant nupdown was
        // always 0, so a spin-polarised ground state silently took the closed-shell branch
        this->nupdown = cal_nupdown_form_occ(this->wg_ks);
        reset_dim_spin2();
    }
    // the Z-vector window: after `reset_dim_spin2`, so nocc/nvirt/openshell are final
#ifdef __MPI
    this->fill_z_window_(paraMat_all_.desc_wfc);
#else
    this->fill_z_window_(nullptr);
#endif

    LR_Util::setup_2d_division(this->paraC_, 1, this->nbasis, this->nbands
#ifdef __MPI
        , paraMat_.blacs_ctxt
#endif
    );

    // clear ks info, new elecstate for excition
    delete this->pelec;   // the ElecStateLCAO allocated above, only needed while reading the KS data
    this->pelec = new elecstate::ElecState();

    // read the ground state charge density and calculate xc kernel
    pgrid().init(this->pw_rho->nx,
        this->pw_rho->ny,
        this->pw_rho->nz,
        this->pw_rho->nplane,
        this->pw_rho->nrxx,
        pw_big->nbz,
        pw_big->bz,
        GlobalV::NPROC);
    Charge chg_gs;
    if (this->inp_->ri_hartree_benchmark == "none") { this->read_ks_chg(chg_gs); }
    this->init_pot(chg_gs);

    // search adjacent atoms and init Gint
    double search_radius = -1.0;
    search_radius = atom_arrange::set_sr_NL(this->ofs_running_,
        this->inp_->out_level,
        orb.get_rcutmax_Phi(),
        ucell.infoNL->get_rcutmax_Beta(),
        PARAM.globalv.gamma_only_local);
    atom_arrange::search(PARAM.globalv.search_pbc,
                         this->ofs_running_,
                         this->gd(),
                         *this->ucell_,
                         search_radius,
                         this->inp_->test_atom_input);
    gint_info_.reset(
        new ModuleGint::GintInfo(
        this->pw_big->nbx,
        this->pw_big->nby,
        this->pw_big->nbz,
        this->pw_rho->nx,
        this->pw_rho->ny,
        this->pw_rho->nz,
        0,
        0,
        this->pw_big->nbzp_start,
        this->pw_big->nbx,
        this->pw_big->nby,
        this->pw_big->nbzp,
        orb.Phi,
        ucell,
        this->gd(),
        this->inp_->nspin,
        PARAM.globalv.gamma_only_local,
        PARAM.globalv.domag,
        this->inp_->device == "gpu",
        this->inp_->nstream));
    ModuleGint::Gint::set_gint_info(gint_info_.get());
    // if EXX from scratch, init 2-center integral and calculate Cs, Vs 
    // when: 
    // 1. EXX xc_kernel
    // 2. cal_force with ground state with EXX functional
#ifdef __EXX
    if (((exx_kernel_list().count(xc_kernel)) && this->inp_->lr_solver != "spectrum")
        || (this->inp_->cal_force && gs_is_hybrid(this->inp_->dft_functional)))
    {
        warn_if_kernel_differs_from_gs(xc_kernel, LR_Util::tolower(this->inp_->dft_functional), this->ofs_running_);
        // `input_conv` already filled `info_ri.coulomb_param` from INPUT.
        exx_info.sync_from_global();
        // populate ABFs/JLE file lists from UnitCell; keep in sync with Exx_NAO::init
        exx_info.info_ri.files_abfs = ucell.abfs_orbital_files;
        exx_info.info_opt_abfs.files_abfs = ucell.abfs_orbital_files;
        exx_info.info_opt_abfs.files_jles = ucell.jle_orbital_files;
        this->exx_lri = std::make_shared<Exx_LRI<T>>(exx_info.info_ri);
        this->exx_lri->init(MPI_COMM_WORLD, ucell,this->kv, orb);
        this->exx_lri->cal_exx_ions(ucell,this->inp_->out_ri_cv);
    }
    // else
#endif
        // ModuleBase::Ylm::set_coefficients() is deprecated
}

template <typename T, typename TR>
void ModuleESolver::ESolver_LR<T, TR>::runner(BaseCell& basecell, const int istep)
{
    basecell.require_kind(BaseCell::Kind::unitcell, __FUNCTION__);
    UnitCell& ucell = static_cast<UnitCell&>(basecell);

    ModuleBase::TITLE("ESolver_LR", "runner");
    ModuleBase::timer::start("ESolver_LR", "runner");

    if (this->ks_)
    {
        // `init_pot_groundstate` leaves the global XC type on the LR kernel, so put it back before
        // the ground-state SCF. `ESolver_KS::runner` is re-entrant: its `before_scf` rebuilds the
        // neighbour lists, the grid tables and the Hamiltonian for the geometry of this step.
        XC_Functional::set_xc_type(ucell.atoms[0].ncpp.xc_func);
        this->ks_->runner(ucell, istep);
        this->etot_gs_ = this->ks_->cal_energy();
        if (!this->ks_initialized_)
        {
            this->initialize_from_ks_(ucell, *this->inp_);
            this->setup_relax_target_();   // needs the final dimensions, i.e. `openshell`
        }
        else { this->refresh_from_ks_(ucell); }
    }

    //allocate 2-particle state and setup 2d division
    this->setup_eigenvectors_X();
    this->pelec->ekb.create(nspin, this->nstates);

    // set once here rather than inside the solver branches: the open-shell `spectrum` branch used
    // to leave it empty, and both `after_all_runners` and the gradient index it
    this->spin_types = this->openshell ? std::vector<std::string>({ "updown" })
                                       : std::vector<std::string>({ "singlet", "triplet" });

    // a relaxation writes these once per ionic step; without the suffix every step would overwrite
    // the last, and the trajectory would be impossible to inspect afterwards
    const std::string step_suffix = this->excited_relax_ ? "_step" + std::to_string(istep) : "";
    auto efile_out = [&](const std::string& label)->std::string {return this->out_dir + "Excitation_Energy_" + label + step_suffix + ".dat";};
    auto vfile_out = [&](const std::string& label)->std::string {return this->out_dir + "Excitation_Amplitude_" + label + step_suffix + "_" + std::to_string(this->my_rank_+1) + ".dat";};
    if (this->inp_->lr_solver == "elpa")
    {
        ModuleBase::WARNING_QUIT("ESolver_LR", "ESolver_LR doesn't support elpa now.");

    }
    else if (this->inp_->lr_solver != "spectrum")
    {
        auto write_states = [&](const std::string& label, const Real<T>* e, const T* v, const int& dim, const int& nst, const int& prec = 8)->void
            {
                if (this->my_rank_ == 0) { assert(nst == LR_Util::write_value(efile_out(label), prec, e, nst)); }
                assert(nst * dim == LR_Util::write_value(vfile_out(label), prec, v, nst, dim));
            };
        std::vector<double> precondition(this->inp_->lr_solver == "lapack" ? 0 : nloc_per_state, 1.0);
        // allocate and initialize A matrix and density matrix
        if (openshell)
        {
            for (int is : {0, 1})
            {
                if (this->inp_->lr_solver != "lapack") {
                    const int offset_is = is * this->paraX_[0].get_local_size();
                    OperatorLRDiag<double> pre_op(this->eig_ks.c + is * nk * (nocc[0] + nvirt[0]), this->paraX_[is], this->nk, this->nocc[is], this->nvirt[is]);
                    pre_op.act(1, offset_is, 1, precondition.data() + offset_is, precondition.data() + offset_is);
                }
            }
            std::cout << "Solving spin-conserving excitation for open-shell system." << std::endl;
            HamiltULR<T> hulr(xc_kernel,
                              nspin,
                              this->nbasis,
                              this->nocc,
                              this->nvirt,
                              *this->ucell_,
                              orb_cutoff_,
                              this->gd(),
                              *this->psi_ks,
                              this->eig_ks,
#ifdef __EXX
                              this->exx_lri,
                              this->exx_info.info_global.hybrid_alpha,
#endif
                              this->pot,
                              this->kv,
                              this->paraX_,
                              this->paraC_,
                              this->paraMat_);
            LR::HSolver::solve(hulr, this->X[0].template data<T>(), nloc_per_state, nstates,
                               this->nk, this->nocc, this->nvirt, this->paraX_,
                               this->pelec->ekb.c, this->inp_->lr_solver,
                               this->inp_->lr_thr, precondition);
            if (this->inp_->out_wfc_lr) { write_states("openshell", this->pelec->ekb.c, this->X[0].template data<T>(), nloc_per_state, nstates); }
        }
        else
        {
            if (this->inp_->lr_solver != "lapack") {
                OperatorLRDiag<double> pre_op(this->eig_ks.c, this->paraX_[0], this->nk, this->nocc[0], this->nvirt[0]);
                pre_op.act(1, nloc_per_state, 1, precondition.data(), precondition.data()); 
            }
            const std::vector<std::string>& spin_types = this->spin_types;
            for (int is = 0;is < nspin;++is)
            {
                std::cout << " Calculating " << spin_types[is] << " excitations" << std::endl;
                HamiltLR<T> hlr(xc_kernel,
                                nspin,
                                this->nbasis,
                                this->nocc,
                                this->nvirt,
                                *this->ucell_,
                                orb_cutoff_,
                                this->gd(),
                                *this->psi_ks,
                                this->eig_ks,
#ifdef __EXX
                                this->exx_lri,
                                this->exx_info.info_global.hybrid_alpha,
#endif
                                this->pot[is],
                                this->kv,
                                this->paraX_,
                                this->paraC_,
                                this->paraMat_,
                                spin_types[is],
                                this->in_dir,
                                this->out_dir,
                                this->inp_->ri_hartree_benchmark,
                                (this->inp_->ri_hartree_benchmark == "aims" ? this->inp_->aims_nbasis : std::vector<int>({})));
                LR::HSolver::solve(hlr, this->X[is].template data<T>(), nloc_per_state, nstates,
                                   this->nk, this->nocc, this->nvirt, this->paraX_,
                                   this->pelec->ekb.c + is * nstates,
                                   this->inp_->lr_solver,
                                   this->inp_->lr_thr,
                                   precondition);
                if (this->inp_->out_wfc_lr) { write_states(spin_types[is], this->pelec->ekb.c + is * nstates, this->X[is].template data<T>(), nloc_per_state, nstates); }
            }
        }
    }
    else    // lr_solver == "spectrum", read the eigenvalues
    {
        auto efile_in = [&](const std::string& label)->std::string {return this->in_dir + "Excitation_Energy_" + label + ".dat";};
        auto vfile_in = [&](const std::string& label)->std::string {return this->in_dir + "Excitation_Amplitude_" + label + "_" + std::to_string(this->my_rank_+1) + ".dat";};
    
        auto read_states = [&](const std::string& label, Real<T>* e, T* v, const int& dim, const int& nst)->void
            {
                if (this->my_rank_ == 0) {
                    assert(nst == LR_Util::read_value(efile_in(label), e, nst));
                    std::cout <<"Rank "<< this->my_rank_ << ": finish reading " << efile_in(label) << std::endl;
                }
#ifdef __MPI
// in velocity gauge, the eigenvalues may be used to calculate the transition dipole, so we'd better broadcast them
                MPI_Bcast(e, nst, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif
                assert(nst * dim == LR_Util::read_value(vfile_in(label), v, nst, dim));
                std::cout <<"Rank "<< this->my_rank_ << ": finish reading " << vfile_in(label) << std::endl;
            };
        std::cout << "reading the excitation states from file: \n";
        if (openshell)
        {
            read_states("openshell", this->pelec->ekb.c, this->X[0].template data<T>(), nloc_per_state, nstates);
        }
        else
        {
            const std::vector<std::string>& spin_types = this->spin_types;
            for (int is = 0;is < nspin;++is) { read_states(spin_types[is], this->pelec->ekb.c + is * nstates, this->X[is].template data<T>(), nloc_per_state, nstates); }
        }
    }
    if (this->excited_relax_)
    {
        // Re-select which root to follow BEFORE taking the gradient, so the force belongs to
        // the same diabatic state as the previous step's.
        this->follow_target_state_(this->ofs_running_);

        // The LR terms only carry the Omega part of the force; the ground-state part is separate
        // and comes straight from the KS solver.
        this->ks_->cal_force(ucell, this->force_gs_);
        this->lr_force_ = this->cal_lr_force_relax_(this->ofs_running_);

        // One line per ionic step with the two halves of the energy and of the gradient. Without
        // it the relaxation only reports a force, and whether E_gs + Omega actually goes down --
        // the thing being minimised -- cannot be read off the log at all.
        const double omega = this->target_omega_();
        auto max_abs = [](const ModuleBase::matrix& m) -> double
            { double v = 0.0; for (int i = 0;i < m.nr * m.nc;++i) { v = std::max(v, std::abs(m.c[i])); } return v; };
        this->ofs_running_ << std::setprecision(8) << std::fixed
            << " EXCITED-STATE RELAX step " << istep
            << ": E_gs = " << this->etot_gs_ * ModuleBase::Ry_to_eV
            << " eV, Omega = " << omega * ModuleBase::Ry_to_eV
            << " eV, E_exc = " << (this->etot_gs_ + omega) * ModuleBase::Ry_to_eV << " eV"
            << std::setprecision(6)
            << " | |F_gs|max = " << max_abs(this->force_gs_) * ModuleBase::Ry_to_eV / ModuleBase::BOHR_TO_A
            << ", |F_Omega|max = " << max_abs(this->lr_force_) * ModuleBase::Ry_to_eV / ModuleBase::BOHR_TO_A
            << " eV/Angstrom" << std::defaultfloat << std::endl;
    }

    ModuleBase::timer::end("ESolver_LR", "runner");
    return;
}

template <typename T, typename TR>
void ModuleESolver::ESolver_LR<T, TR>::after_all_runners(BaseCell& basecell)
{
    basecell.require_kind(BaseCell::Kind::unitcell, __FUNCTION__);
    UnitCell& ucell = static_cast<UnitCell&>(basecell);

    ModuleBase::TITLE("ESolver_LR", "after_all_runners");
    if (this->inp_->ri_hartree_benchmark != "none") { return; } //no need to calculate the spectrum in the benchmark routine

    // cal electron-hole density
    if (this->inp_->out_chg[0])
    {
        LR_Density<T> lr_density(*this->ucell_, kv, gd(), *psi_ks, orb_cutoff_, pgrid(),
            nspin, nocc, nvirt, nbasis,
            paraX_, paraC_, paraMat_, openshell);

        if (openshell)
            for (int is = 0;is < this->nspin;++is)
                lr_density.output_eh_density_all_states(this->X[0].template data<T>(), is, nstates);
        else
            for (int is = 0;is < this->X.size();++is)
                lr_density.output_eh_density_all_states(this->X[is].template data<T>(), is, nstates);
    }

    //cal spectrum
    if (LR_Util::tolower(this->inp_->abs_gauge) == "velocity" )
    {
        const int nspin_tmp = this->inp_->nspin == 2 ? 2 : 1;
        this->velocity_mo = LR_Util::cal_velocity_mo(*this->ucell_, this->gd(), this->tcb(), 
                                                    this->paraMat_, this->paraC_, this->kv, *this->psi_ks, 
                                                    this->nk, nspin_tmp, this->nbasis, this->nocc, this->nvirt);
    }
    std::vector<double> freq(100);
    std::vector<double> abs_wavelen_range({ 20, 200 });//default range
    if (this->inp_->abs_wavelen_range.size() >= 2 && std::abs(this->inp_->abs_wavelen_range[1] - this->inp_->abs_wavelen_range[0]) > 0.02)
    {
        abs_wavelen_range = this->inp_->abs_wavelen_range;
    }
    double lambda_diff = std::abs(abs_wavelen_range[1] - abs_wavelen_range[0]);
    double lambda_min = std::min(abs_wavelen_range[1], abs_wavelen_range[0]);
    for (int i = 0;i < freq.size();++i) { freq[i] = 91.126664 / (lambda_min + 0.01 * static_cast<double>(i + 1) * lambda_diff); }
    // auto spin_types = (nspin == 2 && !openshell) ? std::vector<std::string>({ "singlet", "triplet" }) : std::vector<std::string>({ "updown" });
    // for (int is = 0;is < this->X.size() - 1;++is)
    for (int is = 0;is < this->X.size();++is)
    {
        LR_Spectrum<T> spectrum(nspin, this->nbasis, this->nocc, this->nvirt, *this->pw_rho, *this->psi_ks,
            *this->ucell_, this->kv, this->gd(), this->orb_cutoff_, this->tcb(),
            this->paraX_, this->paraC_, this->paraMat_,
            &this->pelec->ekb.c[is * nstates], this->eig_ks.c, this->X[is].template data<T>(), nstates, openshell,
            LR_Util::tolower(this->inp_->abs_gauge), this->my_rank_, this->out_dir);
        if (LR_Util::tolower(this->inp_->abs_gauge) == "velocity" ) {spectrum.set_vmo(this->velocity_mo.data());}
        spectrum.cal_spectrum();
        spectrum.transition_analysis(spin_types[is]+"_tda");
        if (spin_types[is] != "triplet")        // triplets has no transition dipole and no contribution to the spectrum
        {
            spectrum.optical_absorption_method1(freq, this->inp_->abs_broadening);
            spectrum.write_transition_dipole(this->out_dir +
                "trans_dipole_" + spin_types[is] + "_tda.dat");
            // =============================================== for test ====================================================
            // spectrum.optical_absorption_method2(freq, this->inp_->abs_broadening);
            // if (LR_Util::tolower(this->inp_->abs_gauge) == "velocity")
            // {   // TEST the formula v/omega rather than v/(e_a-e_i)
            //     spectrum.test_transition_dipoles_velocity_omega();
            //     spectrum.write_transition_dipole(this->out_dir +
            //         "trans_dipole_" + spin_types[is] + "_vomega_tda.dat");
            // }
            // =============================================== for test ====================================================
        }
        if (this->inp_->cal_force && !this->excited_relax_) { this->cal_force_and_grad_matrix_(is, this->ofs_running_); }
    }
}
template<typename T, typename TR>
void ModuleESolver::ESolver_LR<T, TR>::set_parallel_orbitals_band(Parallel_Orbitals& pmat, const int nbands_in)
{
#ifdef __MPI
    pmat.set_desc_wfc_Eij(this->nbasis, nbands_in, pmat.get_row_size());
    int err = pmat.set_nloc_wfc_Eij(nbands_in, this->ofs_running_, this->ofs_warning_);
    // Skipped for the aims benchmark: with `aims_nbasis` the per-atom orbital counts behind
    // `iat2iwt` do not match `nbasis`, so the atomic trace would be wrong. The guard came from
    // d4fe3fe84 ("Support different basis number from aims"), and was silently undone by
    // 1e4c1c6af (BSE, #7718) re-adding an unconditional call above it.
    if (this->inp_->ri_hartree_benchmark != "aims")
    {
        pmat.set_atomic_trace(this->ucell_->get_iat2iwt(), this->ucell_->nat, this->nbasis);
    }
#else
    pmat.nrow_bands = this->nbasis;
    pmat.ncol_bands = nbands_in;
#endif
}

template<typename T, typename TR>
void ModuleESolver::ESolver_LR<T, TR>::setup_eigenvectors_X()
{
    ModuleBase::TITLE("ESolver_LR", "setup_eigenvectors_X");
    // this function is called once per `runner`, and `paraX_` is only ever appended to,
    // so without this reset a second ionic step would double its size
    this->paraX_.clear();
    for (int is = 0;is < nspin;++is)
    {
        Parallel_2D px;
        LR_Util::setup_2d_division(px, /*nb2d=*/1, this->nvirt[is], this->nocc[is]
#ifdef __MPI
            , this->paraC_.blacs_ctxt
#endif
        );//nvirt - row, nocc - col 
        this->paraX_.emplace_back(std::move(px));
    }
    this->nloc_per_state = nk * (openshell ? paraX_[0].get_local_size() + paraX_[1].get_local_size() : paraX_[0].get_local_size());

    this->X.resize(openshell ? 1 : nspin, LR_Util::newTensor<T>({ nstates, nloc_per_state }));
    for (auto& x : X) { x.zero(); }

    // if spectrum-only, read the LR-eigenstates from file and return
    if (this->inp_->lr_solver != "spectrum") { set_X_initial_guess(); }
}

template<typename T, typename TR>
void ModuleESolver::ESolver_LR<T, TR>::set_X_initial_guess()
{
    // set the initial guess of X
    for (int is = 0;is < this->nspin;++is)
    {
        const int& no = this->nocc[is];
        const int& nv = this->nvirt[is];
        const int& np = this->npairs[is];
        const Parallel_2D& px = this->paraX_[is];

        // if (E_{lumo}-E_{homo-1} < E_{lumo+1}-E{homo}), mode = 0, else 1(smaller first)
        bool ix_mode = false;   //default
        if (this->eig_ks.nc > no + 1 && no >= 2 && eig_ks(is, no) - eig_ks(is, no - 2) - 1e-5 > eig_ks(is, no + 1) - eig_ks(is, no - 1)) { ix_mode = true; }
        this->ofs_running_ << "setting the initial guess of X of spin" << is << std::endl;
        if (no >= 2 && eig_ks.nc > no) { this->ofs_running_ << "E_{lumo}-E_{homo-1}=" << eig_ks(is, no) - eig_ks(is, no - 2) << std::endl; }
        if (no >= 1 && eig_ks.nc > no + 1) { this->ofs_running_ << "E_{lumo+1}-E{homo}=" << eig_ks(is, no + 1) - eig_ks(is, no - 1) << std::endl; }
        this->ofs_running_ << "mode of X-index: " << ix_mode << std::endl;

        /// global index map between (i,c) and ix
        ModuleBase::matrix ioiv2ix;
        std::vector<std::pair<int, int>> ix2ioiv;
        std::pair<ModuleBase::matrix, std::vector<std::pair<int, int>>> indexmap =
            LR_Util::set_ix_map_diagonal(ix_mode, no, nv);

        ioiv2ix = std::move(std::get<0>(indexmap));
        ix2ioiv = std::move(std::get<1>(indexmap));

        for (int ib = 0; ib < nstates; ++ib)
        {
            const int ipair = ib % np;
            const int occ_global = std::get<0>(ix2ioiv[ipair]);   // occ
            const int virt_global = std::get<1>(ix2ioiv[ipair]);   // virt
            const int ik = ib / np;
            const int xstart_b = ib * nloc_per_state;    //start index of band ib
            const int xstart_bs = (openshell && is == 1) ? xstart_b + nk * paraX_[0].get_local_size() : xstart_b;  // start index of band ib, spin is
            const int is_in_x = openshell ? 0 : is;     // if openshell, spin-up and spin-down are put together
            if (px.in_this_processor(virt_global, occ_global))
            {
                const int xstart_pair = ik * px.get_local_size();
                const int ipair_loc = px.global2local_col(occ_global) * px.get_row_size() + px.global2local_row(virt_global);
                X[is_in_x].data<T>()[xstart_bs + xstart_pair + ipair_loc] = (static_cast<T>(1.0) / static_cast<T>(nk));
            }
        }
    }
}

template<typename T, typename TR>
void ModuleESolver::ESolver_LR<T, TR>::init_pot(const Charge& chg_gs)
{
    using ST = PotHxcLR::SpinType;
    using GX = LR::KernelXC::GxcSpin;
    this->pot.assign(nspin, nullptr);   // assign, not resize: a re-init must drop the previous geometry's
    if (this->inp_->ri_hartree_benchmark != "none") { return; } //no need to initialize potential for Hxc kernel in the RI-benchmark routine

    // The singlet and triplet potentials evaluate the *same* kernel arrays and differ only in which
    // spin combination of them they read, so they share one `KernelXC` to save memory.
    const bool oshell = (nspin == 2) && openshell;
    // $g^{xc}$ (third-order) is only ever needed by the LR gradient, and only for the spin
    // combinations that are actually going to be requested.
    const int gxc_lr = (!this->inp_->cal_force || !LR_Util::has_local_xc(xc_kernel)) ? GX::NoGxc
        : ((nspin == 1) ? GX::Singlet : GX::BothSpins);
    std::shared_ptr<const LR::KernelXC> kernel_lr = PotHxcLR::make_kernel(
        xc_kernel, *this->pw_rho, *this->ucell_, chg_gs, pgrid(), oshell, gxc_lr, this->inp_->lr_init_xc_kernel);
    switch (nspin)
    {
    case 1:
        this->pot[0] = std::make_shared<PotHxcLR>(kernel_lr, xc_kernel, *this->pw_rho, *this->ucell_, chg_gs.nrxx, ST::S1);
        break;
    case 2:
        this->pot[0] = std::make_shared<PotHxcLR>(kernel_lr, xc_kernel, *this->pw_rho, *this->ucell_, chg_gs.nrxx, oshell ? ST::S2_updown : ST::S2_singlet);
        this->pot[1] = std::make_shared<PotHxcLR>(kernel_lr, xc_kernel, *this->pw_rho, *this->ucell_, chg_gs.nrxx, oshell ? ST::S2_updown : ST::S2_triplet);
        break;
    default:
        throw std::invalid_argument("ESolver_LR: nspin must be 1 or 2");
    }
    // ground-state potentials are needed for calculating the excited state force
    if (this->inp_->cal_force)
    {
        this->init_pot_groundstate(chg_gs);
        // `dft_functional == "default"` leaves the raw INPUT string unresolved (it never gets
        // overwritten to the actual functional in use); the functional actually read from the
        // pseudopotential lives in `ucell.atoms[i].ncpp.xc_func` instead. Comparing against the
        // literal "default" string here would always disagree with `xc_kernel`, forcing a
        // separate `kernel_gs` with no g^xc even when the ground state and the LR kernel are the
        // same functional -- and `pot_hxc_gs` (built from that `kernel_gs`) throws the first time
        // `cal_W_from_Z` asks it for g^xc.
        const std::string xc_kernel_gs = (this->inp_->dft_functional == "default")
            ? LR_Util::tolower(this->ucell_->atoms[0].ncpp.xc_func)
            : LR_Util::tolower(this->inp_->dft_functional);
        // `ST::S1` is only correct when nspin=1. `PotHxcLR` builds its `KernelXC` with
        // the input `nspin`, so at nspin=2 the kernel arrays carry 3 spin components per grid point
        // while the S1 integrand indexes them as if there were 1 -- it does not even read a
        // consistent spin combination. Use `ST::S2_gs` there, which is exactly half of S2_singlet,
        // matching the `K_Hxc(singlet) = 2 * pot_hxc_gs` convention of the gradient operators.
        const ST st_gs = (nspin == 1) ? ST::S1_gs : (oshell ? ST::S2_updown : ST::S2_gs);
        // `pot_hxc_gs` supplies the $g^{xc}$ of both $W^c$ and the force term, for either spin.
        // Those two call sites are guarded by `has_local_xc(xc_kernel)` -- the *LR* kernel name --
        // so an `xc_kernel rpa` run never touches them however local `dft_functional` is.
        const int gxc_gs = (!LR_Util::has_local_xc(xc_kernel) || !LR_Util::has_local_xc(xc_kernel_gs)) ? GX::NoGxc
            : ((nspin == 1) ? GX::Singlet : GX::BothSpins);
        // When the LR kernel *is* the ground-state functional -- the usual TDDFT case -- the two
        // `KernelXC` are bit-for-bit identical, so reuse the one just built.
        const bool share_lr = (xc_kernel_gs == xc_kernel) && ((gxc_lr & gxc_gs) == gxc_gs);
        std::shared_ptr<const LR::KernelXC> kernel_gs = share_lr ? kernel_lr
            : PotHxcLR::make_kernel(xc_kernel_gs, *this->pw_rho, *this->ucell_, chg_gs, pgrid(), oshell, gxc_gs, this->inp_->lr_init_xc_kernel);
        this->pot_hxc_gs = std::make_shared<LR::PotHxcLR>(kernel_gs, xc_kernel_gs, *this->pw_rho, *this->ucell_, chg_gs.nrxx, st_gs);
    }
}

template<typename T, typename TR>
void ModuleESolver::ESolver_LR<T, TR>::read_ks_wfc()
{
    assert(this->psi_ks != nullptr);
    this->eig_ks.create(this->kv.get_nks(), this->nbands);
    this->wg_ks.create(this->kv.get_nks(), this->nbands);
    if (this->inp_->ri_hartree_benchmark == "aims")        // for aims benchmark
    {
#ifdef __EXX
        int ncore = 0;
        std::vector<double> eig_ks_vec = RI_Benchmark::read_aims_ebands<double>(this->in_dir + "band_out", nocc_in, nvirt_in, ncore);
        std::cout << "ncore=" << ncore << ", nocc=" << nocc_in << ", nvirt=" << nvirt_in << ", nbands=" << this->nbands << std::endl;
        std::cout << "eig_ks_vec.size()=" << eig_ks_vec.size() << std::endl;
        if(eig_ks_vec.size() != this->nbands) {ModuleBase::WARNING_QUIT("ESolver_LR", "read_aims_ebands failed.");};
        for (int i = 0;i < nbands;++i) { this->eig_ks(0, i) = eig_ks_vec[i]; }
        RI_Benchmark::read_aims_eigenvectors<T>(*this->psi_ks, this->in_dir + "KS_eigenvectors.out", ncore, nbands, nbasis);
#else
        ModuleBase::WARNING_QUIT("ESolver_LR", "RI benchmark is only supported when compile with LibRI.");
#endif
    }
    else if (!ModuleIO::read_wfc_nao(this->in_dir, this->paraMat_, *this->psi_ks,
                this->eig_ks,
                this->wg_ks,
                this->kv.ik2iktot,
                this->kv.get_nkstot(),
                this->inp_->nspin,
                this->inp_->init_wfc_file_format == "binary",
                /*skip_bands=*/this->nocc_max - this->nocc_in)) {
        ModuleBase::WARNING_QUIT("ESolver_LR", "read ground-state wavefunction failed.");
    }

    if (this->inp_->cal_force)
    {    // allocate psi_ks_all and eig_ks_all to read all the bands
        this->psi_ks_all_own_.reset(new psi::Psi<T>(this->kv.get_nks(), paraMat_all_.ncol_bands, paraMat_all_.get_row_size(), this->kv.ngk, true));
        this->psi_ks_all_ = this->psi_ks_all_own_.get();
        this->eig_ks_all.create(this->kv.get_nks(), this->inp_->nbands);
        this->wg_ks_all.create(this->kv.get_nks(), this->inp_->nbands);
        if (!ModuleIO::read_wfc_nao(this->in_dir, paraMat_all_, *this->psi_ks_all_,
                this->eig_ks_all,
                this->wg_ks_all,
                this->kv.ik2iktot,
                this->kv.get_nkstot(),
                this->inp_->nspin,
                this->inp_->init_wfc_file_format == "binary",
                /*skip_bands=*/0))
        {
            this->ofs_running_ << " Read in all the KS wavefunctions for force calculation. " << std::endl;
        }
    }
}

template<typename T, typename TR>
void ModuleESolver::ESolver_LR<T, TR>::fill_z_window_(const int* desc_src)
{
    ModuleBase::TITLE("ESolver_LR", "fill_z_window_");
    if (!this->inp_->cal_force || this->psi_ks_all_ == nullptr) { return; }

    const int start_band = this->nocc_max - *std::max_element(nocc.begin(), nocc.end());
    // Every band the ground state solved, from the window start upward. `eig_ks_all` is the
    // authority on how many there are: on the ks-lr path it is the KS solver's `ekb`, on the
    // file path it was read with `skip_bands = 0`.
    this->nbands_z_ = this->eig_ks_all.nc - start_band;
    if (this->nbands_z_ <= this->nbands)
    {   // nothing to widen: nbands was already at (or below) the X window
        this->nbands_z_ = this->nbands;
    }
    this->nvirt_z_.assign(this->nspin, 0);
    for (int is = 0; is < this->nspin; ++is) { this->nvirt_z_[is] = this->nbands_z_ - this->nocc[is]; }

    LR_Util::setup_2d_division(this->paraC_z_, 1, this->nbasis, this->nbands_z_
#ifdef __MPI
        , this->paraMat_.blacs_ctxt
#endif
    );
    this->paraX_z_.clear();
    for (int is = 0; is < this->nspin; ++is)
    {
        Parallel_2D px;
        LR_Util::setup_2d_division(px, /*nb2d=*/1, this->nvirt_z_[is], this->nocc[is]
#ifdef __MPI
            , this->paraC_z_.blacs_ctxt
#endif
        );
        this->paraX_z_.emplace_back(std::move(px));
    }
    // Open shell solves one eigenproblem whose vector is the concatenation [up | down], so a
    // state's block is as long as both channels together; the closed-shell singlet/triplet
    // algorithm carries a single channel. (Same rule as `nloc_per_state` in
    // `setup_eigenvectors_X`, applied to the widened windows.)
    this->nloc_per_state_z_ = this->nk * (this->openshell
        ? this->paraX_z_[0].get_local_size() + this->paraX_z_[1].get_local_size()
        : this->paraX_z_[0].get_local_size());

#ifdef __MPI
    this->psi_ks_z_.reset(new psi::Psi<T>(this->kv.get_nks(), this->paraC_z_.get_col_size(),
        this->paraC_z_.get_row_size(), this->kv.ngk, true));
#else
    this->psi_ks_z_.reset(new psi::Psi<T>(this->kv.get_nks(), this->nbands_z_, this->nbasis, this->kv.ngk, true));
#endif
    this->eig_ks_z_.create(this->kv.get_nks(), this->nbands_z_);

    for (int ik = 0; ik < this->kv.get_nks(); ++ik)
    {
        // same redistribution `refresh_from_ks_` does for the X window, over more bands
#ifdef __MPI
        Cpxgemr2d(this->nbasis, this->nbands_z_, &(*this->psi_ks_all_)(ik, 0, 0), 1, start_band + 1,
            const_cast<int*>(desc_src), &(*this->psi_ks_z_)(ik, 0, 0), 1, 1,
            this->paraC_z_.desc, this->paraC_z_.blacs_ctxt);
#else
        for (int ib = 0; ib < this->nbands_z_; ++ib)
        {
            const auto* start = &(*this->psi_ks_all_)(ik, start_band + ib, 0);
            std::copy(start, start + this->nbasis, &(*this->psi_ks_z_)(ik, ib, 0));
        }
#endif
        for (int ib = 0; ib < this->nbands_z_; ++ib)
        { this->eig_ks_z_(ik, ib) = this->eig_ks_all(ik, start_band + ib); }
    }
    this->ofs_running_ << "Z-vector window: nbands = " << this->nbands_z_
        << " (X window: " << this->nbands << "), nvirt =";
    for (int is = 0; is < this->nspin; ++is)
    { this->ofs_running_ << " " << this->nvirt_z_[is] << "(X window: " << this->nvirt[is] << ")"; }
    this->ofs_running_ << std::endl;
    if (this->nbands_z_ < this->nbasis)
    {
        // The Z window can only be as wide as the ground state's band count, so a gradient is
        // converged in it only when `nbands` reaches the size of the AO basis. Measured on
        // `08_BeH2/rpa_at_lda` (NLOCAL 17), where Omega = eps_a - eps_i makes the finite
        // difference exact: nbands 11 -> 18% too high, nbands 17 -> 6 digits.
        this->ofs_running_ << " WARNING: the excited-state gradient is not converged with"
            " respect to the Z-vector (CPSCF) space: nbands = " << this->nbands_z_
            << " covers only part of the " << this->nbasis << " AO basis functions (NLOCAL)."
            " Set nbands = " << this->nbasis << " for a converged gradient; nvirt may stay as"
            " it is, since the excitation energies do not depend on this." << std::endl;
    }
}

template<typename T, typename TR>
void ModuleESolver::ESolver_LR<T, TR>::read_ks_chg(Charge& chg_gs)
{
    chg_gs.set_rhopw(this->pw_rho);
    const bool kin_den = XC_Functional::get_ked_flag() || (this->inp_->out_elf[0] > 0); // mohan add 20251202
    chg_gs.allocate(this->nspin, kin_den, XC_Functional::get_ked_flag(), this->inp_->test_charge);
    this->ofs_running_ << " try to read charge from file : ";
    for (int is = 0; is < this->nspin; ++is)
    {
        std::stringstream ssc;
        if (this->nspin == 1) { ssc << this->in_dir << "chg.cube"; }
        else { ssc << this->in_dir << "chgs" << is + 1 << ".cube"; }
        this->ofs_running_ << ssc.str() << std::endl;
        if (ModuleIO::read_vdata_palgrid(pgrid(),
            this->my_rank_,
            this->ofs_running_,
            ssc.str(),
            chg_gs.rho[is],
            this->ucell_->nat))
        this->ofs_running_ << " Read in the charge density: " << ssc.str() << std::endl;
    }
}
template class ModuleESolver::ESolver_LR<double, double>;
template class ModuleESolver::ESolver_LR<std::complex<double>, double>;
