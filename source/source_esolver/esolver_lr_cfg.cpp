// Resolve LR dimensions, validate solver options, and initialize kernel potentials.
#include "esolver_lr_lcao_tddft.h"
#include "source_base/parallel_reduce.h"
#include "source_lcao/module_lr/utils/lr_util.h"
#include "source_lcao/module_lr/utils/lr_io.h"
#include "source_lcao/lcao_nonlocal_info.h"
#include "source_io/module_wf/read_wfc_nao.h"
#include "source_io/module_parameter/parameter.h"
#include <algorithm>
#include <cmath>
#include <set>

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
    double up = 0.0;
    double dn = 0.0;
    for (int ib = 0;ib < wg.nc;++ib) { up += occ_sum_k(0, ib); dn += occ_sum_k(1, ib); }
    // wg is replicated within a pool, but each pool holds different k points.
    if (this->kv.para_k.kpar > 1)
    {
        if (this->kv.para_k.rank_in_pool != 0)
        {
            up = 0.0;
            dn = 0.0;
        }
        Parallel_Reduce::reduce_all(up);
        Parallel_Reduce::reduce_all(dn);
    }
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
}

template<typename T, typename TR>
void ModuleESolver::ESolver_LR<T, TR>::set_dimension()
{
    this->nspin = this->inp_->nspin;
    this->nstates = this->inp_->lr_nstates;
    this->nbasis = PARAM.globalv.nlocal;
    int ks_nbands = this->inp_->nbands;
    if (this->nspin == 2)
    {
        this->nupdown = static_cast<int>(std::lround(this->inp_->nupdown));
        if (this->ks_)
        {
            this->nupdown = cal_nupdown_form_occ(this->ks_->pelec->wg);
        }
        else if (this->inp_->ri_hartree_benchmark != "aims"
                 && this->inp_->ri_hartree_benchmark != "aims-librpa")
        {
            std::vector<double> populations;
            const bool gamma_only = std::is_same<T, double>::value;
            const bool binary = this->inp_->init_wfc_file_format == "binary";
            if (!ModuleIO::read_wfc_nao_spin_populations(this->in_dir, this->kv.get_nkstot(), this->nspin, gamma_only, binary,
                    this->my_rank_, populations))
            {
                ModuleBase::WARNING_QUIT("ESolver_LR", "read complete KS occupations failed");
            }
            this->nupdown = static_cast<int>(std::lround(populations[0]) - std::lround(populations[1]));
        }
    }
    // read_pseudo/ParamUpdater has already resolved nelec and applied nelec_delta.
    this->nocc_max = LR_Util::cal_nocc(this->inp_->nelec, this->nspin, this->nupdown);
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
    this->nocc_in = LR_Util::cal_nocc_window(this->inp_->nocc, this->nocc_max);
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
    this->ofs_running_ << "full occupied bands in largest spin channel: " << nocc_max << std::endl;
    const int skipped_core_bands = nocc_max - nocc_in;
    this->ofs_running_ << "shared core bands omitted from LR window: " << skipped_core_bands << std::endl;
    this->ofs_running_ << "number of occupied bands: " << nocc_in << std::endl;
    this->ofs_running_ << "number of virtual bands: " << nvirt_in << std::endl;
    this->ofs_running_ << "number of Atom orbitals (LCAO-basis size): " << this->nbasis << std::endl;
    // eig_ks is allocated after set_dimension; report the resolved counts here.
    this->ofs_running_ << "number of KS bands: " << ks_nbands << std::endl;
    this->ofs_running_ << "number of KS bands in LR window: " << this->nbands << std::endl;
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
        const int minority_spin = nupdown > 0 ? 1 : 0;
        const int polarization = std::abs(nupdown);
        nocc[minority_spin] -= polarization;
        nvirt[minority_spin] += polarization;
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
        // literal "default" string here would always disagree with `xc_kernel`, even when the
        // ground state and LR use the same functional. Resolve it before comparing the names
        // or constructing a separate GS kernel.
        const std::string xc_kernel_gs = (this->inp_->dft_functional == "default")
            ? LR_Util::tolower(this->ucell_->atoms[0].ncpp.xc_func)
            : LR_Util::tolower(this->inp_->dft_functional);
        // `ST::S1` is only correct when nspin=1. `PotHxcLR` builds its `KernelXC` with
        // the input `nspin`, so at nspin=2 the kernel arrays carry 3 spin components per grid point
        // while the S1 integrand indexes them as if there were 1 -- it does not even read a
        // consistent spin combination. Use `ST::S2_gs` there, which is exactly half of S2_singlet,
        // matching the `K_Hxc(singlet) = 2 * pot_hxc_gs` convention of the gradient operators.
        const ST st_gs = (nspin == 1) ? ST::S1_gs : (oshell ? ST::S2_updown : ST::S2_gs);
        // GS response terms retain this kernel. The g^xc terms originating from the LR K
        // contribution to W and the force use the corresponding LR potential instead.
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

template class ModuleESolver::ESolver_LR<double, double>;
template class ModuleESolver::ESolver_LR<std::complex<double>, double>;
