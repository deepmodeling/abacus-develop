#include "esolver_lr_lcao_tddft.h"
#include "source_lcao/module_lr/utils/lr_util.h"
#include "source_lcao/module_lr/hamilt_casida.h"
#include "source_lcao/module_lr/hamilt_ulr.hpp"
#include "source_lcao/module_lr/potentials/pot_hxc_lrtd.h"
#include "source_lcao/lcao_nonlocal_info.h"
#include "source_lcao/module_lr/hsolver_lrtd.hpp"
#include "source_lcao/module_lr/lr_spectrum.h"
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
#include "source_lcao/module_lr/ri_benchmark/ri_benchmark.h"
#include "source_lcao/module_lr/operator_casida/operator_lr_diag.h" // for precondition
#ifdef __EXX
#include "source_lcao/module_ri/exx_lri_interface.h"
#include "source_hamilt/module_xc/exx_info.h"
#endif

#ifdef __EXX
template<>
void ModuleESolver::ESolver_LR<double>::move_exx_lri(std::shared_ptr<Exx_LRI<double>>& exx_ks)
{
    ModuleBase::TITLE("ESolver_LR<double>", "move_exx_lri");
    this->exx_lri = exx_ks;
    exx_ks = nullptr;
}
template<>
void ModuleESolver::ESolver_LR<std::complex<double>>::move_exx_lri(std::shared_ptr<Exx_LRI<std::complex<double>>>& exx_ks)
{
    ModuleBase::TITLE("ESolver_LR<complex>", "move_exx_lri");
    this->exx_lri = exx_ks;
    exx_ks = nullptr;
}
template<>
void ModuleESolver::ESolver_LR<std::complex<double>>::move_exx_lri(std::shared_ptr<Exx_LRI<double>>& exx_ks)
{
    throw std::runtime_error("ESolver_LR<std::complex<double>>::move_exx_lri: cannot move double to std::complex<double>");
}
template<>
void ModuleESolver::ESolver_LR<double>::move_exx_lri(std::shared_ptr<Exx_LRI<std::complex<double>>>& exx_ks)
{
    throw std::runtime_error("ESolver_LR<double>::move_exx_lri: cannot move std::complex<double> to double");
}
#endif

using namespace LR;

template<typename T, typename TR>
int ModuleESolver::ESolver_LR<T, TR>::cal_nupdown_form_occ(const ModuleBase::matrix& wg)
{   // only for nspin=2
    const int& nk = wg.nr / 2;
    auto occ_sum_k = [&](const int& is, const int& ib)->double { double o = 0.0; for (int ik = 0;ik < nk;++ik) { o += wg(is * nk + ik, ib); } return o;};
    int nupdown = 0;
    for (int ib = 0;ib < wg.nc;++ib)
    {
        const int nu = static_cast<int>(std::lround(occ_sum_k(0, ib)));
        const int nd = static_cast<int>(std::lround(occ_sum_k(1, ib)));
        if ((nu + nd) == 0) { break; }
        nupdown += nu - nd;
    }
    return nupdown;
}

template<typename T, typename TR>
void ModuleESolver::ESolver_LR<T, TR>::setup_2center_table(TwoCenterBundle& two_center_bundle, LCAO_Orbitals& orb, UnitCell& ucell)
{
    // set up 2-center table
#ifdef __FFT_TWO_CENTER
    two_center_bundle.tabulate();
#else
    two_center_bundle.tabulate(inp.lcao_ecut, inp.lcao_dk, inp.lcao_dr, inp.lcao_rmax);
#endif
    if (input.vnl_in_h)
    {
        auto* lcao_nl = new LCAONonlocalInfo();
        lcao_nl->setupNonlocal(ucell.ntype, ucell.atoms, GlobalV::ofs_running, orb,
                               input.basis_type, input.out_element_info,
                               input.lspinorb, input.nspin);
        ucell.infoNL.reset(lcao_nl);
        two_center_bundle.build_beta(ucell.ntype, lcao_nl->get_nonlocal().Beta);
    }
}

template<typename T, typename TR>
void ModuleESolver::ESolver_LR<T, TR>::parameter_check()const
{
    const std::set<std::string> lr_solvers = { "dav", "lapack" , "spectrum", "dav_subspace", "cg", "elpa", "plot" };
    const std::set<std::string> xc_kernels = { "rpa", "lda", "pwlda", "pbe", "hf", "hse", "bse" };
    const std::set<std::string> abs_gauge = { "velocity", "length" };
    if (lr_solvers.find(this->input.lr_solver) == lr_solvers.end()) {
        throw std::invalid_argument("ESolver_LR: unknown type of lr_solver");
    }
    if (xc_kernels.find(this->xc_kernel) == xc_kernels.end()) {
        throw std::invalid_argument("ESolver_LR: unknown type of xc_kernel");
    }
    if (this->nspin != 1 && this->nspin != 2) {
        throw std::invalid_argument("LR-TDDFT only supports nspin = 1 or 2 now");
    }
    if (abs_gauge.find(this->input.abs_gauge) == abs_gauge.end()) {
        throw std::invalid_argument("ESolver_LR: unknown type of abs_gauge");
    }
}

template<typename T, typename TR>
void ModuleESolver::ESolver_LR<T, TR>::set_dimension()
{
    this->nspin = input.nspin;
    this->nstates = input.lr_nstates;
    this->nbasis = PARAM.globalv.nlocal;
    int ks_nbands = input.nbands;
    this->nocc_max = LR_Util::cal_nocc(LR_Util::cal_nelec(*this->ucell_));
    if (input.ri_hartree_benchmark == "aims" || input.ri_hartree_benchmark == "aims-librpa"
        && !input.aims_nbasis.empty())
    {
        // calculate total number of basis funcs, see https://en.cppreference.com/w/cpp/algorithm/inner_product
        this->nbasis = std::inner_product(input.aims_nbasis.begin(), /* iterator1.begin */
                                          input.aims_nbasis.end(),  /* iterator1.end */
                                          this->ucell_->atoms,  /* iterator2.begin */
                                          0,  /* init value */
                                          std::plus<int>(), /* iter op1 */
                                          [](const int& a, const Atom& b) { return a * b.na; }); /* iter op2 */
        std::cout << "nbasis from aims: " << this->nbasis << std::endl;
        for (int it = 0; it < this->ucell_->ntype; ++it)
        {
            this->ucell_->atoms[it].nw = input.aims_nbasis[it];
        }
        const_cast<UnitCell*>(this->ucell_)->set_iat2iwt(1);  // update iat2iwt for aims_nbasis 25-05-23

        int nbands_file = 0;
        int nk_file = 0;
        int nspin_file = 0;
        int nocc_file = 0;
        LR_IO::parse_band_out_file(input.rpa_outdir, nbands_file, nk_file, nspin_file, nocc_file);
        std::cout << "nocc from band_out: " << nocc_file << std::endl;
        ks_nbands = nbands_file;
        this->nocc_max = nocc_file;
    }
    // calculate the number of occupied and unoccupied states
    // which determines the basis size of the excited states    
    this->nocc_in = std::max(1, std::min(input.nocc, this->nocc_max));
    this->nvirt_in = ks_nbands - this->nocc_max;   //nbands-nocc
    if (input.nvirt > this->nvirt_in) { GlobalV::ofs_running << "ESolver_LR: input nvirt is too large to cover by nbands, set nvirt = nbands - nocc = " << this->nvirt_in << std::endl; }
    else if (input.nvirt > 0) { this->nvirt_in = input.nvirt; }
    this->nbands = this->nocc_in + this->nvirt_in;
    this->nk = input.nspin == 2 ? this->kv.get_nks() / 2 : this->kv.get_nks();
    this->nocc.resize(nspin, nocc_in);
    this->nvirt.resize(nspin, nvirt_in);
    if (this->nstates <= 0) { 
        this->nstates = nk * nocc_in * nvirt_in;
        GlobalV::ofs_running << "ESolver_LR: lr_nstates <= 0, set nstates = nk * nocc * nvirt = " << this->nstates << std::endl;
    }
    for (int is = 0;is < nspin;++is) { this->npairs.push_back(nocc[is] * nvirt[is]); }
    GlobalV::ofs_running << "Setting LR-TDDFT parameters: " << std::endl;
    GlobalV::ofs_running << "number of occupied bands: " << nocc_in << std::endl;
    GlobalV::ofs_running << "number of virtual bands: " << nvirt_in << std::endl;
    GlobalV::ofs_running << "number of Atom orbitals (LCAO-basis size): " << this->nbasis << std::endl;
    GlobalV::ofs_running << "number of KS bands: " << this->eig_ks.nc << std::endl;
    GlobalV::ofs_running << "number of excited states to be solved: " << this->nstates << std::endl;
}

template<typename T, typename TR>
void ModuleESolver::ESolver_LR<T, TR>::reset_dim_spin2()
{
	if (nspin != 2) 
	{ 
		return; 
	}
	if (nupdown == 0) 
	{ 
		std::cout << " ** Assuming degenerate spin-up and spin-down states  **" << std::endl; 
	}
	else
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
	if (input.lr_unrestricted) 
	{ 
		this->openshell = true; 
	}
}

template <typename T, typename TR>
ModuleESolver::ESolver_LR<T, TR>::ESolver_LR(const Input_para& inp)
    : input(inp)
#ifdef __EXX
    , exx_info(GlobalC::exx_info)
#endif
{
}

template <typename T, typename TR>
void ModuleESolver::ESolver_LR<T, TR>::before_all_runners(BaseCell& basecell, const Input_para& inp)
{
    basecell.require_kind(BaseCell::Kind::unit_cell, __FUNCTION__);
    UnitCell& ucell = static_cast<UnitCell&>(basecell);
    this->ucell_ = &ucell;
    if (inp.esolver_type == "ks-lr")
    {
        ModuleESolver::ESolver_KS_LCAO<T, TR> ks_solver;
        ks_solver.before_all_runners(basecell, inp);
        ks_solver.runner(basecell, 0);
        this->initialize_from_ks_(std::move(ks_solver), ucell, inp);
    }
    else
    {
        this->initialize_from_unitcell_(ucell, inp);
    }
}

template <typename T, typename TR>
void ModuleESolver::ESolver_LR<T, TR>::initialize_from_ks_(ModuleESolver::ESolver_KS_LCAO<T, TR>&& ks_sol,
                                                 UnitCell& ucell,
                                                 const Input_para& inp)
{
    ModuleBase::TITLE("ESolver_LR", "ESolver_LR(KS)");

	if (this->input.lr_solver == "spectrum") 
	{
		throw std::invalid_argument("when lr_solver==spectrum, esolver_type must be `lr` to skip KS calculation.");
	}

    this->gd = std::move(ks_sol.gd);

    // xc kernel
    this->xc_kernel = LR_Util::tolower(inp.xc_kernel);
    //kv
    this->kv = std::move(ks_sol.kv);

    this->parameter_check();

    this->set_dimension();

    // setup_wd_division is not need to be covered in #ifdef __MPI, see its implementation
    LR_Util::setup_2d_division(this->paraMat_, 1, this->nbasis, this->nbasis);

    this->paraMat_.atom_begin_row = std::move(ks_sol.pv.atom_begin_row);
    this->paraMat_.atom_begin_col = std::move(ks_sol.pv.atom_begin_col);
    this->paraMat_.iat2iwt_ = ucell.get_iat2iwt();

    LR_Util::setup_2d_division(this->paraC_, 1, this->nbasis, this->nbands
#ifdef __MPI
        , this->paraMat_.blacs_ctxt
#endif
    );

    auto move_gs = [&, this]() -> void  // move the ground state info
        {
            this->psi_ks = ks_sol.psi;
            ks_sol.psi = nullptr;
            //only need the eigenvalues. the 'elecstates' of excited states is different from ground state.
            this->eig_ks = std::move(ks_sol.pelec->ekb);
        };
#ifdef __MPI
	if (this->nbands == input.nbands)
	{ 
		move_gs(); 
	}
    else    // copy the part of ground state info according to paraC_
    {
        this->psi_ks = new psi::Psi<T>(this->kv.get_nks(), 
                                       this->paraC_.get_col_size(), 
                                       this->paraC_.get_row_size(),
                                       this->kv.ngk,
                                       true);
        this->eig_ks.create(this->kv.get_nks(), this->nbands);
        const int start_band = this->nocc_max - *std::max_element(nocc.begin(), nocc.end());
        for (int ik = 0;ik < this->kv.get_nks();++ik)
        {
            Cpxgemr2d(this->nbasis, this->nbands, &(*ks_sol.psi)(ik, 0, 0), 1, start_band + 1, ks_sol.pv.desc_wfc,
                &(*this->psi_ks)(ik, 0, 0), 1, 1, this->paraC_.desc, this->paraC_.blacs_ctxt);
            for (int ib = 0;ib < this->nbands;++ib) { this->eig_ks(ik, ib) = ks_sol.pelec->ekb(ik, start_band + ib); }
        }
    }
#else
    move_gs();
#endif
    if (nspin == 2)
    {
        this->nupdown = cal_nupdown_form_occ(ks_sol.pelec->wg);
        reset_dim_spin2();
    }
    this->gint_info_ = std::move(ks_sol.gint_info_);
    // move pw basis
    if (this->pw_rho_flag)
    {
        this->pw_rho_flag = true;
        delete this->pw_rho;    // newed in ESolver_FP::ESolver_FP
    }
    this->pw_rho = ks_sol.pw_rho;
    ks_sol.pw_rho = nullptr;
    //init potential and calculate kernels using ground state charge
    init_pot(*ks_sol.pelec->charge);

#ifdef __EXX
    if (xc_kernel == "hf" || xc_kernel == "hse")
    {
        // if the same kernel is calculated in the esolver_ks, move it
        std::string dft_functional = LR_Util::tolower(input.dft_functional);
        if (ks_sol.exx_nao.exd && std::is_same<T, double>::value && xc_kernel == dft_functional) {
            this->move_exx_lri(ks_sol.exx_nao.exd->exx_ptr);
        } else if (ks_sol.exx_nao.exc && std::is_same<T, std::complex<double>>::value && xc_kernel == dft_functional) {
            this->move_exx_lri(ks_sol.exx_nao.exc->exx_ptr);
        } else    // construct C, V from scratch
        {
            // set ccp_type according to the xc_kernel
            if (xc_kernel == "hf") { exx_info.info_global.ccp_type = Conv_Coulomb_Pot_K::Ccp_Type::Hf; }
            else if (xc_kernel == "hse") { exx_info.info_global.ccp_type = Conv_Coulomb_Pot_K::Ccp_Type::Erfc; }
            exx_info.sync_from_global();
            this->exx_lri = std::make_shared<Exx_LRI<T>>(exx_info.info_ri);
            this->exx_lri->init(MPI_COMM_WORLD, ucell,this->kv, ks_sol.orb_);
            this->exx_lri->cal_exx_ions(ucell,input.out_ri_cv);
        }
    }
#endif
    this->pelec = new elecstate::ElecStateLCAO<T>();
    orb_cutoff_ = ks_sol.orb_.cutoffs();
    if (LR_Util::tolower(input.abs_gauge) == "velocity")
    {
        this->two_center_bundle_ = std::move(ks_sol.two_center_bundle_);
    }
}

template <typename T, typename TR>
void ModuleESolver::ESolver_LR<T, TR>::initialize_from_unitcell_(UnitCell& ucell, const Input_para& inp)
{
    ModuleBase::TITLE("ESolver_LR", "ESolver_LR(from scratch)");
    // xc kernel
    this->xc_kernel = LR_Util::tolower(inp.xc_kernel);

    // necessary steps in ESolver_FP
    ESolver_FP::before_all_runners(ucell, inp);
    this->pelec = new elecstate::ElecStateLCAO<T>();

    // necessary steps in ESolver_KS::before_all_runners : symmetry and k-points
    if (ModuleSymmetry::Symmetry::symm_flag == 1)
    {
        const int cal_symm_repr[2] = {input.cal_symm_repr[0], input.cal_symm_repr[1]};
        ucell.symm.analy_sys(ucell.lat, ucell.st, ucell.atoms, GlobalV::ofs_running,
                             input.symmetry_prec, input.nspin, input.calculation, cal_symm_repr);
        ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "SYMMETRY");
    }
    const bool use_ibz = false;
    const std::string global_out_dir = PARAM.globalv.global_out_dir;
    const bool gamma_only_local = PARAM.globalv.gamma_only_local;
    const double kspacing[3] = {input.kspacing[0], input.kspacing[1], input.kspacing[2]};
    const std::string kmesh_type = input.kmesh_type;
    const double koffset[3] = {input.koffset[0], input.koffset[1], input.koffset[2]};
    this->kv.set(ucell, ucell.symm, input.kpoint_file, input.nspin, ucell.G, ucell.latvec, GlobalV::ofs_running, use_ibz, global_out_dir, gamma_only_local, kspacing, kmesh_type, koffset);
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "INIT K-POINTS");
    ModuleIO::print_parameters(ucell, this->kv, inp);

    this->parameter_check();

    /// read orbitals and build the interpolation table
    two_center_bundle_.build_orb(ucell.ntype, ucell.orbital_fn.data(), inp.orbital_dir);

    LCAO_Orbitals orb;
    two_center_bundle_.to_LCAO_Orbitals(orb, inp.lcao_ecut, inp.lcao_dk, inp.lcao_dr, inp.lcao_rmax,
                                        inp.out_element_info, inp.cal_force);
    orb_cutoff_ = orb.cutoffs();
    if (LR_Util::tolower(input.abs_gauge) == "velocity")
    {
        setup_2center_table(this->two_center_bundle_, orb, ucell);
    }

    this->set_dimension();
    //  setup 2d-block distribution for AO-matrix and KS wfc
    LR_Util::setup_2d_division(this->paraMat_, 1, this->nbasis, this->nbasis);
#ifdef __MPI
    this->paraMat_.set_desc_wfc_Eij(this->nbasis, this->nbands, paraMat_.get_row_size());
    int err = this->paraMat_.set_nloc_wfc_Eij(this->nbands, GlobalV::ofs_running, GlobalV::ofs_warning);
    this->paraMat_.set_atomic_trace(ucell.get_iat2iwt(), ucell.nat, this->nbasis);
    if (input.ri_hartree_benchmark != "aims") { this->paraMat_.set_atomic_trace(ucell.get_iat2iwt(), ucell.nat, this->nbasis); }
#else
    this->paraMat_.nrow_bands = this->nbasis;
    this->paraMat_.ncol_bands = this->nbands;
#endif

    // read the ground state info
    // now ModuleIO::read_wfc_nao needs `Parallel_Orbitals` and can only read all the bands
    // it need improvement to read only the bands needed
    this->psi_ks = new psi::Psi<T>(this->kv.get_nks(),
                                   this->paraMat_.ncol_bands,
                                   this->paraMat_.get_row_size(), 
                                   this->kv.ngk,
                                   true);
    this->read_ks_wfc();
    if (nspin == 2)
    {
        this->nupdown = cal_nupdown_form_occ(this->pelec->wg);
        reset_dim_spin2();
    }

    LR_Util::setup_2d_division(this->paraC_, 1, this->nbasis, this->nbands
#ifdef __MPI
        , paraMat_.blacs_ctxt
#endif
    );

    // clear ks info, new elecstate for excition
    this->pelec = new elecstate::ElecState();

    // read the ground state charge density and calculate xc kernel
    Pgrid.init(this->pw_rho->nx,
        this->pw_rho->ny,
        this->pw_rho->nz,
        this->pw_rho->nplane,
        this->pw_rho->nrxx,
        pw_big->nbz,
        pw_big->bz,
        GlobalV::NPROC);
    Charge chg_gs;
    if (input.ri_hartree_benchmark == "none") { this->read_ks_chg(chg_gs); }
    this->init_pot(chg_gs);

    // search adjacent atoms and init Gint
    double search_radius = -1.0;
    search_radius = atom_arrange::set_sr_NL(GlobalV::ofs_running,
        input.out_level,
        orb.get_rcutmax_Phi(),
        ucell.infoNL->get_rcutmax_Beta(),
        PARAM.globalv.gamma_only_local);
    atom_arrange::search(PARAM.globalv.search_pbc,
                         GlobalV::ofs_running,
                         this->gd,
                         *this->ucell_,
                         search_radius,
                         input.test_atom_input);
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
        this->gd));
    ModuleGint::Gint::set_gint_info(gint_info_.get());
    // if EXX from scratch, init 2-center integral and calculate Cs, Vs 
#ifdef __EXX
    if ((xc_kernel == "hf" || xc_kernel == "hse") && this->input.lr_solver != "spectrum")
    {
        // set ccp_type according to the xc_kernel
        if (xc_kernel == "hf") { exx_info.info_global.ccp_type = Conv_Coulomb_Pot_K::Ccp_Type::Hf; }
        else if (xc_kernel == "hse") { exx_info.info_global.ccp_type = Conv_Coulomb_Pot_K::Ccp_Type::Erfc; }
        this->exx_lri = std::make_shared<Exx_LRI<T>>(exx_info.info_ri);
        this->exx_lri->init(MPI_COMM_WORLD, ucell,this->kv, orb);
        this->exx_lri->cal_exx_ions(ucell,input.out_ri_cv);
    }
    // else
#endif
        // ModuleBase::Ylm::set_coefficients() is deprecated
}

template <typename T, typename TR>
void ModuleESolver::ESolver_LR<T, TR>::runner(BaseCell& basecell, const int istep)
{
    basecell.require_kind(BaseCell::Kind::unit_cell, __FUNCTION__);
    UnitCell& ucell = static_cast<UnitCell&>(basecell);

    ModuleBase::TITLE("ESolver_LR", "runner");
    ModuleBase::timer::start("ESolver_LR", "runner");
    //allocate 2-particle state and setup 2d division
    this->setup_eigenvectors_X();
    this->pelec->ekb.create(nspin, this->nstates);

    auto efile_out = [&](const std::string& label)->std::string {return PARAM.globalv.global_out_dir + "Excitation_Energy_" + label + ".dat";};
    auto vfile_out = [&](const std::string& label)->std::string {return PARAM.globalv.global_out_dir + "Excitation_Amplitude_" + label + "_" + std::to_string(GlobalV::MY_RANK) + ".dat";};
    if (this->input.lr_solver == "elpa")
    {
        ModuleBase::WARNING_QUIT("ESolver_LR", "ESolver_LR doesn't support elpa now.");

    }
    else if (this->input.lr_solver != "spectrum")
    {
        auto write_states = [&](const std::string& label, const Real<T>* e, const T* v, const int& dim, const int& nst, const int& prec = 8)->void
            {
                if (GlobalV::MY_RANK == 0) { assert(nst == LR_Util::write_value(efile_out(label), prec, e, nst)); }
                assert(nst * dim == LR_Util::write_value(vfile_out(label), prec, v, nst, dim));
            };
        std::vector<double> precondition(this->input.lr_solver == "lapack" ? 0 : nloc_per_state, 1.0);
        // allocate and initialize A matrix and density matrix
        if (openshell)
        {
            for (int is : {0, 1})
            {
                if (input.lr_solver != "lapack") {
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
                              this->gd,
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
                               this->pelec->ekb.c, this->input.lr_solver,
                               this->input.lr_thr, precondition);
            if (input.out_wfc_lr) { write_states("openshell", this->pelec->ekb.c, this->X[0].template data<T>(), nloc_per_state, nstates); }
        }
        else
        {
            if (input.lr_solver != "lapack") {
                OperatorLRDiag<double> pre_op(this->eig_ks.c, this->paraX_[0], this->nk, this->nocc[0], this->nvirt[0]);
                pre_op.act(1, nloc_per_state, 1, precondition.data(), precondition.data()); 
            }
            auto spin_types = std::vector<std::string>({ "singlet", "triplet" });
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
                                this->gd,
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
                                input.ri_hartree_benchmark,
                                (input.ri_hartree_benchmark == "aims" ? input.aims_nbasis : std::vector<int>({})));
                LR::HSolver::solve(hlr, this->X[is].template data<T>(), nloc_per_state, nstates,
                                   this->nk, this->nocc, this->nvirt, this->paraX_,
                                   this->pelec->ekb.c + is * nstates,
                                   this->input.lr_solver,
                                   this->input.lr_thr,
                                   precondition);
                if (input.out_wfc_lr) { write_states(spin_types[is], this->pelec->ekb.c + is * nstates, this->X[is].template data<T>(), nloc_per_state, nstates); }
            }
        }
    }
    else    // lr_solver == "spectrum", read the eigenvalues
    {
        auto efile_in = [&](const std::string& label)->std::string {return PARAM.globalv.global_readin_dir + "Excitation_Energy_" + label + ".dat";};
        auto vfile_in = [&](const std::string& label)->std::string {return PARAM.globalv.global_readin_dir + "Excitation_Amplitude_" + label + "_" + std::to_string(GlobalV::MY_RANK) + ".dat";};
    
        auto read_states = [&](const std::string& label, Real<T>* e, T* v, const int& dim, const int& nst)->void
            {
                if (GlobalV::MY_RANK == 0) {
                    assert(nst == LR_Util::read_value(efile_in(label), e, nst));
                    std::cout <<"Rank "<< GlobalV::MY_RANK << ": finish reading " << efile_in(label) << std::endl;
                }
#ifdef __MPI
// in velocity gauge, the eigenvalues may be used to calculate the transition dipole, so we'd better broadcast them
                MPI_Bcast(e, nst, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif
                assert(nst * dim == LR_Util::read_value(vfile_in(label), v, nst, dim));
                std::cout <<"Rank "<< GlobalV::MY_RANK << ": finish reading " << vfile_in(label) << std::endl;
            };
        std::cout << "reading the excitation states from file: \n";
        if (openshell)
        {
            read_states("openshell", this->pelec->ekb.c, this->X[0].template data<T>(), nloc_per_state, nstates);
        }
        else
        {
            auto spin_types = std::vector<std::string>({ "singlet", "triplet" });
            for (int is = 0;is < nspin;++is) { read_states(spin_types[is], this->pelec->ekb.c + is * nstates, this->X[is].template data<T>(), nloc_per_state, nstates); }
        }
    }
    ModuleBase::timer::end("ESolver_LR", "runner");
    return;
}

template <typename T, typename TR>
void ModuleESolver::ESolver_LR<T, TR>::after_all_runners(BaseCell& basecell)
{
    basecell.require_kind(BaseCell::Kind::unit_cell, __FUNCTION__);
    UnitCell& ucell = static_cast<UnitCell&>(basecell);

    ModuleBase::TITLE("ESolver_LR", "after_all_runners");
    if (input.ri_hartree_benchmark != "none") { return; } //no need to calculate the spectrum in the benchmark routine
    //cal spectrum
    if (LR_Util::tolower(this->input.abs_gauge) == "velocity" )
    {
        const int nspin_tmp = input.nspin == 2 ? 2 : 1;
        this->velocity_mo = LR_Util::cal_velocity_mo(*this->ucell_, this->gd, this->two_center_bundle_, 
                                                    this->paraMat_, this->paraC_, this->kv, *this->psi_ks, 
                                                    this->nk, nspin_tmp, this->nbasis, this->nocc, this->nvirt);
    }
    std::vector<double> freq(100);
    std::vector<double> abs_wavelen_range({ 20, 200 });//default range
    if (input.abs_wavelen_range.size() >= 2 && std::abs(input.abs_wavelen_range[1] - input.abs_wavelen_range[0]) > 0.02)
    {
        abs_wavelen_range = input.abs_wavelen_range;
    }
    double lambda_diff = std::abs(abs_wavelen_range[1] - abs_wavelen_range[0]);
    double lambda_min = std::min(abs_wavelen_range[1], abs_wavelen_range[0]);
    for (int i = 0;i < freq.size();++i) { freq[i] = 91.126664 / (lambda_min + 0.01 * static_cast<double>(i + 1) * lambda_diff); }
    auto spin_types = (nspin == 2 && !openshell) ? std::vector<std::string>({ "singlet", "triplet" }) : std::vector<std::string>({ "updown" });
    for (int is = 0;is < this->X.size();++is)
    {
        LR_Spectrum<T> spectrum(nspin, this->nbasis, this->nocc, this->nvirt, *this->pw_rho, *this->psi_ks,
            *this->ucell_, this->kv, this->gd, this->orb_cutoff_, this->two_center_bundle_,
            this->paraX_, this->paraC_, this->paraMat_,
            &this->pelec->ekb.c[is * nstates], this->eig_ks.c, this->X[is].template data<T>(), nstates, openshell,
            LR_Util::tolower(input.abs_gauge));
        if (LR_Util::tolower(this->input.abs_gauge) == "velocity" ) {spectrum.set_vmo(this->velocity_mo.data());}
        spectrum.cal_spectrum();
        spectrum.transition_analysis(spin_types[is]+"_tda");
        if (spin_types[is] != "triplet")        // triplets has no transition dipole and no contribution to the spectrum
        {
            spectrum.optical_absorption_method1(freq, input.abs_broadening);
            spectrum.write_transition_dipole(PARAM.globalv.global_out_dir + 
                "trans_dipole_" + spin_types[is] + "_tda.dat");
            // =============================================== for test ====================================================
            // spectrum.optical_absorption_method2(freq, input.abs_broadening);
            // if (LR_Util::tolower(input.abs_gauge) == "velocity")
            // {   // TEST the formula v/omega rather than v/(e_a-e_i)
            //     spectrum.test_transition_dipoles_velocity_omega();
            //     spectrum.write_transition_dipole(PARAM.globalv.global_out_dir + 
            //         "trans_dipole_" + spin_types[is] + "_vomega_tda.dat");
            // }
            // =============================================== for test ====================================================
        }
    }
}

template<typename T, typename TR>
void ModuleESolver::ESolver_LR<T, TR>::setup_eigenvectors_X()
{
    ModuleBase::TITLE("ESolver_LR", "setup_eigenvectors_X");
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

    auto spin_types = (nspin == 2 && !openshell) ? std::vector<std::string>({ "singlet", "triplet" }) : std::vector<std::string>({ "updown" });
    // if spectrum-only, read the LR-eigenstates from file and return
    if (this->input.lr_solver != "spectrum") { set_X_initial_guess(); }
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
        GlobalV::ofs_running << "setting the initial guess of X of spin" << is << std::endl;
        if (no >= 2 && eig_ks.nc > no) { GlobalV::ofs_running << "E_{lumo}-E_{homo-1}=" << eig_ks(is, no) - eig_ks(is, no - 2) << std::endl; }
        if (no >= 1 && eig_ks.nc > no + 1) { GlobalV::ofs_running << "E_{lumo+1}-E{homo}=" << eig_ks(is, no + 1) - eig_ks(is, no - 1) << std::endl; }
        GlobalV::ofs_running << "mode of X-index: " << ix_mode << std::endl;

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
    this->pot.resize(nspin, nullptr);
    if (this->input.ri_hartree_benchmark != "none") { return; } //no need to initialize potential for Hxc kernel in the RI-benchmark routine
    switch (nspin)
    {
        using ST = PotHxcLR::SpinType;
    case 1:
        this->pot[0] = std::make_shared<PotHxcLR>(xc_kernel, *this->pw_rho, *this->ucell_, chg_gs, Pgrid, ST::S1, input.lr_init_xc_kernel);
        break;
    case 2:
        this->pot[0] = std::make_shared<PotHxcLR>(xc_kernel, *this->pw_rho, *this->ucell_, chg_gs, Pgrid, openshell ? ST::S2_updown : ST::S2_singlet, input.lr_init_xc_kernel);
        this->pot[1] = std::make_shared<PotHxcLR>(xc_kernel, *this->pw_rho, *this->ucell_, chg_gs, Pgrid, openshell ? ST::S2_updown : ST::S2_triplet, input.lr_init_xc_kernel);
        break;
    default:
        throw std::invalid_argument("ESolver_LR: nspin must be 1 or 2");
    }
}

template<typename T, typename TR>
void ModuleESolver::ESolver_LR<T, TR>::read_ks_wfc()
{
    assert(this->psi_ks != nullptr);
    this->pelec->ekb.create(this->kv.get_nks(), this->nbands);
    this->pelec->wg.create(this->kv.get_nks(), this->nbands);

    if (input.ri_hartree_benchmark == "aims")        // for aims benchmark
    {
#ifdef __EXX
        int ncore = 0;
        std::vector<double> eig_ks_vec = RI_Benchmark::read_aims_ebands<double>(PARAM.globalv.global_readin_dir + "band_out", nocc_in, nvirt_in, ncore);
        std::cout << "ncore=" << ncore << ", nocc=" << nocc_in << ", nvirt=" << nvirt_in << ", nbands=" << this->nbands << std::endl;
        std::cout << "eig_ks_vec.size()=" << eig_ks_vec.size() << std::endl;
        if(eig_ks_vec.size() != this->nbands) {ModuleBase::WARNING_QUIT("ESolver_LR", "read_aims_ebands failed.");};
        for (int i = 0;i < nbands;++i) { this->pelec->ekb(0, i) = eig_ks_vec[i]; }
        RI_Benchmark::read_aims_eigenvectors<T>(*this->psi_ks, PARAM.globalv.global_readin_dir + "KS_eigenvectors.out", ncore, nbands, nbasis);
#else
        ModuleBase::WARNING_QUIT("ESolver_LR", "RI benchmark is only supported when compile with LibRI.");
#endif
    }
	else if (!ModuleIO::read_wfc_nao(PARAM.globalv.global_readin_dir, this->paraMat_, *this->psi_ks, 
				this->pelec->ekb,
                this->pelec->wg,
				this->kv.ik2iktot,
				this->kv.get_nkstot(),
                input.nspin,
				/*skip_bands=*/this->nocc_max - this->nocc_in)) {
        ModuleBase::WARNING_QUIT("ESolver_LR", "read ground-state wavefunction failed.");
    }
    this->eig_ks = std::move(this->pelec->ekb);
}

template<typename T, typename TR>
void ModuleESolver::ESolver_LR<T, TR>::read_ks_chg(Charge& chg_gs)
{
    chg_gs.set_rhopw(this->pw_rho);
    const bool kin_den = chg_gs.kin_density(); // mohan add 20251202
    chg_gs.allocate(this->nspin, kin_den);
    GlobalV::ofs_running << " try to read charge from file : ";
    for (int is = 0; is < this->nspin; ++is)
    {
        std::stringstream ssc;
        ssc << PARAM.globalv.global_readin_dir << "chgs" << is + 1 << ".cube";
        GlobalV::ofs_running << ssc.str() << std::endl;
        if (ModuleIO::read_vdata_palgrid(Pgrid,
            GlobalV::MY_RANK,
            GlobalV::ofs_running,
            ssc.str(),
            chg_gs.rho[is],
            this->ucell_->nat)) {
            GlobalV::ofs_running << " Read in the charge density: " << ssc.str() << std::endl;
        } else {    // prenspin for nspin=4 is not supported currently
            ModuleBase::WARNING_QUIT(
                "init_rho",
                "!!! Couldn't find the charge file !!! The default directory \n of " + ssc.str() +" is OUT.suffix, "
                "or you must set read_file_dir \n to a specific directory. ");
        }
    }
}
template class ModuleESolver::ESolver_LR<double, double>;
template class ModuleESolver::ESolver_LR<std::complex<double>, double>;
