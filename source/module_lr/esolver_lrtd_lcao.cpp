#include "esolver_lrtd_lcao.h"
#include "utils/gint_move.hpp"
#include "utils/lr_util.h"
#include "hamilt_casida.h"
#include "module_lr/potentials/pot_hxc_lrtd.h"
#include "module_lr/hsolver_lrtd.h"
#include "module_lr/lr_spectrum.h"
#include <memory>
#include "module_hamilt_lcao/hamilt_lcaodft/hamilt_lcao.h"
#include "module_io/read_wfc_nao.h"
#include "module_io/cube_io.h"
#include "module_io/print_info.h"
#include "module_cell/module_neighbor/sltk_atom_arrange.h"
#include "module_lr/utils/lr_util_print.h"
#include "module_base/scalapack_connector.h"
#include "module_parameter/parameter.h"
#include "module_lr/ri_benchmark/ri_benchmark.h"

#ifdef __EXX
template<>
void LR::ESolver_LR<double>::move_exx_lri(std::shared_ptr<Exx_LRI<double>>& exx_ks)
{
    ModuleBase::TITLE("ESolver_LR<double>", "move_exx_lri");
    this->exx_lri = exx_ks;
    exx_ks = nullptr;
}
template<>
void LR::ESolver_LR<std::complex<double>>::move_exx_lri(std::shared_ptr<Exx_LRI<std::complex<double>>>& exx_ks)
{
    ModuleBase::TITLE("ESolver_LR<complex>", "move_exx_lri");
    this->exx_lri = exx_ks;
    exx_ks = nullptr;
}
template<>
void LR::ESolver_LR<std::complex<double>>::move_exx_lri(std::shared_ptr<Exx_LRI<double>>& exx_ks)
{
    throw std::runtime_error("ESolver_LR<std::complex<double>>::move_exx_lri: cannot move double to complex<double>");
}
template<>
void LR::ESolver_LR<double>::move_exx_lri(std::shared_ptr<Exx_LRI<std::complex<double>>>& exx_ks)
{
    throw std::runtime_error("ESolver_LR<double>::move_exx_lri: cannot move complex<double> to double");
}
#endif
template<>void LR::ESolver_LR<double>::set_gint() { this->gint_ = &this->gint_g_;this->gint_g_.gridt = &this->gt_; }
template<>void LR::ESolver_LR<std::complex<double>>::set_gint() { this->gint_ = &this->gint_k_; this->gint_k_.gridt = &this->gt_; }

inline double getreal(std::complex<double> x) { return x.real(); }
inline double getreal(double x) { return x; }

inline void redirect_log(const bool& out_alllog)
{
    GlobalV::ofs_running.close();
    std::stringstream   ss;
    if (out_alllog)
    {
        ss << PARAM.globalv.global_out_dir << "running_lr_" << GlobalV::MY_RANK + 1 << ".log";
        GlobalV::ofs_running.open(ss.str());
    }
    else
    {
        if (GlobalV::MY_RANK == 0)
        {
            ss << PARAM.globalv.global_out_dir << "running_lr.log";
            GlobalV::ofs_running.open(ss.str());
        }
    }
}

template<typename T, typename TR>
void LR::ESolver_LR<T, TR>::parameter_check()const
{
    std::set<std::string> lr_solvers = { "dav", "lapack" , "spectrum", "dav_subspace" };
    std::set<std::string> xc_kernels = { "rpa", "lda", "pbe", "hf" , "hse" };
    if (lr_solvers.find(this->input.lr_solver) == lr_solvers.end()) {
        throw std::invalid_argument("ESolver_LR: unknown type of lr_solver");
}
    if (xc_kernels.find(this->xc_kernel) == xc_kernels.end()) {
        throw std::invalid_argument("ESolver_LR: unknown type of xc_kernel");
}
    if (this->nspin != 1 && this->nspin != 2) {
        throw std::invalid_argument("LR-TDDFT only supports nspin = 1 or 2 now");
}
}

template<typename T, typename TR>
void LR::ESolver_LR<T, TR>::set_dimension()
{
    this->nspin = PARAM.inp.nspin;
    if (nspin == 2) { std::cout << "** Assuming the spin-up and spin-down states are degenerate. **" << std::endl;
}
    this->nstates = input.lr_nstates;
    this->nbasis = PARAM.globalv.nlocal;
    // calculate the number of occupied and unoccupied states
    // which determines the basis size of the excited states
    this->nocc_max = LR_Util::cal_nocc(LR_Util::cal_nelec(ucell));
    this->nocc = std::max(1, std::min(input.nocc, this->nocc_max));
    this->nvirt = PARAM.inp.nbands - this->nocc_max;   //nbands-nocc
    if (input.nvirt > this->nvirt) {
        GlobalV::ofs_warning << "ESolver_LR: input nvirt is too large to cover by nbands, set nvirt = nbands - nocc = " << this->nvirt << std::endl;
    } else if (input.nvirt > 0) { this->nvirt = input.nvirt;
}
    this->nbands = this->nocc + this->nvirt;
    this->npairs = this->nocc * this->nvirt;
    this->nk = this->kv.get_nks() / this->nspin;
    if (this->nstates > this->nocc * this->nvirt * this->nk) {
        throw std::invalid_argument("ESolver_LR: nstates > nocc*nvirt*nk");
}

    GlobalV::ofs_running << "Setting LR-TDDFT parameters: " << std::endl;
    GlobalV::ofs_running << "number of occupied bands: " << this->nocc << std::endl;
    GlobalV::ofs_running << "number of virtual bands: " << this->nvirt << std::endl;
    GlobalV::ofs_running << "number of Atom orbitals (LCAO-basis size): " << this->nbasis << std::endl;
    GlobalV::ofs_running << "number of KS bands: " << this->eig_ks.nc << std::endl;
    GlobalV::ofs_running << "number of electron-hole pairs (2-particle basis size): " << this->npairs << std::endl;
    GlobalV::ofs_running << "number of excited states to be solved: " << this->nstates << std::endl;
    if (input.ri_hartree_benchmark == "aims" && !input.aims_nbasis.empty())
    {
        this->nbasis = [&]() -> int { int nbas = 0; for (int it = 0;it < ucell.ntype;++it) { nbas += ucell.atoms[it].na * input.aims_nbasis[it]; };return nbas;}();
        std::cout << "nbasis from aims: " << this->nbasis << std::endl;
    }
}

template <typename T, typename TR>
LR::ESolver_LR<T, TR>::ESolver_LR(ModuleESolver::ESolver_KS_LCAO<T, TR>&& ks_sol,
    const Input_para& inp, UnitCell& ucell)
    : input(inp), ucell(ucell)
#ifdef __EXX
    , exx_info(GlobalC::exx_info)
#endif
{
    redirect_log(inp.out_alllog);
    ModuleBase::TITLE("ESolver_LR", "ESolver_LR(KS)");

    if (this->input.lr_solver == "spectrum") {
        throw std::invalid_argument("when lr_solver==spectrum, esolver_type must be set to `lr` to skip the KS calculation.");
}

    // xc kernel
    this->xc_kernel = inp.xc_kernel;
    std::transform(xc_kernel.begin(), xc_kernel.end(), xc_kernel.begin(), tolower);
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
    if (this->nbands == PARAM.inp.nbands) { move_gs(); }
    else    // copy the part of ground state info according to paraC_
    {
        this->psi_ks = new psi::Psi<T>(this->kv.get_nks(), this->paraC_.get_col_size(), this->paraC_.get_row_size());
        this->eig_ks.create(this->kv.get_nks(), this->nbands);
        const int start_band = this->nocc_max - this->nocc;
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

    //grid integration
    this->gt_ = std::move(ks_sol.GridT);
    if (std::is_same<T, double>::value) { this->gint_g_ = std::move(ks_sol.GG); }
    else { this->gint_k_ = std::move(ks_sol.GK); }
    this->set_gint();

    // move pw basis
    delete this->pw_rho;    // newed in ESolver_FP::ESolver_FP
    this->pw_rho = ks_sol.pw_rho;
    ks_sol.pw_rho = nullptr;
    //init potential and calculate kernels using ground state charge
    init_pot(*ks_sol.pelec->charge);

#ifdef __EXX
    if (xc_kernel == "hf" || xc_kernel == "hse")
    {
        // if the same kernel is calculated in the esolver_ks, move it
        std::string dft_functional = input.dft_functional;
        std::transform(dft_functional.begin(), dft_functional.end(), dft_functional.begin(), tolower);
        if (ks_sol.exx_lri_double && std::is_same<T, double>::value && xc_kernel == dft_functional) {
            this->move_exx_lri(ks_sol.exx_lri_double);
        } else if (ks_sol.exx_lri_complex && std::is_same<T, std::complex<double>>::value && xc_kernel == dft_functional) {
            this->move_exx_lri(ks_sol.exx_lri_complex);
        } else    // construct C, V from scratch
        {
            // set ccp_type according to the xc_kernel
            if (xc_kernel == "hf") { exx_info.info_global.ccp_type = Conv_Coulomb_Pot_K::Ccp_Type::Hf; }
            else if (xc_kernel == "hse") { exx_info.info_global.ccp_type = Conv_Coulomb_Pot_K::Ccp_Type::Hse; }
            this->exx_lri = std::make_shared<Exx_LRI<T>>(exx_info.info_ri);
            this->exx_lri->init(MPI_COMM_WORLD, this->kv, ks_sol.orb_);
            this->exx_lri->cal_exx_ions(input.out_ri_cv);
        }
    }
#endif
    this->pelec = new elecstate::ElecStateLCAO<T>();
    orb_cutoff_ = ks_sol.orb_.cutoffs();
}

template <typename T, typename TR>
LR::ESolver_LR<T, TR>::ESolver_LR(const Input_para& inp, UnitCell& ucell) : input(inp), ucell(ucell)
#ifdef __EXX
, exx_info(GlobalC::exx_info)
#endif
{
    redirect_log(inp.out_alllog);
    ModuleBase::TITLE("ESolver_LR", "ESolver_LR(from scratch)");
    // xc kernel
    this->xc_kernel = inp.xc_kernel;
    std::transform(xc_kernel.begin(), xc_kernel.end(), xc_kernel.begin(), tolower);

    // necessary steps in ESolver_FP
    ESolver_FP::before_all_runners(inp, ucell);
    this->pelec = new elecstate::ElecStateLCAO<T>();

    // necessary steps in ESolver_KS::before_all_runners : symmetry and k-points
    if (ModuleSymmetry::Symmetry::symm_flag == 1)
    {
        GlobalC::ucell.symm.analy_sys(ucell.lat, ucell.st, ucell.atoms, GlobalV::ofs_running);
        ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "SYMMETRY");
    }
    this->kv.set(ucell.symm, PARAM.inp.kpoint_file, PARAM.inp.nspin, ucell.G, ucell.latvec, GlobalV::ofs_running);
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "INIT K-POINTS");
    ModuleIO::setup_parameters(ucell, this->kv);

    this->parameter_check();

    /// read orbitals and build the interpolation table
    two_center_bundle_.build_orb(ucell.ntype, ucell.orbital_fn);

    LCAO_Orbitals orb;
    two_center_bundle_.to_LCAO_Orbitals(orb, inp.lcao_ecut, inp.lcao_dk, inp.lcao_dr, inp.lcao_rmax);
    orb_cutoff_ = orb.cutoffs();

    this->set_dimension();
    //  setup 2d-block distribution for AO-matrix and KS wfc
    LR_Util::setup_2d_division(this->paraMat_, 1, this->nbasis, this->nbasis);
#ifdef __MPI
    this->paraMat_.set_desc_wfc_Eij(this->nbasis, this->nbands, paraMat_.get_row_size());
    int err = this->paraMat_.set_nloc_wfc_Eij(this->nbands, GlobalV::ofs_running, GlobalV::ofs_warning);
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
        this->paraMat_.get_row_size());
    this->read_ks_wfc();

    LR_Util::setup_2d_division(this->paraC_, 1, this->nbasis, this->nbands
#ifdef __MPI
        , paraMat_.blacs_ctxt
#endif
    );

    //allocate 2-particle state and setup 2d division
    this->pelec = new elecstate::ElecState();

    // read the ground state charge density and calculate xc kernel
    GlobalC::Pgrid.init(this->pw_rho->nx,
        this->pw_rho->ny,
        this->pw_rho->nz,
        this->pw_rho->nplane,
        this->pw_rho->nrxx,
        pw_big->nbz,
        pw_big->bz);
    Charge chg_gs;
    if (input.ri_hartree_benchmark != "aims") { this->read_ks_chg(chg_gs); }
    this->init_pot(chg_gs);

    // search adjacent atoms and init Gint
    std::cout << "ucell.infoNL.get_rcutmax_Beta(): " << GlobalC::ucell.infoNL.get_rcutmax_Beta() << std::endl;
    double search_radius = -1.0;
    search_radius = atom_arrange::set_sr_NL(GlobalV::ofs_running,
        PARAM.inp.out_level,
        orb.get_rcutmax_Phi(),
        GlobalC::ucell.infoNL.get_rcutmax_Beta(),
        PARAM.globalv.gamma_only_local);
    atom_arrange::search(PARAM.inp.search_pbc,
        GlobalV::ofs_running,
        GlobalC::GridD,
        this->ucell,
        search_radius,
        PARAM.inp.test_atom_input);
    this->set_gint();
    this->gint_->gridt = &this->gt_;
    this->gint_->reset_DMRGint(1);

    // (3) Periodic condition search for each grid.
    double dr_uniform = 0.001;
    std::vector<double> rcuts;
    std::vector<std::vector<double>> psi_u;
    std::vector<std::vector<double>> dpsi_u;
    std::vector<std::vector<double>> d2psi_u;

    Gint_Tools::init_orb(dr_uniform, rcuts, GlobalC::ucell, orb, psi_u, dpsi_u, d2psi_u);
    this->gt_.set_pbc_grid(this->pw_rho->nx,
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
        PARAM.inp.nstream);
    psi_u.clear();
    psi_u.shrink_to_fit();
    dpsi_u.clear();
    dpsi_u.shrink_to_fit();
    d2psi_u.clear();
    d2psi_u.shrink_to_fit();

    if (std::is_same<T, std::complex<double>>::value)
    {
        this->gt_.cal_nnrg(&this->paraMat_, orb.cutoffs());
        this->gint_k_.allocate_pvpR();   // uses gt_.nnrg
    }
    this->gint_->prep_grid(this->gt_,
        this->pw_big->nbx,
        this->pw_big->nby,
        this->pw_big->nbzp,
        this->pw_big->nbzp_start,
        this->pw_rho->nxyz,
        this->pw_big->bx,
        this->pw_big->by,
        this->pw_big->bz,
        this->pw_big->bxyz,
        this->pw_big->nbxx,
        this->pw_rho->ny,
        this->pw_rho->nplane,
        this->pw_rho->startz_current,
        &ucell,
        &orb);
    this->gint_->initialize_pvpR(ucell, &GlobalC::GridD, 1);    // always use nspin=1 for transition density

    // if EXX from scratch, init 2-center integral and calculate Cs, Vs 
#ifdef __EXX
    if ((xc_kernel == "hf" || xc_kernel == "hse") && this->input.lr_solver != "spectrum")
    {
        // set ccp_type according to the xc_kernel
        if (xc_kernel == "hf") { exx_info.info_global.ccp_type = Conv_Coulomb_Pot_K::Ccp_Type::Hf; }
        else if (xc_kernel == "hse") { exx_info.info_global.ccp_type = Conv_Coulomb_Pot_K::Ccp_Type::Hse; }
        this->exx_lri = std::make_shared<Exx_LRI<T>>(exx_info.info_ri);
        this->exx_lri->init(MPI_COMM_WORLD, this->kv, orb);
        this->exx_lri->cal_exx_ions(input.out_ri_cv);
    }
    // else
#endif
        // ModuleBase::Ylm::set_coefficients() is deprecated
}

template <typename T, typename TR>
void LR::ESolver_LR<T, TR>::runner(int istep, UnitCell& cell)
{
    ModuleBase::TITLE("ESolver_LR", "runner");
    //allocate 2-particle state and setup 2d division
    this->setup_eigenvectors_X();
    this->pelec->ekb.create(nspin, this->nstates);

    if (this->input.lr_solver != "spectrum")
    {
        // allocate and initialize A matrix and density matrix
        std::vector<std::string> spin_type = { "Spin Singlet", "Spin Triplet" };
        for (int is = 0;is < nspin;++is)
        {
            if (nspin == 2) { std::cout << "Calculating " << spin_type[is] << " excitations" << std::endl; }
            hamilt::Hamilt<T>* phamilt = new HamiltCasidaLR<T>(xc_kernel, nspin, this->nbasis, this->nocc, this->nvirt, this->ucell, orb_cutoff_, GlobalC::GridD, this->psi_ks, this->eig_ks,
#ifdef __EXX
                this->exx_lri, this->exx_info.info_global.hybrid_alpha,
#endif
                this->gint_, this->pot[is], this->kv, & this->paraX_, & this->paraC_, & this->paraMat_, 
                spin_type[is], input.ri_hartree_benchmark, (input.ri_hartree_benchmark == "aims" ? input.aims_nbasis : std::vector<int>({})));
            // solve the Casida equation
            HSolverLR<T> hsol(nk, this->npairs, is, this->input.out_wfc_lr);
            hsol.set_diagethr(hsol.diag_ethr, 0, 0, std::max(1e-13, this->input.lr_thr));
            hsol.solve(phamilt, *this->X[is], this->pelec, this->input.lr_solver/*,
                !std::set<std::string>({ "hf", "hse" }).count(this->xc_kernel)*/);  //whether the kernel is Hermitian
            delete phamilt;
        }
    }
    else    // read the eigenvalues
    {
        std::ifstream ifs(PARAM.globalv.global_readin_dir + "Excitation_Energy.dat");
        std::cout << "reading the excitation energies from file: \n";
        for (int is = 0;is < nspin;++is)
        {
            for (int i = 0;i < this->nstates;++i) { ifs >> this->pelec->ekb(is, i); }
            for (int i = 0;i < this->nstates;++i) { std::cout << this->pelec->ekb(is, i) << " "; }
        }
    }
    return;
}

template <typename T, typename TR>
void LR::ESolver_LR<T, TR>::after_all_runners()
{
    ModuleBase::TITLE("ESolver_LR", "after_all_runners");
    if (input.ri_hartree_benchmark != "none") { return; } //no need to calculate the spectrum in the benchmark routine
    //cal spectrum
    std::vector<double> freq(100);
    std::vector<double> abs_wavelen_range({ 20, 200 });//default range
    if (input.abs_wavelen_range.size() == 2 && std::abs(input.abs_wavelen_range[1] - input.abs_wavelen_range[0]) > 0.02)
    {
        abs_wavelen_range = input.abs_wavelen_range;
    }
    double lambda_diff = std::abs(abs_wavelen_range[1] - abs_wavelen_range[0]);
    double lambda_min = std::min(abs_wavelen_range[1], abs_wavelen_range[0]);
    for (int i = 0;i < freq.size();++i) { freq[i] = 91.126664 / (lambda_min + 0.01 * static_cast<double>(i + 1) * lambda_diff); }
    for (int is = 0;is < this->nspin;++is)
    {
        LR_Spectrum<T> spectrum(&this->pelec->ekb.c[is * this->nstates], *this->X[is],
            this->nspin, this->nbasis, this->nocc, this->nvirt, this->gint_, *this->pw_rho, *this->psi_ks,
            this->ucell, this->kv, this->paraX_, this->paraC_, this->paraMat_);
        spectrum.oscillator_strength();
        spectrum.transition_analysis(is);
        spectrum.optical_absorption(freq, input.abs_broadening, is);
    }
}


template<typename T, typename TR>
void LR::ESolver_LR<T, TR>::setup_eigenvectors_X()
{
    ModuleBase::TITLE("ESolver_LR", "setup_eigenvectors_X");
    LR_Util::setup_2d_division(this->paraX_, 1, this->nvirt, this->nocc
#ifdef __MPI
        , this->paraC_.blacs_ctxt
#endif
    );//nvirt - row, nocc - col 
    this->X.resize(this->nspin);
    const std::vector<std::string> spin_types = { "Spin Singlet", "Spin Triplet" };
    // if spectrum-only, read the LR-eigenstates from file and return
    if (this->input.lr_solver == "spectrum")
    {
        std::cout << "reading the excitation amplitudes from file: \n";
        for (int is = 0; is < this->nspin; ++is)
        {
            this->X[is] = std::make_shared<psi::Psi<T>>(LR_Util::read_psi_bandfirst<T>(
                PARAM.globalv.global_readin_dir + "Excitation_Amplitude_" + spin_types[is], GlobalV::MY_RANK));
        }
    }
    else
    {
        for (int is = 0; is < this->nspin; ++is)
        {
            this->X[is] = std::make_shared<psi::Psi<T>>(this->nk,
                this->nstates,
                this->paraX_.get_local_size(),
                nullptr,
                false); // band(state)-first
            this->X[is]->zero_out();
        }
        set_X_initial_guess();
    }
}

template<typename T, typename TR>
void LR::ESolver_LR<T, TR>::set_X_initial_guess()
{
    // set the initial guess of X
  // if (E_{lumo}-E_{homo-1} < E_{lumo+1}-E{homo}), mode = 0, else 1(smaller first)
    bool ix_mode = false;   //default
    if (this->eig_ks.nc > nocc + 1 && nocc >= 2 && eig_ks(0, nocc) - eig_ks(0, nocc - 2) - 1e-5 > eig_ks(0, nocc + 1) - eig_ks(0, nocc - 1)) { ix_mode = true; }
    GlobalV::ofs_running << "setting the initial guess of X: " << std::endl;
    if (nocc >= 2 && eig_ks.nc > nocc) { GlobalV::ofs_running << "E_{lumo}-E_{homo-1}=" << eig_ks(0, nocc) - eig_ks(0, nocc - 2) << std::endl; }
    if (nocc >= 1 && eig_ks.nc > nocc + 1) { GlobalV::ofs_running << "E_{lumo+1}-E{homo}=" << eig_ks(0, nocc + 1) - eig_ks(0, nocc - 1) << std::endl; }
    GlobalV::ofs_running << "mode of X-index: " << ix_mode << std::endl;

    /// global index map between (i,c) and ix
    ModuleBase::matrix ioiv2ix;
    std::vector<std::pair<int, int>> ix2ioiv;
    std::pair<ModuleBase::matrix, std::vector<std::pair<int, int>>> indexmap =
        LR_Util::set_ix_map_diagonal(ix_mode, nocc, nvirt);

    ioiv2ix = std::move(std::get<0>(indexmap));
    ix2ioiv = std::move(std::get<1>(indexmap));

    // use unit vectors as the initial guess
    // for (int i = 0; i < std::min(this->nstates * PARAM.inp.pw_diag_ndim, nocc * nvirt); i++)
    for (int is = 0;is < this->nspin;++is)
    {
        for (int s = 0; s < nstates; ++s)
        {
            this->X[is]->fix_b(s);
            int ipair = s % (npairs);
            int occ_global = std::get<0>(ix2ioiv[ipair]);   // occ
            int virt_global = std::get<1>(ix2ioiv[ipair]);   // virt
            int ik = s / (npairs);
            if (this->paraX_.in_this_processor(virt_global, occ_global))
                (*X[is])(ik, this->paraX_.global2local_col(occ_global) * this->paraX_.get_row_size()
                    + this->paraX_.global2local_row(virt_global))
                = (static_cast<T>(1.0) / static_cast<T>(nk));
        }
        this->X[is]->fix_b(0);  //recover the pointer
    }
}

template<typename T, typename TR>
void LR::ESolver_LR<T, TR>::init_pot(const Charge& chg_gs)
{
    this->pot.resize(nspin, nullptr);
    if (this->input.ri_hartree_benchmark != "none") { return; } //no need to initialize potential for Hxc kernel in the RI-benchmark routine
    switch (nspin)
    {
    case 1:
        this->pot[0] = std::make_shared<PotHxcLR>(xc_kernel, this->pw_rho, &ucell, &chg_gs, PotHxcLR::SpinType::S1);
        break;
    case 2:
        this->pot[0] = std::make_shared<PotHxcLR>(xc_kernel, this->pw_rho, &ucell, &chg_gs, PotHxcLR::SpinType::S2_singlet);
        this->pot[1] = std::make_shared<PotHxcLR>(xc_kernel, this->pw_rho, &ucell, &chg_gs, PotHxcLR::SpinType::S2_triplet);
        break;
    default:
        throw std::invalid_argument("ESolver_LR: nspin must be 1 or 2");
    }
}

template<typename T, typename TR>
void LR::ESolver_LR<T, TR>::read_ks_wfc()
{
    assert(this->psi_ks != nullptr);
    this->pelec->ekb.create(this->kv.get_nks(), this->nbands);
    this->pelec->wg.create(this->kv.get_nks(), this->nbands);

    if (input.ri_hartree_benchmark == "aims")        // for aims benchmark
    {
#ifdef __EXX
        int ncore = 0;
        std::vector<double> eig_ks_vec = RI_Benchmark::read_aims_ebands<double>(PARAM.globalv.global_readin_dir + "band_out", nocc, nvirt, ncore);
        std::cout << "ncore=" << ncore << ", nocc=" << nocc << ", nvirt=" << nvirt << ", nbands=" << this->nbands << std::endl;
        std::cout << "eig_ks_vec.size()=" << eig_ks_vec.size() << std::endl;
        if(eig_ks_vec.size() != this->nbands) {ModuleBase::WARNING_QUIT("ESolver_LR", "read_aims_ebands failed.");};
        for (int i = 0;i < nbands;++i) { this->pelec->ekb(0, i) = eig_ks_vec[i]; }
        RI_Benchmark::read_aims_eigenvectors<T>(*this->psi_ks, PARAM.globalv.global_readin_dir + "KS_eigenvectors.out", ncore, nbands, nbasis);
#else
        ModuleBase::WARNING_QUIT("ESolver_LR", "RI benchmark is only supported when compile with LibRI.");
#endif
    }
    else if (!ModuleIO::read_wfc_nao(PARAM.globalv.global_readin_dir, this->paraMat_, *this->psi_ks, this->pelec,
        /*skip_bands=*/this->nocc_max - this->nocc)) {
        ModuleBase::WARNING_QUIT("ESolver_LR", "read ground-state wavefunction failed.");
    }
    this->eig_ks = std::move(this->pelec->ekb);
}

template<typename T, typename TR>
void LR::ESolver_LR<T, TR>::read_ks_chg(Charge& chg_gs)
{
    chg_gs.set_rhopw(this->pw_rho);
    chg_gs.allocate(this->nspin);
    GlobalV::ofs_running << " try to read charge from file : ";
    for (int is = 0; is < this->nspin; ++is)
    {
        std::stringstream ssc;
        ssc << PARAM.globalv.global_readin_dir << "SPIN" << is + 1 << "_CHG.cube";
        GlobalV::ofs_running << ssc.str() << std::endl;
        double ef;
        if (ModuleIO::read_cube(
#ifdef __MPI
            & (GlobalC::Pgrid),
#endif
            GlobalV::MY_RANK,
            is,
            GlobalV::ofs_running,
            this->nspin,
            ssc.str(),
            chg_gs.rho[is],
            this->pw_rho->nx,
            this->pw_rho->ny,
            this->pw_rho->nz,
            ef,
            &(GlobalC::ucell),
            chg_gs.prenspin)) {
            GlobalV::ofs_running << " Read in the charge density: " << ssc.str() << std::endl;
        } else {    // prenspin for nspin=4 is not supported currently
            ModuleBase::WARNING_QUIT(
                "init_rho",
                "!!! Couldn't find the charge file !!! The default directory \n of SPIN1_CHG.cube is OUT.suffix, "
                "or you must set read_file_dir \n to a specific directory. ");
}
    }
}
template class LR::ESolver_LR<double, double>;
template class LR::ESolver_LR<std::complex<double>, double>;
