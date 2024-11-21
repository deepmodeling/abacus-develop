#include "esolver_lrtd_lcao.h"
#include "utils/gint_move.hpp"
#include "utils/lr_util.h"
#include "hamilt_casida.h"
#include "hamilt_ulr.hpp"
#include "module_lr/potentials/pot_hxc_lrtd.h"
#include "module_lr/hsolver_lrtd.hpp"
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
#include "module_lr/operator_casida/operator_lr_diag.h" // for precondition

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

inline int cal_nupdown_form_occ(const ModuleBase::matrix& wg)
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
void LR::ESolver_LR<T, TR>::parameter_check()const
{
    const std::set<std::string> lr_solvers = { "dav", "lapack" , "spectrum", "dav_subspace", "cg" };
    const std::set<std::string> xc_kernels = { "rpa", "lda", "pwlda", "pbe", "hf" , "hse" };
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
    this->nstates = input.lr_nstates;
    this->nbasis = PARAM.globalv.nlocal;
    // calculate the number of occupied and unoccupied states
    // which determines the basis size of the excited states
    this->nocc_max = LR_Util::cal_nocc(LR_Util::cal_nelec(ucell));
    this->nocc_in = std::max(1, std::min(input.nocc, this->nocc_max));
    this->nvirt_in = PARAM.inp.nbands - this->nocc_max;   //nbands-nocc
    if (input.nvirt > this->nvirt_in) { GlobalV::ofs_warning << "ESolver_LR: input nvirt is too large to cover by nbands, set nvirt = nbands - nocc = " << this->nvirt_in << std::endl; }
    else if (input.nvirt > 0) { this->nvirt_in = input.nvirt; }
    this->nbands = this->nocc_in + this->nvirt_in;
    this->nk = this->kv.get_nks() / this->nspin;
    this->nocc.resize(nspin, nocc_in);
    this->nvirt.resize(nspin, nvirt_in);
    for (int is = 0;is < nspin;++is) { this->npairs.push_back(nocc[is] * nvirt[is]); }
    GlobalV::ofs_running << "Setting LR-TDDFT parameters: " << std::endl;
    GlobalV::ofs_running << "number of occupied bands: " << nocc_in << std::endl;
    GlobalV::ofs_running << "number of virtual bands: " << nvirt_in << std::endl;
    GlobalV::ofs_running << "number of Atom orbitals (LCAO-basis size): " << this->nbasis << std::endl;
    GlobalV::ofs_running << "number of KS bands: " << this->eig_ks.nc << std::endl;
    GlobalV::ofs_running << "number of excited states to be solved: " << this->nstates << std::endl;
    if (input.ri_hartree_benchmark == "aims" && !input.aims_nbasis.empty())
    {
        // calculate total number of basis funcs, see https://en.cppreference.com/w/cpp/algorithm/inner_product
        this->nbasis = std::inner_product(input.aims_nbasis.begin(), /* iterator1.begin */
                                          input.aims_nbasis.end(),  /* iterator1.end */
                                          ucell.atoms,  /* iterator2.begin */
                                          0,  /* init value */
                                          std::plus<int>(), /* iter op1 */
                                          [](const int& a, const Atom& b) { return a * b.na; }); /* iter op2 */
        std::cout << "nbasis from aims: " << this->nbasis << std::endl;
    }
}

template<typename T, typename TR>
void LR::ESolver_LR<T, TR>::reset_dim_spin2()
{
    if (nspin != 2) { return; }
    if (nupdown == 0) { std::cout << "** Assuming the spin-up and spin-down states are degenerate. **" << std::endl; }
    else
    {
        this->openshell = true;
        nupdown > 0 ? ((nocc[1] -= nupdown) && (nvirt[1] += nupdown)) : ((nocc[0] += nupdown) && (nvirt[0] -= nupdown));
        npairs = { nocc[0] * nvirt[0], nocc[1] * nvirt[1] };
        std::cout << "** Solve the spin-up and spin-down states separately for open-shell system. **" << std::endl;
    }
    for (int is : {0, 1}) { if (npairs[is] <= 0) { throw std::invalid_argument(std::string("ESolver_LR: npairs (nocc*nvirt) <= 0 for spin") + std::string(is == 0 ? "up" : "down")); } }
    if (nstates > (npairs[0] + npairs[1]) * nk) { throw std::invalid_argument("ESolver_LR: nstates > nocc*nvirt*nk"); }
    if (input.lr_unrestricted) { this->openshell = true; }
}

template <typename T, typename TR>
LR::ESolver_LR<T, TR>::ESolver_LR(ModuleESolver::ESolver_KS_LCAO<T, TR>&& ks_sol,
    const Input_para& inp, UnitCell& ucell)
    : input(inp), ucell(ucell)
#ifdef __EXX
    , exx_info(GlobalC::exx_info)
#endif
{
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
        const int start_band = this->nocc_max - std::max(nocc[0], nocc[1]);
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

    //grid integration
    this->gt_ = std::move(ks_sol.GridT);
    if (std::is_same<T, double>::value) { this->gint_g_ = std::move(ks_sol.GG); }
    else { this->gint_k_ = std::move(ks_sol.GK); }
    this->set_gint();
    this->gint_->reset_DMRGint(1);

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
        GlobalC::GridD,
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
    ModuleBase::timer::tick("ESolver_LR", "runner");
    //allocate 2-particle state and setup 2d division
    this->setup_eigenvectors_X();
    this->pelec->ekb.create(nspin, this->nstates);

    auto efile = [&](const std::string& label)->std::string {return PARAM.globalv.global_out_dir + "Excitation_Energy_" + label + ".dat";};
    auto vfile = [&](const std::string& label)->std::string {return PARAM.globalv.global_out_dir + "Excitation_Amplitude_" + label + "_" + std::to_string(GlobalV::MY_RANK) + ".dat";};
    if (this->input.lr_solver != "spectrum")
    {
        auto write_states = [&](const std::string& label, const Real<T>* e, const T* v, const int& dim, const int& nst, const int& prec = 8)->void
            {
                if (GlobalV::MY_RANK == 0) { assert(nst == LR_Util::write_value(efile(label), prec, e, nst)); }
                assert(nst * dim == LR_Util::write_value(vfile(label), prec, v, nst, dim));
            };
        std::vector<double> precondition(this->input.lr_solver == "lapack" ? 0 : nloc_per_band, 1.0);
        // allocate and initialize A matrix and density matrix
        if (openshell)
        {
            for (int is : {0, 1})
            {
                const int offset_is = is * this->paraX_[0].get_local_size();
                OperatorLRDiag<double> pre_op(this->eig_ks.c + is * nk * (nocc[0] + nvirt[0]), this->paraX_[is], this->nk, this->nocc[is], this->nvirt[is]);
                if (input.lr_solver != "lapack") { pre_op.act(1, offset_is, 1, precondition.data() + offset_is, precondition.data() + offset_is); }
            }
            std::cout << "Solving spin-conserving excitation for open-shell system." << std::endl;
            HamiltULR<T> hulr(xc_kernel, nspin, this->nbasis, this->nocc, this->nvirt, this->ucell, orb_cutoff_, GlobalC::GridD, *this->psi_ks, this->eig_ks,
#ifdef __EXX
                this->exx_lri, this->exx_info.info_global.hybrid_alpha,
#endif
                this->gint_, this->pot, this->kv, this->paraX_, this->paraC_, this->paraMat_);
            LR::HSolver::solve(hulr, this->X[0].template data<T>(), nloc_per_band, nstates, this->pelec->ekb.c, this->input.lr_solver, this->input.lr_thr, precondition);
            if (input.out_wfc_lr) { write_states("openshell", this->pelec->ekb.c, this->X[0].template data<T>(), nloc_per_band, nstates); }
        }
        else
        {
            OperatorLRDiag<double> pre_op(this->eig_ks.c, this->paraX_[0], this->nk, this->nocc[0], this->nvirt[0]);
            if (input.lr_solver != "lapack") { pre_op.act(1, nloc_per_band, 1, precondition.data(), precondition.data()); }
            auto spin_types = std::vector<std::string>({ "singlet", "triplet" });
            for (int is = 0;is < nspin;++is)
            {
                std::cout << "Calculating " << spin_types[is] << " excitations" << std::endl;
                HamiltLR<T> hlr(xc_kernel, nspin, this->nbasis, this->nocc, this->nvirt, this->ucell, orb_cutoff_, GlobalC::GridD, *this->psi_ks, this->eig_ks,
#ifdef __EXX
                    this->exx_lri, this->exx_info.info_global.hybrid_alpha,
#endif
                    this->gint_, this->pot[is], this->kv, this->paraX_, this->paraC_, this->paraMat_,
                    spin_types[is], input.ri_hartree_benchmark, (input.ri_hartree_benchmark == "aims" ? input.aims_nbasis : std::vector<int>({})));
                // solve the Casida equation
                LR::HSolver::solve(hlr, this->X[is].template data<T>(), nloc_per_band, nstates,
                    this->pelec->ekb.c + is * nstates, this->input.lr_solver, this->input.lr_thr, precondition/*,
                        !std::set<std::string>({ "hf", "hse" }).count(this->xc_kernel)*/);  //whether the kernel is Hermitian
                if (input.out_wfc_lr) { write_states(spin_types[is], this->pelec->ekb.c + is * nstates, this->X[is].template data<T>(), nloc_per_band, nstates); }
            }
        }
    }
    else    // read the eigenvalues
    {
        auto read_states = [&](const std::string& label, Real<T>* e, T* v, const int& dim, const int& nst)->void
            {
                if (GlobalV::MY_RANK == 0) { assert(nst == LR_Util::read_value(efile(label), e, nst)); }
                assert(nst * dim == LR_Util::read_value(vfile(label), v, nst, dim));
            };
        std::cout << "reading the excitation amplitudes from file: \n";
        if (openshell)
        {
            read_states("openshell", this->pelec->ekb.c, this->X[0].template data<T>(), nloc_per_band, nstates);
        }
        else
        {
            auto spin_types = std::vector<std::string>({ "singlet", "triplet" });
            for (int is = 0;is < nspin;++is) { read_states(spin_types[is], this->pelec->ekb.c + is * nstates, this->X[is].template data<T>(), nloc_per_band, nstates); }
        }
    }
    ModuleBase::timer::tick("ESolver_LR", "runner");
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
        LR_Spectrum<T> spectrum(nspin, this->nbasis, this->nocc, this->nvirt, this->gint_, *this->pw_rho, *this->psi_ks,
            this->ucell, this->kv, GlobalC::GridD, this->orb_cutoff_,
            this->paraX_, this->paraC_, this->paraMat_,
            &this->pelec->ekb.c[is * nstates], this->X[is].template data<T>(), nstates, openshell);
        spectrum.transition_analysis(spin_types[is]);
        spectrum.optical_absorption(freq, input.abs_broadening, spin_types[is]);
    }
}

template<typename T, typename TR>
void LR::ESolver_LR<T, TR>::setup_eigenvectors_X()
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
    this->nloc_per_band = nk * (openshell ? paraX_[0].get_local_size() + paraX_[1].get_local_size() : paraX_[0].get_local_size());

    this->X.resize(openshell ? 1 : nspin, LR_Util::newTensor<T>({ nstates, nloc_per_band }));
    for (auto& x : X) { x.zero(); }

    auto spin_types = (nspin == 2 && !openshell) ? std::vector<std::string>({ "singlet", "triplet" }) : std::vector<std::string>({ "updown" });
    // if spectrum-only, read the LR-eigenstates from file and return
    if (this->input.lr_solver != "spectrum") { set_X_initial_guess(); }
}

template<typename T, typename TR>
void LR::ESolver_LR<T, TR>::set_X_initial_guess()
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
            const int xstart_b = ib * nloc_per_band;    //start index of band ib
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
void LR::ESolver_LR<T, TR>::init_pot(const Charge& chg_gs)
{
    this->pot.resize(nspin, nullptr);
    if (this->input.ri_hartree_benchmark != "none") { return; } //no need to initialize potential for Hxc kernel in the RI-benchmark routine
    switch (nspin)
    {
        using ST = PotHxcLR::SpinType;
    case 1:
        this->pot[0] = std::make_shared<PotHxcLR>(xc_kernel, *this->pw_rho, ucell, chg_gs, GlobalC::Pgrid, ST::S1, input.lr_init_xc_kernel);
        break;
    case 2:
        this->pot[0] = std::make_shared<PotHxcLR>(xc_kernel, *this->pw_rho, ucell, chg_gs, GlobalC::Pgrid, openshell ? ST::S2_updown : ST::S2_singlet, input.lr_init_xc_kernel);
        this->pot[1] = std::make_shared<PotHxcLR>(xc_kernel, *this->pw_rho, ucell, chg_gs, GlobalC::Pgrid, openshell ? ST::S2_updown : ST::S2_triplet, input.lr_init_xc_kernel);
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
    else if (!ModuleIO::read_wfc_nao(PARAM.globalv.global_readin_dir, this->paraMat_, *this->psi_ks, this->pelec,
        /*skip_bands=*/this->nocc_max - this->nocc_in)) {
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
        if (ModuleIO::read_vdata_palgrid(GlobalC::Pgrid,
            GlobalV::MY_RANK,
            GlobalV::ofs_running,
            ssc.str(),
            chg_gs.rho[is],
            ucell.nat)) {
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
