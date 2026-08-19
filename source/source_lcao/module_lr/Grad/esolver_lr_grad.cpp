#include "source_esolver/esolver_lr_lcao_tddft.h"
#include "source_lcao/module_lr/Grad/multipliers/zeq_solver.h"
#include "source_lcao/module_lr/Grad/multipliers/cal_edm_from_multipliers.h"
#include "source_lcao/module_lr/Grad/force/lr_force.h"
#include "source_estate/module_dm/cal_dm_psi.h"
#include "source_io/module_output/output_log.h"

using namespace LR;

template <typename Tstream>
inline void print_force(const std::vector<ModuleBase::matrix>& force, Tstream& ofs)
{
    const int nstate = force.size();
    ofs << "Gradients of each excited state: (eV/Angstrom)" << std::endl;
    ofs << std::setprecision(6) << std::setw(6) << "state" << std::setw(6) << "atom"
        << std::setw(15) << "x" << std::setw(15) << "y" << std::setw(15) << "z" << std::endl;
    const double fac = ModuleBase::Ry_to_eV / ModuleBase::BOHR_TO_A;
    for (int i = 0;i < nstate;++i)
    {
        for (int iat = 0;iat < force[i].nr;++iat)
        {
            std::string istate = iat == 0 ? std::to_string(i) : " ";
            ofs << std::setw(6) << istate << std::setw(6) << iat << std::setw(6) << "force";
            for (int ixyz = 0;ixyz < 3;++ixyz) { ofs << std::setw(15) << force[i](iat, ixyz) * fac; }
            ofs << std::endl;
        }
    }
}

// check C_uaC_va-C_uiC_vi of lumo-homo, nocc=1,  nk=1
template<typename T>
inline void test_dm_diff_H2(const T* dm, const psi::Psi<T>& c, const int nbasis)
{
    std::cout << "difference dm cal: " << std::endl;
    LR_Util::print_value(dm, nbasis, nbasis);
    std::cout << "difference dm ref: " << std::endl;
    for (int i = 0;i < nbasis;++i)
    {
        for (int j = 0;j < nbasis;++j)
        {
            std::cout << c(0, 1, i) * c(0, 1, j) - c(0, 0, i) * c(0, 0, j) << " ";
        }
        std::cout << std::endl;
    }
}

// check e_aC_uaC_va-e_iC_uiC_vi of lumo-homo, nocc=1,  nk=1
template<typename T>
inline void test_edm_H2(const T* const edm, const double* const eig_ks, const psi::Psi<T>& c, const int nbasis)
{
    std::cout << "edm cal: " << std::endl;
    LR_Util::print_value(edm, nbasis, nbasis);
    std::cout << "edm ref: " << std::endl;
    for (int i = 0;i < nbasis;++i)
    {
        for (int j = 0;j < nbasis;++j)
        {
            std::cout << eig_ks[1] * c(0, 1, i) * c(0, 1, j) - eig_ks[0] * c(0, 0, i) * c(0, 0, j) << " ";
        }
        std::cout << std::endl;
    }
}

template<typename T, typename TR>
void ModuleESolver::ESolver_LR<T, TR>::init_pot_groundstate(const Charge& chg_gs)
{
    ModuleBase::TITLE("ESolver_LR", "init_pot_gs");
    std::vector<std::string> pot_register;
    if (PARAM.inp.vl_in_h)
    {
        //! 11) calculate the structure factor
        this->sf.setup(&(*this->ucell_), Pgrid, this->pw_rhod);
        this->locpp.init_vloc((*this->ucell_), this->pw_rho);
        pot_register.push_back("local");
    }
    if(PARAM.inp.vh_in_h)
    {
        pot_register.push_back("hartree");
    }
    pot_register.push_back("xc");
    
    // initialize the ground state potential
    this->pot_gs = LR_Util::make_unique<elecstate::Potential>(this->pw_rhod, this->pw_rho,
        &(*this->ucell_), &this->locpp.vloc, &this->sf, &this->solvent,
       &this->etxc_gs, &this->vtxc_gs);
    this->pot_gs.get()->pot_register(pot_register);
    XC_Functional::set_xc_type((*this->ucell_).atoms[0].ncpp.xc_func);    // set XC type of the ground state
    this->pot_gs->init_pot(&chg_gs);  // call update_from_charge inside
    if (LR_Util::has_local_xc(this->xc_kernel))
    {
        XC_Functional::set_xc_type(this->xc_kernel);    // recover the excited state xc kernel type
    }
    if (PARAM.inp.test_force)
    {
        this->pot_gs_hartree = LR_Util::make_unique<elecstate::Potential>(this->pw_rhod, this->pw_rho,
            &(*this->ucell_), &this->locpp.vloc, &this->sf, &this->solvent,
            &this->etxc_gs, &this->vtxc_gs);
        this->pot_gs_hartree->pot_register({ "hartree" });
        this->pot_gs_hartree->init_pot(&chg_gs);  // call update_from_charge inside
    }
}

template<typename T, typename TR>
ct::Tensor ModuleESolver::ESolver_LR<T, TR>::solve_zvector_eqation(const int ispin)
{
    ModuleBase::TITLE("ESolver_LR", "cal_force");
    ModuleBase::timer::start("ESolver_LR", "solve_zvector_eqation");
    ct::Tensor Z = LR_Util::newTensor<T>({ this->nstates, this->nloc_per_state });
    // construct and solve the Z-vector equation
    Z_vector_equation(this->X[ispin].template data<T>(), Z.template data<T>(),
        this->xc_kernel, this->nstates, this->nspin, this->nbasis, this->nocc, this->nvirt,
        (*this->ucell_), orb_cutoff_, this->gd, *this->psi_ks, this->eig_ks,
#ifdef __EXX    
        std::weak_ptr<Exx_LRI<T>>(this->exx_lri), this->exx_info.info_global.hybrid_alpha,
#endif
        std::weak_ptr<PotHxcLR>(this->pot[ispin]), std::weak_ptr<PotHxcLR>(this->pot_hxc_gs),
        this->kv, this->paraX_, this->paraC_, this->paraMat_, this->spin_types[ispin]);
    ModuleBase::timer::end("ESolver_LR", "solve_zvector_eqation");
    return Z;
}

template<typename T, typename TR>
std::vector<ModuleBase::matrix> ModuleESolver::ESolver_LR<T, TR>::cal_force(const int ispin)
{
    if (PARAM.inp.test_force && ispin == 0) { this->test_force(); }

    const ct::Tensor& Z = this->solve_zvector_eqation(ispin);

    ModuleBase::TITLE("ESolver_LR", "cal_force");
    ModuleBase::timer::start("ESolver_LR", "cal_force");
    const auto& c = LR_Util::get_psi_spin(*this->psi_ks, ispin, this->nk);   // wavefunction coefficients of ground state

    // calculate the force (the partial gradient of Lagrangian)
    LR_Force<T> lr_force((*this->ucell_), this->kv.kvec_d, this->paraMat_,
        *this->pw_rhod, *this->pw_rho, this->locpp, this->sf, this->gd, this->two_center_bundle_
#ifdef __EXX
        , std::weak_ptr<Exx_LRI<T>>(this->exx_lri), this->exx_info.info_global.hybrid_alpha
#endif
    );
    GlobalV::ofs_running << "Start to calculate excited-state force of " << this->spin_types[ispin] << std::endl;
    // ground state dm for currrent spin (only for test the correctness of the force)
    // elecstate::DensityMatrix<T, T> dm_gs(this->paraMat_, 1, this->kv.kvec_d, this->nk);

    // for each state, calculate dm_trans, dm_relaxed_diff, edm and force
    std::vector<ModuleBase::matrix> forces(this->nstates);
    for (int istate = 0;istate < this->nstates;++istate)
    {
        const int offset = istate * this->nloc_per_state;
        // The imag part will be cancelled in the force calculation, so we use double DM(R) to calculate force. 
        // But complex transition DM(R) is still used in energy density matrix calculation.
        const auto& dm_trans_k = cal_dm_trans_pblas(this->X[ispin].template data<T>() + offset, this->paraX_[ispin], c, this->paraC_, this->nbasis, this->nocc[ispin], this->nvirt[ispin], this->paraMat_);
        // D(X) complex, for the EXX (LibRI) force. Built FIRST and left UN-symmetrized:
        // the exchange kernel (mu kappa | nu lambda) puts the two indices of one D^X into
        // different electron coordinates, so Tr[D^X D^X K_exx] = (aa|ii) requires the full
        // non-symmetric D^X. Symmetrizing would give 1/2[(aa|ii)+(ai|ia)], which is wrong.
        // (`cal_force_exx_dm_trans` feeds the same tensor to both slots, so it is consistent.)
        auto dm_trans =   // D(X) complex
            LR_Util::build_dm_from_dmk<T, T>(dm_trans_k,
                this->paraMat_, this->nk, this->kv.kvec_d, (*this->ucell_), this->gd, this->orb_cutoff_);
        LR_Util::transpose_DMR(dm_trans, (*this->ucell_).nat);
        // D(X) real, for the grid Hxc force. The Coulomb kernel (mu nu | kappa lambda) is
        // symmetric within each index pair, so it only ever sees the symmetric part of D^X.
        // In `PulayForceStress::cal_pulay_fs`, `cal_gint_rho` (which builds v) symmetrizes
        // implicitly, while `cal_gint_fvl`'s internal factor 2 assumes D_{mu nu} = D_{nu mu}.
        // Passing an un-symmetrized D^X makes the two slots of the bilinear form disagree.
        // NOTE: `build_dm_from_dmk` symmetrizes `dm_trans_k` IN PLACE, hence the ordering.
        auto dm_trans_real =   // D(X), double (FIXME: not enough for periodic system!)
            LR_Util::build_dm_from_dmk<T, double>(dm_trans_k,
                this->paraMat_, this->nk, this->kv.kvec_d, (*this->ucell_), this->gd, this->orb_cutoff_,
                /*symmetrize=*/true);
        LR_Util::transpose_DMR(dm_trans_real, (*this->ucell_).nat);
        // LR_Util::print_DMR(dm_trans, "dm_trans of istate " + std::to_string(istate));
        // difference density matrix 
        std::vector<ct::Tensor> dm_diff_k = cal_dm_diff_pblas(this->X[ispin].template data<T>() + offset, this->paraX_[ispin], c, this->paraC_, this->nbasis, this->nocc[ispin], this->nvirt[ispin], this->paraMat_);
        std::cout << "dm_diff_k T(k) before symmetrization, istate " + std::to_string(istate) << std::endl;
        LR_Util::print_value(dm_diff_k[0].data<T>(), this->paraMat_.get_col_size(), this->paraMat_.get_row_size());
        // for (auto& d : dm_diff_k) { LR_Util::matsym(d.data<T>(), this->nbasis, this->paraMat_); }   // symmetrize
        // std::cout << "dm_diff_k T(k) after symmetrization, istate " + std::to_string(istate) << std::endl;
        // LR_Util::print_value(dm_diff_k[0].data<T>(), this->paraMat_.get_col_size(), this->paraMat_.get_row_size());

        const std::vector<ct::Tensor>& dm_relaxed_k = cal_dm_trans_pblas(Z.template data<T>() + offset, this->paraX_[ispin], c, this->paraC_, this->nbasis, this->nocc[ispin], this->nvirt[ispin], this->paraMat_);
        std::cout << "dm_relaxed_k Z(k) before symmetrization, istate " + std::to_string(istate) << std::endl;
        LR_Util::print_value(dm_relaxed_k[0].data<T>(), this->paraMat_.get_col_size(), this->paraMat_.get_row_size());
        for (auto& d : dm_relaxed_k) { LR_Util::matsym(d.data<T>(), this->nbasis, this->paraMat_); }    // symmetrize
        std::cout << "dm_relaxed_k Z(k) after symmetrization, istate " + std::to_string(istate) << std::endl;
        LR_Util::print_value(dm_relaxed_k[0].data<T>(), this->paraMat_.get_col_size(), this->paraMat_.get_row_size());
        // relaxed difference density matrix
        const std::vector<ct::Tensor>& relaxed_diff_dm_k = dm_diff_k + dm_relaxed_k;
        const elecstate::DensityMatrix<T, T>& diff_dm =
            LR_Util::build_dm_from_dmk<T, T>(dm_diff_k,
                this->paraMat_, this->nk, this->kv.kvec_d, (*this->ucell_), this->gd, this->orb_cutoff_);
        const elecstate::DensityMatrix<T, T>& relaxed_diff_dm =
            LR_Util::build_dm_from_dmk<T, T>(relaxed_diff_dm_k,
                this->paraMat_, this->nk, this->kv.kvec_d, (*this->ucell_), this->gd, this->orb_cutoff_);
        // LR_Util::print_DMR(relaxed_diff_dm, "relaxed_diff_dm T+Z (Z symmetrized) of istate " + std::to_string(istate));

        // elecstate::DensityMatrix<T, T> relaxed_diff_dm =    // T+D(Z), (R) can be complex
        //     LR_Util::build_dm_from_dmk<T, T>(
        //         // LR_Util::operator+(
        //             cal_dm_diff_pb las(this->X[ispin].template data<T>() + offset, this->paraX_[ispin], c, this->paraC_, this->nbasis, this->nocc[ispin], this->nvirt[ispin], this->paraMat_)
        //             + cal_dm_trans_pblas(Z.template data<T>() + offset, this->paraX_[ispin], c, this->paraC_, this->nbasis, this->nocc[ispin], this->nvirt[ispin], this->paraMat_)
        //             ,// ),
        //         this->paraMat_, this->nk, this->kv.kvec_d, (*this->ucell_), this->gd, this->orb_cutoff_);
        // LR_Util::print_DMR(relaxed_diff_dm, "relaxed_diff_dm of istate " + std::to_string(istate));
        elecstate::DensityMatrix<T, double> relaxed_diff_dm_real(&this->paraMat_, 1, this->kv.kvec_d, this->nk);
        LR_Util::initialize_DMR(relaxed_diff_dm_real, this->paraMat_, (*this->ucell_), this->gd, this->orb_cutoff_);
        LR_Util::get_DMR_real_imag_part(relaxed_diff_dm, relaxed_diff_dm_real, 'R');

        // get edm of type DensityMatrix
        // weak_ptr here is to avoid "could not match 'weak_ptr' against 'shared_ptr'"
        // but why there're no bug in the previous code (HamiltLR and HamiltULF)? 
        // seems because those two are classes having constructors 
        // but `cal_edm_from_XZ_istate` here is a functions
        std::weak_ptr<PotHxcLR> pot_weak = this->pot[ispin];
        std::weak_ptr<PotHxcLR> pot_hxc_gs_weak = this->pot_hxc_gs;
#ifdef __EXX
        std::weak_ptr<Exx_LRI<T>> exx_lri_weak = this->exx_lri;
#endif
        const std::vector<ct::Tensor>& edm_k =
            cal_edm_from_XZ_istate(this->X[ispin].template data<T>() + offset,
                Z.template data<T>() + offset,
                this->pelec->ekb.c[ ispin * nstates + istate],
                // pack the following as a struct or use parameter package
                this->eig_ks.c, dm_trans,
                c, this->nspin, this->nbasis, this->nocc, this->nvirt, (*this->ucell_), this->orb_cutoff_,
#ifdef __EXX
                exx_lri_weak, this->exx_info.info_global.hybrid_alpha,
#endif
                pot_weak, pot_hxc_gs_weak,
                this->kv, this->gd, this->paraX_, this->paraC_, this->paraMat_,
                this->xc_kernel, this->spin_types[ispin]);
        if (PARAM.inp.test_force && nocc[0] == 1 && nvirt[0] == 1)
        {
            const std::vector<ct::Tensor>& dm_diff = cal_dm_diff_pblas(this->X[0].template data<T>() + offset, this->paraX_[0], c, this->paraC_, this->nbasis, this->nocc[0], this->nvirt[0], this->paraMat_);
            // test_dm_diff_H2<T>(relaxed_diff_dm.get_DMK_pointer(0), c, this->nbasis);
            test_dm_diff_H2<T>(dm_diff[0].data<T>(), c, this->nbasis);
            test_edm_H2<T>(edm_k[0].data<T>(), this->eig_ks.c, c, this->nbasis);
        }
        elecstate::DensityMatrix<T, double> edm_real = LR_Util::build_dm_from_dmk<T, double>(edm_k,
            this->paraMat_, this->nk, this->kv.kvec_d, (*this->ucell_), this->gd, this->orb_cutoff_, /*symmetrize=*/true);
        // print edm_real (R)
        if (PARAM.inp.test_force)
        {
            LR_Util::save_DMR(edm_real, "data-EDMR-sparse", this->paraMat_);
            // LR_Util::print_DMR(edm_real, "edm_real (R) of istate " + std::to_string(istate));
        }

        ModuleBase::matrix force_hxc_dmtrans = lr_force.cal_force_hxc_dmtrans(dm_trans_real, *this->pot[ispin]);
        if (PARAM.inp.test_force)
            ModuleIO::print_force(GlobalV::ofs_running, (*this->ucell_), "HXC DMTRANS FORCE (eV/Angstrom)", force_hxc_dmtrans, false);

        const elecstate::DensityMatrix<T, double>& dm_gs = this->cal_dm_gs();
        ModuleBase::matrix force_hamiltgs_relaxed_diff = lr_force.cal_force_hamilt_gs_dm_relaxed_diff(relaxed_diff_dm_real, dm_gs, false, this->pot_hxc_gs.get());
        if (PARAM.inp.test_force)
            ModuleIO::print_force(GlobalV::ofs_running, (*this->ucell_), "H_GS-(T+Z) FORCE (without EXX) (eV/Angstrom)", force_hamiltgs_relaxed_diff, false);

        ModuleBase::matrix force_overlap_edm = lr_force.cal_force_overlap_edm(edm_real);    // "-" sign has been included in the force factor
        if (PARAM.inp.test_force)
            ModuleIO::print_force(GlobalV::ofs_running, (*this->ucell_), "OVERLAP-EDM FORCE (eV/Angstrom)", force_overlap_edm, false);

        if (PARAM.inp.test_force)
        {
            // test H[T] force (Z=0), non-EXX part
            elecstate::DensityMatrix<T, double> diff_dm_real(&this->paraMat_, 1, this->kv.kvec_d, this->nk);
            LR_Util::initialize_DMR(diff_dm_real, this->paraMat_, (*this->ucell_), this->gd, this->orb_cutoff_);
            LR_Util::get_DMR_real_imag_part(diff_dm, diff_dm_real, 'R');

            GlobalV::ofs_running << "========== [TEST H_GS-(T) force (Z=0), non-EXX part] ===========" << std::endl;
            ModuleBase::matrix force_hamiltgs_diff = lr_force.cal_force_hamilt_gs_dm_relaxed_diff(diff_dm_real, dm_gs);
            ModuleIO::print_force(GlobalV::ofs_running, (*this->ucell_), "H_GS-T FORCE (without EXX) (eV/Angstrom)", force_hamiltgs_diff, false);
            GlobalV::ofs_running << "========== [\\TEST H_GS-(T) force (Z=0), non-EXX part] ===========" << std::endl;
        }


#ifdef __EXX
        const double& alpha = this->exx_info.info_global.hybrid_alpha;

        if (LR::exx_kernel_list().count(xc_kernel))
        {
            const auto& Ds_trans = LR_Util::get_exx_Ds_spin1(dm_trans, (*this->ucell_), this->kv, this->paraMat_);
            ModuleBase::matrix force_exx_dmtrans = lr_force.cal_force_exx_dm_trans(Ds_trans, alpha * 4.0);  // cancel the two 0.5s in Ds
            if (PARAM.inp.test_force)
                ModuleIO::print_force(GlobalV::ofs_running, (*this->ucell_), "EXX DMTRANS FORCE (eV/Angstrom)", force_exx_dmtrans, false);
            force_hxc_dmtrans += force_exx_dmtrans;

        }

        if (LR::exx_kernel_list().count(PARAM.inp.dft_functional))
        {
            const auto& Ds_gs = LR_Util::get_exx_Ds_spin1(dm_gs, (*this->ucell_), this->kv, this->paraMat_);    // returns 0.5*D[0]
            const auto& Ds_relaxed_diff = LR_Util::get_exx_Ds_spin1(relaxed_diff_dm, (*this->ucell_), this->kv, this->paraMat_);   // returns 0.5*D[0]
            // LR_Util::print_CV(Ds_relaxed_diff, "Ds_relaxed_diff for EXX force");
            ModuleBase::matrix force_exx_gs_relaxed_diff = lr_force.cal_force_exx_gs_dm_relaxed_diff(Ds_gs, Ds_relaxed_diff, alpha * 4.0);  // cancel the two 0.5s in Ds
            if (PARAM.inp.test_force)
                ModuleIO::print_force(GlobalV::ofs_running, (*this->ucell_), "EXX GS-(T+Z) FORCE (eV/Angstrom)", force_exx_gs_relaxed_diff, false);
            force_hamiltgs_relaxed_diff += force_exx_gs_relaxed_diff;

            if (PARAM.inp.test_force)
            {
                // test H[T] force (Z=0), EXX part
                const auto& Ds_diff = LR_Util::get_exx_Ds_spin1(diff_dm, (*this->ucell_), this->kv, this->paraMat_);   // returns 0.5*D[0]
                GlobalV::ofs_running << "========== [TEST H_GS-(T) force (Z=0), EXX part] ===========" << std::endl;
                ModuleBase::matrix force_exx_gs_diff = lr_force.cal_force_exx_gs_dm_relaxed_diff(Ds_gs, Ds_diff, alpha * 4.0);  // cancel the two 0.5s in Ds
                ModuleIO::print_force(GlobalV::ofs_running, (*this->ucell_), "H_GS-T EXX FORCE (Z=0) (eV/Angstrom)", force_exx_gs_diff, false);
                GlobalV::ofs_running << "========== [\\TEST H_GS-(T) force (Z=0), EXX part] ===========" << std::endl;
            }
        }
#endif
        forces[istate] = force_hxc_dmtrans + force_hamiltgs_relaxed_diff + force_overlap_edm;
    }
    ModuleBase::timer::end("ESolver_LR", "cal_force");
    // total force
    print_force(forces, std::cout);
    print_force(forces, GlobalV::ofs_running);
    return forces;
}

template<typename T, typename TR>
elecstate::DensityMatrix<T, double> ModuleESolver::ESolver_LR<T, TR>::cal_dm_gs()
{
    elecstate::DensityMatrix<T, double> dm_gs(&this->paraMat_, this->nspin, this->kv.kvec_d, this->nk);
    elecstate::cal_dm_psi(&this->paraMat_all_, this->wg_ks_all, *this->psi_ks_all, dm_gs);   // nbands is important here
    LR_Util::initialize_DMR(dm_gs, this->paraMat_, (*this->ucell_), this->gd, this->orb_cutoff_);   // nbands is not important here
    dm_gs.cal_DMR();
    return dm_gs;
}

template<typename T, typename TR>
void ModuleESolver::ESolver_LR<T, TR>::test_force()
{
    LR_Force<T> lr_force((*this->ucell_), this->kv.kvec_d, this->paraMat_, *this->pw_rhod, *this->pw_rho,
        this->locpp, this->sf, this->gd, this->two_center_bundle_
#ifdef __EXX
        , std::weak_ptr<Exx_LRI<T>>(this->exx_lri), this->exx_info.info_global.hybrid_alpha
#endif
    );

    const elecstate::DensityMatrix<T, double>& dm_gs = this->cal_dm_gs();
    // LR_Util::print_DMR(dm_gs, "DM(R) of ground state");
    ///========================== test 1: reproduce the force of ground state =========================
    // energy density matrix of the ground state
    elecstate::DensityMatrix<T, double> edm_gs(&this->paraMat_, this->nspin, this->kv.kvec_d, this->nk);   //DX
    ModuleBase::matrix wg_ekb_ks_all(nspin, PARAM.inp.nbands);
    std::transform(this->wg_ks_all.c, this->wg_ks_all.c + nspin * PARAM.inp.nbands,
        this->eig_ks_all.c, wg_ekb_ks_all.c, std::multiplies<double>());
    elecstate::cal_dm_psi(&this->paraMat_all_, wg_ekb_ks_all, *this->psi_ks_all, edm_gs);
    LR_Util::initialize_DMR(edm_gs, this->paraMat_, (*this->ucell_), this->gd, this->orb_cutoff_);
    edm_gs.cal_DMR();
    // ground-state force
    ModuleBase::matrix force_gs = lr_force.reproduce_force_gs(kv, dm_gs, edm_gs);
    ModuleIO::print_force(GlobalV::ofs_running, (*this->ucell_), "Ground State FORCE (eV/Angstrom)", force_gs, false);
    /// ======================================= END test 1 =========================================
    ///========================== test 2: reproduce the DX Hartree term =========================
    ModuleBase::matrix f_hxc_potgs = lr_force.reproduce_force_gs_loc(dm_gs, *this->pot_gs_hartree);
    ModuleIO::print_force(GlobalV::ofs_running, (*this->ucell_), "GS Hartree force calculated by 'cal_pulay_fs' from potential (eV/Angstrom)", f_hxc_potgs, false);
    ModuleBase::matrix f_hxc_potlr = lr_force.cal_force_hxc_dmtrans(dm_gs, *this->pot[0]);
    ModuleIO::print_force(GlobalV::ofs_running, (*this->ucell_), "GS Hxc force calculated by 'LR_Force' from kernel (eV/Angstrom)", f_hxc_potlr, false);
    // `cal_force_hxc_dmtrans` now includes the Pulay -> Pulay+Hellmann-Feynman factor 2 itself,
    // so this must match the ground-state Hartree force directly (dm_gs is already symmetric).
    /// ======================================= END test 2 =========================================
    ///========================== test 3: H2 SZ 4-center gradients =========================
    if (this->nbasis == 2 && (*this->ucell_).nat == 2)
    {
        // lr_force.cal_H2_sz_center2_deriv(orb_cutoff_, kv);  // for gradient
        // lr_force.cal_H2_sz_center4(orb_cutoff_, kv, /*is_grad=*/false);  // for 4-center integrals
        // lr_force.cal_H2_sz_center4(orb_cutoff_, kv, /*is_grad=*/true);   // for gradient
        // exit(0);
    }
}

template class ModuleESolver::ESolver_LR<double, double>;
template class ModuleESolver::ESolver_LR<std::complex<double>, double>;
