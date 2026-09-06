#include "esolver_lr_lcao_bse.h"
#include <array>
#include "source_cell/module_neighbor/sltk_atom_arrange.h"
#include "source_io/module_parameter/parameter.h"
#include "source_io/module_output/print_info.h"
#include "source_hamilt/module_gint/gint.h"
#include "source_lcao/module_bse/hamilt_bse.h"
#include "source_lcao/module_lr/lr_spectrum.h"
#include "source_lcao/module_lr/utils/exciton_plotter.h"
#include "source_lcao/module_lr/utils/lr_io.h"

namespace ModuleESolver
{

template <typename T, typename TR>
void ESolver_BSE<T, TR>::before_all_runners(BaseCell& basecell, const Input_para& inp)
{
    basecell.require_kind(BaseCell::Kind::unitcell, __FUNCTION__);
    UnitCell& ucell = static_cast<UnitCell&>(basecell);

    ModuleBase::TITLE("ESolver_BSE", "before_all_runners");
    ModuleBase::timer::start("ESolver_BSE", "before_all_runners");
    this->ucell_ = &ucell;
    // xc kernel
    this->xc_kernel = LR_Util::tolower(inp.xc_kernel);

    // necessary steps in ESolver_FP
    ModuleESolver::ESolver_FP::before_all_runners(basecell, inp);
    this->pelec = new elecstate::ElecStateLCAO<T>();

    this->kRlist = LR_IO::RI_kRlist(*this->ucell_,
                                    &this->kv,
                                    this->rpa_dir,
                                    inp.bse_use_fine_kgrid,
                                    this->out_dir);
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "Set K-POINTS and R-list for RI");
    ModuleIO::print_parameters(ucell, this->kv, inp);

    this->parameter_check();

    /// read orbitals and build the interpolation table
    this->two_center_bundle_.build_orb(ucell.ntype, ucell.orbital_fn.data(), inp.orbital_dir);

    this->two_center_bundle_.to_LCAO_Orbitals(this->orb_, inp.lcao_ecut, inp.lcao_dk, inp.lcao_dr, inp.lcao_rmax,
                                              inp.out_element_info, inp.cal_force);
    this->orb_cutoff_ = this->orb_.cutoffs();
    if (LR_Util::tolower(this->inp_->abs_gauge) == "velocity")
    {
        this->setup_2center_table(this->two_center_bundle_, this->orb_, ucell);
    }

    this->set_dimension();
    //  setup 2d-block distribution for AO-matrix and KS wfc
    LR_Util::setup_2d_division(this->paraMat_, 1, this->nbasis, this->nbasis);
#ifdef __MPI
    this->paraMat_.set_desc_wfc_Eij(this->nbasis, this->nbands, this->paraMat_.get_row_size());
    int err = this->paraMat_.set_nloc_wfc_Eij(this->nbands, GlobalV::ofs_running, GlobalV::ofs_warning);
    this->paraMat_.set_atomic_trace(ucell.get_iat2iwt(), ucell.nat, this->nbasis);
#else
    this->paraMat_.nrow_bands = this->nbasis;
    this->paraMat_.ncol_bands = this->nbands;
#endif

    this->psi_ks = new psi::Psi<T>(this->kv.get_nks(),
                                   this->paraMat_.ncol_bands,
                                   this->paraMat_.get_row_size(),
                                   this->kv.ngk,
                                   true);
    this->psi_ks_global = new psi::Psi<T>(this->kv.get_nks(),
                                          this->nbands,
                                          this->nbasis,
                                          this->kv.ngk,
                                          true);
    this->read_ks_wfc();
    // NOTE: openshell is not implemented in BSE
    if (this->nspin == 2)
    {
        this->nupdown = this->cal_nupdown_form_occ(this->pelec->wg);
        this->reset_dim_spin2();
    }

    LR_Util::setup_2d_division(this->paraC_, this->paraMat_.get_block_size(), this->nbasis, this->nbands
#ifdef __MPI
            , this->paraMat_.blacs_ctxt
#endif
        );

    this->Pgrid.init(this->pw_rho->nx,
               this->pw_rho->ny,
               this->pw_rho->nz,
               this->pw_rho->nplane,
               this->pw_rho->nrxx,
               this->pw_big->nbz,
               this->pw_big->bz,
               GlobalV::NPROC);

    // search adjacent atoms and init Gint
    double search_radius = -1.0;
    search_radius = atom_arrange::set_sr_NL(GlobalV::ofs_running,
                                            inp.out_level,
                                            this->orb_.get_rcutmax_Phi(),
                                            ucell.infoNL->get_rcutmax_Beta(),
                                            PARAM.globalv.gamma_only_local);
    atom_arrange::search(PARAM.globalv.search_pbc,
                         GlobalV::ofs_running,
                         this->gd,
                         *this->ucell_,
                         search_radius,
                         inp.test_atom_input);

        this->gint_info_.reset(new ModuleGint::GintInfo(
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
            this->orb_.Phi,
            ucell,
            this->gd));
            ModuleGint::Gint::set_gint_info(this->gint_info_.get());

    this->pot.resize(this->nspin, nullptr);
    if (this->inp_->lr_solver != "spectrum" && this->inp_->lr_solver != "plot")
    {
        this->mo_lri = LR_Util::make_unique<BSE::MolecularLRI<T>>(*this->ucell_,
            this->nk,
            this->kRlist,
            this->nocc[0],
            this->nvirt[0],
            *this->psi_ks_global,
            this->inp_->bse_q_approx_mode,
            this->inp_->bse_q_approx_threshold,
            this->inp_->out_ri_cv,
            this->out_dir,
            GlobalV::MY_RANK,
            GlobalV::NPROC);

        if (!this->inp_->bse_ri_hartree && this->inp_->ri_hartree_benchmark == "none")
        {
            Charge chg_gs;
            this->read_ks_chg(chg_gs);
            this->init_pot(chg_gs);
        }
    }
    ModuleBase::timer::end("ESolver_BSE", "before_all_runners");
}

template <typename T, typename TR>
void ESolver_BSE<T, TR>::runner(BaseCell& basecell, const int istep)
{
    basecell.require_kind(BaseCell::Kind::unitcell, __FUNCTION__);
    UnitCell& ucell = static_cast<UnitCell&>(basecell);

    ModuleBase::TITLE("ESolver_BSE", "runner");
    ModuleBase::timer::start("ESolver_BSE", "runner");
    //allocate 2-particle state and setup 2d division
    this->allocate_eigen_infos();

    auto efile_out = [&](const std::string& label)->std::string {
        return this->out_dir + "Excitation_Energy_" + label + ".dat";};
    auto vfile_out = [&](const std::string& label)->std::string {
        return this->out_dir + "Excitation_Amplitude_" + label + "_" + std::to_string(GlobalV::MY_RANK) + ".dat";};
    auto efile_in = [&](const std::string& label)->std::string {
        return this->in_dir + "Excitation_Energy_" + label + ".dat";};
    auto vfile_in = [&](const std::string& label)->std::string {
        return this->in_dir + "Excitation_Amplitude_" + label + "_" + std::to_string(GlobalV::MY_RANK) + ".dat";};

    if (this->inp_->lr_solver == "elpa")
    {
        if (this->inp_->bse_spin_types == std::vector<std::string>{"ipa"})
        {
            this->ipa_solver();
        }
        else
        {
            std::cout << "Calculating Casida/BSE matrix directly." << std::endl;
            assert(this->xc_kernel == "bse");
            this->lri_init();
            BSE::HamiltBSE<T> bse_matrix(this->nspin, this->nbasis, this->nocc, this->nvirt, *this->ucell_,
                                    this->orb_cutoff_, this->gd, *this->psi_ks, *this->psi_ks_global, this->eig_gw,
                                    *this->mo_lri,
                                    this->pot[0], this->kv, this->paraX_, this->paraC_, this->paraMat_,
                                    this->inp_->bse_spin_types,
                                    this->inp_->bse_tda,
                                    this->inp_->bse_ri_hartree,
                                    this->inp_->bse_mem_save,
                                    this->inp_->bse_continue,
                                    this->inp_->out_bse_ab,
                                    this->out_dir,
                                    this->in_dir,
                                    GlobalV::MY_RANK,
                                    GlobalV::NPROC,
                                    this->inp_->ri_hartree_benchmark);

            auto write_tda_states = [&](const std::string& label,
                                        const Real<T>* e,
                                        const T* v,
                                        const int& dim,
                                        const int& nst,
                                        const int& prec = 8) -> void {
                if (GlobalV::MY_RANK == 0)
                {
                    assert(nst == LR_Util::write_value(efile_out(label), prec, e, nst));
                }
                assert(nst * dim == LR_Util::write_value(vfile_out(label), prec, v, nst, dim));
                ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "write tda states " + label);
            };
            auto write_full_states = [&](const std::string& label,
                                         const Real<T>* e,
                                         const T* X,
                                         const T* Y,
                                         const int& dim,
                                         const int& nst,
                                         const int& prec = 8) -> void {
                if (GlobalV::MY_RANK == 0)
                {
                    assert(nst == LR_Util::write_value(efile_out("full_"+label), prec, e, nst));
                }
                assert(nst * dim == LR_Util::write_value(vfile_out("full_X_"+label), prec, X, nst, dim));
                assert(nst * dim == LR_Util::write_value(vfile_out("full_Y_"+label), prec, Y, nst, dim));
                ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "write full states " + label);
            };

            if ((this->inp_->bse_tda == "both" || this->inp_->bse_tda == "tda"))
            {
                for (int is = 0; is < this->inp_->bse_spin_types.size(); ++is)
                {
                    bse_matrix.tda_solver(is,
                                          this->nstates,
                                          &this->tda_ene[is * this->nstates],
                                          this->X[is].template data<T>());

                    std::cout << "eigenvalues: (Ry)" << std::endl;
                    int write_nstates = std::min(this->nstates, 20);
                    LR_Util::print_value(&this->tda_ene[is * this->nstates], write_nstates);
                    std::cout << "eigenvalues: (eV)" << std::endl;
                    for (int i = 0; i < write_nstates; ++i)
                    {
                        std::cout << this->tda_ene[is * this->nstates + i] * ModuleBase::Ry_to_eV << " ";
                    }
                    std::cout << std::endl;
                    std::cout << "Excition binding energies (eV):"
                              << (direct_gap - tda_ene[is * this->nstates]) * ModuleBase::Ry_to_eV << std::endl;

                    if (this->inp_->out_wfc_lr)
                    {
                        write_tda_states(this->inp_->bse_spin_types[is],
                                         &this->tda_ene[is * this->nstates],
                                         this->X[is].template data<T>(),
                                         this->nloc_per_state,
                                         this->nstates);
                    }
                    malloc_trim(0);
                }
            }
            if ((this->inp_->bse_tda == "both" || this->inp_->bse_tda == "full"))
            {
                for (int is = 0; is < this->inp_->bse_spin_types.size(); ++is)
                {
                    bse_matrix.full_solver(is,
                                           this->nstates,
                                           &this->full_ene[is * this->nstates],
                                           this->full_X[is].template data<T>(),
                                           this->full_Y[is].template data<T>());

                    std::cout << "eigenvalues: (Ry)" << std::endl;
                    int write_nstates = std::min(this->nstates, 20);
                    LR_Util::print_value(&this->full_ene[is * this->nstates], write_nstates);
                    for (int i = 0; i < write_nstates; ++i)
                    {
                        std::cout << this->full_ene[is * this->nstates + i] * ModuleBase::Ry_to_eV << " ";
                    }
                    std::cout << std::endl;
                    std::cout << "Excition binding energies (eV):"
                              << (direct_gap - full_ene[is * this->nstates]) * ModuleBase::Ry_to_eV << std::endl;

                    if (this->inp_->out_wfc_lr)
                    {
                        write_full_states(this->inp_->bse_spin_types[is],
                                          &this->full_ene[is * this->nstates],
                                          this->full_X[is].template data<T>(),
                                          this->full_Y[is].template data<T>(),
                                          this->nloc_per_state,
                                          this->nstates);
                    }
                    malloc_trim(0);
                }
            }
        }
    }
    else if (this->inp_->lr_solver == "spectrum" || this->inp_->lr_solver == "plot")
    {
        std::cout << "Reading BSE excitation states from file." << std::endl;
        auto read_tda_states = [&](const std::string& label, Real<T>* e, T* v, const int& dim, const int& nst)->void
        {
            if (GlobalV::MY_RANK == 0) {
                assert(nst == LR_Util::read_value(efile_in(label), e, nst));
                ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "finish reading " + efile_in(label));
            }
#ifdef __MPI
            MPI_Bcast(e, nst, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif
            assert(nst * dim == LR_Util::read_value(vfile_in(label), v, nst, dim));
            ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "finish reading " + vfile_in(label));
        };
        auto read_full_states = [&](const std::string& label, Real<T>* e, T* X, T* Y, const int& dim, const int& nst)->void
        {
            if (GlobalV::MY_RANK == 0) {
                assert(nst == LR_Util::read_value(efile_in("full_"+label), e, nst));
                ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "finish reading " + efile_in("full_"+label));
            }
#ifdef __MPI
            MPI_Bcast(e, nst, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif
            assert(nst * dim == LR_Util::read_value(vfile_in("full_X_"+label), X, nst, dim));
            ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "finish reading " + vfile_in("full_X_"+label));
            assert(nst * dim == LR_Util::read_value(vfile_in("full_Y_"+label), Y, nst, dim));
            ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "finish reading " + vfile_in("full_Y_"+label));
        };

        if (this->inp_->bse_tda == "both" || this->inp_->bse_tda == "tda")
        {
            for (int is = 0; is < this->inp_->bse_spin_types.size(); ++is)
            {
                read_tda_states(this->inp_->bse_spin_types[is],
                                &this->tda_ene[is * this->nstates],
                                this->X[is].template data<T>(),
                                this->nloc_per_state,
                                this->nstates);
            }
        }
        if (this->inp_->bse_tda == "both" || this->inp_->bse_tda == "full")
        {
            for (int is = 0; is < this->inp_->bse_spin_types.size(); ++is)
            {
                read_full_states(this->inp_->bse_spin_types[is],
                                 &this->full_ene[is * this->nstates],
                                 this->full_X[is].template data<T>(),
                                 this->full_Y[is].template data<T>(),
                                 this->nloc_per_state,
                                 this->nstates);
                // check whether |X|^2 - |Y|^2 = 1
                for (int i = 0; i < this->nstates; ++i) 
                {
                    double norm_xy = 0.0;
                    for (int j = 0; j < this->nloc_per_state; ++j) {
                        norm_xy += std::norm(this->full_X[is].template data<T>()[i * this->nloc_per_state + j])
                                   - std::norm(this->full_Y[is].template data<T>()[i * this->nloc_per_state + j]);
                    }
                    Parallel_Reduce::reduce_all(norm_xy);
                    if (std::abs(norm_xy - 1.0) > 1e-6){
                        std::cout << "| CHECK WARNING: for full excitation " << i
                                  << ", |X|^2 - |Y|^2 = " << std::setprecision(10) << norm_xy << std::endl;
                    }
                }
            }
        }
    }
    else
    {
        ModuleBase::WARNING_QUIT("ESolver_BSE", "lr_solver must be elpa, plot or spectrum");
    }
    ModuleBase::timer::end("ESolver_BSE", "runner");
    return;
}

template <typename T, typename TR>
void ESolver_BSE<T, TR>::after_all_runners(BaseCell& basecell)
{
    basecell.require_kind(BaseCell::Kind::unitcell, __FUNCTION__);
    UnitCell& ucell = static_cast<UnitCell&>(basecell);

    ModuleBase::TITLE("ESolver_BSE", "after_all_runners");
    ModuleBase::timer::start("ESolver_BSE", "after_all_runners");
    const std::string& output_dir = this->out_dir;
    const std::set<std::string> benchmarks = {"abacus-librpa", "abacus", "none" };
    if (benchmarks.find(this->inp_->ri_hartree_benchmark) == benchmarks.end())
    {
        return;
    } // no need to calculate the spectrum

    if (this->inp_->lr_solver == "plot")
    {
        for (int is = 0; is < this->X.size(); ++is)
        {
            std::cout << "plot BSE exciton wavefunction for state: " << this->inp_->plot_istate
                << ", spin type: " << this->inp_->bse_spin_types[is] << std::endl;
            LR_Util::ExcitonPlotter<T> eplot(this->nspin, this->nbasis, this->nocc, this->nvirt, *this->psi_ks,
                                            *this->ucell_, this->kv, this->gd, this->orb_cutoff_, this->Pgrid, *this->pw_rho,
                                            this->paraX_, this->paraC_, this->paraMat_,
                                            output_dir,
                                            &this->tda_ene[is * this->nstates], this->X[is].template data<T>(),
                                            false/*openshell*/, &this->orb_);
            const std::string plot_type = LR_Util::tolower(this->inp_->exciton_plot_type);
            const std::string plot_format = LR_Util::tolower(this->inp_->exciton_plot_format);
            const bool write_slice = (plot_format == "slice" || plot_format == "both");
            const bool write_cube = (plot_format == "cube" || plot_format == "both");

            if (plot_format != "cube" && plot_format != "slice" && plot_format != "both")
            {
                ModuleBase::WARNING_QUIT("ESolver_BSE", "exciton_plot_format must be cube, slice, or both");
            }

            if (plot_type == "conditional")
            {
                if (plot_format != "slice")
                {
                    ModuleBase::WARNING_QUIT(
                        "ESolver_BSE",
                        "conditional exciton density only supports exciton_plot_format = slice");
                }
                if (this->inp_->exciton_fixed_coordinate.size() != 6)
                {
                    ModuleBase::WARNING_QUIT(
                        "ESolver_BSE",
                        "exciton_fixed_coordinate must contain six values: hole x y z followed by electron x y z");
                }
                const std::array<double, 3> r_h_fix = {this->inp_->exciton_fixed_coordinate[0],
                                                       this->inp_->exciton_fixed_coordinate[1],
                                                       this->inp_->exciton_fixed_coordinate[2]};
                const std::array<double, 3> r_e_fix = {this->inp_->exciton_fixed_coordinate[3],
                                                       this->inp_->exciton_fixed_coordinate[4],
                                                       this->inp_->exciton_fixed_coordinate[5]};
                eplot.plot_cond_slice(this->inp_->plot_istate, r_h_fix,
                    this->inp_->exciton_slice_plane,
                    this->inp_->exciton_slice_pos,
                    this->inp_->exciton_slice_npoints,
                    this->inp_->exciton_slice_range, "elec");
                eplot.plot_cond_slice(this->inp_->plot_istate, r_e_fix,
                    this->inp_->exciton_slice_plane,
                    this->inp_->exciton_slice_pos,
                    this->inp_->exciton_slice_npoints,
                    this->inp_->exciton_slice_range, "hole");
            }
            else if (plot_type == "average")
            {
                if (write_cube)
                {
                    // Average hole density: integrates out the electron coordinate
                    eplot.plot_average_density(this->inp_->plot_istate, "hole");
                    // Average electron density: integrates out the hole coordinate
                    eplot.plot_average_density(this->inp_->plot_istate, "elec");
                }
                if (write_slice)
                {
                    eplot.plot_average_slice(this->inp_->plot_istate, "hole",
                        this->inp_->exciton_slice_plane,
                        this->inp_->exciton_slice_pos,
                        this->inp_->exciton_slice_npoints,
                        this->inp_->exciton_slice_range);
                    eplot.plot_average_slice(this->inp_->plot_istate, "elec",
                        this->inp_->exciton_slice_plane,
                        this->inp_->exciton_slice_pos,
                        this->inp_->exciton_slice_npoints,
                        this->inp_->exciton_slice_range);
                }
            }
            else
            {
                ModuleBase::WARNING_QUIT("ESolver_BSE", "exciton_plot_type must be average or conditional");
            }
        }
    }
    if (this->inp_->lr_solver == "spectrum" || this->inp_->lr_solver == "elpa")
    {
        std::cout << "Calculating BSE optical absorption spectrum." << std::endl;
        if (LR_Util::tolower(this->inp_->abs_gauge) == "velocity" )
        {
            const int nspin_tmp = this->inp_->nspin == 2 ? 2 : 1;
            this->velocity_mo = LR_Util::cal_velocity_mo(*this->ucell_, this->gd, this->two_center_bundle_, 
                                                        this->paraMat_, this->paraC_, this->kv, *this->psi_ks, 
                                                        this->nk, nspin_tmp, this->nbasis, this->nocc, this->nvirt);
        }
        if (this->inp_->bse_tda == "both" || this->inp_->bse_tda == "tda")
        {
            for (int is = 0; is < this->X.size(); ++is)
            {
                LR::LR_Spectrum<T> spectrum(this->nspin, this->nbasis, this->nocc, this->nvirt, *this->pw_rho, *this->psi_ks,
                                            *this->ucell_, this->kv, this->gd, this->orb_cutoff_, this->two_center_bundle_,
                                            this->paraX_, this->paraC_, this->paraMat_,
                                            &this->tda_ene[is * this->nstates], this->eig_ks.c,
                                            this->X[is].template data<T>(), this->nstates, false/*openshell*/,
                                            LR_Util::tolower(this->inp_->abs_gauge), GlobalV::MY_RANK, output_dir);
                if (LR_Util::tolower(this->inp_->abs_gauge) == "velocity")
                {
                    spectrum.set_vmo(this->velocity_mo.data());
                }
                spectrum.cal_spectrum();
                spectrum.transition_analysis(this->inp_->bse_spin_types[is]+"_tda");              
                if (this->inp_->bse_spin_types[is] != "triplet")        // triplets has no transition dipole and no contribution to the spectrum
                {
                    spectrum.write_transition_dipole(output_dir +
                        "trans_dipole_" + this->inp_->bse_spin_types[is] + "_tda.dat");
                    // ============================== for test ==============================
                    if (LR_Util::tolower(this->inp_->abs_gauge) == "velocity")
                    {   //// TEST the formula v/omega rather than v/(e_a-e_i)
                        // spectrum.test_transition_dipoles_velocity_omega();
                        // spectrum.write_transition_dipole(out_dir + 
                        //     "trans_dipole_" + spin_types[is] + "_vomega_tda.dat");
                    }
                    // ============================== for test ==============================
                }
            }        
        }
        if (this->inp_->bse_tda == "both" || this->inp_->bse_tda == "full")
        {
            for (int is = 0;is < this->full_X.size();++is)
            {
                LR::LR_Spectrum<T> spectrum(this->nspin, this->nbasis, this->nocc, this->nvirt, *this->pw_rho, *this->psi_ks,
                                            *this->ucell_, this->kv, this->gd, this->orb_cutoff_, this->two_center_bundle_,
                                            this->paraX_, this->paraC_, this->paraMat_,
                                            &this->full_ene[is * this->nstates], this->eig_ks.c,
                                            this->full_X[is].template data<T>(), this->nstates, false/*openshell*/,
                                            LR_Util::tolower(this->inp_->abs_gauge), GlobalV::MY_RANK, output_dir);
                if (LR_Util::tolower(this->inp_->abs_gauge) == "velocity")
                {
                    spectrum.set_vmo(this->velocity_mo.data());
                }
                spectrum.set_Y(this->full_Y[is].template data<T>());
                spectrum.set_full(true);
                spectrum.cal_spectrum();
                spectrum.transition_analysis(this->inp_->bse_spin_types[is]+"_full");
                if (this->inp_->bse_spin_types[is] != "triplet")        // triplets has no transition dipole and no contribution to the spectrum
                {
                    spectrum.write_transition_dipole(output_dir +
                        "trans_dipole_" + this->inp_->bse_spin_types[is] + "_full.dat");
                }
            }
        }
    }
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "ESolver_BSE::after_all_runners");
    ModuleBase::timer::end("ESolver_BSE", "after_all_runners");
}

template<typename T, typename TR>
void ESolver_BSE<T, TR>::ipa_solver()
{// if ipa, assign X as identity matrix directly
    ModuleBase::TITLE("ESolver_BSE", "ipa_solver");
    ModuleBase::timer::start("ESolver_BSE", "ipa_solver");
    std::cout << "Independent particle approximation is used, assign X as identity matrix directly." << std::endl;
    assert(this->inp_->bse_tda == "tda");
    std::vector<double> ev(this->nk * this->nocc[0] * this->nvirt[0], 0.0);
    for (int ik = 0; ik < this->nk; ++ik)
    {
        for (int i = 0; i < this->nocc[0]; ++i)
        {
            for (int a = 0; a < this->nvirt[0]; ++a)
            {
                int index = ik * this->nocc[0] * this->nvirt[0] + i * this->nvirt[0] + a;
                ev[index] = this->eig_gw(ik, this->nocc[0] + a) - this->eig_gw(ik, i);
            }
        }
    }

    std::vector<int> indices(ev.size());
    std::iota(indices.begin(), indices.end(), 0); // [0, 1, 2, ..., size-1]

    std::sort(indices.begin(), indices.end(), [&](int lhs, int rhs) {
        return ev[lhs] < ev[rhs];
    });
    std::sort(ev.begin(), ev.end());
    std::copy_n(ev.data(), this->nstates, this->tda_ene.data());
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (std::size_t istate = 0; istate < this->nstates; ++istate)
    {
        int sorted_index = indices[istate];
        int ik = sorted_index / (this->nocc[0] * this->nvirt[0]);
        int loffset_X = (istate * this->nk + ik) * this->paraX_[0].get_local_size();
        int i = (sorted_index / this->nvirt[0]) % this->nocc[0];
        int a = sorted_index % this->nvirt[0];
        int col_loc = this->paraX_[0].global2local_col(i);
        int row_loc = this->paraX_[0].global2local_row(a);
        if (col_loc == -1 || row_loc == -1) continue;
        this->X[0].template data<T>()[loffset_X + col_loc * this->paraX_[0].get_row_size() + row_loc] = 1.0;
    }
    ModuleBase::timer::end("ESolver_BSE", "ipa_solver");
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "IPA solver");
}

template<typename T, typename TR>
void ESolver_BSE<T, TR>::lri_init()
{
    ModuleBase::TITLE("ESolver_BSE", "LRI init");
    using TA = int;
    using TC = std::array<int, 3>;
    using TAC = std::pair<TA, TC>;
    // start reading Ws and Cs
    std::map<TA, std::map<TAC, RI::Tensor<T>>> Cs_in;
    std::map<TA, std::map<TAC, RI::Tensor<T>>> Vs_in;
    std::map<TA, std::map<TAC, RI::Tensor<T>>> Ws_in;
    // if (GlobalV::MY_RANK == 0) // comment to read from all processes to avoid communication
    // {
        Cs_in = LRI_CV_Tools::read_Cs_ao_all<T>(this->rpa_dir);
        if (this->inp_->ri_hartree_benchmark == "aims-librpa" )
        {
            Vs_in = LR_IO::read_coulomb_mat_general_k<T, T>(this->rpa_dir, Cs_in, this->kRlist);
        }
        else if (this->inp_->ri_hartree_benchmark == "none" || this->inp_->ri_hartree_benchmark == "abacus-librpa" )
        {
            Vs_in = LR_IO::read_coulomb_mat_k<T, T>(this->rpa_dir, Cs_in, this->kRlist);
        }
        Ws_in = LR_IO::read_Ws<T, T>(Vs_in, this->kRlist.Rlist);
    // }
#ifdef __MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    this->mo_lri->init(Cs_in, Vs_in, Ws_in, this->exx_info.info_ri);
    malloc_trim(0);
    ModuleBase::TITLE("ESolver_BSE", "Finish LRI init");
}

template<typename T, typename TR>
void ESolver_BSE<T, TR>::read_ks_wfc()
{
    assert(this->psi_ks != nullptr);
    this->pelec->ekb.create(this->kv.get_nks(), this->nbands);
    this->pelec->wg.create(this->kv.get_nks(), this->nbands);
    this->eig_gw.create(this->kv.get_nks(), this->nbands);

    int ncore = 0; // skip core bands
    int nbands_file = 0;
    int nk_file = 0;
    int nspin_file = 0;
    int nocc_file = 0;
    int nspin_tmp = this->inp_->nspin == 2 ? 2 : 1;
    LR_IO::parse_band_out_file(this->rpa_dir, nbands_file, nk_file, nspin_file, nocc_file);
    if (nk_file != this->nk) {
        ModuleBase::WARNING_QUIT("ESolver_BSE", "Inconsistence: The nk in `band_out` is " + std::to_string(nk_file)
             + ", while BSE::nk is " + std::to_string(this->nk));
    }
    std::vector<double> eig_gw_info;
    if (this->inp_->bse_use_fine_kgrid)
    {
        eig_gw_info = LR_IO::read_energy_qp_from_band_files(this->kv, this->nocc[0], this->nvirt[0], ncore,
                                        this->rpa_dir, this->nk, nspin_tmp, nspin_file);
        LR_IO::read_librpa_eigenvectors_from_band_files<T>(*this->psi_ks, *this->psi_ks_global,
            this->rpa_dir, ncore, nbands_file, nspin_tmp, nspin_file, GlobalV::MY_RANK, this->paraMat_);
    }
    else
    {
        eig_gw_info = LR_IO::read_energy_qp(this->nocc[0], this->nvirt[0],
                                        this->rpa_dir, ncore, this->nk, nspin_tmp, nspin_file);
        LR_IO::read_librpa_eigenvectors<T>(*this->psi_ks, *this->psi_ks_global,
            this->rpa_dir, ncore, nbands_file, nspin_tmp, nspin_file, GlobalV::MY_RANK, this->paraMat_);
    }
    int cbm_k(0), vbm_k(0), direct_k(0);
    for (int iks = 0; iks < this->kv.get_nks(); ++iks) {
        for (int ib = 0; ib < this->nbands; ++ib) {
            this->pelec->wg(iks, ib) = eig_gw_info[iks * this->nbands *3 + ib * 3 + 0];
            this->pelec->ekb(iks, ib) = eig_gw_info[iks * this->nbands *3 + ib * 3 + 1];
            this->eig_gw(iks, ib) = eig_gw_info[iks * this->nbands *3 + ib * 3 + 2];
        }
        double cbm = this->eig_gw(iks, this->nocc[0]);
        for (int ib = this->nocc[0]; ib < this->nbands; ++ib) { // in case of non-ordered bands
            double e = this->eig_gw(iks, ib);
            if (e < cbm) cbm = e;
        }
        double vbm = this->eig_gw(iks, this->nocc[0]-1);
        for (int ib = 0; ib < this->nocc[0]-1; ++ib) {
            double e = this->eig_gw(iks, ib);
            if (e > vbm) vbm = e;
        }
        if (iks == 0) {
            this->cbm_energy = cbm;
            this->vbm_energy = vbm;
            this->direct_gap = cbm - vbm;
        }
        else {
            if (this->cbm_energy > cbm) {
                this->cbm_energy = cbm;
                cbm_k = iks;
            }
            if (this->vbm_energy < vbm) {
                this->vbm_energy = vbm;
                vbm_k = iks;
            }
            if (this->direct_gap > cbm - vbm) {
                this->direct_gap = cbm - vbm;
                direct_k = iks;
            }
        }
    }
    std::cout << "VBM energy (eV): " << this->vbm_energy * ModuleBase::Ry_to_eV << " at k " << vbm_k << std::endl;
    std::cout << "CBM energy (eV): " << this->cbm_energy * ModuleBase::Ry_to_eV << " at k " << cbm_k << std::endl;
    std::cout << "Indirect gap (eV): " << (this->cbm_energy - this->vbm_energy) * ModuleBase::Ry_to_eV << std::endl;
    std::cout << "Direct gap (eV): " << this->direct_gap * ModuleBase::Ry_to_eV << " at k " << direct_k << std::endl;

    this->eig_ks = std::move(this->pelec->ekb);
}

template<typename T, typename TR>
void ESolver_BSE<T, TR>::init_pot(const Charge& chg_gs)
{
    switch (this->nspin)
    {
        using ST = LR::PotHxcLR::SpinType;
    case 1: case 2:
        this->pot[0] = std::make_shared<LR::PotHxcLR>(this->xc_kernel, *this->pw_rho, *this->ucell_, chg_gs, this->Pgrid,
            ST::S1, this->inp_->lr_init_xc_kernel);
        break;
    // case 2:
    //     this->pot[0] = std::make_shared<PotHxcLR>(xc_kernel, *this->pw_rho, ucell, chg_gs, Pgrid, openshell ? ST::S2_updown : ST::S2_singlet, this->inp_->lr_init_xc_kernel);
    //     this->pot[1] = std::make_shared<PotHxcLR>(xc_kernel, *this->pw_rho, ucell, chg_gs, Pgrid, openshell ? ST::S2_updown : ST::S2_triplet, this->inp_->lr_init_xc_kernel);
    //     break;
    default:
        throw std::invalid_argument("ESolver_BSE: nspin must be 1 or 2");
    }
}

template<typename T, typename TR>
void ESolver_BSE<T, TR>::allocate_eigen_infos()
{
    ModuleBase::TITLE("ESolver_BSE", "allocate_eigen_infos");

    for (int is = 0; is < this->nspin; ++is)
    {
        Parallel_2D px;
        LR_Util::setup_2d_division(px, /*nb2d=*/1, this->nvirt[is], this->nocc[is]
#ifdef __MPI
            , this->paraC_.blacs_ctxt
#endif
        );
        this->paraX_.emplace_back(std::move(px));
    }
    this->nloc_per_state = this->nk
                           * (this->openshell ? this->paraX_[0].get_local_size() + this->paraX_[1].get_local_size()
                                              : this->paraX_[0].get_local_size());

    int n_spin_types = this->inp_->bse_spin_types.size();
    if (this->inp_->bse_tda == "both" || this->inp_->bse_tda == "tda") {
        BSE_Util::print_mem_estimate("TDA BSE eigen states",
                                     n_spin_types * static_cast<std::size_t>(this->nstates)
                                         * (1 + this->nloc_per_state),
                                     sizeof(T));
        this->tda_ene.resize(n_spin_types * this->nstates);
        this->X.resize(n_spin_types, LR_Util::newTensor<T>({ this->nstates, this->nloc_per_state }));
        for (auto& x : this->X) { x.zero(); }
    }
    if (this->inp_->bse_tda == "both" || this->inp_->bse_tda == "full") {
        BSE_Util::print_mem_estimate("full BSE eigen states",
                                     n_spin_types * static_cast<std::size_t>(this->nstates)
                                         * (1 + 2 * this->nloc_per_state),
                                     sizeof(T));
        this->full_ene.resize(n_spin_types * this->nstates);
        this->full_X.resize(n_spin_types, LR_Util::newTensor<T>({ this->nstates, this->nloc_per_state }));
        this->full_Y.resize(n_spin_types, LR_Util::newTensor<T>({ this->nstates, this->nloc_per_state }));
        for (auto& x : this->full_X) { x.zero(); }
        for (auto& y : this->full_Y) { y.zero(); }
    }
}

template class ESolver_BSE<double, double>;
template class ESolver_BSE<std::complex<double>, double>;
} // namespace ModuleESolver
