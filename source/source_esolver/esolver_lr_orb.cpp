// LR orbital distributions, excitation amplitudes, and ground-state file input.
#include "esolver_lr_lcao_tddft.h"
#include "source_lcao/module_lr/utils/lr_util.h"
#include "source_lcao/module_lr/ri_benchmark/ri_benchmark.h"
#include "source_io/module_wf/read_wfc_nao.h"
#include "source_io/module_output/cube_io.h"
#include "source_base/module_external/scalapack_connector.h"
#include <algorithm>

using namespace LR;

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
    const int block_size = this->paraC_.get_block_size();
    for (int is = 0;is < nspin;++is)
    {
        Parallel_2D px;
        LR_Util::setup_2d_division(px, block_size, this->nvirt[is], this->nocc[is]
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
            ModuleBase::WARNING_QUIT("ESolver_LR", "read all ground-state wavefunctions for force calculation failed.");
        }
        this->ofs_running_ << " Read in all the KS wavefunctions for force calculation. " << std::endl;
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

    const int block_size = this->paraC_.get_block_size();
    LR_Util::setup_2d_division(this->paraC_z_, block_size, this->nbasis, this->nbands_z_
#ifdef __MPI
        , this->paraMat_.blacs_ctxt
#endif
    );
    this->paraX_z_.clear();
    for (int is = 0; is < this->nspin; ++is)
    {
        Parallel_2D px;
        LR_Util::setup_2d_division(px, block_size, this->nvirt_z_[is], this->nocc[is]
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
