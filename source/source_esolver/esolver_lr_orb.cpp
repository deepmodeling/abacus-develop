// LR orbital distributions, excitation amplitudes, and ground-state file input.
#include "esolver_lr_lcao_tddft.h"
#include "source_base/tool_quit.h"
#include "source_base/tool_title.h"
#include "source_lcao/module_lr/lr_gs.h"
#include "source_lcao/module_lr/utils/lr_util.h"

#include <cassert>
#include <complex>
#include <memory>

template <typename T, typename TR>
void ModuleESolver::ESolver_LR<T, TR>::set_parallel_orbitals_band(Parallel_Orbitals &pmat, const int nbands_in)
{
    LR::set_band_layout(pmat, this->nbasis, nbands_in, *this->ucell_, this->inp_->ri_hartree_benchmark == "aims");
}

template <typename T, typename TR>
void ModuleESolver::ESolver_LR<T, TR>::setup_eigenvectors_X()
{
    ModuleBase::TITLE("ESolver_LR", "setup_eigenvectors_X");
    LR::AmpLayout layout = LR::make_amp_layout(this->paraC_, this->nocc, this->nvirt, this->nk, this->openshell);
    this->paraX_ = std::move(layout.para_x);
    this->nloc_per_state = layout.nloc;

    this->X.resize(this->openshell ? 1 : this->nspin, LR_Util::newTensor<T>({this->nstates, this->nloc_per_state}));
    for (container::Tensor& x : this->X)
    {
        x.zero();
    }

    // if spectrum-only, read the LR-eigenstates from file and return
    if (this->inp_->lr_solver != "spectrum")
    {
        LR::set_x_init<T>(this->X, this->eig_ks, this->nocc, this->nvirt, this->npairs, this->paraX_, this->nstates,
                          this->nloc_per_state, this->nk, this->openshell, this->ofs_running_);
    }
}

template <typename T, typename TR>
void ModuleESolver::ESolver_LR<T, TR>::read_ks_wfc()
{
    assert(this->psi_ks != nullptr);
    this->eig_ks.create(this->kv.get_nks(), this->nbands);
    this->wg_ks.create(this->kv.get_nks(), this->nbands);
    if (this->inp_->ri_hartree_benchmark == "aims") // for aims benchmark
    {
        LR::read_wfc_aims<T>(*this->psi_ks, this->eig_ks, this->in_dir, this->nocc_in, this->nvirt_in, this->nbands,
                             this->nbasis);
    }
    else if (!LR::read_wfc_win<T>(this->in_dir, this->paraMat_, *this->psi_ks, this->eig_ks, this->wg_ks, this->kv,
                                  this->inp_->nspin, this->inp_->init_wfc_file_format == "binary",
                                  /*skip_bands=*/this->nocc_max - this->nocc_in))
    {
        ModuleBase::WARNING_QUIT("ESolver_LR", "read ground-state wavefunction failed.");
    }

    if (this->inp_->cal_force)
    { // allocate psi_ks_all and eig_ks_all to read all the bands
        this->psi_ks_all_own_.reset(new psi::Psi<T>(
            this->kv.get_nks(), paraMat_all_.ncol_bands, paraMat_all_.para_2d().nrloc(), this->kv.ngk, true));
        this->psi_ks_all_ = this->psi_ks_all_own_.get();
        this->eig_ks_all.create(this->kv.get_nks(), this->inp_->nbands);
        this->wg_ks_all.create(this->kv.get_nks(), this->inp_->nbands);
        if (!LR::read_wfc_all<T>(this->in_dir, paraMat_all_, *this->psi_ks_all_, this->eig_ks_all, this->wg_ks_all,
                                 this->kv, this->inp_->nspin, this->inp_->init_wfc_file_format == "binary"))
        {
            ModuleBase::WARNING_QUIT("ESolver_LR", "read all ground-state wavefunctions for force calculation failed.");
        }
        this->ofs_running_ << " Read in all the KS wavefunctions for force calculation. " << std::endl;
    }
}

template <typename T, typename TR>
void ModuleESolver::ESolver_LR<T, TR>::fill_z_window_(const int *desc_src)
{
    ModuleBase::TITLE("ESolver_LR", "fill_z_window_");
    LR::ZWin<T> z = LR::build_z_win<T>(this->psi_ks_all_, this->eig_ks_all, this->nocc, this->nvirt, this->nocc_max,
                                       this->nbands, this->nbasis, this->nk, this->openshell, this->kv, this->paraC_,
                                       this->paraMat_, desc_src, this->inp_->cal_force, !this->ks_initialized_,
                                       this->ofs_running_);
    this->psi_ks_z_ = std::move(z.psi);
    this->eig_ks_z_ = std::move(z.eig);
    this->paraC_z_ = std::move(z.para_c);
    this->paraX_z_ = std::move(z.para_x);
    this->nvirt_z_ = std::move(z.nvirt);
    this->nbands_z_ = z.nbands;
    this->nloc_per_state_z_ = z.nloc;
}

template <typename T, typename TR>
void ModuleESolver::ESolver_LR<T, TR>::read_ks_chg(Charge &chg_gs)
{
    LR::read_gs_chg(this->in_dir, *this->inp_, this->nspin, *this->ucell_, pgrid(), this->pw_rho, chg_gs,
                    this->my_rank_, this->ofs_running_);
}

template class ModuleESolver::ESolver_LR<double, double>;
template class ModuleESolver::ESolver_LR<std::complex<double>, double>;
