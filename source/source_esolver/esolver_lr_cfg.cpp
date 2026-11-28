#include "source_lcao/module_lr/kernel_pots.h"
// Resolve LR dimensions and initialize kernel potentials.
#include "esolver_lr_lcao_tddft.h"
#include "source_io/module_parameter/parameter.h"
#include "source_lcao/module_lr/lr_dim.h"
#include <type_traits>
#include <vector>

template <typename T, typename TR>
void ModuleESolver::ESolver_LR<T, TR>::set_dim()
{
    LR::DimSrc src;
    src.inp = this->inp_;
    src.nlocal = PARAM.globalv.nlocal;
    src.nks = this->kv.get_nks();
    src.nkstot = this->kv.get_nkstot();
    src.rank = this->my_rank_;
    src.gamma_only = std::is_same<T, double>::value;
    src.wg_ks = this->ks_ ? &this->ks_->pelec->wg : nullptr;
    src.kpar = this->kv.para_k.kpar;
    src.rank_in_pool = this->kv.para_k.rank_in_pool;
    src.in_dir = this->in_dir;
    src.ucell = const_cast<UnitCell*>(this->ucell_);
    src.log = &this->ofs_running_;

    const LR::Dim dim = LR::resolve_dim(src);

    this->nspin = dim.nspin;
    this->nstates = dim.nstates;
    this->nbasis = dim.nbasis;
    this->nbands = dim.nbands;
    this->nk = dim.nk;
    this->nupdown = dim.nupdown;
    this->nocc_in = dim.nocc_in;
    this->nocc_max = dim.nocc_max;
    this->nvirt_in = dim.nvirt_in;
    this->openshell = dim.openshell;
    this->nocc = dim.nocc;
    this->nvirt = dim.nvirt;
    this->npairs = dim.npairs;
}

template <typename T, typename TR>
void ModuleESolver::ESolver_LR<T, TR>::reset_dim_spin()
{
    LR::Dim dim;
    dim.nspin = this->nspin;
    dim.nstates = this->nstates;
    dim.nk = this->nk;
    dim.nupdown = this->nupdown;
    dim.openshell = this->openshell;
    dim.nocc = this->nocc;
    dim.nvirt = this->nvirt;
    dim.npairs = this->npairs;

    LR::apply_spin_pol(dim, *this->inp_);

    this->nocc = dim.nocc;
    this->nvirt = dim.nvirt;
    this->npairs = dim.npairs;
    this->openshell = dim.openshell;
}

template <typename T, typename TR>
void ModuleESolver::ESolver_LR<T, TR>::init_pot(const Charge &chg_gs)
{
    this->pot.assign(nspin, nullptr);
    if (this->inp_->ri_hartree_benchmark != "none")
    {
        return;
    }
    const auto &grid = pgrid();
    auto prepared =
        LR::prepare_lr_potentials(*this->inp_, xc_kernel, *this->pw_rho, *this->ucell_, chg_gs, grid, nspin, openshell);
    this->pot = std::move(prepared.pot);
    if (this->inp_->cal_force)
    {
        this->init_pot_groundstate(chg_gs);
        this->pot_hxc_gs = LR::prepare_gs_potential(
            *this->inp_, xc_kernel, *this->pw_rho, *this->ucell_, chg_gs, grid, nspin, openshell, prepared);
    }
}

template class ModuleESolver::ESolver_LR<double, double>;
template class ModuleESolver::ESolver_LR<std::complex<double>, double>;
