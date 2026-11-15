#include "esolver_lr_lcao_tddft.h"

#include "source_base/timer.h"
#include "source_lcao/module_lr/gradient_inputs.h"
#include "source_lcao/module_lr/lr_force.h"
#include "source_lcao/module_lr/nacv_cs.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <stdexcept>

namespace
{
/** @brief Time the complete NAC orchestration, including exception paths. */
class NacOrchestrationTimer
{
  public:
    NacOrchestrationTimer() { ModuleBase::timer::start("ESolver_LR", "cal_nac"); }
    ~NacOrchestrationTimer() { ModuleBase::timer::end("ESolver_LR", "cal_nac"); }
    NacOrchestrationTimer(const NacOrchestrationTimer &) = delete;
    NacOrchestrationTimer &operator=(const NacOrchestrationTimer &) = delete;
};
} // namespace

namespace ModuleESolver
{
template <typename T, typename TR> LR::NacvResult ESolver_LR<T, TR>::cal_nac(int, int, int)
{
    throw std::logic_error("Analytic auxiliary TDA NAC is implemented for real Gamma orbitals only");
}

template <> LR::NacvResult ESolver_LR<double, double>::cal_nac(int bra_state, int ket_state, int spin_channel)
{
    const NacOrchestrationTimer timer;
    if (!this->inp_ || this->inp_->esolver_type != "ks-lr" || !this->inp_->gamma_only || !this->inp_->cal_force ||
        !this->psi_ks_z_ || !this->pelec || !this->ucell_ || !this->conv_esolver)
    {
        throw std::logic_error("Analytic NAC requires converged KS/LR data and the force-response window");
    }
    if (this->nk != 1 || this->inp_->kpar != 1 || spin_channel < 0 || spin_channel >= this->nspin)
    {
        throw std::logic_error("Analytic NAC requires a real Gamma LR channel and kpar=1");
    }
    if (this->openshell && (spin_channel != 0 || this->inp_->lr_nac_spin != "updown"))
    {
        throw std::invalid_argument("Unrestricted analytic NAC requires lr_nac_spin=updown");
    }
    if (!this->openshell && this->inp_->lr_nac_spin == "updown")
    {
        throw std::invalid_argument("lr_nac_spin=updown requires an unrestricted LR reference");
    }
    const double absolute_magnetization = this->ucell_->magnet.abs_mag;
    if (!std::isfinite(absolute_magnetization) || absolute_magnetization < 0.0)
    {
        throw std::invalid_argument("Analytic NAC requires finite reference magnetization");
    }
    if (!this->openshell && absolute_magnetization > 1e-6)
    {
        throw std::invalid_argument(
            "Analytic NAC requires a closed-shell reference with absolute magnetization <= 1e-6");
    }
    if (bra_state < 0 || ket_state < 0 || bra_state > this->nstates || ket_state > this->nstates ||
        bra_state == ket_state)
    {
        throw std::out_of_range("Analytic NAC needs two distinct surfaces inside the solved state window");
    }
    const bool reverse = bra_state > ket_state;
    const int lower = std::min(bra_state, ket_state);
    const int upper = std::max(bra_state, ket_state);
    const bool ground_bra = lower == 0;
    if (this->openshell)
    {
        for (int spin = 0; spin < 2; ++spin)
        {
            LR::validate_nac_window(this->nocc[spin], this->nvirt[spin], this->nbasis);
        }
        LR::validate_nac_occupations(this->wg_ks_all, this->nocc, this->nbasis);
    }
    else
    {
        LR::validate_nac_window(this->nocc[spin_channel], this->nvirt[spin_channel], this->nbasis);
        LR::validate_nac_occupations(
            this->wg_ks_all,
            this->nocc[spin_channel],
            this->nbasis,
            this->nspin);
    }
    const std::vector<double> &residuals = this->lr_residuals_.at(spin_channel);
    if (residuals.size() != static_cast<std::size_t>(this->nstates))
    {
        throw std::logic_error("Analytic NAC requires certified LR eigenpair residuals");
    }
    const int ket_index = upper - 1;
    const ct::Tensor ket = this->pad_X_to_z_(spin_channel, ket_index, 1);
    ct::Tensor bra;
    const double *bra_pointer = nullptr;
    double bra_energy = 0.0;
    double bra_residual = 0.0;
    if (!ground_bra)
    {
        const int bra_index = lower - 1;
        bra = this->pad_X_to_z_(spin_channel, bra_index, 1);
        bra_pointer = bra.data<double>();
        bra_energy = this->pelec->ekb(spin_channel, bra_index);
        bra_residual = residuals.at(bra_index);
    }
    const double ket_energy = this->pelec->ekb(spin_channel, ket_index);
    const double gap_resolution = LR::nacv_gap_resolution(
        bra_energy,
        ket_energy,
        bra_residual,
        residuals.at(ket_index));
    const LR::GradientInputs<double> inputs = this->gradient_inputs_();
    const module_dm::DensityMatrix<double, double> ground_density = this->cal_dm_gs();
    const LR::ForceSettings settings{this->nspin,
        this->inp_->test_force,
        this->inp_->vl_in_h,
        this->inp_->vh_in_h,
        this->inp_->dft_functional};
    LR::LR_Force<double> force_terms(
        *this->ucell_,
        this->kv.kvec_d,
        this->paraMat_,
        *this->pw_rhod,
        *this->pw_rho,
        this->vloc(),
        this->sfac(),
        this->gd(),
        this->tcb(),
        settings,
        this->ofs_running_
#ifdef __EXX
        ,
        inputs.exx_lri,
        inputs.hybrid_alpha
#endif
    );
    LR::NacvResult result = LR::evaluate_aux_tda_nac(
        inputs,
        force_terms,
        ground_density,
        this->tcb(),
        spin_channel,
        ground_bra,
        bra_pointer,
        ket.data<double>(),
        bra_energy,
        ket_energy,
        gap_resolution,
        this->inp_->lr_grad_solver);
    if (reverse)
    {
        LR::reverse_nac_pair(result);
    }
    return result;
}

template LR::NacvResult ESolver_LR<std::complex<double>, double>::cal_nac(int, int, int);
} // namespace ModuleESolver
