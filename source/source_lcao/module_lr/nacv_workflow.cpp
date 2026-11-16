#include "nacv_workflow.h"

#include "gradient_inputs.h"
#include "hamilt_casida.h"
#include "hamilt_ulr.h"
#include "hsolver_lrtd.h"
#include "lr_amp.h"
#include "nacv_cs.h"
#include "source_base/timer.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <stdexcept>

namespace LR
{
namespace
{
/** @brief Scope the full geometry-local NAC batch, including exception paths. */
class NacWorkflowTimer
{
  public:
    NacWorkflowTimer() { ModuleBase::timer::start("LR", "run_nac_requests"); }
    ~NacWorkflowTimer() { ModuleBase::timer::end("LR", "run_nac_requests"); }
    NacWorkflowTimer(const NacWorkflowTimer &) = delete;
    NacWorkflowTimer &operator=(const NacWorkflowTimer &) = delete;
};

/**
 * @brief Widen one root using the existing shared amplitude-distribution helper.
 * @param[in] data Borrowed original and response-window distributions.
 * @param[in] channel Spin adaptation index.
 * @param[in] state Zero-based root index.
 * @return Owned response-window amplitude row in the original orbital gauge.
 */
ct::Tensor pad_nac_state(const NacvInputs<double> &data, int channel, int state)
{
    const auto &inputs = data.gradient;
    return pad_amplitudes<double>(
        data.amplitudes,
        data.open_shell,
        data.amplitude_layouts,
        inputs.px,
        inputs.nocc,
        data.virtuals,
        inputs.nk,
        data.local_amplitude_size,
        inputs.nloc,
        channel,
        state,
        1);
}

/**
 * @brief Preserve the existing single-pair preparation and force contractions.
 * @param[in] data Borrowed geometry-local batch data.
 * @param[in] bra_state Surface index; zero denotes the ground reference.
 * @param[in] ket_state Distinct surface index inside the solved window.
 * @param[in] spin_channel Adaptation index, or zero for unrestricted LR.
 * @return Owned raw/ETF vectors in bohr^-1 and numerators in Ry/bohr.
 */
NacvResult evaluate_nac_pair(
    const NacvInputs<double> &data,
    int bra_state,
    int ket_state,
    int spin_channel)
{

    const double absolute_magnetization = data.gradient.ucell.magnet.abs_mag;
    if (!std::isfinite(absolute_magnetization) || absolute_magnetization < 0.0)
    {
        throw std::invalid_argument("Analytic NAC requires finite reference magnetization");
    }
    if (!data.open_shell && absolute_magnetization > 1e-6)
    {
        throw std::invalid_argument(
            "Analytic NAC requires a closed-shell reference with absolute magnetization <= 1e-6");
    }
    if (bra_state < 0 || ket_state < 0 || bra_state > data.energies.nc || ket_state > data.energies.nc ||
        bra_state == ket_state)
    {
        throw std::out_of_range("Analytic NAC needs two distinct surfaces inside the solved state window");
    }
    const bool reverse = bra_state > ket_state;
    const int lower = std::min(bra_state, ket_state);
    const int upper = std::max(bra_state, ket_state);
    const bool ground_bra = lower == 0;
    const std::vector<double> &residuals = data.residuals.at(spin_channel);
    if (residuals.size() != static_cast<std::size_t>(data.energies.nc))
    {
        throw std::logic_error("Analytic NAC requires certified LR eigenpair residuals");
    }
    const int ket_index = upper - 1;
    const ct::Tensor ket = pad_nac_state(data, spin_channel, ket_index);
    ct::Tensor bra;
    const double *bra_pointer = nullptr;
    double bra_energy = 0.0;
    double bra_residual = 0.0;
    if (!ground_bra)
    {
        const int bra_index = lower - 1;
        bra = pad_nac_state(data, spin_channel, bra_index);
        bra_pointer = bra.data<double>();
        bra_energy = data.energies(spin_channel, bra_index);
        bra_residual = residuals.at(bra_index);
    }
    const double ket_energy = data.energies(spin_channel, ket_index);
    const double gap_resolution = LR::nacv_gap_resolution(
        bra_energy,
        ket_energy,
        bra_residual,
        residuals.at(ket_index));
    const LR::GradientInputs<double> &inputs = data.gradient;
    const auto &ground_density = data.ground_density;
    const LR::ForceSettings settings{data.gradient.nspin,
        data.gradient.test_force,
        data.request.local_potential,
        data.request.hartree_potential,
        data.gradient.dft_functional};
    LR::LR_Force<double> force_terms(
        data.gradient.ucell,
        data.gradient.kv.kvec_d,
        data.gradient.pmat,
        data.density_basis,
        data.smooth_basis,
        data.local_potential,
        data.structure_factor,
        data.gradient.gd,
        data.two_center_bundle,
        settings,
        data.gradient.ofs
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
        data.two_center_bundle,
        spin_channel,
        ground_bra,
        bra_pointer,
        ket.data<double>(),
        bra_energy,
        ket_energy,
        gap_resolution,
        data.request.response_solver);
    if (reverse)
    {
        LR::reverse_nac_pair(result);
    }
    return result;
}
} // namespace

template <typename T, typename THamilt>
std::vector<double> collect_nac_eigenpairs(
    const THamilt &hamiltonian,
    const T *amplitudes,
    const double *energies,
    int local_size,
    int state_count,
    bool requested,
    int ionic_step,
    const std::string &spin_name,
    int rank,
    std::ostream &output)
{
    if (!requested)
    {
        return {};
    }
    std::vector<double> residuals = HSolver::eigenpair_residuals(
        hamiltonian,
        amplitudes,
        energies,
        local_size,
        state_count);
    if (rank == 0)
    {
        print_nac_eigenpairs(
            output,
            ionic_step,
            spin_name,
            energies,
            residuals);
    }
    return residuals;
}

template <typename T> void run_nac_requests(const NacvInputs<T> &data, int)
{
    if (!data.request.pairs.empty())
    {
        throw std::logic_error("Analytic auxiliary TDA NAC is implemented for real Gamma orbitals only");
    }
}

template <> void run_nac_requests<double>(const NacvInputs<double> &data, int ionic_step)
{
    const NacWorkflowTimer timer;
    const auto &request = data.request;
    const int channel = request.spin_name == "triplet" ? 1 : 0;
    if (data.gradient.nk != 1 || channel < 0 || channel >= data.gradient.nspin)
    {
        throw std::logic_error("Analytic NAC requires a real Gamma LR channel and kpar=1");
    }
    if (data.open_shell && (channel != 0 || request.spin_name != "updown"))
    {
        throw std::invalid_argument("Unrestricted analytic NAC requires lr_nac_spin=updown");
    }
    if (!data.open_shell && request.spin_name == "updown")
    {
        throw std::invalid_argument("lr_nac_spin=updown requires an unrestricted LR reference");
    }
    if (data.open_shell)
    {
        for (int spin = 0; spin < 2; ++spin)
        {
            validate_nac_window(data.gradient.nocc[spin], data.virtuals[spin], data.gradient.nbasis);
        }
    }
    else
    {
        validate_nac_window(data.gradient.nocc[channel], data.virtuals[channel], data.gradient.nbasis);
    }
    for (std::size_t pair = 0; pair < request.pairs.size(); pair += 2)
    {
        const int bra = request.pairs[pair];
        const int ket = request.pairs[pair + 1];
        const NacvResult coupling = evaluate_nac_pair(
            data,
            bra,
            ket,
            channel);
        if (data.gradient.my_rank == 0)
        {
            print_nac_result(
                data.gradient.ofs,
                coupling,
                ionic_step,
                bra,
                ket,
                request.spin_name);
        }
    }
}

template void run_nac_requests<std::complex<double>>(const NacvInputs<std::complex<double>> &, int);

template std::vector<double> collect_nac_eigenpairs<double, HamiltLR<double>>(
    const HamiltLR<double> &,
    const double *,
    const double *,
    int,
    int,
    bool,
    int,
    const std::string &,
    int,
    std::ostream &);

template std::vector<double> collect_nac_eigenpairs<double, HamiltULR<double>>(
    const HamiltULR<double> &,
    const double *,
    const double *,
    int,
    int,
    bool,
    int,
    const std::string &,
    int,
    std::ostream &);

template std::vector<double> collect_nac_eigenpairs<std::complex<double>, HamiltLR<std::complex<double>>>(
    const HamiltLR<std::complex<double>> &,
    const std::complex<double> *,
    const double *,
    int,
    int,
    bool,
    int,
    const std::string &,
    int,
    std::ostream &);

template std::vector<double> collect_nac_eigenpairs<std::complex<double>, HamiltULR<std::complex<double>>>(
    const HamiltULR<std::complex<double>> &,
    const std::complex<double> *,
    const double *,
    int,
    int,
    bool,
    int,
    const std::string &,
    int,
    std::ostream &);

} // namespace LR
