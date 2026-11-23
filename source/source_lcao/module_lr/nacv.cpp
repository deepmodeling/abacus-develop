/** @file nacv.cpp
 * @brief Evaluate geometry-local NAC batches with reusable response and force operators.
 */
#include "nacv.h"
#include "nacv_ao.h"
#include "nacv_response.h"
#include "cal_edm.h"
#include "gradient_inputs.h"
#include "hamilt_casida.h"
#include "hamilt_ulr.h"
#include "lr_amp.h"
#include "lr_force.h"
#include "nacv_tools.h"
#include "source_base/parallel_reduce.h"
#include "source_base/timer.h"
#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdint>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace
{
/** @brief Measure the complete analytic evaluation, including spin-zero and rejected-input paths. */
class NacvTimer
{
  public:
    NacvTimer() { ModuleBase::timer::start("LR", "evaluate_aux_tda_nac"); }
    ~NacvTimer() { ModuleBase::timer::end("LR", "evaluate_aux_tda_nac"); }
};

/**
 * @brief Identify the explicit collinear unrestricted LR channel.
 * @param[in] inputs Borrowed LR channel labels and spin layout.
 * @return True for the single concatenated up/down amplitude channel.
 */
bool unrestricted_nac(const LR::GradientInputs<double> &inputs)
{
    return inputs.spin_types.size() == 1 && inputs.spin_types[0] == "updown";
}

/**
 * @brief Solve the GS A+B adjoint for a reference-to-excited auxiliary NAC.
 * @param[in] inputs Validated singlet geometry, orbitals, kernels and distributions.
 * @param[in] ket Local TDA excitation amplitudes in the full response window.
 * @param[in] energy Positive ket excitation energy in Ry.
 * @param[in,out] context Geometry-local orbital-response operators and direct factors.
 * @param[in] open_shell True for the concatenated alpha/beta orbital response.
 * @return Owned dimensionless local adjoint tensor of shape1-by-inputs.nloc.
 * @details The RHS is energy*ket for unrestricted spin orbitals and
 * sqrt(2)*energy*ket for a closed-shell singlet, in Ry. The existing
 * response solver checks its residual before the contraction uses the result.
 */
ct::Tensor ground_adjoint(
    const LR::GradientInputs<double> &inputs,
    const double *ket,
    double energy,
    bool open_shell,
    LR::NacvResponse &context)
{
    const int ld = inputs.nloc;
    ct::Tensor rhs = LR_Util::newTensor<double>({1, ld});
    ct::Tensor solution = LR_Util::newTensor<double>({1, ld});
    double scale = std::sqrt(2.0) * energy;
    if (open_shell)
    {
        scale = energy;
    }
    for (int index = 0; index < ld; ++index)
    {
        rhs.data<double>()[index] = scale * ket[index];
    }
    context.solve(rhs.data<double>(), solution.data<double>(), 1);
    return solution;
}

} // namespace

namespace LR
{
/**
 * @brief Evaluate a real-Gamma auxiliary-TDA analytic derivative coupling.
 * @param[in] inputs Explicit ground-state geometry, orbital layouts, kernel and log.
 * @param[in,out] force_terms Existing analytic force contractions for this geometry.
 * @param[in] ground_density Ground density in the force convention, two channels for unrestricted orbitals.
 * @param[in] spin_channel Closed-shell adaptation index, or zero for the concatenated up/down channel.
 * @param[in] bra_ground True for S0-to-excited coupling; otherwise an excited pair.
 * @param[in] bra Local widened bra amplitudes; ignored for S0. May be null only
 * when bra_ground or inputs.nloc is zero.
 * @param[in] ket Local widened ket amplitudes in the same canonical MO gauge.
 * @param[in] bra_energy_ry Excitation energy of the bra, zero for S0, in Rydberg.
 * @param[in] ket_energy_ry Positive ket excitation energy in Rydberg.
 * @param[in] gap_resolution_ry Nonnegative eigenpair-resolution bound in Rydberg.
 * @param[in,out] response_context Geometry-local response operators and direct factors.
 * @pre The geometry-local NAC entry has validated the reference and selected roots.
 * @return Owned raw/ETF vectors in bohr^-1, numerators in Ry/bohr and diagnostics.
 * @throw std::invalid_argument For unsupported geometry, spin, energies or layouts.
 * @details Implements direct differentiation of the auxiliary determinant expansion
 * (Ou et al., JCP 141, 024114, 2014). EE uses common-average-energy quadratic
 * polarization and a separate antisymmetric AO connection; GE uses the GS A+B
 * adjoint plus its symmetric overlap correction. Independent overlap-FD validation
 * is required before treating a newly supported functional as production-ready.
 */
NacvResult evaluate_aux_tda_nac(
    const GradientInputs<double> &inputs,
    LR_Force<double> &force_terms,
    const module_dm::DensityMatrix<double, double> &ground_density,
    int spin_channel,
    bool bra_ground,
    const double *bra,
    const double *ket,
    double bra_energy_ry,
    double ket_energy_ry,
    double gap_resolution_ry,
    NacvResponse &response_context)
{
    const NacvTimer timer;
    NacvResult result;
    result.gap_rydberg = ket_energy_ry - bra_energy_ry;
    result.raw.create(inputs.ucell.nat, 3);
    result.etf.create(inputs.ucell.nat, 3);
    result.ao_connection.create(inputs.ucell.nat, 3);
    result.numerator_raw.create(inputs.ucell.nat, 3);
    result.numerator_etf.create(inputs.ucell.nat, 3);
    const bool open_shell = unrestricted_nac(inputs);
    if (bra_ground && !open_shell && inputs.spin_types[spin_channel] == "triplet")
    {
        result.gap_resolved = std::abs(result.gap_rydberg) > gap_resolution_ry;
        return result;
    }
    const bool emit_state_forces = false;
    ModuleBase::matrix pair_force(inputs.ucell.nat, 3);
    ct::Tensor density;
    if (bra_ground)
    {
        density = nac_density_ao(
            inputs,
            true,
            open_shell,
            bra,
            ket,
            spin_channel);
        const ct::Tensor response = ground_adjoint(
            inputs,
            ket,
            ket_energy_ry,
            open_shell,
            response_context);
        force_terms.prepare_force_batch(ground_density);
        ct::Tensor zero_amplitudes = LR_Util::newTensor<double>({1, inputs.nloc});
        zero_amplitudes.zero();
        const std::vector<double> energy{ket_energy_ry};
        std::vector<ModuleBase::matrix> linear_force;
        if (open_shell)
        {
            linear_force = evaluate_open_shell_force(
                inputs,
                force_terms,
                ground_density,
                zero_amplitudes,
                response,
                energy,
                0,
                emit_state_forces);
        }
        else
        {
            linear_force = evaluate_closed_shell_force(
                inputs,
                force_terms,
                ground_density,
                zero_amplitudes,
                response,
                energy,
                0,
                spin_channel,
                emit_state_forces);
        }
        pair_force = linear_force[0];
        ct::Tensor correction =
            LR_Util::newTensor<double>({inputs.pmat.para_2d().ncloc(), inputs.pmat.para_2d().nrloc()});
        const double correction_scale = 0.5 * ket_energy_ry;
        for (std::int64_t index = 0; index < inputs.pmat.para_2d().size(); ++index)
        {
            correction.data<double>()[index] = correction_scale * density.data<double>()[index];
        }
        const module_dm::DensityMatrix<double, double> edm = sparse_density(inputs, correction, true);
        pair_force += force_terms.cal_force_overlap_edm(edm);
    }
    else
    {
        density = nac_density_ao(
            inputs,
            false,
            open_shell,
            bra,
            ket,
            spin_channel);
        ct::Tensor combinations = LR_Util::newTensor<double>({2, inputs.nloc});
        const double inverse_sqrt_two = 1.0 / std::sqrt(2.0);
        for (int index = 0; index < inputs.nloc; ++index)
        {
            combinations.data<double>()[index] = (bra[index] + ket[index]) * inverse_sqrt_two;
            combinations.data<double>()[inputs.nloc + index] = (bra[index] - ket[index]) * inverse_sqrt_two;
        }
        ct::Tensor response = LR_Util::newTensor<double>({2, inputs.nloc});
        const int pair_count = 2;
        response_context.solve_polarized(combinations.data<double>(), response.data<double>(), pair_count);
        force_terms.prepare_force_batch(ground_density);
        // Both polarized amplitudes have the same Rayleigh quotient. Using it in
        // both forces gives the symmetric metric term at the mean energy.
        // The negative force half-difference yields the ETF numerator; the AO connection is added separately.
        const double mean_energy = 0.5 * (bra_energy_ry + ket_energy_ry);
        const std::vector<double> energies{mean_energy, mean_energy};
        std::vector<ModuleBase::matrix> forces;
        if (open_shell)
        {
            forces = evaluate_open_shell_force(
                inputs,
                force_terms,
                ground_density,
                combinations,
                response,
                energies,
                0,
                emit_state_forces);
        }
        else
        {
            forces = evaluate_closed_shell_force(
                inputs,
                force_terms,
                ground_density,
                combinations,
                response,
                energies,
                0,
                spin_channel,
                emit_state_forces);
        }
        pair_force = (forces[0] - forces[1]) * 0.5;
    }
    const auto &overlap_derivatives = force_terms.batch_overlap_derivatives();
    result.ao_connection = antisymmetric_ao_connection(inputs, density, overlap_derivatives);
    assemble_nacv_result(result, pair_force, gap_resolution_ry);
    return result;
}

} // namespace LR

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
template <typename T, typename THamilt>
std::vector<double> nac_eigenpair_residuals(
    const THamilt &,
    const T *,
    const double *,
    int,
    int);
NacvResult evaluate_nac_pair(
    const NacvInputs<double> &,
    int,
    int,
    int,
    LR_Force<double> &,
    NacvResponse &);
void validate_nac_data(const NacvInputs<double> &);
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
    std::vector<double> residuals = nac_eigenpair_residuals(
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
    validate_nac_data(data);
    const auto &request = data.request;
    const int channel = request.spin_name == "triplet" ? 1 : 0;
    const ForceSettings settings{data.gradient.nspin,
        data.gradient.test_force,
        request.local_potential,
        request.hartree_potential,
        data.gradient.dft_functional};
    LR_Force<double> force_terms(
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
        data.gradient.exx_lri,
        data.gradient.hybrid_alpha
#endif
    );
    const NacvResponseControls controls{request.response_tolerance, request.response_max_iterations};
    NacvResponse response_context(
        data.gradient,
        channel,
        data.open_shell,
        request.response_solver,
        controls);

    for (std::size_t pair = 0; pair < request.pairs.size(); pair += 2)
    {
        const int bra = request.pairs[pair];
        const int ket = request.pairs[pair + 1];
        const NacvResult coupling = evaluate_nac_pair(
            data,
            bra,
            ket,
            channel,
            force_terms,
            response_context);
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

namespace
{

/**
 * @brief Reapply the LR operator and measure global eigenpair residual norms.
 * @param[in] hm Distributed operator for the current geometry.
 * @param[in] amplitudes State-major local eigenvectors.
 * @param[in] energies Excitation energies in Ry.
 * @param[in] local_size Local entries per root, possibly zero on empty ranks.
 * @param[in] state_count Number of solved roots.
 * @return One global normalized Euclidean residual per root in Ry.
 */
template <typename T, typename THamilt>
std::vector<double> nac_eigenpair_residuals(
    const THamilt &hm,
    const T *amplitudes,
    const double *energies,
    int local_size,
    int state_count)
{
    ModuleBase::timer::start("LR", "nac_eigenpair_residuals");
    const std::size_t total_size = static_cast<std::size_t>(local_size) * state_count;
    const std::size_t storage_size = std::max<std::size_t>(1, total_size);
    const T zero = T(0);
    std::vector<T> applied(storage_size, zero);
    hm.hPsi(
        amplitudes,
        applied.data(),
        local_size,
        state_count);
    std::vector<double> residuals(state_count);
    for (int state = 0; state < state_count; ++state)
    {
        double residual_square = 0.0;
        double norm_square = 0.0;
        for (int row = 0; row < local_size; ++row)
        {
            const std::size_t index = static_cast<std::size_t>(state) * local_size + row;
            const T response_residual = applied[index] - energies[state] * amplitudes[index];
            const double residual = std::abs(response_residual);
            const double magnitude = std::abs(amplitudes[index]);
            residual_square += residual * residual;
            norm_square += magnitude * magnitude;
        }
        Parallel_Reduce::reduce_all(residual_square);
        Parallel_Reduce::reduce_all(norm_square);
        if (!std::isfinite(residual_square) || !std::isfinite(norm_square) || norm_square <= 0.0)
        {
            ModuleBase::timer::end("LR", "nac_eigenpair_residuals");
            throw std::runtime_error("Invalid LR eigenpair residual for analytic NAC");
        }
        const double normalized_square = residual_square / norm_square;
        residuals[state] = std::sqrt(normalized_square);
    }
    ModuleBase::timer::end("LR", "nac_eigenpair_residuals");
    return residuals;
}

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
    int spin_channel,
    LR_Force<double> &force_terms,
    NacvResponse &response_context)
{
    const bool reverse = bra_state > ket_state;
    const int lower = std::min(bra_state, ket_state);
    const int upper = std::max(bra_state, ket_state);
    const bool ground_bra = lower == 0;
    const std::vector<double> &residuals = data.residuals.at(spin_channel);
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
    LR::NacvResult result = LR::evaluate_aux_tda_nac(
        inputs,
        force_terms,
        ground_density,
        spin_channel,
        ground_bra,
        bra_pointer,
        ket.data<double>(),
        bra_energy,
        ket_energy,
        gap_resolution,
        response_context);
    if (reverse)
    {
        LR::reverse_nac_pair(result);
    }
    return result;
}

/**
 * @brief Validate the computed reference once at the NAC batch boundary.
 * @param[in] data Borrowed, solver-owned geometry and converged LR state data.
 * @details INPUT checked mode flags, pair structure and CG controls. The upper
 * surface bound is deferred until the solver resolves its root count.
 * Occupations, complete orbital windows and positive requested excitation
 * energies can only be checked after the actual KS/LR reference is available.
 */
void validate_nac_data(const NacvInputs<double> &data)
{
    const auto &inputs = data.gradient;
    const auto &request = data.request;
    const int channel = request.spin_name == "triplet" ? 1 : 0;
    validate_nac_surfaces(request.pairs, data.energies.nc);
    if (data.open_shell)
    {
        if (request.spin_name != "updown")
        {
            throw std::invalid_argument("Unrestricted analytic NAC requires lr_nac_spin=updown");
        }
        for (int spin = 0; spin < 2; ++spin)
        {
            validate_nac_window(inputs.nocc[spin], data.virtuals[spin], inputs.nbasis);
        }
        validate_nac_occupations(data.occupations, inputs.nocc, inputs.nbasis);
    }
    else
    {
        if (request.spin_name == "updown" || inputs.ucell.magnet.abs_mag > 1e-6)
        {
            throw std::invalid_argument("Spin-adapted analytic NAC requires a closed-shell reference");
        }
        validate_nac_window(inputs.nocc[channel], data.virtuals[channel], inputs.nbasis);
        validate_nac_occupations(
            data.occupations,
            inputs.nocc[channel],
            inputs.nbasis,
            inputs.nspin);
    }
    std::string reference_xc = inputs.dft_functional;
    if (reference_xc == "default")
    {
        reference_xc = inputs.ucell.atoms[0].ncpp.xc_func;
    }
    validate_nac_functionals(reference_xc, inputs.xc_kernel);
    for (const int state : request.pairs)
    {
        if (state > 0 && !(data.energies(channel, state - 1) > 0.0))
        {
            throw std::invalid_argument("Analytic NAC requires positive requested excitation energies");
        }
    }
}
} // namespace

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
