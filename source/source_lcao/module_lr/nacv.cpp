/** @file nacv.cpp
 * @brief Evaluate geometry-local NAC batches with reusable response and force operators.
 */
#include "nacv.h"
#include "cal_edm.h"
#include "cal_hs_grad.h"
#include "gradient_inputs.h"
#include "hamilt_casida.h"
#include "hamilt_ulr.h"
#include "hamilt_zeq_l.h"
#include "hamilt_zeq_r.h"
#include "lr_amp.h"
#include "lr_force.h"
#include "nacv_tools.h"
#include "pulay_hc.h"
#include "source_base/parallel_reduce.h"
#include "source_base/timer.h"
#include "zeq_solver.h"
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

namespace LR
{
template <typename T> struct GradientInputs;
/** @brief INPUT-validated CG controls, with an absolute residual tolerance in Ry. */
struct NacvResponseControls
{
    double tolerance;
    int max_iterations;
};
/**
 * @brief Reuse orbital-response operators and direct factors within one NAC geometry.
 * @details The object is owned by a single run_nac_requests call and borrows
 * its complete orbital window and kernels. GS response and LR RHS kernels keep
 * their distinct XC purposes; mutable operator density buffers are used sequentially.
 */
class NacvResponse
{
  public:
    /**
     * @param[in] inputs Borrowed converged reference; must outlive this object.
     * @param[in] channel Selected singlet/triplet or unrestricted channel.
     * @param[in] open_shell Whether amplitudes concatenate alpha and beta blocks.
     * @param[in] solver Existing NAC response solver name.
     * @param[in] controls Existing accuracy settings, unchanged.
     */
    NacvResponse(
        const GradientInputs<double> &inputs,
        int channel,
        bool open_shell,
        const std::string &solver,
        const NacvResponseControls &controls);
    ~NacvResponse();
    /**
     * @param[in,out] rhs Packed local RHS in Ry, preserved by the solver.
     * @param[out] solution Dimensionless local Z, with count rows of inputs.nloc.
     * @param[in] count Number of RHS rows.
     */
    void solve(double *rhs, double *solution, int count);
    /**
     * @param[in] amplitudes Local state-major polarized TDA vectors.
     * @param[out] solution Corresponding Z rows in the same layout.
     * @param[in] count Number of polarized vectors.
     */
    void solve_polarized(const double *amplitudes, double *solution, int count);

  private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};
} // namespace LR

namespace LR
{
template <typename T> struct NacvInputs;
/**
 * @brief Own geometry-fixed force intermediates for one NAC request batch.
 * @details The reference and parallel layouts are borrowed for the batch lifetime.
 * Potentials and derivative matrices are built lazily, including on empty MPI ranks.
 */
class NacvForce final : public LR_Force<double>
{
  public:
    /**
     * @param[in] data Geometry-local reference, density and integration dependencies.
     * @param[in] settings Resolved force switches, with no global configuration reads.
     */
    NacvForce(const NacvInputs<double> &data, const ForceSettings &settings);
    /** @return Immutable geometry intermediates, prepared on their first use. */
    const ForceGeometry &geometry();

  private:
    const NacvInputs<double> &data_;
    std::vector<hamilt::HContainer<double>> overlap_;
    std::vector<hamilt::HContainer<double>> kinetic_;
    std::unique_ptr<elecstate::Potential> local_;
    std::unique_ptr<elecstate::Potential> ground_;
    std::unique_ptr<ForceGeometry> geometry_;
    const ForceGeometry *geometry_terms() override;
    void prepare_geometry();
};
} // namespace LR

namespace LR
{
NacvForce::NacvForce(const NacvInputs<double> &data, const ForceSettings &settings)
    : LR_Force<double>(
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
          ),
      data_(data)
{
}
const ForceGeometry &NacvForce::geometry() { return *geometry_terms(); }
const ForceGeometry *NacvForce::geometry_terms()
{
    if (!geometry_)
    {
        prepare_geometry();
    }
    return geometry_.get();
}
void NacvForce::prepare_geometry()
{
    const NacvTimerScope timer("NacvForce", "prepare_geometry");
    const auto &g = data_.gradient;
    overlap_ = cal_hs_grad(
        'S',
        g.ucell,
        g.pmat,
        g.gd,
        data_.two_center_bundle);
    kinetic_ = cal_hs_grad(
        'T',
        g.ucell,
        g.pmat,
        g.gd,
        data_.two_center_bundle);
    local_ = local_potential();
    ground_ = dm_to_hxc_potential(data_.ground_density);
    geometry_.reset(new ForceGeometry{overlap_,
        kinetic_,
        *local_,
        *ground_});
}
} // namespace LR

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
 * @brief Transform the particle-minus-hole interstate density to the AO layout.
 * @param[in] inputs Borrowed geometry, orbitals and distributed matrix layouts.
 * @param[in] spin_channel Channel supplying the occupied and virtual windows.
 * @param[in] orbital_spin Spin slice supplying the canonical MO coefficients.
 * @param[in] bra Local bra TDA amplitudes in inputs.px[spin_channel].
 * @param[in] ket Local ket amplitudes in the same canonical MO gauge and layout.
 * @return Owned AO coefficient tensor in inputs.pmat's local layout.
 * @pre The two states are distinct orthogonal roots; the reference-overlap
 * density is omitted. Inputs have passed the geometry-local NAC batch entry.
 * @details Uses D(p,q)=<bra|a_q^dagger a_p|ket>, with the spin sum implicit.
 */
ct::Tensor pair_density_ao(
    const LR::GradientInputs<double> &inputs,
    int spin_channel,
    int orbital_spin,
    const double *bra,
    const double *ket)
{
    const Parallel_Blkcyclic &px = inputs.px[spin_channel];
    const int nocc = inputs.nocc[spin_channel];
    const int nvirt = inputs.nvirt[spin_channel];
    const int naos = inputs.nbasis;
    const psi::Psi<double> orbitals = LR_Util::get_psi_spin(inputs.psi_ks, orbital_spin, inputs.nk);
    orbitals.fix_k(0);
    Parallel_Blkcyclic particle_layout;
    Parallel_Blkcyclic hole_layout;
    LR_Util::setup_2d_division(
        particle_layout,
        px.nb(),
        naos,
        nocc
#ifdef __MPI
        ,
        px.ictxt()
#endif
    );
    LR_Util::setup_2d_division(
        hole_layout,
        px.nb(),
        naos,
        nvirt
#ifdef __MPI
        ,
        px.ictxt()
#endif
    );
    const std::size_t particle_size = std::max<std::int64_t>(1, particle_layout.size());
    const std::size_t hole_size = std::max<std::int64_t>(1, hole_layout.size());
    const std::size_t ao_size = std::max<std::int64_t>(1, inputs.pmat.para_2d().size());
    std::vector<double> particle_bra(particle_size);
    std::vector<double> particle_ket(particle_size);
    std::vector<double> hole_bra(hole_size);
    std::vector<double> hole_ket(hole_size);
    std::vector<double> ao(ao_size);
    std::vector<double> hole(ao_size);
    double empty = 0.0;
    const double *bra_data = bra;
    const double *ket_data = ket;
    if (px.size() == 0)
    {
        bra_data = &empty;
        ket_data = &empty;
    }
#ifdef __MPI
    LR::CvX(
        orbitals.get_pointer(),
        inputs.pc,
        bra_data,
        px,
        naos,
        nocc,
        nvirt,
        particle_bra.data(),
        particle_layout);
    LR::CvX(
        orbitals.get_pointer(),
        inputs.pc,
        ket_data,
        px,
        naos,
        nocc,
        nvirt,
        particle_ket.data(),
        particle_layout);
    LR::CoXT(
        orbitals.get_pointer(),
        inputs.pc,
        bra_data,
        px,
        naos,
        nocc,
        nvirt,
        hole_bra.data(),
        hole_layout);
    LR::CoXT(
        orbitals.get_pointer(),
        inputs.pc,
        ket_data,
        px,
        naos,
        nocc,
        nvirt,
        hole_ket.data(),
        hole_layout);
#else
    LR::CvX(
        orbitals.get_pointer(),
        bra_data,
        naos,
        nocc,
        nvirt,
        particle_bra.data());
    LR::CvX(
        orbitals.get_pointer(),
        ket_data,
        naos,
        nocc,
        nvirt,
        particle_ket.data());
    LR::CoXT(
        orbitals.get_pointer(),
        bra_data,
        naos,
        nocc,
        nvirt,
        hole_bra.data());
    LR::CoXT(
        orbitals.get_pointer(),
        ket_data,
        naos,
        nocc,
        nvirt,
        hole_ket.data());
#endif
    // D(p,q)=<bra|a_q^dagger a_p|ket>: particle=ket*bra^T, hole=-bra*ket^T.
    LR::matdot(
        particle_ket.data(),
        particle_bra.data(),
        particle_layout,
        ao.data(),
        inputs.pmat.para_2d());
    LR::matdot(
        hole_bra.data(),
        hole_ket.data(),
        hole_layout,
        hole.data(),
        inputs.pmat.para_2d());
    ct::Tensor density = LR_Util::newTensor<double>({inputs.pmat.para_2d().ncloc(), inputs.pmat.para_2d().nrloc()});
    for (std::int64_t index = 0; index < inputs.pmat.para_2d().size(); ++index)
    {
        density.data<double>()[index] = ao[index] - hole[index];
    }
    return density;
}

/**
 * @brief Transform one spin contribution to the reference-to-ket AO density.
 * @param[in] inputs Validated real-Gamma orbital context.
 * @param[in] spin_channel Window and distributed amplitude layout.
 * @param[in] orbital_spin Spin slice of the canonical orbitals.
 * @param[in] spin_factor One for a spin orbital; sqrt(2) for a closed-shell singlet.
 * @param[in] ket Local amplitudes in the selected occupied-major spin block.
 * @return Owned local AO coefficient tensor including the explicit spin factor.
 * @pre All ranks enter with the same validated spin context.
 */
ct::Tensor ground_density_ao(
    const LR::GradientInputs<double> &inputs,
    int spin_channel,
    int orbital_spin,
    const double *ket,
    double spin_factor)
{
    const psi::Psi<double> orbitals = LR_Util::get_psi_spin(inputs.psi_ks, orbital_spin, inputs.nk);
    const int nocc = inputs.nocc[spin_channel];
    const int nvirt = inputs.nvirt[spin_channel];
    double empty = 0.0;
    const double *amplitudes = ket;
    if (inputs.px[spin_channel].size() == 0)
    {
        amplitudes = &empty;
    }
#ifdef __MPI
    std::vector<ct::Tensor> density = LR::cal_dm_trans_pblas(
        amplitudes,
        inputs.px[spin_channel],
        orbitals,
        inputs.pc,
        inputs.nbasis,
        nocc,
        nvirt,
        inputs.pmat.para_2d());
#else
    std::vector<ct::Tensor> density = LR::cal_dm_trans_blas(
        amplitudes,
        orbitals,
        nocc,
        nvirt);
#endif
    for (std::int64_t index = 0; index < inputs.pmat.para_2d().size(); ++index)
    {
        density[0].data<double>()[index] *= spin_factor;
    }
    return density[0];
}

/**
 * @brief Sum the spin-resolved interstate AO densities for an auxiliary CIS pair.
 * @param[in] inputs Validated orbitals and complete spin windows.
 * @param[in] ground True for a reference-to-excited transition.
 * @param[in] bra Packed bra amplitudes; unused for the reference.
 * @param[in] ket Packed ket amplitudes, up block followed by down for UKS.
 * @param[in] spin_channel Closed-shell LR adaptation, or zero for UKS.
 * @return Owned spin-summed AO coefficients for the orbital connection.
 */
ct::Tensor nac_density_ao(
    const LR::GradientInputs<double> &inputs,
    bool ground,
    const double *bra,
    const double *ket,
    int spin_channel)
{
    const bool open_shell = unrestricted_nac(inputs);
    const int spin_count = open_shell ? 2 : 1;
    ct::Tensor total = LR_Util::newTensor<double>({inputs.pmat.para_2d().ncloc(), inputs.pmat.para_2d().nrloc()});
    total.zero();
    for (int spin = 0; spin < spin_count; ++spin)
    {
        int window = spin_channel;
        int orbital_spin = 0;
        std::int64_t offset = 0;
        if (open_shell)
        {
            window = spin;
            orbital_spin = spin;
            if (spin == 1)
            {
                offset = inputs.px[0].size();
            }
        }
        const double *bra_spin = nullptr;
        const double *ket_spin = nullptr;
        if (inputs.px[window].size() > 0)
        {
            ket_spin = ket + offset;
            if (!ground)
            {
                bra_spin = bra + offset;
            }
        }
        ct::Tensor contribution;
        if (ground)
        {
            double factor = 1.0;
            if (!open_shell)
            {
                factor = std::sqrt(2.0);
            }
            contribution = ground_density_ao(
                inputs,
                window,
                orbital_spin,
                ket_spin,
                factor);
        }
        else
        {
            contribution = pair_density_ao(
                inputs,
                window,
                orbital_spin,
                bra_spin,
                ket_spin);
        }
        for (std::int64_t index = 0; index < inputs.pmat.para_2d().size(); ++index)
        {
            total.data<double>()[index] += contribution.data<double>()[index];
        }
    }
    return total;
}

/**
 * @brief Construct sparse real-space density blocks for existing force contractions.
 * @param[in] inputs Borrowed geometry, neighbor lists and parallel AO layout.
 * @param[in] matrix Local AO coefficient tensor in inputs.pmat's layout.
 * @param[in] symmetrize True to apply the existing symmetric-density convention.
 * @return Density matrix with owned values and borrowed parallel-layout context.
 * @pre The referenced context remains valid while the returned matrix is used.
 */
module_dm::DensityMatrix<double, double> sparse_density(
    const LR::GradientInputs<double> &inputs, const ct::Tensor &matrix, bool symmetrize)
{
    const std::vector<ct::Tensor> blocks{matrix};
    return LR_Util::build_dm_from_dmk<double, double>(
        blocks,
        inputs.pmat,
        inputs.nk,
        inputs.kv.kvec_d,
        inputs.ucell,
        inputs.gd,
        inputs.orb_cutoff,
        symmetrize);
}

/**
 * @brief Contract the antisymmetric moving-AO connection for an auxiliary state pair.
 * @param[in] inputs Borrowed current geometry and validated parallel AO layout.
 * @param[in] density Local interstate AO coefficient tensor, including spin factors.
 * @param[in] geometry Borrowed overlap derivatives for this geometry.
 * @return Owned atom-major nat-by-3 connection in bohr^-1.
 * @details Contracts the antisymmetric density with the bra-center overlap
 * derivative; the symmetric metric-response term is handled separately.
 */
ModuleBase::matrix antisymmetric_ao_connection(
    const LR::GradientInputs<double> &inputs, const ct::Tensor &density, const LR::ForceGeometry &geometry)
{
    const std::int64_t size = inputs.pmat.para_2d().size();
    const std::int64_t storage_size = std::max<std::int64_t>(1, size);
    std::vector<double> transposed(storage_size, 0.0);
    if (size > 0)
    {
        std::copy_n(density.data<double>(), size, transposed.data());
    }
#ifdef __MPI
    LR_Util::mattrans(transposed.data(), inputs.nbasis, inputs.pmat.para_2d());
#else
    for (int column = 0; column < inputs.nbasis; ++column)
    {
        for (int row = 0; row < inputs.nbasis; ++row)
        {
            const std::size_t output = static_cast<std::size_t>(column) * inputs.nbasis + row;
            const std::size_t input = static_cast<std::size_t>(row) * inputs.nbasis + column;
            transposed[output] = density.data<double>()[input];
        }
    }
#endif
    ct::Tensor difference = LR_Util::newTensor<double>({inputs.pmat.para_2d().ncloc(), inputs.pmat.para_2d().nrloc()});
    for (std::int64_t index = 0; index < size; ++index)
    {
        difference.data<double>()[index] = density.data<double>()[index] - transposed[index];
    }
    const module_dm::DensityMatrix<double, double> dm = sparse_density(inputs, difference, false);
    // cal_hs_grad stores the separation derivative, minus the bra-center derivative.
    // The Pulay contraction contributes 2*factor*D*dHS. Thus -1/4 gives
    // Gamma : (B-B^T)/2 from D-D^T, with Gamma=D^T. No symmetric metric term is added here.
    const double contraction_factor = -0.25;
    return PulayForceStress::cal_pulay_fs(
        dm,
        inputs.ucell,
        geometry.overlap,
        contraction_factor);
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
    NacvForce &force_terms,
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
            bra,
            ket,
            spin_channel);
        const ct::Tensor response = ground_adjoint(
            inputs,
            ket,
            ket_energy_ry,
            open_shell,
            response_context);
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
    result.ao_connection = antisymmetric_ao_connection(inputs, density, force_terms.geometry());
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
    NacvForce &,
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
    NacvForce force_terms(data, settings);
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
    NacvForce &force_terms,
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
 * @details INPUT already checked mode flags, pair indices and CG controls.
 * Occupations, complete orbital windows and positive requested excitation
 * energies can only be checked after the actual KS/LR reference is available.
 */
void validate_nac_data(const NacvInputs<double> &data)
{
    const auto &inputs = data.gradient;
    const auto &request = data.request;
    const int channel = request.spin_name == "triplet" ? 1 : 0;
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

namespace LR
{
struct NacvResponse::Impl
{
    const GradientInputs<double> &inputs;
    const int channel;
    const bool open_shell;
    const std::string solver;
    const NacvResponseControls controls;
    const std::string input_directory;
    std::unique_ptr<Z_vector_L<double>> closed_left;
    std::unique_ptr<Z_vector_UL<double>> open_left;
    std::unique_ptr<Z_vector_R<double>> closed_right;
    std::unique_ptr<Z_vector_UR<double>> open_right;
    Parallel_Blkcyclic matrix_layout;
    std::unique_ptr<NacvFactorization> factorization;

    Impl(
        const GradientInputs<double> &data,
        int selected_channel,
        bool unrestricted,
        const std::string &selected_solver,
        const NacvResponseControls &settings)
        : inputs(data), channel(selected_channel), open_shell(unrestricted), solver(selected_solver), controls(settings)
    {
    }

    void prepare_left()
    {
        const NacvTimerScope timer("NacvResponse", "prepare_left");
        const auto &g = inputs;
        if (open_shell)
        {
            if (open_left)
            {
                return;
            }
            const int dimension = g.nk * (g.nocc[0] * g.nvirt[0] + g.nocc[1] * g.nvirt[1]);
            open_left.reset(new Z_vector_UL<double>(
                g.xc_kernel,
                g.nspin,
                g.nbasis,
                g.nocc,
                g.nvirt,
                g.ucell,
                g.orb_cutoff,
                g.gd,
                g.psi_ks,
                g.eig_ks,
#ifdef __EXX
                g.exx_lri,
                g.hybrid_alpha,
#endif
                g.pot_hxc_gs,
                g.kv,
                g.px,
                g.pc,
                g.pmat,
                g.dft_functional,
                g.nk,
                g.nloc,
                dimension));
        }
        else
        {
            if (closed_left)
            {
                return;
            }
            const auto &spin = g.spin_types[channel];
            closed_left.reset(new Z_vector_L<double>(
                g.xc_kernel,
                g.nspin,
                g.nbasis,
                g.nocc,
                g.nvirt,
                g.ucell,
                g.orb_cutoff,
                g.gd,
                g.psi_ks,
                g.eig_ks,
#ifdef __EXX
                g.exx_lri,
                g.hybrid_alpha,
#endif
                g.pot_hxc_gs,
                g.kv,
                g.px,
                g.pc,
                g.pmat,
                spin,
                input_directory,
                g.out_dir,
                g.dft_functional));
        }
    }

    void prepare_right()
    {
        const NacvTimerScope timer("NacvResponse", "prepare_right");
        const auto &g = inputs;
        const auto &spin = g.spin_types[channel];
        const std::weak_ptr<PotHxcLR> potential = g.pot[channel];
        if (open_shell)
        {
            if (open_right)
            {
                return;
            }
            const int dimension = g.nk * (g.nocc[0] * g.nvirt[0] + g.nocc[1] * g.nvirt[1]);
            open_right.reset(new Z_vector_UR<double>(
                g.xc_kernel,
                g.nspin,
                g.nbasis,
                g.nocc,
                g.nvirt,
                g.ucell,
                g.orb_cutoff,
                g.gd,
                g.psi_ks,
                g.eig_ks,
#ifdef __EXX
                g.exx_lri,
                g.hybrid_alpha,
#endif
                potential,
                g.pot_hxc_gs,
                g.kv,
                g.px,
                g.pc,
                g.pmat,
                g.ks_solver,
                g.dft_functional,
                g.nk,
                g.nloc,
                dimension));
        }
        else
        {
            if (closed_right)
            {
                return;
            }
            closed_right.reset(new Z_vector_R<double>(
                g.xc_kernel,
                g.nspin,
                g.nbasis,
                g.nocc,
                g.nvirt,
                g.ucell,
                g.orb_cutoff,
                g.gd,
                g.psi_ks,
                g.eig_ks,
#ifdef __EXX
                g.exx_lri,
                g.hybrid_alpha,
#endif
                potential,
                g.pot_hxc_gs,
                g.kv,
                g.px,
                g.pc,
                g.pmat,
                input_directory,
                g.out_dir,
                g.dft_functional,
                spin));
        }
    }

    template <typename Operator> void prepare_matrix(Operator &left, int spin_count)
    {
        const NacvTimerScope timer("NacvResponse", "prepare_matrix");
        const int dimension = zvec_global_dim(left, spin_count);
        if (solver == "lapack")
        {
            auto matrix = left.matrix();
            factorization.reset(new NacvFactorization(std::move(matrix), dimension));
        }
#ifdef __MPI
        else
        {
            const auto &px = inputs.px[0];
            const int process_dimension = std::max(px.nrproc(), px.ncproc());
            const int block = std::max(1, std::min(32, dimension / process_dimension));
            LR_Util::setup_2d_division(
                matrix_layout,
                block,
                dimension,
                dimension,
                px.ictxt());
            auto matrix = zvec_hessian_2d<double>(
                left,
                spin_count,
                inputs.nloc,
                matrix_layout);
            const bool cholesky = solver == "scalapack_chol";
            factorization.reset(new NacvFactorization(std::move(matrix), matrix_layout, cholesky));
        }
#else
        else
        {
            throw std::logic_error("ScaLAPACK NAC response requires an MPI build");
        }
#endif
    }

    template <typename Operator>
    void solve_with(
        Operator &left,
        double *rhs,
        double *solution,
        int count,
        int spin_count)
    {
        const int ld = inputs.nloc;
        if (solver == "cg")
        {
            const auto apply = [&left, ld, count](const double *in, double *out)
            {
                left.hPsi(
                    in,
                    out,
                    ld,
                    count);
            };
            solve_nac_cg(
                solution,
                rhs,
                ld,
                count,
                apply,
                inputs.test_force,
                controls.tolerance,
                controls.max_iterations);
            return;
        }
        if (!factorization)
        {
            prepare_matrix(left, spin_count);
        }
        const int dimension = zvec_global_dim(left, spin_count);
        const std::size_t full_size = static_cast<std::size_t>(dimension) * count;
        std::vector<double> full(full_size, 0.0);
        for (int state = 0; state < count; ++state)
        {
            double *target = full.data() + static_cast<std::size_t>(state) * dimension;
            const double *source = rhs + state * ld;
#ifdef __MPI
            zvec_local_to_full(
                left,
                spin_count,
                source,
                target);
#else
            std::copy_n(source, dimension, target);
#endif
        }
        if (solver == "lapack")
        {
            factorization->solve(full, count);
        }
#ifdef __MPI
        else
        {
            Parallel_Blkcyclic rhs_layout;
            LR_Util::setup_2d_division(
                rhs_layout,
                matrix_layout.nb(),
                dimension,
                count,
                matrix_layout.ictxt());
            const std::size_t storage = std::max<std::size_t>(1, rhs_layout.size());
            std::vector<double> local(storage, 0.0);
            for (int col = 0; col < rhs_layout.ncloc(); ++col)
            {
                for (int row = 0; row < rhs_layout.nrloc(); ++row)
                {
                    local[col * rhs_layout.nrloc() + row] =
                        full[rhs_layout.icglb(col) * dimension + rhs_layout.irglb(row)];
                }
            }
            factorization->solve(local, rhs_layout);
            std::fill(full.begin(), full.end(), 0.0);
            LR_Util::gather_2d_to_full(
                rhs_layout,
                local.data(),
                full.data(),
                false,
                dimension,
                count);
        }
#endif
        for (int state = 0; state < count; ++state)
        {
            const double *source = full.data() + static_cast<std::size_t>(state) * dimension;
            double *target = solution + state * ld;
#ifdef __MPI
            zvec_full_to_local(
                left,
                spin_count,
                source,
                target);
#else
            std::copy_n(source, dimension, target);
#endif
        }
        if (inputs.test_force)
        {
            LR_Util::print_value(solution, count, ld);
        }
    }
};

NacvResponse::NacvResponse(
    const GradientInputs<double> &inputs,
    int channel,
    bool open_shell,
    const std::string &solver,
    const NacvResponseControls &controls)
{
    impl_.reset(new Impl(
        inputs,
        channel,
        open_shell,
        solver,
        controls));
}
NacvResponse::~NacvResponse() = default;

void NacvResponse::solve(double *rhs, double *solution, int count)
{
    const NacvTimerScope timer("NacvResponse", "solve");
    if (!impl_->closed_left && !impl_->open_left)
    {
        impl_->prepare_left();
    }
    if (impl_->open_shell)
    {
        impl_->solve_with(
            *impl_->open_left,
            rhs,
            solution,
            count,
            2);
    }
    else
    {
        impl_->solve_with(
            *impl_->closed_left,
            rhs,
            solution,
            count,
            1);
    }
}
void NacvResponse::solve_polarized(const double *amplitudes, double *solution, int count)
{
    const NacvTimerScope timer("NacvResponse", "solve_polarized");
    if (!impl_->closed_right && !impl_->open_right)
    {
        impl_->prepare_right();
    }
    const int ld = impl_->inputs.nloc;
    ct::Tensor rhs = LR_Util::newTensor<double>({count, ld});
    rhs.zero();
    if (impl_->open_shell)
    {
        impl_->open_right->hPsi(
            amplitudes,
            rhs.data<double>(),
            ld,
            count);
    }
    else
    {
        impl_->closed_right->hPsi(
            amplitudes,
            rhs.data<double>(),
            ld,
            count);
    }
    solve(rhs.data<double>(), solution, count);
}
} // namespace LR
