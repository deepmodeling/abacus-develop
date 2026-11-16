#include "nacv_cs.h"

#include "cal_edm.h"
#include "cal_hs_grad.h"
#include "gradient_inputs.h"
#include "hamilt_zeq_l.h"
#include "nacv_response.h"
#include "pulay_hc.h"
#include "source_base/timer.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
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
 * @param[in] integrals Current two-center integrals used for overlap derivatives.
 * @return Owned atom-major nat-by-3 connection in bohr^-1.
 * @details Contracts the antisymmetric density with the bra-center overlap
 * derivative; the symmetric metric-response term is handled separately.
 */
ModuleBase::matrix antisymmetric_ao_connection(
    const LR::GradientInputs<double> &inputs, const ct::Tensor &density, const TwoCenterBundle &integrals)
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
    const std::vector<hamilt::HContainer<double>> overlap_derivative = cal_hs_grad(
        'S',
        inputs.ucell,
        inputs.pmat,
        inputs.gd,
        integrals);
    // cal_hs_grad stores the separation derivative, minus the bra-center derivative.
    // The Pulay contraction contributes 2*factor*D*dHS. Thus -1/4 gives
    // Gamma : (B-B^T)/2 from D-D^T, with Gamma=D^T. No symmetric metric term is added here.
    const double contraction_factor = -0.25;
    return PulayForceStress::cal_pulay_fs(
        dm,
        inputs.ucell,
        overlap_derivative,
        contraction_factor);
}

/**
 * @brief Solve the GS A+B adjoint for a reference-to-excited auxiliary NAC.
 * @param[in] inputs Validated singlet geometry, orbitals, kernels and distributions.
 * @param[in] ket Local TDA excitation amplitudes in the full response window.
 * @param[in] energy Positive ket excitation energy in Ry.
 * @param[in] solver Explicit existing Z-vector solver selection.
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
    const std::string &solver,
    bool open_shell,
    const LR::NacvResponseControls &controls)
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
    if (open_shell)
    {
        const int global_dimension = inputs.nocc[0] * inputs.nvirt[0] + inputs.nocc[1] * inputs.nvirt[1];
        LR::Z_vector_UL<double> hessian(
            inputs.xc_kernel,
            inputs.nspin,
            inputs.nbasis,
            inputs.nocc,
            inputs.nvirt,
            inputs.ucell,
            inputs.orb_cutoff,
            inputs.gd,
            inputs.psi_ks,
            inputs.eig_ks,
#ifdef __EXX
            inputs.exx_lri,
            inputs.hybrid_alpha,
#endif
            inputs.pot_hxc_gs,
            inputs.kv,
            inputs.px,
            inputs.pc,
            inputs.pmat,
            inputs.dft_functional,
            inputs.nk,
            ld,
            global_dimension);
        LR::solve_nac_zeq_with(
            solution.data<double>(),
            rhs.data<double>(),
            ld,
            1,
            hessian,
            2,
            solver,
            false,
            controls);
        return solution;
    }
    const std::string spin = "singlet";
    const std::string input_directory;
    LR::Z_vector_L<double> hessian(
        inputs.xc_kernel,
        inputs.nspin,
        inputs.nbasis,
        inputs.nocc,
        inputs.nvirt,
        inputs.ucell,
        inputs.orb_cutoff,
        inputs.gd,
        inputs.psi_ks,
        inputs.eig_ks,
#ifdef __EXX
        inputs.exx_lri,
        inputs.hybrid_alpha,
#endif
        inputs.pot_hxc_gs,
        inputs.kv,
        inputs.px,
        inputs.pc,
        inputs.pmat,
        spin,
        input_directory,
        inputs.out_dir,
        inputs.dft_functional);
    LR::solve_nac_zeq_with(
        solution.data<double>(),
        rhs.data<double>(),
        ld,
        1,
        hessian,
        1,
        solver,
        false,
        controls);
    return solution;
}

} // namespace

namespace LR
{
NacvResult evaluate_aux_tda_nac(
    const GradientInputs<double> &inputs,
    LR_Force<double> &force_terms,
    const module_dm::DensityMatrix<double, double> &ground_density,
    const TwoCenterBundle &two_center_bundle,
    int spin_channel,
    bool bra_ground,
    const double *bra,
    const double *ket,
    double bra_energy_ry,
    double ket_energy_ry,
    double gap_resolution_ry,
    const std::string &response_solver,
    const NacvResponseControls &controls)
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
            response_solver,
            open_shell,
            controls);
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
        const std::string input_directory;
        const std::string &spin = inputs.spin_types[spin_channel];
        std::weak_ptr<PotHxcLR> potential = inputs.pot[spin_channel];
        const int pair_count = 2;
        solve_nac_zvector(
            inputs,
            combinations.data<double>(),
            response.data<double>(),
            pair_count,
            potential,
            spin,
            input_directory,
            open_shell,
            response_solver,
            controls);
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
    result.ao_connection = antisymmetric_ao_connection(inputs, density, two_center_bundle);
    assemble_nacv_result(result, pair_force, gap_resolution_ry);
    return result;
}

} // namespace LR
