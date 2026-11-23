/** @file nacv_ao.cpp
 * @brief Contract auxiliary TDA interstate densities with moving AO derivatives.
 */
#include "nacv_ao.h"
#include "cal_edm.h"
#include "gradient_inputs.h"
#include "pulay_hc.h"
#include <algorithm>
#include <cmath>
#include <cstdint>

namespace LR
{
namespace
{
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

} // namespace

ct::Tensor nac_density_ao(
    const LR::GradientInputs<double> &inputs,
    bool ground,
    bool open_shell,
    const double *bra,
    const double *ket,
    int spin_channel)
{
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

ModuleBase::matrix antisymmetric_ao_connection(
    const GradientInputs<double> &inputs,
    const ct::Tensor &density,
    const std::vector<hamilt::HContainer<double>> &overlap_derivatives)
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
        overlap_derivatives,
        contraction_factor);
}
} // namespace LR
