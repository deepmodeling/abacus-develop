#include "elecstate_tools.h"

#include "occupy.h"
#include "source_base/parallel_comm.h"
#include "source_base/parallel_reduce.h"

#include <algorithm>
#include <numeric>
#include <vector>

namespace
{
// Return the global index of the first band stored by this band group.
int get_band_offset(const int local_nbands, const int global_nbands)
{
    if (local_nbands < 0 || global_nbands < 0)
    {
        ModuleBase::WARNING_QUIT("get_band_offset", "band counts cannot be negative");
    }
#ifndef __MPI
    // A serial process must own the complete band set.
    if (local_nbands != global_nbands)
    {
        ModuleBase::WARNING_QUIT("get_band_offset", "distributed bands require an MPI build");
    }
    return 0;
#else
    int band_group = 0;
    int band_groups = 0;
    MPI_Comm_rank(BP_WORLD, &band_group);
    MPI_Comm_size(BP_WORLD, &band_groups);

    // The collected counts identify whether bands are replicated or partitioned
    // across band groups without relying on a particular eigensolver name.
    std::vector<int> band_counts(band_groups);
    MPI_Allgather(&local_nbands, 1, MPI_INT, band_counts.data(), 1, MPI_INT, BP_WORLD);

    // Non-BPCG SDFT replicates the complete deterministic band set in every
    // band group. Only complementary BPCG shards need a prefix offset.
    const bool bands_are_replicated
        = std::all_of(band_counts.begin(), band_counts.end(), [global_nbands](const int count) { return count == global_nbands; });
    if (bands_are_replicated)
    {
        return 0;
    }

    const int gathered_nbands = std::accumulate(band_counts.begin(), band_counts.end(), 0);
    if (gathered_nbands != global_nbands)
    {
        ModuleBase::WARNING_QUIT("get_band_offset", "local band counts do not match global nbands");
    }
    // BPCG assigns contiguous global band ranges in band-group order.
    return std::accumulate(band_counts.begin(), band_counts.begin() + band_group, 0);
#endif
}
} // namespace

namespace elecstate
{
void calEBand(const ModuleBase::matrix& ekb, const ModuleBase::matrix& wg, fenergy& f_en)
{
    ModuleBase::TITLE("ElecState", "calEBand");
    // Each rank accumulates the contribution from its local k-points and bands.
    double eband = 0.0;
#ifdef _OPENMP
#pragma omp parallel for collapse(2) reduction(+ : eband)
#endif
    for (int ik = 0; ik < ekb.nr; ++ik)
    {
        for (int ibnd = 0; ibnd < ekb.nc; ibnd++)
        {
            eband += ekb(ik, ibnd) * wg(ik, ibnd);
        }
    }
    f_en.eband = eband;

#ifdef __MPI
    // Combine contributions distributed by both KPAR and BNDPAR.
    const int npool = GlobalV::KPAR * PARAM.inp.bndpar;
    Parallel_Reduce::reduce_double_allpool(npool, GlobalV::NPROC_IN_POOL, f_en.eband);
#endif
    return;
}

void calculate_weights(const ModuleBase::matrix& ekb,
                       ModuleBase::matrix& wg,
                       const K_Vectors* klist,
                       Efermi& eferm,
                       fenergy& f_en,
                       std::vector<double>& nelec_spin,
                       const int global_nbands,
                       const bool skip_weights)
{
    ModuleBase::TITLE("ElecState", "calculate_weights");
    // Preserve occupations supplied through ocp_set instead of recomputing them.
    if (skip_weights == true)
        return;

    const int nbands = ekb.nc;
    const int nks = ekb.nr;
    if (!(Occupy::use_gaussian_broadening || Occupy::fixed_occupations))
    {
        // Taoni fix smearing_method=fixed for BPCG on 2026-08-21
        // Integer occupations use global band indices even when ekb is a local
        // contiguous BPCG shard.
        const int band_offset = get_band_offset(nbands, global_nbands);
        if (PARAM.globalv.two_fermi)
        {
            Occupy::iweights(nks, klist->wk, nbands, band_offset, nelec_spin[0], ekb, eferm.ef_up, wg, 0, klist->isk);
            Occupy::iweights(nks, klist->wk, nbands, band_offset, nelec_spin[1], ekb, eferm.ef_dw, wg, 1, klist->isk);
            // ef = ( ef_up + ef_dw ) / 2.0_dp need??? mohan add 2012-04-16
            // Keep independent Fermi levels for the two spin channels.
        }
        else
        {
            // A spin selector of -1 requests the combined-spin occupation path.
            Occupy::iweights(nks, klist->wk, nbands, band_offset, PARAM.inp.nelec, ekb, eferm.ef, wg, -1, klist->isk);
        }
    }
    else if (Occupy::use_gaussian_broadening)
    {
        if (PARAM.globalv.two_fermi)
        {
            double demet_up = 0.0;
            double demet_dw = 0.0;
            Occupy::gweights(nks,
                             klist->wk,
                             nbands,
                             nelec_spin[0],
                             Occupy::gaussian_parameter,
                             Occupy::gaussian_type,
                             ekb,
                             eferm.ef_up,
                             demet_up,
                             wg,
                             0,
                             klist->isk);
            Occupy::gweights(nks,
                             klist->wk,
                             nbands,
                             nelec_spin[1],
                             Occupy::gaussian_parameter,
                             Occupy::gaussian_type,
                             ekb,
                             eferm.ef_dw,
                             demet_dw,
                             wg,
                             1,
                             klist->isk);
            f_en.demet = demet_up + demet_dw;
        }
        else
        {
            // A spin selector of -1 requests the combined-spin occupation path.
            Occupy::gweights(nks,
                             klist->wk,
                             nbands,
                             PARAM.inp.nelec,
                             Occupy::gaussian_parameter,
                             Occupy::gaussian_type,
                             ekb,
                             eferm.ef,
                             f_en.demet,
                             wg,
                             -1,
                             klist->isk);
        }
#ifdef __MPI
        // demet is accumulated independently on every k-point and band partition.
        const int npool = GlobalV::KPAR * PARAM.inp.bndpar;
        Parallel_Reduce::reduce_double_allpool(npool, GlobalV::NPROC_IN_POOL, f_en.demet);
#endif
    }
    else if (Occupy::fixed_occupations)
    {
        ModuleBase::WARNING_QUIT("calculate_weights", "other occupations, not implemented");
    }

    return;
}

void fixed_weights(const std::vector<double>& ocp_kb,
                   const int& nbands,
                   const double& nelec,
                   const K_Vectors* klist,
                   ModuleBase::matrix& wg,
                   bool& skip_weights)
{
    assert(nbands > 0);
    assert(nelec > 0.0);

    const double ne_thr = 1.0e-5;

    // ocp_kb is a global k-major, band-minor array supplied by the user.
    const std::size_t expected_size = static_cast<std::size_t>(klist->get_nkstot()) * static_cast<std::size_t>(nbands);
    if (expected_size != ocp_kb.size())
    {
        ModuleBase::WARNING_QUIT("ElecState::fixed_weights", "size of occupation array is wrong , please check ocp_set");
    }

    double num_elec = 0.0;
    for (std::size_t i = 0; i < ocp_kb.size(); ++i)
    {
        num_elec += ocp_kb[i];
    }

    if (std::abs(num_elec - nelec) > ne_thr)
    {
        ModuleBase::WARNING_QUIT("ElecState::fixed_weights", "total number of occupations is wrong , please check ocp_set");
    }

    // Translate this rank's local k-point and band indices into ocp_kb indices.
    const int band_offset = get_band_offset(wg.nc, nbands);
    if (klist->ik2iktot.size() < static_cast<std::size_t>(wg.nr) || band_offset < 0 || band_offset + wg.nc > nbands)
    {
        ModuleBase::WARNING_QUIT("ElecState::fixed_weights", "invalid distributed occupation layout");
    }

    for (int ik = 0; ik < wg.nr; ++ik)
    {
        const int global_k = klist->ik2iktot[ik];
        if (global_k < 0 || global_k >= klist->get_nkstot())
        {
            ModuleBase::WARNING_QUIT("ElecState::fixed_weights", "invalid global k-point index");
        }
        for (int ib = 0; ib < wg.nc; ++ib)
        {
            // Taoni fix ocp_set under kpar and bndpar on 2026-08-21
            const int global_band = band_offset + ib;
            const std::size_t occupation_index
                = static_cast<std::size_t>(global_k) * static_cast<std::size_t>(nbands) + static_cast<std::size_t>(global_band);
            wg(ik, ib) = ocp_kb[occupation_index];
        }
    }
    // Subsequent eigensolver calls must not overwrite these user occupations.
    skip_weights = true;
}
} // namespace elecstate
