#include "source_lcao/module_extrap/wf_orthonormalize_lcao.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

#ifdef __MPI
#include <mpi.h>
#endif

namespace ModuleExtrap
{
namespace
{

inline std::size_t idx(const int row, const int col, const int ncol) noexcept
{
    return static_cast<std::size_t>(row) * static_cast<std::size_t>(ncol) + static_cast<std::size_t>(col);
}

int active_bands_for_state(const ModuleBase::matrix& occupations,
                           const int istate,
                           const int nbands,
                           const double occupation_threshold)
{
    if (occupations.c == nullptr || occupations.nr <= istate || occupations.nc <= 0)
    {
        return nbands;
    }

    // In ABACUS, Psi may contain more local eigenvectors than the occupation
    // matrix explicitly stores.  This is already supported by cal_dm_psi(),
    // which simply ignores local bands without a corresponding global weight.
    const int nband_with_occupation = std::min(nbands, occupations.nc);
    int nactive = 0;
    for (int ib = 0; ib < nband_with_occupation; ++ib)
    {
        if (std::abs(occupations(istate, ib)) > occupation_threshold)
        {
            nactive = ib + 1;
        }
    }
    return nactive;
}

void build_overlap_in_band_space(const double* overlap,
                                 const double* coeff,
                                 const int nactive_bands,
                                 const int nbasis,
                                 std::vector<double>& metric)
{
    metric.assign(static_cast<std::size_t>(nactive_bands) * static_cast<std::size_t>(nactive_bands), 0.0);
    std::vector<double> coeff_s(static_cast<std::size_t>(nactive_bands) * static_cast<std::size_t>(nbasis), 0.0);

    // Coefficients are stored as a row-major band-by-basis matrix C.  Build the
    // band-space metric explicitly as O = C S C^T to avoid ambiguity about the
    // row-major layout when using Fortran BLAS wrappers.
    for (int ib = 0; ib < nactive_bands; ++ib)
    {
        for (int mu = 0; mu < nbasis; ++mu)
        {
            double value = 0.0;
            for (int nu = 0; nu < nbasis; ++nu)
            {
                value += coeff[idx(ib, nu, nbasis)] * overlap[idx(nu, mu, nbasis)];
            }
            coeff_s[idx(ib, mu, nbasis)] = value;
        }
    }

    for (int ib = 0; ib < nactive_bands; ++ib)
    {
        for (int jb = 0; jb <= ib; ++jb)
        {
            double value = 0.0;
            for (int mu = 0; mu < nbasis; ++mu)
            {
                value += coeff_s[idx(ib, mu, nbasis)] * coeff[idx(jb, mu, nbasis)];
            }
            metric[idx(ib, jb, nactive_bands)] = value;
            metric[idx(jb, ib, nactive_bands)] = value;
        }
    }
}

void fill_metric_diagnostics(const std::vector<double>& metric,
                             const int n,
                             WfOrthonormalizeResult& result)
{
    if (n <= 0 || metric.empty())
    {
        return;
    }

    result.min_metric_diag = metric[idx(0, 0, n)];
    result.max_metric_diag = metric[idx(0, 0, n)];
    result.max_metric_abs = 0.0;
    result.max_metric_asymmetry = 0.0;

    for (int i = 0; i < n; ++i)
    {
        const double diag = metric[idx(i, i, n)];
        result.min_metric_diag = std::min(result.min_metric_diag, diag);
        result.max_metric_diag = std::max(result.max_metric_diag, diag);
        for (int j = 0; j < n; ++j)
        {
            result.max_metric_abs = std::max(result.max_metric_abs, std::abs(metric[idx(i, j, n)]));
            result.max_metric_asymmetry = std::max(result.max_metric_asymmetry,
                                                   std::abs(metric[idx(i, j, n)] - metric[idx(j, i, n)]));
        }
    }
}

bool cholesky_lower_in_place(std::vector<double>& metric,
                             const int n,
                             const double pivot_threshold,
                             WfOrthonormalizeResult& result)

{
    for (int j = 0; j < n; ++j)
    {
        double diag = metric[idx(j, j, n)];
        for (int k = 0; k < j; ++k)
        {
            const double ljk = metric[idx(j, k, n)];
            diag -= ljk * ljk;
        }

        if (!(diag > pivot_threshold) || !std::isfinite(diag))
        {
            result.failed_pivot_index = j;
            result.failed_pivot = diag;
            return false;
        }

        const double ljj = std::sqrt(diag);
        metric[idx(j, j, n)] = ljj;

        for (int i = j + 1; i < n; ++i)
        {
            double value = metric[idx(i, j, n)];
            for (int k = 0; k < j; ++k)
            {
                value -= metric[idx(i, k, n)] * metric[idx(j, k, n)];
            }
            metric[idx(i, j, n)] = value / ljj;
        }

        for (int k = j + 1; k < n; ++k)
        {
            metric[idx(j, k, n)] = 0.0;
        }
    }
    return true;
}

void apply_inverse_cholesky_left(const std::vector<double>& lower,
                                 double* coeff,
                                 const int nactive_bands,
                                 const int nbasis)
{
    std::vector<double> column(static_cast<std::size_t>(nactive_bands), 0.0);

    // Solve L * C_new(:,mu) = C_old(:,mu) for every basis function mu.  Only
    // occupied/fractionally occupied bands are transformed; unoccupied bands do
    // not contribute to the density matrix and are left unchanged.
    for (int mu = 0; mu < nbasis; ++mu)
    {
        for (int ib = 0; ib < nactive_bands; ++ib)
        {
            double value = coeff[idx(ib, mu, nbasis)];
            for (int jb = 0; jb < ib; ++jb)
            {
                value -= lower[idx(ib, jb, nactive_bands)] * column[jb];
            }
            column[ib] = value / lower[idx(ib, ib, nactive_bands)];
        }

        for (int ib = 0; ib < nactive_bands; ++ib)
        {
            coeff[idx(ib, mu, nbasis)] = column[ib];
        }
    }
}

double max_orthonormality_deviation(const double* overlap,
                                    const double* coeff,
                                    const int nactive_bands,
                                    const int nbasis)
{
    std::vector<double> metric;
    build_overlap_in_band_space(overlap, coeff, nactive_bands, nbasis, metric);

    double max_dev = 0.0;
    for (int i = 0; i < nactive_bands; ++i)
    {
        for (int j = 0; j < nactive_bands; ++j)
        {
            const double expected = (i == j) ? 1.0 : 0.0;
            max_dev = std::max(max_dev, std::abs(metric[idx(i, j, nactive_bands)] - expected));
        }
    }
    return max_dev;
}


WfcExtrapStatus allreduce_dense_vector(std::vector<double>& values, const Parallel_Orbitals& pv)
{
#ifdef __MPI
    if (!pv.is_serial)
    {
        const MPI_Comm comm = pv.comm();
        if (comm == MPI_COMM_NULL)
        {
            return WfcExtrapStatus::InvalidInput;
        }
        if (values.size() > static_cast<std::size_t>(std::numeric_limits<int>::max()))
        {
            return WfcExtrapStatus::InvalidInput;
        }
        MPI_Allreduce(MPI_IN_PLACE,
                      values.data(),
                      static_cast<int>(values.size()),
                      MPI_DOUBLE,
                      MPI_SUM,
                      comm);
    }
#endif
    return WfcExtrapStatus::Success;
}

WfcExtrapStatus assemble_global_overlap(const double* overlap,
                                        const Parallel_Orbitals& pv,
                                        std::vector<double>& overlap_global)
{
    const int nbasis_global = pv.get_global_row_size();
    if (overlap == nullptr || nbasis_global <= 0 || pv.get_global_col_size() != nbasis_global)
    {
        return WfcExtrapStatus::DimensionMismatch;
    }

    overlap_global.assign(static_cast<std::size_t>(nbasis_global) * static_cast<std::size_t>(nbasis_global), 0.0);

    const int nrow_local = pv.get_row_size();
    const int ncol_local = pv.get_col_size();
    for (int icol = 0; icol < ncol_local; ++icol)
    {
        const int icol_global = pv.local2global_col(icol);
        for (int irow = 0; irow < nrow_local; ++irow)
        {
            const int irow_global = pv.local2global_row(irow);
            // ABACUS/ScaLAPACK local matrices are stored column-major:
            // local(irow, icol) -> overlap[icol * nrow_local + irow].
            overlap_global[idx(irow_global, icol_global, nbasis_global)]
                = overlap[static_cast<std::size_t>(icol) * static_cast<std::size_t>(nrow_local)
                          + static_cast<std::size_t>(irow)];
        }
    }

    return allreduce_dense_vector(overlap_global, pv);
}

int local_band_to_global(const int ib_local, const Parallel_Orbitals& pv)
{
    return (ib_local / pv.nb * pv.dim1 + pv.coord[1]) * pv.nb + ib_local % pv.nb;
}

int local_wfc_col_to_global(const int ib_local, const Parallel_Orbitals& pv, const bool coeff_columns_are_basis)
{
    // Gamma-only direct diagonalization solvers store Psi columns with the same
    // distribution as the AO/basis columns, i.e. local Psi column -> global basis
    // column.  Other NAO paths may store only the requested band columns; in that
    // case the band-column distribution follows the ScaLAPACK formula used by
    // Parallel_Orbitals::set_nloc_wfc_Eij().
    return coeff_columns_are_basis ? pv.local2global_col(ib_local) : local_band_to_global(ib_local, pv);
}

WfcExtrapStatus assemble_global_coeff_active(const double* coeff_local,
                                             const Parallel_Orbitals& pv,
                                             const int ncoeff_local,
                                             const int nbasis_local,
                                             const int nactive_bands,
                                             const bool coeff_columns_are_basis,
                                             std::vector<double>& coeff_global)
{
    const int nbasis_global = pv.get_global_row_size();
    if (coeff_local == nullptr || nbasis_global <= 0 || ncoeff_local <= 0 || nbasis_local <= 0
        || nactive_bands < 0 || pv.get_row_size() != nbasis_local)
    {
        return WfcExtrapStatus::DimensionMismatch;
    }

    coeff_global.assign(static_cast<std::size_t>(nactive_bands) * static_cast<std::size_t>(nbasis_global), 0.0);

    for (int ib_local = 0; ib_local < ncoeff_local; ++ib_local)
    {
        const int ib_global = local_wfc_col_to_global(ib_local, pv, coeff_columns_are_basis);
        if (ib_global >= nactive_bands)
        {
            continue;
        }
        for (int mu_local = 0; mu_local < nbasis_local; ++mu_local)
        {
            const int mu_global = pv.local2global_row(mu_local);
            coeff_global[idx(ib_global, mu_global, nbasis_global)]
                = coeff_local[static_cast<std::size_t>(ib_local) * static_cast<std::size_t>(nbasis_local)
                              + static_cast<std::size_t>(mu_local)];
        }
    }

    return allreduce_dense_vector(coeff_global, pv);
}

void scatter_global_coeff_active(const std::vector<double>& coeff_global,
                                 const Parallel_Orbitals& pv,
                                 const int ncoeff_local,
                                 const int nbasis_local,
                                 const int nactive_bands,
                                 const bool coeff_columns_are_basis,
                                 double* coeff_local)
{
    const int nbasis_global = pv.get_global_row_size();
    for (int ib_local = 0; ib_local < ncoeff_local; ++ib_local)
    {
        const int ib_global = local_wfc_col_to_global(ib_local, pv, coeff_columns_are_basis);
        if (ib_global >= nactive_bands)
        {
            continue;
        }
        for (int mu_local = 0; mu_local < nbasis_local; ++mu_local)
        {
            const int mu_global = pv.local2global_row(mu_local);
            coeff_local[static_cast<std::size_t>(ib_local) * static_cast<std::size_t>(nbasis_local)
                        + static_cast<std::size_t>(mu_local)]
                = coeff_global[idx(ib_global, mu_global, nbasis_global)];
        }
    }
}

} // namespace

WfOrthonormalizeResult reorthonormalize_gamma_lcao(const double* overlap,
                                                   const Parallel_Orbitals& pv,
                                                   psi::Psi<double>& wfc,
                                                   const ModuleBase::matrix& occupations,
                                                   const double occupation_threshold,
                                                   const double pivot_threshold,
                                                   const double check_tolerance)
{
    WfOrthonormalizeResult result;

    if (overlap == nullptr || !(occupation_threshold >= 0.0) || !(pivot_threshold >= 0.0) || !(check_tolerance > 0.0))
    {
        result.status = WfcExtrapStatus::InvalidInput;
        return result;
    }

    const int nstate = wfc.get_nk();
    const int ncoeff_local = wfc.get_nbands();
    const int nbasis_local = wfc.get_nbasis();
    const int nbasis_global = pv.get_global_row_size();
    const int ncoeff_global_basis = pv.get_global_col_size();
    const int ncoeff_global_bands = pv.get_wfc_global_nbands();
    const bool coeff_columns_are_basis = (ncoeff_local == pv.get_col_size());
    const bool coeff_columns_are_bands = (ncoeff_local == pv.ncol_bands);
    const int ncoeff_global = coeff_columns_are_basis ? ncoeff_global_basis : ncoeff_global_bands;
    result.nstate = nstate;
    result.nbands = ncoeff_global;
    result.nbasis = nbasis_global;

    if (nstate <= 0 || ncoeff_local <= 0 || nbasis_local <= 0 || nbasis_global <= 0 || ncoeff_global <= 0
        || ncoeff_global > nbasis_global || pv.get_row_size() != nbasis_local
        || pv.get_global_col_size() != nbasis_global || (!coeff_columns_are_basis && !coeff_columns_are_bands))
    {
        result.status = WfcExtrapStatus::DimensionMismatch;
        return result;
    }

    if (!wfc.get_k_first())
    {
        result.status = WfcExtrapStatus::Unsupported;
        return result;
    }

    if (occupations.c == nullptr || occupations.nr < nstate || occupations.nc <= 0)
    {
        result.status = WfcExtrapStatus::DimensionMismatch;
        return result;
    }

    std::vector<double> overlap_global;
    result.status = assemble_global_overlap(overlap, pv, overlap_global);
    if (!result.ok())
    {
        return result;
    }

    // Operate on a full owned copy first.  The caller's Psi is overwritten only after
    // all state/spin channels pass the positive-definiteness and residual checks.
    std::vector<double> trial(wfc.size(), 0.0);
    const double* coeff_begin = wfc.get_pointer() - wfc.get_psi_bias();
    std::copy(coeff_begin, coeff_begin + wfc.size(), trial.begin());

    double max_deviation_all = 0.0;
    int max_active_bands = 0;
    for (int istate = 0; istate < nstate; ++istate)
    {
        double* coeff_local = trial.data()
                            + static_cast<std::size_t>(istate) * static_cast<std::size_t>(ncoeff_local)
                                  * static_cast<std::size_t>(nbasis_local);

        const int nactive_bands = active_bands_for_state(occupations, istate, ncoeff_global, occupation_threshold);
        max_active_bands = std::max(max_active_bands, nactive_bands);
        result.nactive_bands = nactive_bands;

        // Empty spin channels may occur in edge cases.  They do not contribute to
        // the density matrix, so there is no occupied subspace to reorthonormalize.
        if (nactive_bands == 0)
        {
            continue;
        }

        std::vector<double> coeff_global;
        result.status = assemble_global_coeff_active(coeff_local,
                                                     pv,
                                                     ncoeff_local,
                                                     nbasis_local,
                                                     nactive_bands,
                                                     coeff_columns_are_basis,
                                                     coeff_global);
        if (!result.ok())
        {
            result.failed_state = istate;
            return result;
        }

        std::vector<double> metric;
        build_overlap_in_band_space(overlap_global.data(), coeff_global.data(), nactive_bands, nbasis_global, metric);
        fill_metric_diagnostics(metric, nactive_bands, result);

        if (!cholesky_lower_in_place(metric, nactive_bands, pivot_threshold, result))
        {
            result.status = WfcExtrapStatus::OrthogonalizationFailed;
            result.failed_state = istate;
            return result;
        }

        apply_inverse_cholesky_left(metric, coeff_global.data(), nactive_bands, nbasis_global);

        const double max_dev = max_orthonormality_deviation(overlap_global.data(),
                                                            coeff_global.data(),
                                                            nactive_bands,
                                                            nbasis_global);
        if (!std::isfinite(max_dev) || max_dev > check_tolerance)
        {
            result.status = WfcExtrapStatus::OrthogonalizationFailed;
            result.failed_state = istate;
            result.max_deviation = max_dev;
            return result;
        }
        max_deviation_all = std::max(max_deviation_all, max_dev);

        scatter_global_coeff_active(coeff_global,
                                    pv,
                                    ncoeff_local,
                                    nbasis_local,
                                    nactive_bands,
                                    coeff_columns_are_basis,
                                    coeff_local);
    }

    wfc.set_all_psi(trial.data(), trial.size());
    result.status = WfcExtrapStatus::Success;
    result.nactive_bands = max_active_bands;
    result.max_deviation = max_deviation_all;
    return result;
}

} // namespace ModuleExtrap
