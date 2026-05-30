#include "source_lcao/module_extrap/wf_orthonormalize_lcao.h"

#include "source_base/module_external/blas_connector.h"

#include <algorithm>
#include <cmath>
#include <vector>

namespace ModuleExtrap
{
namespace
{

inline std::size_t idx(const int row, const int col, const int ncol) noexcept
{
    return static_cast<std::size_t>(row) * static_cast<std::size_t>(ncol) + static_cast<std::size_t>(col);
}

void build_overlap_in_band_space(const double* overlap,
                                 const double* coeff,
                                 const int nbands,
                                 const int nbasis,
                                 std::vector<double>& metric)
{
    metric.assign(static_cast<std::size_t>(nbands) * static_cast<std::size_t>(nbands), 0.0);
    std::vector<double> coeff_s(static_cast<std::size_t>(nbands) * static_cast<std::size_t>(nbasis), 0.0);

    const double one = 1.0;
    const double zero = 0.0;

    // Coefficients are stored as a row-major band-by-basis matrix X.  With the
    // current row-major overlap S, the band-space metric is O = X S X^T.
    BlasConnector::gemm('N',
                        'N',
                        nbands,
                        nbasis,
                        nbasis,
                        one,
                        coeff,
                        nbasis,
                        overlap,
                        nbasis,
                        zero,
                        coeff_s.data(),
                        nbasis);

    BlasConnector::gemm('N',
                        'T',
                        nbands,
                        nbands,
                        nbasis,
                        one,
                        coeff_s.data(),
                        nbasis,
                        coeff,
                        nbasis,
                        zero,
                        metric.data(),
                        nbands);

    // Remove tiny asymmetry introduced by BLAS roundoff before the Cholesky step.
    for (int i = 0; i < nbands; ++i)
    {
        for (int j = 0; j < i; ++j)
        {
            const double value = 0.5 * (metric[idx(i, j, nbands)] + metric[idx(j, i, nbands)]);
            metric[idx(i, j, nbands)] = value;
            metric[idx(j, i, nbands)] = value;
        }
    }
}

bool cholesky_lower_in_place(std::vector<double>& metric, const int n, const double pivot_threshold)
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
                                 const int nbands,
                                 const int nbasis)
{
    std::vector<double> column(static_cast<std::size_t>(nbands), 0.0);

    // Solve L * C_new(:,mu) = C_old(:,mu) for every basis function mu.
    for (int mu = 0; mu < nbasis; ++mu)
    {
        for (int ib = 0; ib < nbands; ++ib)
        {
            double value = coeff[idx(ib, mu, nbasis)];
            for (int jb = 0; jb < ib; ++jb)
            {
                value -= lower[idx(ib, jb, nbands)] * column[jb];
            }
            column[ib] = value / lower[idx(ib, ib, nbands)];
        }

        for (int ib = 0; ib < nbands; ++ib)
        {
            coeff[idx(ib, mu, nbasis)] = column[ib];
        }
    }
}

double max_orthonormality_deviation(const double* overlap,
                                    const double* coeff,
                                    const int nbands,
                                    const int nbasis)
{
    std::vector<double> metric;
    build_overlap_in_band_space(overlap, coeff, nbands, nbasis, metric);

    double max_dev = 0.0;
    for (int i = 0; i < nbands; ++i)
    {
        for (int j = 0; j < nbands; ++j)
        {
            const double expected = (i == j) ? 1.0 : 0.0;
            max_dev = std::max(max_dev, std::abs(metric[idx(i, j, nbands)] - expected));
        }
    }
    return max_dev;
}

} // namespace

WfOrthonormalizeResult reorthonormalize_gamma_lcao(const double* overlap,
                                                   psi::Psi<double>& wfc,
                                                   const double pivot_threshold,
                                                   const double check_tolerance)
{
    WfOrthonormalizeResult result;

    if (overlap == nullptr || !(pivot_threshold >= 0.0) || !(check_tolerance > 0.0))
    {
        result.status = WfcExtrapStatus::InvalidInput;
        return result;
    }

    const int nstate = wfc.get_nk();
    const int nbands = wfc.get_nbands();
    const int nbasis = wfc.get_nbasis();

    if (nstate <= 0 || nbands <= 0 || nbasis <= 0 || nbands > nbasis)
    {
        result.status = WfcExtrapStatus::DimensionMismatch;
        return result;
    }

    if (!wfc.get_k_first())
    {
        result.status = WfcExtrapStatus::Unsupported;
        return result;
    }

    // Operate on a full owned copy first.  The caller's Psi is overwritten only after
    // all state/spin channels pass the positive-definiteness and residual checks.
    std::vector<double> trial(wfc.size(), 0.0);
    const double* coeff_begin = wfc.get_pointer() - wfc.get_psi_bias();
    std::copy(coeff_begin, coeff_begin + wfc.size(), trial.begin());

    double max_deviation_all = 0.0;
    for (int istate = 0; istate < nstate; ++istate)
    {
        double* coeff = trial.data()
                      + static_cast<std::size_t>(istate) * static_cast<std::size_t>(nbands)
                            * static_cast<std::size_t>(nbasis);

        std::vector<double> metric;
        build_overlap_in_band_space(overlap, coeff, nbands, nbasis, metric);

        if (!cholesky_lower_in_place(metric, nbands, pivot_threshold))
        {
            result.status = WfcExtrapStatus::OrthogonalizationFailed;
            result.failed_state = istate;
            return result;
        }

        apply_inverse_cholesky_left(metric, coeff, nbands, nbasis);

        const double max_dev = max_orthonormality_deviation(overlap, coeff, nbands, nbasis);
        if (!std::isfinite(max_dev) || max_dev > check_tolerance)
        {
            result.status = WfcExtrapStatus::OrthogonalizationFailed;
            result.failed_state = istate;
            result.max_deviation = max_dev;
            return result;
        }
        max_deviation_all = std::max(max_deviation_all, max_dev);
    }

    wfc.set_all_psi(trial.data(), trial.size());
    result.status = WfcExtrapStatus::Success;
    result.max_deviation = max_deviation_all;
    return result;
}

} // namespace ModuleExtrap
