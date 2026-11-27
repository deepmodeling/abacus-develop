#include "source_base/memory_budget.h"

#include <iomanip>
#include <limits>
#include <sstream>
#include <string>

namespace ModuleBase
{

namespace
{

std::size_t saturating_mul(const std::size_t a, const std::size_t b)
{
    if (a == 0 || b == 0)
    {
        return 0;
    }
    if (a > std::numeric_limits<std::size_t>::max() / b)
    {
        return std::numeric_limits<std::size_t>::max();
    }
    return a * b;
}

std::size_t saturating_add(const std::size_t a, const std::size_t b)
{
    if (a > std::numeric_limits<std::size_t>::max() - b)
    {
        return std::numeric_limits<std::size_t>::max();
    }
    return a + b;
}

} // namespace

LcaoDenseEstimate estimate_lcao_dense(const std::size_t nrloc,
                                      const std::size_t ncloc,
                                      const std::size_t nbands,
                                      const std::size_t nks,
                                      const int nspin,
                                      const bool gamma_only,
                                      const bool gamma_direct,
                                      const bool use_elpa,
                                      const std::size_t element_size)
{
    LcaoDenseEstimate estimate;
    const std::size_t local_elements = saturating_mul(nrloc, ncloc);
    const std::size_t matrix_bytes = saturating_mul(local_elements, element_size);

    estimate.hk_bytes = matrix_bytes;
    estimate.sk_bytes = matrix_bytes;
    estimate.dmk_bytes = saturating_mul(matrix_bytes, nks);

    std::size_t psi_elements = 0;
    if (gamma_only && gamma_direct)
    {
        // Gamma-only direct solvers keep ncol = ncloc, so Psi is a full
        // local N x N block rather than nbands x nrloc.
        psi_elements = saturating_mul(local_elements, static_cast<std::size_t>(nspin > 0 ? nspin : 1));
    }
    else
    {
        psi_elements = saturating_mul(saturating_mul(nks, nbands), nrloc);
    }
    estimate.psi_bytes = saturating_mul(psi_elements, element_size);

    if (use_elpa)
    {
        // ELPA allocates one additional local work matrix of the same shape.
        estimate.elpa_work_bytes = matrix_bytes;
    }

    const std::size_t hs_bytes = saturating_add(estimate.hk_bytes, estimate.sk_bytes);
    const std::size_t hsdmk_bytes = saturating_add(hs_bytes, estimate.dmk_bytes);
    const std::size_t psi_elpa_bytes = saturating_add(estimate.psi_bytes, estimate.elpa_work_bytes);
    estimate.per_rank_total_bytes = saturating_add(hsdmk_bytes, psi_elpa_bytes);
    return estimate;
}

bool is_gamma_direct_ks_solver(const std::string& ks_solver)
{
    return ks_solver == "genelpa" || ks_solver == "elpa" || ks_solver == "lapack"
           || ks_solver == "pexsi" || ks_solver == "cusolver" || ks_solver == "cusolvermp";
}

bool is_elpa_ks_solver(const std::string& ks_solver)
{
    return ks_solver == "genelpa" || ks_solver == "elpa";
}

std::string format_memory_bytes(const std::size_t bytes)
{
    static const char* units[] = {"B", "KB", "MB", "GB", "TB", "PB"};
    double value = static_cast<double>(bytes);
    int unit = 0;
    while (value >= 1024.0 && unit < 5)
    {
        value /= 1024.0;
        ++unit;
    }
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(2) << value << " " << units[unit];
    return oss.str();
}

} // namespace ModuleBase
