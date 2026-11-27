#include "source_lcao/lcao_memory_check.h"

#include <complex>
#include <iostream>
#include <sstream>

#include "source_base/memory_budget.h"

namespace LCAO_domain
{

template <typename TK>
void check_lcao_dense_memory(const Parallel_Orbitals& pv,
                             const K_Vectors& kv,
                             const Input_para& inp,
                             bool gamma_only,
                             int my_rank,
                             std::ostream& log)
{
    // Only rank 0 prints; the estimate is informational only and never aborts.
    if (my_rank != 0)
    {
        return;
    }

    const std::size_t nbasis_global = static_cast<std::size_t>(pv.get_wfc_global_nbasis());
    if (nbasis_global == 0)
    {
        return;
    }

    // Estimate the memory of one full (non-distributed) set of the dense LCAO
    // matrices by feeding the global NLOCAL x NLOCAL dimensions into the
    // estimator.  This is independent of the 2D block-cyclic distribution and
    // of any per-rank budget accounting.
    const std::size_t nbands_global = static_cast<std::size_t>(pv.get_wfc_global_nbands());
    const std::size_t nks_s = static_cast<std::size_t>(kv.get_nks());
    const bool gamma_direct = gamma_only && ModuleBase::is_gamma_direct_ks_solver(inp.ks_solver);
    const bool use_elpa = ModuleBase::is_elpa_ks_solver(inp.ks_solver);
    const std::size_t element_size = sizeof(TK);

    const ModuleBase::LcaoDenseEstimate est = ModuleBase::estimate_lcao_dense(
        nbasis_global, nbasis_global, nbands_global, nks_s, inp.nspin,
        gamma_only, gamma_direct, use_elpa, element_size);

    std::ostringstream oss;
    oss << "\n ABACUS: LCAO dense-matrix memory estimate (lower bound)\n";
    oss << "   NLOCAL=" << pv.get_wfc_global_nbasis()
        << "  nbands=" << pv.get_wfc_global_nbands()
        << "  nspin=" << inp.nspin
        << "  nks=" << kv.get_nks()
        << "  gamma_only=" << (gamma_only ? "true" : "false")
        << "  solver=" << inp.ks_solver << "\n";
    oss << "   required (full matrices, one copy): H "
        << ModuleBase::format_memory_bytes(est.hk_bytes)
        << " + S " << ModuleBase::format_memory_bytes(est.sk_bytes)
        << " + DMK " << ModuleBase::format_memory_bytes(est.dmk_bytes)
        << " + Psi " << ModuleBase::format_memory_bytes(est.psi_bytes);
    if (use_elpa)
    {
        oss << " + ELPA-work " << ModuleBase::format_memory_bytes(est.elpa_work_bytes);
    }
    oss << "\n";
    oss << "     = " << ModuleBase::format_memory_bytes(est.per_rank_total_bytes)
        << " total\n";
    // Write to the running log only, not stdout, so the estimate block does
    // not clutter the screen.
    log << oss.str() << std::endl;
}

template void check_lcao_dense_memory<double>(const Parallel_Orbitals&,
                                              const K_Vectors&,
                                              const Input_para&,
                                              bool,
                                              int,
                                              std::ostream&);
template void check_lcao_dense_memory<std::complex<double>>(const Parallel_Orbitals&,
                                                            const K_Vectors&,
                                                            const Input_para&,
                                                            bool,
                                                            int,
                                                            std::ostream&);

} // namespace LCAO_domain
