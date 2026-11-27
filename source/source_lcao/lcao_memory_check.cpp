#include "source_lcao/lcao_memory_check.h"

#include <complex>
#include <iostream>
#include <limits>
#include <sstream>

#include "source_base/memory_budget.h"
#include "source_base/tool_quit.h"

namespace LCAO_domain
{

namespace
{
const char* decision_to_string(const ModuleBase::MemoryDecision d)
{
    switch (d)
    {
    case ModuleBase::MemoryDecision::Ok:
        return "OK";
    case ModuleBase::MemoryDecision::Warn:
        return "WARN";
    case ModuleBase::MemoryDecision::Abort:
        return "ABORT";
    }
    return "UNKNOWN";
}
} // namespace

template <typename TK>
void check_lcao_dense_memory(const Parallel_Orbitals& pv,
                             const K_Vectors& kv,
                             const Input_para& inp,
                             bool gamma_only,
                             int my_rank,
                             int ranks_per_node,
                             std::ostream& log)
{
    // Local block dimensions are valid only after init_basis_lcao has run;
    // skip the check when the 2D division is not set up yet.
    const int nrloc = pv.para_2d().nrloc();
    const int ncloc = pv.para_2d().ncloc();
    if (nrloc <= 0 || ncloc <= 0)
    {
        return;
    }

    const std::size_t nrloc_s = static_cast<std::size_t>(nrloc);
    const std::size_t ncloc_s = static_cast<std::size_t>(ncloc);
    // Match the column count Setup_Psi::allocate_psi uses to size Psi
    // (source_psi/setup_psi.cpp).  Under MPI it takes pv.ncol_bands -- the
    // local band block, which differs from the global inp.nbands -- so using
    // inp.nbands here would mis-estimate the Psi allocation.  Gamma-only
    // direct solvers instead keep a full ncloc-wide block, a case handled by
    // the gamma_direct branch of estimate_lcao_dense.
    const int nbands = pv.ncol_bands;
    const std::size_t nbands_s = static_cast<std::size_t>(nbands);
    const std::size_t nks_s = static_cast<std::size_t>(kv.get_nks());
    const bool gamma_direct = gamma_only && ModuleBase::is_gamma_direct_ks_solver(inp.ks_solver);
    const bool use_elpa = ModuleBase::is_elpa_ks_solver(inp.ks_solver);
    const std::size_t element_size = sizeof(TK);

    const ModuleBase::LcaoDenseEstimate est = ModuleBase::estimate_lcao_dense(
        nrloc_s, ncloc_s, nbands_s, nks_s, inp.nspin, gamma_only, gamma_direct, use_elpa, element_size);

    const ModuleBase::MemoryStatus status = ModuleBase::get_memory_status();

    // ranks_per_node is injected by the caller (Parallel_Global::mpi_number);
    // guard against a non-positive value (serial build reports 1).
    const std::size_t ranks_per_node_sz = ranks_per_node > 0 ? static_cast<std::size_t>(ranks_per_node) : 1;
    // A /proc/meminfo reading is host-wide and shared by every rank on the
    // node, so divide it by ranks_per_node.  A cgroup limit, however, may
    // already be the rank's own domain (e.g. one container per rank); dividing
    // it again would under-count by the rank count and spuriously abort, so
    // use the cgroup value directly as the per-rank budget.
    const std::size_t budget = status.from_cgroup
                                   ? status.available_bytes
                                   : ModuleBase::per_rank_memory_budget(status.available_bytes, ranks_per_node_sz);

    // The projection (Warn branch) is a future refinement; pass the lower
    // bound as both arguments so decide_memory returns Abort as soon as the
    // lower bound does not fit.
    const ModuleBase::MemoryDecision decision = ModuleBase::decide_memory(
        est.per_rank_total_bytes, est.per_rank_total_bytes, budget);

    const bool is_rank0 = (my_rank == 0);
    // Only an unreadable source (neither cgroup nor /proc/meminfo) makes the
    // budget unknown.  A readable cgroup that computes to zero available --
    // usage at its limit with no reclaimable cache -- is a real, enforceable
    // zero and must still drive the Abort decision below.
    const bool memory_unknown = !status.known;

    // Node-level total is informational only (the decision compares per-rank
    // vs per-rank budget); guard the multiplication against overflow.
    std::size_t node_total_bytes = 0;
    if (est.per_rank_total_bytes <= std::numeric_limits<std::size_t>::max() / ranks_per_node_sz)
    {
        node_total_bytes = est.per_rank_total_bytes * ranks_per_node_sz;
    }
    else
    {
        node_total_bytes = std::numeric_limits<std::size_t>::max();
    }

    if (is_rank0)
    {
        std::ostringstream oss;
        oss << "\n ABACUS: LCAO dense-matrix memory estimate (lower bound)\n";
        oss << "   NLOCAL=" << pv.get_wfc_global_nbasis()
            << "  ncol_bands(local)=" << nbands
            << "  nspin=" << inp.nspin
            << "  nks=" << kv.get_nks()
            << "  gamma_only=" << (gamma_only ? "true" : "false")
            << "  solver=" << inp.ks_solver << "\n";
        oss << "   per-rank local block: nrloc=" << nrloc << " x ncloc=" << ncloc << "\n";
        oss << "   required (per rank): H " << ModuleBase::format_memory_bytes(est.hk_bytes)
            << " + S " << ModuleBase::format_memory_bytes(est.sk_bytes)
            << " + DMK " << ModuleBase::format_memory_bytes(est.dmk_bytes)
            << " + Psi " << ModuleBase::format_memory_bytes(est.psi_bytes);
        if (use_elpa)
        {
            oss << " + ELPA-work " << ModuleBase::format_memory_bytes(est.elpa_work_bytes);
        }
        oss << "\n";
        oss << "     = " << ModuleBase::format_memory_bytes(est.per_rank_total_bytes)
            << " per rank";
        if (!status.from_cgroup)
        {
            // Host memory is shared by the node, so report the node-level total.
            oss << "   (ranks/node=" << ranks_per_node_sz
                << ", node total ~ " << ModuleBase::format_memory_bytes(node_total_bytes) << ")";
        }
        oss << "\n";
        if (memory_unknown)
        {
            oss << "   available: unknown (could not read cgroup or /proc/meminfo)\n";
            oss << "   decision: print-only (available memory unknown; skip abort)\n";
        }
        else
        {
            oss << "   available: " << ModuleBase::format_memory_bytes(budget)
                << " per rank   (source: " << (status.from_cgroup ? "cgroup" : "/proc/meminfo") << ")\n";
            oss << "   decision: " << decision_to_string(decision) << "\n";
        }
        // Write to the running log only, not stdout, so the estimate block
        // does not clutter the screen (WARNING_QUIT still reports an abort).
        log << oss.str() << std::endl;
    }

    if (decision == ModuleBase::MemoryDecision::Abort && !memory_unknown)
    {
        // Every rank must quit; otherwise the non-aborting ranks hang on the
        // next MPI call while rank 0 has already exited.
        ModuleBase::WARNING_QUIT("check_lcao_dense_memory",
            "insufficient memory for LCAO dense matrices; see running log for the memory estimate");
    }
}

template void check_lcao_dense_memory<double>(const Parallel_Orbitals&,
                                              const K_Vectors&,
                                              const Input_para&,
                                              bool,
                                              int,
                                              int,
                                              std::ostream&);
template void check_lcao_dense_memory<std::complex<double>>(const Parallel_Orbitals&,
                                                            const K_Vectors&,
                                                            const Input_para&,
                                                            bool,
                                                            int,
                                                            int,
                                                            std::ostream&);

} // namespace LCAO_domain
