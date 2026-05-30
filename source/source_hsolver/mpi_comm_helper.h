#ifndef MPI_COMM_HELPER_H
#define MPI_COMM_HELPER_H

/**
 * @file mpi_comm_helper.h
 * @brief Blocking MPI communication helpers for eigenvalue solvers.
 *
 * Provides type-safe wrappers for common MPI communication patterns:
 *   - reduce_pool: MPI_Allreduce with MPI_SUM
 *   - bcast: MPI_Bcast
 *
 * Also includes CommStrategy enum for adaptive communication strategy
 * selection based on problem size.
 *
 * All MPI operations are guarded by #ifdef __MPI with no-op fallbacks.
 */

#ifdef __MPI
#include <mpi.h>
#include <vector>
#include <cassert>
#endif

#include <complex>
#include <type_traits>

namespace hsolver {

/**
 * @brief Tracks outstanding non-blocking MPI requests and waits for completion.
 *
 * Accumulates MPI_Request handles from non-blocking operations and provides
 * a single wait_all() call to synchronize.
 */
class MPIRequestTracker {
public:
#ifdef __MPI
    /// Add a request to the tracker
    void add(MPI_Request req) { requests_.push_back(req); }

    /// Wait for all outstanding requests to complete
    void wait_all() {
        if (!requests_.empty()) {
            MPI_Waitall(static_cast<int>(requests_.size()),
                        requests_.data(),
                        MPI_STATUSES_IGNORE);
            requests_.clear();
        }
    }

    /// Check if any requests are pending
    bool has_pending() const { return !requests_.empty(); }

    /// Get number of pending requests
    int pending_count() const { return static_cast<int>(requests_.size()); }

    /// Reset the tracker (cancel all pending requests)
    void reset() {
        for (auto& req : requests_) {
            MPI_Cancel(&req);
            MPI_Request_free(&req);
        }
        requests_.clear();
    }

    ~MPIRequestTracker() { reset(); }

private:
    std::vector<MPI_Request> requests_;
#else
    // No-op implementations for serial builds
    void wait_all() {}
    bool has_pending() const { return false; }
    int pending_count() const { return 0; }
    void reset() {}
#endif
};

/**
 * @brief Non-blocking MPI communication operations.
 *
 * All functions are safe to call in serial mode (they become no-ops).
 */
namespace MPICommHelper {

// =========================================================================
// Blocking reduce / broadcast — type-dispatching via mpi_type trait
// =========================================================================

/// Type trait mapping C++ types to MPI_Datatype.
template <typename T> struct mpi_type {
    static constexpr MPI_Datatype value = MPI_BYTE;
};
template <> struct mpi_type<double> {
    static constexpr MPI_Datatype value = MPI_DOUBLE;
};
template <> struct mpi_type<std::complex<double>> {
    static constexpr MPI_Datatype value = MPI_DOUBLE_COMPLEX;
};
template <> struct mpi_type<std::complex<float>> {
    static constexpr MPI_Datatype value = MPI_C_FLOAT_COMPLEX;
};
template <> struct mpi_type<int> {
    static constexpr MPI_Datatype value = MPI_INT;
};

#ifdef __MPI

/**
 * @brief Pool reduce (MPI_SUM). Uses blocking MPI_Allreduce.
 */
template <typename T>
inline void reduce_pool(T* buffer, int count, MPI_Comm comm) {
    MPI_Allreduce(MPI_IN_PLACE, buffer, count, mpi_type<T>::value, MPI_SUM, comm);
}

/**
 * @brief Broadcast. Uses blocking MPI_Bcast.
 */
template <typename T>
inline void bcast(T* buffer, int count, int root, MPI_Comm comm) {
    MPI_Bcast(buffer, count, mpi_type<T>::value, root, comm);
}

#else // !__MPI — serial no-ops

template <typename T>
inline void reduce_pool(T*, int, int) {}
template <typename T>
inline void bcast(T*, int, int, int) {}

#endif

} // namespace MPICommHelper

// =========================================================================
// Communication strategy selection.
// Kept as a simple enum + helper function rather than a separate header
// to avoid over-engineering. Use the resolve() function to select a
// strategy based on problem size.
// =========================================================================

/// Communication strategy for MPI operations.
enum class CommStrategy : int {
    kBlocking    = 0,  ///< Original blocking MPI calls (safe, no extra memory)
    kNonBlocking = 1,  ///< Non-blocking MPI with overlap (default)
    kPipelined   = 2,  ///< Double-buffered pipeline (best for large problems)
    kAdaptive    = 3   ///< Automatic selection based on problem size
};

/// Resolve the effective strategy. If kAdaptive, picks based on problem size:
/// dimensions larger than 100000 use kPipelined, otherwise kNonBlocking.
inline CommStrategy resolve_comm_strategy(CommStrategy strategy,
                                          int dim, int nband) {
    if (strategy != CommStrategy::kAdaptive) {
        return strategy;
    }
    return (dim * nband > 100000) ? CommStrategy::kPipelined
                                  : CommStrategy::kNonBlocking;
}

} // namespace hsolver

#endif // MPI_COMM_HELPER_H
