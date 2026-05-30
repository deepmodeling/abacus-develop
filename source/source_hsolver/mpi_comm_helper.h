#ifndef MPI_COMM_HELPER_H
#define MPI_COMM_HELPER_H

/**
 * @file mpi_comm_helper.h
 * @brief Non-blocking MPI communication helpers for eigenvalue solver optimization.
 *
 * This module provides non-blocking versions of common MPI communication patterns
 * used in the diagonalization module. It enables:
 *   - Non-blocking broadcast (MPI_Ibcast wrapper)
 *   - Non-blocking reduce-to-all (MPI_Iallreduce wrapper)
 *   - Pipelined communication with request tracking
 *
 * All operations are guarded by #ifdef __MPI. When MPI is not available,
 * all functions become no-ops.
 *
 * Usage example:
 * @code
 *   MPIRequestTracker tracker;
 *   tracker.nbcast(vcc, nbase * nband, MPI_DOUBLE_COMPLEX, 0, comm);
 *   // ... do local work while broadcast proceeds ...
 *   tracker.wait_all();
 * @endcode
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
 * Each function posts a non-blocking operation and adds the MPI_Request
 * to the provided tracker. Call tracker.wait_all() to synchronize.
 *
 * All functions are safe to call in serial mode (they become no-ops).
 */
namespace MPICommHelper {

// =========================================================================
// Non-blocking broadcast
// =========================================================================

#ifdef __MPI
/**
 * @brief Non-blocking broadcast (like MPI_Ibcast).
 *
 * @tparam T Element type (must match the MPI_Datatype)
 * @param buffer Pointer to data buffer
 * @param count  Number of elements
 * @param datatype MPI datatype for the elements
 * @param root   Root rank for broadcast
 * @param comm   MPI communicator
 * @param tracker Request tracker to hold the MPI_Request
 */
template <typename T>
inline void nbcast(T* buffer, int count, MPI_Datatype datatype,
                   int root, MPI_Comm comm, MPIRequestTracker& tracker) {
    MPI_Request req;
    MPI_Ibcast(buffer, count, datatype, root, comm, &req);
    tracker.add(req);
}

// Convenience: keep nallreduce for internal use
template <typename T>
inline void nallreduce(T* buffer, int count, MPI_Datatype datatype,
                       MPI_Op op, MPI_Comm comm, MPIRequestTracker& tracker) {
    MPI_Request req;
    MPI_Iallreduce(MPI_IN_PLACE, buffer, count, datatype, op, comm, &req);
    tracker.add(req);
}

// =========================================================================
// Non-blocking reduce / broadcast — type-dispatching via mpi_type trait
// =========================================================================

/// Type trait mapping C++ types to MPI_Datatype.
template <typename T> struct mpi_type {
    static constexpr MPI_Datatype value = MPI_BYTE; // fallback, should not be used
};
template <> struct mpi_type<double> {
    static constexpr MPI_Datatype value = MPI_DOUBLE;
};
template <> struct mpi_type<float> {
    static constexpr MPI_Datatype value = MPI_FLOAT;
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

/**
 * @brief Non-blocking pool reduce (MPI_SUM, non-blocking).
 *
 * Works for double, std::complex<double>, std::complex<float> via mpi_type.
 */
template <typename T>
inline void nreduce_pool(T* buffer, int count,
                         MPI_Comm comm, MPIRequestTracker& tracker) {
    nallreduce(buffer, count, mpi_type<T>::value, MPI_SUM, comm, tracker);
}

/**
 * @brief Non-blocking broadcast (MPI_Ibcast).
 *
 * Works for double, std::complex<double>, std::complex<float> via mpi_type.
 */
template <typename T>
inline void nbcast(T* buffer, int count, int root,
                   MPI_Comm comm, MPIRequestTracker& tracker) {
    MPI_Request req;
    MPI_Ibcast(buffer, count, mpi_type<T>::value, root, comm, &req);
    tracker.add(req);
}

// =========================================================================
// Non-blocking point-to-point (for PLinearTransform optimization)
// =========================================================================

/**
 * @brief Post non-blocking send.
 */
template <typename T>
inline void nsend(const T* buffer, int count, MPI_Datatype datatype,
                  int dest, int tag, MPI_Comm comm, MPIRequestTracker& tracker) {
    MPI_Request req;
    MPI_Issend(buffer, count, datatype, dest, tag, comm, &req);
    tracker.add(req);
}

/**
 * @brief Post non-blocking receive.
 */
template <typename T>
inline void nrecv(T* buffer, int count, MPI_Datatype datatype,
                  int source, int tag, MPI_Comm comm, MPIRequestTracker& tracker) {
    MPI_Request req;
    MPI_Irecv(buffer, count, datatype, source, tag, comm, &req);
    tracker.add(req);
}

#endif // __MPI

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
