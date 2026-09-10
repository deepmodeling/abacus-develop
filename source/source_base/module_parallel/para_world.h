#ifndef PARA_WORLD_H
#define PARA_WORLD_H

#include <cstring>
#include <memory>
#include <string>

#ifdef __MPI
#include "mpi.h"
#endif

namespace Parallel
{

/**
 * @brief Value type describing one MPI communication domain.
 *
 * A ParaWorld couples a domain tag (a short string constant, see
 * para_tag.h) with the communicator, rank and size of the current
 * process inside that domain. It replaces loose globals such as
 * GlobalV::RANK_IN_POOL / POOL_WORLD by an object that functions
 * receive explicitly.
 *
 * The communicator is stored as an opaque handle so that the class
 * layout is identical in serial and MPI builds. Binaries that mix
 * translation units compiled with different __MPI settings (e.g. unit
 * tests linked against the MPI-compiled base library) would otherwise
 * be an ODR violation with undefined behavior. In serial builds rank()
 * always returns 0 and size() always returns 1, so call sites compile
 * unchanged in both serial and MPI builds.
 */
class ParaWorld
{
public:
    virtual ~ParaWorld() = default;

    /// Domain tag string.
    const std::string& tag() const
    {
        return tag_;
    }

    /// Rank of this process inside the domain (0 in serial builds).
    int rank() const
    {
        return rank_;
    }

    /// Number of processes in the domain (1 in serial builds).
    int size() const
    {
        return size_;
    }

    /**
     * @brief True if this process belongs to the domain.
     *
     * In MPI builds this means comm() != MPI_COMM_NULL; in serial
     * builds a default/empty domain is invalid, everything else valid.
     */
    bool valid() const;

#ifdef __MPI
    /// Underlying MPI communicator (MPI builds only).
    MPI_Comm comm() const
    {
        MPI_Comm comm = MPI_COMM_NULL;
        static_assert(sizeof(MPI_Comm) <= sizeof(comm_),
                      "MPI_Comm does not fit into the opaque handle");
        std::memcpy(&comm, &comm_, sizeof(MPI_Comm));
        return comm;
    }
#endif

    /**
     * @brief Build a serial (size=1, rank=0) domain for the given tag.
     *
     * Safe degradation used by tests and by ParaCollection when a tag
     * is not found.
     */
    static ParaWorld serial(const std::string& tag)
    {
        return ParaWorld(tag);
    }

    /**
     * @brief Build a serial domain as a heap-allocated unique_ptr.
     *
     * Convenience factory for ParaCollection::add().
     */
    static std::unique_ptr<ParaWorld> make_serial(const std::string& tag)
    {
        return std::unique_ptr<ParaWorld>(new ParaWorld(tag));
    }

#ifdef __MPI
    /**
     * @brief Build a domain wrapping an MPI communicator.
     *
     * Factory for setup functions that need to create ParaWorld objects
     * from split communicators.
     */
    static ParaWorld make_mpi(const std::string& tag, const MPI_Comm& comm)
    {
        return ParaWorld(tag, comm);
    }

    static std::unique_ptr<ParaWorld> make_mpi_ptr(const std::string& tag, const MPI_Comm& comm)
    {
        return std::unique_ptr<ParaWorld>(new ParaWorld(tag, comm));
    }
#endif

protected:
    /**
     * @brief Construct a serial (single-process) domain.
     *
     * Usable in both serial and MPI builds; in MPI builds the
     * communicator is set to MPI_COMM_SELF. Mainly intended for
     * tests and safe fall-back behavior.
     *
     * @param[in] tag  domain tag string (must be non-empty for a
     *                 meaningful domain; empty tag marks "no domain")
     */
    explicit ParaWorld(const std::string& tag);

#ifdef __MPI
    /**
     * @brief Construct a domain wrapping an existing MPI communicator.
     *
     * @param[in] tag   domain tag string
     * @param[in] comm  MPI communicator (may be MPI_COMM_NULL, which
     *                  yields an invalid domain on this rank)
     */
    ParaWorld(const std::string& tag, const MPI_Comm& comm);
#endif

private:
#ifdef __MPI
    /// Wrap an MPI communicator into the opaque handle storage.
    static void* handle_from_comm(const MPI_Comm& comm)
    {
        void* handle = nullptr;
        std::memcpy(&handle, &comm, sizeof(MPI_Comm));
        return handle;
    }
#endif

    std::string tag_;   ///< domain tag
    int rank_;          ///< rank inside domain
    int size_;          ///< number of processes in domain
    // Opaque communicator handle, present in both serial and MPI builds
    // so that the class layout never depends on the __MPI macro (see the
    // class comment). Never owned/freed here.
    void* comm_;
};

} // namespace Parallel

#endif // PARA_WORLD_H
