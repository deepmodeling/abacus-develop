#ifndef PARA_WORLD_H
#define PARA_WORLD_H

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
 * In serial builds (no __MPI) the communicator member does not
 * exist; rank() always returns 0 and size() always returns 1, so
 * call sites compile unchanged in both serial and MPI builds.
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
        return comm_;
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
    std::string tag_;   ///< domain tag
    int rank_;          ///< rank inside domain
    int size_;          ///< number of processes in domain
#ifdef __MPI
    MPI_Comm comm_;     ///< wrapped communicator (never owned/freed here)
#endif
};

} // namespace Parallel

#endif // PARA_WORLD_H
