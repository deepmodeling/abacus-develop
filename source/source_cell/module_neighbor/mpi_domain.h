#ifndef MODULE_NEIGHBOR_MPI_DOMAIN_H
#define MODULE_NEIGHBOR_MPI_DOMAIN_H

#include <array>
#include <vector>

#ifdef __MPI
#include <mpi.h>
#endif

namespace ModuleNeighbor
{

#ifdef __MPI
using NeighborMpiComm = MPI_Comm;
#else
// Keep the serial interface buildable without defining global MPI symbols. Some
// non-MPI unit tests include mpi.h after this header, so declaring MPI_Comm here
// would conflict with the real OpenMPI typedef.
using NeighborMpiComm = int;
const NeighborMpiComm kSerialMpiCommWorld = 0;
const NeighborMpiComm kSerialMpiCommNull = -1;
#endif

struct DomainBounds
{
    std::array<double, 3> lower;
    std::array<double, 3> upper;
};

struct MpiAtomRecord
{
    double x;
    double y;
    double z;
    int global_index;
    std::array<int, 3> pbc_shift;
    bool is_ghost;

    MpiAtomRecord();
    MpiAtomRecord(const double x_in,
                  const double y_in,
                  const double z_in,
                  const int global_index_in,
                  const std::array<int, 3>& pbc_shift_in = std::array<int, 3>{{0, 0, 0}},
                  const bool is_ghost_in = false);
};

class MpiDomain
{
  public:
    MpiDomain();
    MpiDomain(const MpiDomain&) = delete;
    MpiDomain& operator=(const MpiDomain&) = delete;
    MpiDomain(MpiDomain&& other);
    MpiDomain& operator=(MpiDomain&& other);
    ~MpiDomain();

    void initialize(NeighborMpiComm parent_comm,
                    const std::array<double, 3>& global_lower,
                    const std::array<double, 3>& global_upper,
                    const double ghost_cutoff,
                    const bool pbc);

    bool initialized() const;
    int rank() const;
    int size() const;
    NeighborMpiComm cart_comm() const;
    const std::array<int, 3>& dims() const;
    const std::array<int, 3>& coords() const;
    const std::array<int, 3>& periods() const;
    const DomainBounds& local_bounds() const;
    double ghost_cutoff() const;

    bool owns(const double x, const double y, const double z) const;
    std::vector<int> select_local_atoms(const std::vector<MpiAtomRecord>& atoms) const;

    // Baseline ghost exchange: each rank publishes its owned atoms, then every
    // rank keeps the remote periodic images that touch its local ghost shell.
    // This is intentionally simple and correctness-first; neighbor-only MPI
    // exchange can replace the internals without changing this interface.
    std::vector<MpiAtomRecord> exchange_ghost_atoms(const std::vector<MpiAtomRecord>& local_atoms) const;

  private:
    void reset();
    void compute_local_bounds();
    bool inside_local(const std::array<double, 3>& point) const;
    bool inside_expanded(const std::array<double, 3>& point) const;
    bool append_ghost_images(const MpiAtomRecord& atom, std::vector<MpiAtomRecord>& ghosts) const;

    NeighborMpiComm cart_comm_;
    bool owns_comm_;
    bool initialized_;
    int rank_;
    int size_;
    std::array<int, 3> dims_;
    std::array<int, 3> coords_;
    std::array<int, 3> periods_;
    std::array<double, 3> global_lower_;
    std::array<double, 3> global_upper_;
    std::array<double, 3> global_length_;
    DomainBounds local_bounds_;
    double ghost_cutoff_;
};

} // namespace ModuleNeighbor

#endif
