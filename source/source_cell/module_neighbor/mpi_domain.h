#ifndef MODULE_NEIGHBOR_MPI_DOMAIN_H
#define MODULE_NEIGHBOR_MPI_DOMAIN_H

#include <array>
#include <map>
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
    // Fractional half-open interval [lower, upper) owned by one Cartesian rank.
    std::array<double, 3> lower;
    std::array<double, 3> upper;
};

struct MpiAtomRecord
{
    // Cartesian position in the same units as the lattice vectors supplied to
    // MpiDomain::initialize_lattice().
    double x;
    double y;
    double z;

    // Identity in the undistributed UnitCell. type and natom may remain -1
    // when only global_index is required by the caller.
    int global_index;
    int type;
    int natom;

    // Integer lattice translation applied to this record. Ghost exchange adds
    // the target-rank periodic translation to any shift already stored here.
    std::array<int, 3> pbc_shift;

    // True only for records received from another spatial owner.
    bool is_ghost;

    MpiAtomRecord();
    MpiAtomRecord(const double x_in,
                  const double y_in,
                  const double z_in,
                  const int global_index_in,
                  const std::array<int, 3>& pbc_shift_in = std::array<int, 3>{{0, 0, 0}},
                  const bool is_ghost_in = false,
                  const int type_in = -1,
                  const int natom_in = -1);
};

struct MpiGhostExchangeStats
{
    // Number of distinct ranks in the fractional ghost communication stencil.
    int neighbor_rank_count = 0;
    // Number of stencil ranks that receive at least one record from this rank.
    int nonempty_send_rank_count = 0;
    // Atom records sent across all target ranks. The same local atom is counted
    // more than once when it overlaps multiple rank ghost shells.
    int sent_atom_count = 0;
    // Ghost records retained after receiving and deduplicating payloads.
    int received_ghost_count = 0;
    // Number of packed double elements sent, including all fields of every
    // atom record. This is the MPI payload volume, not an atom count.
    int sent_payload_count = 0;
    // Number of packed double elements received before atom unpacking and
    // duplicate removal.
    int received_payload_count = 0;
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

    // Initialize an orthogonal domain whose lattice vectors are the three
    // global box extents. The object duplicates parent_comm through a Cartesian
    // communicator and owns that duplicate until reset or destruction.
    void initialize(NeighborMpiComm parent_comm,
                    const std::array<double, 3>& global_lower,
                    const std::array<double, 3>& global_upper,
                    const double ghost_cutoff,
                    const bool pbc);

    // Initialize a general lattice domain. origin, lattice_vectors and
    // ghost_cutoff must use one consistent Cartesian length unit. Spatial
    // ownership is evaluated in fractional coordinates, including triclinic
    // cells; pbc controls both rank topology and coordinate wrapping.
    void initialize_lattice(NeighborMpiComm parent_comm,
                            const std::array<double, 3>& origin,
                            const std::array<std::array<double, 3>, 3>& lattice_vectors,
                            const double ghost_cutoff,
                            const bool pbc);

    bool initialized() const;
    int rank() const;
    int size() const;
    NeighborMpiComm cart_comm() const;
    const std::array<int, 3>& dims() const;
    const std::array<int, 3>& coords() const;
    const std::array<int, 3>& periods() const;
    // Bounds are stored in fractional coordinates. This representation remains
    // valid for both orthogonal and triclinic lattice domains.
    const DomainBounds& local_bounds() const;
    const std::array<double, 3>& fractional_ghost_padding() const;
    double ghost_cutoff() const;

    std::array<double, 3> cartesian_to_fractional(const double x, const double y, const double z) const;
    std::array<double, 3> fractional_to_cartesian(const double fx, const double fy, const double fz) const;

    // Return whether the point belongs to this rank's fractional half-open
    // interval. Periodic coordinates are normalized before the ownership test.
    bool owns(const double x, const double y, const double z) const;

    // Return indices of input records owned by this rank, preserving input order.
    std::vector<int> select_local_atoms(const std::vector<MpiAtomRecord>& atoms) const;

    // Send local atom images that overlap the required Cartesian-rank ghost
    // stencil and return deduplicated ghost records received by this rank.
    // local_atoms must contain only records owned by the calling rank.
    std::vector<MpiAtomRecord> exchange_ghost_atoms(const std::vector<MpiAtomRecord>& local_atoms) const;

    // The optional statistics describe this rank's exchange only; callers must
    // perform MPI reductions when global totals or maxima are required.
    std::vector<MpiAtomRecord> exchange_ghost_atoms(const std::vector<MpiAtomRecord>& local_atoms,
                                                    MpiGhostExchangeStats* stats) const;

  private:
    void reset();
    void compute_local_bounds();
    std::array<double, 3> normalize_fractional(const std::array<double, 3>& fractional) const;
    bool inside_local(const std::array<double, 3>& fractional) const;
    DomainBounds bounds_for_coords(const std::array<int, 3>& coords) const;
    bool inside_local_for_bounds(const DomainBounds& bounds,
                                 const std::array<int, 3>& coords,
                                 const std::array<double, 3>& fractional) const;
    bool inside_expanded_for_bounds(const DomainBounds& bounds,
                                    const std::array<double, 3>& fractional) const;
    bool neighbor_from_direction(const std::array<int, 3>& direction,
                                 int& neighbor_rank,
                                 std::array<int, 3>& target_coords,
                                 std::array<int, 3>& image_shift) const;

    // Number of rank layers needed on one axis to cover the fractional ghost
    // padding. This may exceed one when the cutoff is wider than a subdomain.
    int neighbor_offset_limit(int axis) const;
    std::vector<int> neighbor_ranks() const;
    std::map<int, std::vector<MpiAtomRecord>>
    build_neighbor_send_records(const std::vector<MpiAtomRecord>& local_atoms) const;

    NeighborMpiComm cart_comm_;
    bool owns_comm_;
    bool initialized_;
    int rank_;
    int size_;
    std::array<int, 3> dims_;
    std::array<int, 3> coords_;
    std::array<int, 3> periods_;
    std::array<double, 3> origin_;
    std::array<std::array<double, 3>, 3> lattice_vectors_;
    // reciprocal_rows_[i] maps a Cartesian displacement to fractional axis i.
    std::array<std::array<double, 3>, 3> reciprocal_rows_;
    std::array<double, 3> fractional_ghost_padding_;
    DomainBounds local_bounds_;
    double ghost_cutoff_;
};

} // namespace ModuleNeighbor

#endif
