#ifndef MODULE_NEIGHBOR_REBUILD_STATE_H
#define MODULE_NEIGHBOR_REBUILD_STATE_H

#include "mpi_domain.h"
#include "source_base/matrix3.h"
#include "source_base/vector3.h"

#include <vector>

class UnitCell;

namespace ModuleNeighbor
{

enum class NeighborRebuildReason
{
    // The cached candidate list remains valid for the current configuration.
    Reuse = 0,
    // No reference configuration has been recorded.
    FirstBuild,
    // A zero skin requests a complete rebuild at every ionic step.
    SkinDisabled,
    // The configured Verlet skin differs from the stored reference.
    SkinChanged,
    // The periodic-boundary mode changed.
    BoundaryChanged,
    // Atom types, counts or ordering changed.
    AtomLayoutChanged,
    // The physical neighbor cutoff changed.
    SearchRadiusChanged,
    // lat0 or the lattice vectors changed.
    LatticeChanged,
    // The maximum atomic displacement exceeded half of the skin.
    DisplacementExceeded
};

// Return a stable diagnostic label for logs and tests.
const char* neighbor_rebuild_reason_name(NeighborRebuildReason reason);

struct NeighborRebuildStats
{
    // Counts accumulated since the most recent reset().
    int rebuild_count = 0;
    int reuse_count = 0;

    // Global maximum displacement for the latest decision, in Bohr.
    double last_max_displacement_bohr = 0.0;
};

/**
 * Track the reference configuration used by a Verlet-style neighbor list.
 *
 * The stored coordinates are fractional for periodic systems and Cartesian
 * (in lat0 units) otherwise. A changed cell, atom layout, boundary mode or
 * physical cutoff invalidates the reference immediately.
 */
class NeighborRebuildState
{
  public:
    explicit NeighborRebuildState(double skin_bohr = 3.0);

    // Discard the reference configuration and accumulated statistics while
    // preserving the currently configured skin.
    void reset();

    // Set the non-negative Verlet skin in Bohr. A zero skin disables reuse.
    void set_skin_bohr(double skin_bohr);

    // Evaluate the local configuration and update last_reason()/displacement
    // statistics. The reference itself changes only after mark_rebuilt().
    bool needs_rebuild(const UnitCell& ucell, double physical_radius_bohr, bool pbc);
#ifdef __MPI
    // Collective counterpart of needs_rebuild(). All ranks receive one rebuild
    // decision using the highest-priority reason and global maximum displacement.
    bool needs_rebuild_mpi(const UnitCell& ucell,
                           double physical_radius_bohr,
                           bool pbc,
                           NeighborMpiComm communicator);
#endif

    // Store the current configuration as the new Verlet reference and record
    // one completed rebuild.
    void mark_rebuilt(const UnitCell& ucell, double physical_radius_bohr, bool pbc);

    // Record that the previously built candidate list was reused.
    void mark_reused();

    bool initialized() const;
    double skin_bohr() const;
    double rebuild_threshold_bohr() const;
    NeighborRebuildReason last_reason() const;
    const NeighborRebuildStats& stats() const;

  private:
    NeighborRebuildReason rebuild_reason(const UnitCell& ucell,
                                         double physical_radius_bohr,
                                         bool pbc) const;
    double max_displacement_bohr(const UnitCell& ucell) const;
    std::vector<ModuleBase::Vector3<double>> collect_reference_coordinates(const UnitCell& ucell,
                                                                           bool pbc) const;

    double skin_bohr_;
    double reference_skin_bohr_;
    bool initialized_;
    bool pbc_;
    int ntype_;
    int nat_;
    double physical_radius_bohr_;
    double reference_lat0_;
    ModuleBase::Matrix3 reference_latvec_;
    std::vector<int> reference_atom_counts_;
    std::vector<ModuleBase::Vector3<double>> reference_coordinates_;
    NeighborRebuildStats stats_;
    NeighborRebuildReason last_reason_;
};

} // namespace ModuleNeighbor

#endif
