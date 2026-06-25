#include "neighbor_rebuild_state.h"

#include "source_cell/unitcell.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace
{
constexpr double metadata_tolerance = 1.0e-12;

bool changed(const double lhs, const double rhs)
{
    return std::abs(lhs - rhs) > metadata_tolerance;
}

bool lattice_changed(const ModuleBase::Matrix3& lhs, const ModuleBase::Matrix3& rhs)
{
    return changed(lhs.e11, rhs.e11) || changed(lhs.e12, rhs.e12) || changed(lhs.e13, rhs.e13)
           || changed(lhs.e21, rhs.e21) || changed(lhs.e22, rhs.e22) || changed(lhs.e23, rhs.e23)
           || changed(lhs.e31, rhs.e31) || changed(lhs.e32, rhs.e32) || changed(lhs.e33, rhs.e33);
}

double norm(const ModuleBase::Vector3<double>& value)
{
    return std::sqrt(value.x * value.x + value.y * value.y + value.z * value.z);
}

int rebuild_reason_priority(const ModuleNeighbor::NeighborRebuildReason reason)
{
    switch (reason)
    {
        case ModuleNeighbor::NeighborRebuildReason::Reuse:
            return 0;
        case ModuleNeighbor::NeighborRebuildReason::LatticeChanged:
            return 1;
        case ModuleNeighbor::NeighborRebuildReason::SearchRadiusChanged:
            return 2;
        case ModuleNeighbor::NeighborRebuildReason::AtomLayoutChanged:
            return 3;
        case ModuleNeighbor::NeighborRebuildReason::BoundaryChanged:
            return 4;
        case ModuleNeighbor::NeighborRebuildReason::SkinChanged:
            return 5;
        case ModuleNeighbor::NeighborRebuildReason::SkinDisabled:
            return 6;
        case ModuleNeighbor::NeighborRebuildReason::FirstBuild:
            return 7;
        case ModuleNeighbor::NeighborRebuildReason::DisplacementExceeded:
            return 8;
    }
    throw std::runtime_error("Unknown neighbor rebuild reason.");
}

ModuleNeighbor::NeighborRebuildReason rebuild_reason_from_priority(const int priority)
{
    switch (priority)
    {
        case 0:
            return ModuleNeighbor::NeighborRebuildReason::Reuse;
        case 1:
            return ModuleNeighbor::NeighborRebuildReason::LatticeChanged;
        case 2:
            return ModuleNeighbor::NeighborRebuildReason::SearchRadiusChanged;
        case 3:
            return ModuleNeighbor::NeighborRebuildReason::AtomLayoutChanged;
        case 4:
            return ModuleNeighbor::NeighborRebuildReason::BoundaryChanged;
        case 5:
            return ModuleNeighbor::NeighborRebuildReason::SkinChanged;
        case 6:
            return ModuleNeighbor::NeighborRebuildReason::SkinDisabled;
        case 7:
            return ModuleNeighbor::NeighborRebuildReason::FirstBuild;
        case 8:
            return ModuleNeighbor::NeighborRebuildReason::DisplacementExceeded;
    }
    throw std::runtime_error("Invalid collective neighbor rebuild reason.");
}
} // namespace

namespace ModuleNeighbor
{

const char* neighbor_rebuild_reason_name(const NeighborRebuildReason reason)
{
    switch (reason)
    {
        case NeighborRebuildReason::Reuse:
            return "reuse";
        case NeighborRebuildReason::FirstBuild:
            return "first_build";
        case NeighborRebuildReason::SkinDisabled:
            return "skin_disabled";
        case NeighborRebuildReason::SkinChanged:
            return "skin_changed";
        case NeighborRebuildReason::BoundaryChanged:
            return "boundary_changed";
        case NeighborRebuildReason::AtomLayoutChanged:
            return "atom_layout_changed";
        case NeighborRebuildReason::SearchRadiusChanged:
            return "search_radius_changed";
        case NeighborRebuildReason::LatticeChanged:
            return "lattice_changed";
        case NeighborRebuildReason::DisplacementExceeded:
            return "displacement_exceeded";
    }
    return "unknown";
}

NeighborRebuildState::NeighborRebuildState(const double skin_bohr)
    : skin_bohr_(skin_bohr),
      reference_skin_bohr_(skin_bohr),
      initialized_(false),
      pbc_(false),
      ntype_(0),
      nat_(0),
      physical_radius_bohr_(0.0),
      reference_lat0_(0.0),
      last_reason_(NeighborRebuildReason::FirstBuild)
{
    if (!std::isfinite(skin_bohr) || skin_bohr < 0.0)
    {
        throw std::invalid_argument("Neighbor-list skin must be finite and non-negative.");
    }
}

void NeighborRebuildState::reset()
{
    initialized_ = false;
    pbc_ = false;
    ntype_ = 0;
    nat_ = 0;
    physical_radius_bohr_ = 0.0;
    reference_lat0_ = 0.0;
    reference_atom_counts_.clear();
    reference_coordinates_.clear();
    stats_ = NeighborRebuildStats();
    last_reason_ = NeighborRebuildReason::FirstBuild;
}

void NeighborRebuildState::set_skin_bohr(const double skin_bohr)
{
    if (!std::isfinite(skin_bohr) || skin_bohr < 0.0)
    {
        throw std::invalid_argument("Neighbor-list skin must be finite and non-negative.");
    }
    skin_bohr_ = skin_bohr;
}

bool NeighborRebuildState::needs_rebuild(const UnitCell& ucell,
                                         const double physical_radius_bohr,
                                         const bool pbc)
{
    last_reason_ = rebuild_reason(ucell, physical_radius_bohr, pbc);
    if (last_reason_ != NeighborRebuildReason::Reuse)
    {
        stats_.last_max_displacement_bohr = 0.0;
        return true;
    }

    stats_.last_max_displacement_bohr = max_displacement_bohr(ucell);
    if (stats_.last_max_displacement_bohr > rebuild_threshold_bohr() + metadata_tolerance)
    {
        last_reason_ = NeighborRebuildReason::DisplacementExceeded;
        return true;
    }
    return false;
}

#ifdef __MPI
bool NeighborRebuildState::needs_rebuild_mpi(const UnitCell& ucell,
                                             const double physical_radius_bohr,
                                             const bool pbc,
                                             const NeighborMpiComm communicator)
{
    const NeighborRebuildReason local_reason = rebuild_reason(ucell, physical_radius_bohr, pbc);
    int local_reason_priority = rebuild_reason_priority(local_reason);
    int global_reason_priority = local_reason_priority;
    if (communicator != MPI_COMM_NULL)
    {
        MPI_Allreduce(&local_reason_priority,
                      &global_reason_priority,
                      1,
                      MPI_INT,
                      MPI_MAX,
                      communicator);
    }
    last_reason_ = rebuild_reason_from_priority(global_reason_priority);
    if (last_reason_ != NeighborRebuildReason::Reuse)
    {
        stats_.last_max_displacement_bohr = 0.0;
        return true;
    }

    const double local_max = max_displacement_bohr(ucell);
    double global_max = local_max;
    if (communicator != MPI_COMM_NULL)
    {
        MPI_Allreduce(&local_max, &global_max, 1, MPI_DOUBLE, MPI_MAX, communicator);
    }
    stats_.last_max_displacement_bohr = global_max;
    if (global_max > rebuild_threshold_bohr() + metadata_tolerance)
    {
        last_reason_ = NeighborRebuildReason::DisplacementExceeded;
        return true;
    }
    return false;
}
#endif

void NeighborRebuildState::mark_rebuilt(const UnitCell& ucell,
                                        const double physical_radius_bohr,
                                        const bool pbc)
{
    initialized_ = true;
    pbc_ = pbc;
    ntype_ = ucell.ntype;
    nat_ = ucell.nat;
    physical_radius_bohr_ = physical_radius_bohr;
    reference_skin_bohr_ = skin_bohr_;
    reference_lat0_ = ucell.lat0;
    reference_latvec_ = ucell.latvec;
    reference_atom_counts_.resize(ucell.ntype);
    for (int type = 0; type < ucell.ntype; ++type)
    {
        reference_atom_counts_[type] = ucell.atoms[type].na;
    }
    reference_coordinates_ = collect_reference_coordinates(ucell, pbc);
    stats_.last_max_displacement_bohr = 0.0;
    ++stats_.rebuild_count;
}

void NeighborRebuildState::mark_reused()
{
    if (!initialized_)
    {
        throw std::runtime_error("Cannot reuse an uninitialized neighbor list.");
    }
    ++stats_.reuse_count;
    last_reason_ = NeighborRebuildReason::Reuse;
}

bool NeighborRebuildState::initialized() const
{
    return initialized_;
}

double NeighborRebuildState::skin_bohr() const
{
    return skin_bohr_;
}

double NeighborRebuildState::rebuild_threshold_bohr() const
{
    return 0.5 * skin_bohr_;
}

NeighborRebuildReason NeighborRebuildState::last_reason() const
{
    return last_reason_;
}

const NeighborRebuildStats& NeighborRebuildState::stats() const
{
    return stats_;
}

NeighborRebuildReason NeighborRebuildState::rebuild_reason(const UnitCell& ucell,
                                                           const double physical_radius_bohr,
                                                           const bool pbc) const
{
    if (!std::isfinite(physical_radius_bohr) || physical_radius_bohr < 0.0)
    {
        throw std::invalid_argument("Neighbor search radius must be finite and non-negative.");
    }
    if (!initialized_)
    {
        return NeighborRebuildReason::FirstBuild;
    }
    if (skin_bohr_ <= metadata_tolerance)
    {
        return NeighborRebuildReason::SkinDisabled;
    }
    if (changed(skin_bohr_, reference_skin_bohr_))
    {
        return NeighborRebuildReason::SkinChanged;
    }
    if (pbc != pbc_)
    {
        return NeighborRebuildReason::BoundaryChanged;
    }
    if (ucell.ntype != ntype_ || ucell.nat != nat_
        || static_cast<int>(reference_atom_counts_.size()) != ucell.ntype)
    {
        return NeighborRebuildReason::AtomLayoutChanged;
    }
    if (changed(physical_radius_bohr, physical_radius_bohr_))
    {
        return NeighborRebuildReason::SearchRadiusChanged;
    }
    if (changed(ucell.lat0, reference_lat0_) || lattice_changed(ucell.latvec, reference_latvec_))
    {
        return NeighborRebuildReason::LatticeChanged;
    }
    for (int type = 0; type < ucell.ntype; ++type)
    {
        if (reference_atom_counts_[type] != ucell.atoms[type].na)
        {
            return NeighborRebuildReason::AtomLayoutChanged;
        }
    }
    return NeighborRebuildReason::Reuse;
}

double NeighborRebuildState::max_displacement_bohr(const UnitCell& ucell) const
{
    const std::vector<ModuleBase::Vector3<double>> current = collect_reference_coordinates(ucell, pbc_);
    if (current.size() != reference_coordinates_.size())
    {
        return std::numeric_limits<double>::infinity();
    }

    double maximum = 0.0;
    for (std::size_t index = 0; index < current.size(); ++index)
    {
        ModuleBase::Vector3<double> delta = current[index] - reference_coordinates_[index];
        if (pbc_)
        {
            delta.x -= std::round(delta.x);
            delta.y -= std::round(delta.y);
            delta.z -= std::round(delta.z);
            delta = delta * ucell.latvec;
        }
        maximum = std::max(maximum, norm(delta) * ucell.lat0);
    }
    return maximum;
}

std::vector<ModuleBase::Vector3<double>>
NeighborRebuildState::collect_reference_coordinates(const UnitCell& ucell, const bool pbc) const
{
    std::vector<ModuleBase::Vector3<double>> coordinates;
    coordinates.reserve(ucell.nat);
    const ModuleBase::Matrix3 inverse_lattice = pbc ? ucell.latvec.Inverse() : ModuleBase::Matrix3();
    for (int type = 0; type < ucell.ntype; ++type)
    {
        for (int natom = 0; natom < ucell.atoms[type].na; ++natom)
        {
            const ModuleBase::Vector3<double>& tau = ucell.atoms[type].tau[natom];
            coordinates.push_back(pbc ? tau * inverse_lattice : tau);
        }
    }
    return coordinates;
}

} // namespace ModuleNeighbor
