#include "source_cell/module_neighlist/domain_decomposition.h"
#include "source_base/parallel_cell.h"
#include "source_cell/unitcell.h"
#include "source_cell/mdcell.h"
#include "source_cell/module_neighlist/neighbor_search.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstring>
#include <limits>
#include <map>
#include <stdexcept>
#include <utility>

namespace
{
const int tag_migration_count_base = 100;
const int tag_migration_atoms_base = 200;
const int tag_ghost_count = 9100;
const int tag_ghost_atoms = 9101;
const int tag_ghost_positions = 9110;
const int tag_ghost_force_size = 9200;
const int tag_ghost_forces = 9201;
}

// MPI-dependent transport and communicator management are confined to helpers.

DomainDecomposition::DomainDecomposition()
    : rank_(0), size_(1), dims_{{1, 1, 1}}, coords_{{0, 0, 0}},
      margin_{{0.0, 0.0, 0.0}}, lat0_(1.0), cutoff_(0.0), skin_(0.0)
{
}

DomainDecomposition::~DomainDecomposition()
{
#ifdef __MPI
    if (owns_cart_comm_ && cart_comm_ != MPI_COMM_NULL)
    {
        MPI_Comm_free(&cart_comm_);
    }
    owns_cart_comm_ = false;
#endif
}

void DomainDecomposition::init(const ModuleBase::CommunicationDomain& comm_domain,
                               const ModuleBase::Matrix3& latvec,
                               double lat0,
                               double cutoff,
                               double skin)
{
    initialized_ = false;
    rank_ = comm_domain.rank();
    size_ = comm_domain.size();
    dims_ = {{1, 1, 1}};
    coords_ = {{0, 0, 0}};
#ifdef __MPI
    if (owns_cart_comm_ && cart_comm_ != MPI_COMM_NULL)
    {
        MPI_Comm_free(&cart_comm_);
    }
    owns_cart_comm_ = false;
    comm_ = comm_domain.communicator();
    if (comm_ == MPI_COMM_NULL)
    {
        throw std::runtime_error("DomainDecomposition requires a valid communication domain.");
    }
    int dims[3] = {0, 0, 0};
    MPI_Dims_create(size_, 3, dims);
    int periods[3] = {1, 1, 1};
    MPI_Cart_create(comm_, 3, dims, periods, 0, &cart_comm_);
    owns_cart_comm_ = cart_comm_ != MPI_COMM_NULL;
    MPI_Comm_rank(cart_comm_, &rank_);
    int coords[3] = {0, 0, 0};
    MPI_Cart_coords(cart_comm_, rank_, 3, coords);
    std::copy(dims, dims + 3, dims_.begin());
    std::copy(coords, coords + 3, coords_.begin());
#endif
    update_geometry_(latvec, lat0, cutoff, skin);
}

void DomainDecomposition::migrate_owned_atoms(MDCell& cell)
{
    synchronize_geometry_(cell);
    std::vector<LocalAtom>& owned_atoms = cell.owned_atoms_;
    const int direction_count = 6;
    const int axis[direction_count] = {0, 0, 1, 1, 2, 2};
    const int step[direction_count] = {-1, 1, -1, 1, -1, 1};
    std::array<int, direction_count> neighbors;
    for (int idir = 0; idir < direction_count; ++idir)
    {
        std::array<int, 3> neighbor_coords = coords_;
        neighbor_coords[axis[idir]] = positive_mod(neighbor_coords[axis[idir]] + step[idir], dims_[axis[idir]]);
        neighbors[idir] = rank_from_coords(neighbor_coords);
    }

    std::vector<LocalAtom> pending_atoms;
    pending_atoms.swap(owned_atoms);
    std::vector<LocalAtom> retained_atoms;
    retained_atoms.reserve(pending_atoms.size());
    const std::array<int, 3> no_shift = {{0, 0, 0}};

    long long global_outgoing = 0;
    do
    {
        std::array<std::vector<PackedAtom>, direction_count> send_atoms;
        for (std::size_t i = 0; i < pending_atoms.size(); ++i)
        {
            LocalAtom atom = std::move(pending_atoms[i]);
            atom.frac = wrapped_frac_from_cart(atom.cart);
            atom.cart = atom.frac * latvec_;

            std::array<int, 3> owner_coords;
            const double frac[3] = {atom.frac.x, atom.frac.y, atom.frac.z};
            for (int idim = 0; idim < 3; ++idim)
            {
                owner_coords[idim] = std::min(static_cast<int>(std::floor(frac[idim] * dims_[idim])), dims_[idim] - 1);
            }
            atom.owner_rank = rank_from_coords(owner_coords);
            if (atom.owner_rank == rank_)
            {
                retained_atoms.push_back(std::move(atom));
                continue;
            }

            int direction = -1;
            for (int idim = 0; idim < 3 && direction < 0; ++idim)
            {
                int delta = owner_coords[idim] - coords_[idim];
                if (delta > dims_[idim] / 2) delta -= dims_[idim];
                if (delta < -dims_[idim] / 2) delta += dims_[idim];
                if (delta != 0) direction = 2 * idim + (delta > 0 ? 1 : 0);
            }
            assert(direction >= 0);
            send_atoms[direction].push_back(pack_atom(atom, no_shift));
        }
        pending_atoms.clear();

        std::array<int, direction_count> send_counts;
        long long local_outgoing = 0;
        for (int idir = 0; idir < direction_count; ++idir)
        {
            const std::size_t bytes = send_atoms[idir].size() * sizeof(PackedAtom);
            if (bytes > static_cast<std::size_t>(std::numeric_limits<int>::max()))
            {
                throw std::overflow_error("DomainDecomposition migration send count exceeds int range.");
            }
            send_counts[idir] = static_cast<int>(bytes);
            local_outgoing += static_cast<long long>(send_atoms[idir].size());
        }
        global_outgoing = local_outgoing;
#ifdef __MPI
        MPI_Allreduce(MPI_IN_PLACE, &global_outgoing, 1, MPI_LONG_LONG, MPI_SUM, comm_);
#endif
        if (global_outgoing == 0) break;

        std::array<std::vector<PackedAtom>, direction_count> recv_atoms;
        exchange_migration_(neighbors, send_atoms, send_counts, recv_atoms);

        for (int idir = 0; idir < direction_count; ++idir)
        {
            for (std::size_t i = 0; i < recv_atoms[idir].size(); ++i)
            {
                pending_atoms.push_back(unpack_owned_atom(recv_atoms[idir][i]));
            }
        }
    } while (global_outgoing > 0);

    owned_atoms.swap(retained_atoms);
    exchange_ghost_atoms(cell);
}

void DomainDecomposition::exchange_ghost_atoms(MDCell& cell)
{
    synchronize_geometry_(cell);
    const std::vector<LocalAtom>& owned_atoms = cell.owned_atoms_;
    std::vector<LocalAtom>& ghost_atoms = cell.ghost_atoms_;
    ghost_atoms.clear();
    ghost_layout_valid_ = false;
    ghost_slots_.clear();
    build_ghost_exchange_slots(ghost_slots_);
    std::vector<GhostExchangeSlot>& slots = ghost_slots_;

    const int nlayer[3] = {neighbor_layer(0), neighbor_layer(1), neighbor_layer(2)};
    const int span_y = 2 * nlayer[1] + 1;
    const int span_z = 2 * nlayer[2] + 1;
    const int lookup_size = (2 * nlayer[0] + 1) * span_y * span_z;
    std::vector<int> slot_lookup(static_cast<std::size_t>(lookup_size), -1);
    for (std::size_t islot = 0; islot < slots.size(); ++islot)
    {
        const std::array<int, 3>& offset = slots[islot].offset;
        const int index = (offset[0] + nlayer[0]) * span_y * span_z
                          + (offset[1] + nlayer[1]) * span_z
                          + (offset[2] + nlayer[2]);
        slot_lookup[static_cast<std::size_t>(index)] = static_cast<int>(islot);
    }

    const auto collect_offsets = [&](const LocalAtom& atom, const int dim, std::vector<int>& offsets) {
        offsets.clear();
        for (int delta = -nlayer[dim]; delta <= nlayer[dim]; ++delta)
        {
            if (delta == 0)
            {
                offsets.push_back(0);
                continue;
            }

            const std::array<int, 3> offset = {{dim == 0 ? delta : 0,
                                                dim == 1 ? delta : 0,
                                                dim == 2 ? delta : 0}};
            std::array<int, 3> target_coords;
            std::array<int, 3> image_shift;
            target_for_offset(offset, target_coords, image_shift);

            const double frac_values[3] = {
                atom.frac.x + image_shift[0],
                atom.frac.y + image_shift[1],
                atom.frac.z + image_shift[2]
            };
            const double lo = static_cast<double>(target_coords[dim]) / dims_[dim];
            const double hi = static_cast<double>(target_coords[dim] + 1) / dims_[dim];
            if (frac_values[dim] >= lo - margin_[dim] &&
                frac_values[dim] < hi + margin_[dim])
            {
                offsets.push_back(delta);
            }
        }
    };

    std::vector<std::vector<PackedAtom>> send_buffers(slots.size());
    std::vector<int> x_offsets;
    std::vector<int> y_offsets;
    std::vector<int> z_offsets;
    for (size_t iat = 0; iat < owned_atoms.size(); ++iat)
    {
        const LocalAtom& atom = owned_atoms[iat];
        collect_offsets(atom, 0, x_offsets);
        collect_offsets(atom, 1, y_offsets);
        collect_offsets(atom, 2, z_offsets);

        for (const int dx : x_offsets)
        {
            for (const int dy : y_offsets)
            {
                for (const int dz : z_offsets)
                {
                    if (dx == 0 && dy == 0 && dz == 0)
                    {
                        continue;
                    }
                    const int lookup_index = (dx + nlayer[0]) * span_y * span_z
                                             + (dy + nlayer[1]) * span_z
                                             + (dz + nlayer[2]);
                    const int slot_index = slot_lookup[static_cast<std::size_t>(lookup_index)];
                    assert(slot_index >= 0);
                    const GhostExchangeSlot& slot = slots[static_cast<std::size_t>(slot_index)];
                    send_buffers[static_cast<std::size_t>(slot_index)].push_back(pack_atom(atom, slot.image_shift));
                    slots[static_cast<std::size_t>(slot_index)].send_atom_indices.push_back(static_cast<int>(iat));
                }
            }
        }
    }

    for (std::size_t islot = 0; islot < slots.size(); ++islot)
    {
        GhostExchangeSlot& slot = slots[islot];
        const std::vector<PackedAtom>& send_atoms = send_buffers[islot];
        slot.ghost_begin = ghost_atoms.size();

        if (send_atoms.size() > static_cast<std::size_t>(std::numeric_limits<int>::max()))
        {
            throw std::overflow_error("DomainDecomposition ghost send count exceeds int range.");
        }

        int send_count = static_cast<int>(send_atoms.size());
        int recv_count = 0;
        exchange_bytes_(&send_count, sizeof(send_count), &recv_count, sizeof(recv_count),
                        slot.send_rank, slot.recv_rank, tag_ghost_count);
        if (recv_count < 0)
        {
            throw std::runtime_error("Invalid ghost receive count.");
        }

        std::vector<PackedAtom> recv_atoms(static_cast<size_t>(recv_count));
        const std::size_t send_bytes_size = send_atoms.size() * sizeof(PackedAtom);
        const std::size_t recv_bytes_size = recv_atoms.size() * sizeof(PackedAtom);
        if (send_bytes_size > static_cast<std::size_t>(std::numeric_limits<int>::max()) ||
            recv_bytes_size > static_cast<std::size_t>(std::numeric_limits<int>::max()))
        {
            throw std::overflow_error("DomainDecomposition ghost message exceeds MPI int byte count range.");
        }
        const int send_bytes = static_cast<int>(send_bytes_size);
        const int recv_bytes = static_cast<int>(recv_bytes_size);

        exchange_bytes_(send_atoms.data(), send_bytes, recv_atoms.data(), recv_bytes,
                        slot.send_rank, slot.recv_rank, tag_ghost_atoms);

        for (size_t i = 0; i < recv_atoms.size(); ++i)
        {
            ghost_atoms.push_back(unpack_ghost_atom(recv_atoms[i]));
        }
        slot.ghost_count = recv_count;
    }
    ghost_layout_valid_ = true;
    cell.clear_forces_(ghost_atoms);
    cell.neighbor_layout_valid_ = false;
}

void DomainDecomposition::update_ghost_atom_positions(MDCell& cell)
{
    synchronize_geometry_(cell);
    if (!ghost_layout_valid_)
    {
        exchange_ghost_atoms(cell);
        return;
    }

    const std::vector<LocalAtom>& owned_atoms = cell.owned_atoms_;
    std::vector<LocalAtom>& ghost_atoms = cell.ghost_atoms_;
    for (std::size_t islot = 0; islot < ghost_slots_.size(); ++islot)
    {
        const GhostExchangeSlot& slot = ghost_slots_[islot];
        std::vector<double> send_frac(3 * slot.send_atom_indices.size(), 0.0);
        for (std::size_t i = 0; i < slot.send_atom_indices.size(); ++i)
        {
            const LocalAtom& atom = owned_atoms[static_cast<std::size_t>(slot.send_atom_indices[i])];
            send_frac[3 * i] = atom.frac.x;
            send_frac[3 * i + 1] = atom.frac.y;
            send_frac[3 * i + 2] = atom.frac.z;
        }
        std::vector<double> recv_frac(3 * static_cast<std::size_t>(slot.ghost_count), 0.0);
        exchange_bytes_(send_frac.data(), send_frac.size() * sizeof(double),
                        recv_frac.data(), recv_frac.size() * sizeof(double),
                        slot.send_rank, slot.recv_rank, tag_ghost_positions);
        for (int i = 0; i < slot.ghost_count; ++i)
        {
            LocalAtom& ghost = ghost_atoms[slot.ghost_begin + static_cast<std::size_t>(i)];
            ghost.frac.set(recv_frac[3 * i], recv_frac[3 * i + 1], recv_frac[3 * i + 2]);
            const std::array<int, 3>& image_shift = slot.send_rank == rank_ && slot.recv_rank == rank_
                                                         ? slot.image_shift
                                                         : slot.recv_image_shift;
            const ModuleBase::Vector3<double> image_frac(ghost.frac.x + image_shift[0],
                                                          ghost.frac.y + image_shift[1],
                                                          ghost.frac.z + image_shift[2]);
            ghost.cart = image_frac * latvec_;
            ghost.force.set(0.0, 0.0, 0.0);
        }
    }
    cell.clear_forces_(ghost_atoms);
}

void DomainDecomposition::accumulate_ghost_forces(MDCell& cell)
{
    // Do not recreate geometry or ghost mappings between force evaluation and return.
    if (!initialized_ || !ghost_layout_valid_)
    {
        throw std::runtime_error("Ghost forces require an initialized exchange layout.");
    }

    std::vector<LocalAtom>& owned_atoms = cell.owned_atoms_;
    std::vector<LocalAtom>& ghost_atoms = cell.ghost_atoms_;
    std::map<std::pair<int, std::int64_t>, std::size_t> owned_lookup;
    for (std::size_t iat = 0; iat < owned_atoms.size(); ++iat)
    {
        const LocalAtom& atom = owned_atoms[iat];
        owned_lookup[std::make_pair(atom.type, atom.type_index)] = iat;
    }

    std::vector<int> peer_ranks;
    peer_ranks.reserve(ghost_slots_.size() * 2);
    for (const GhostExchangeSlot& slot : ghost_slots_)
    {
        if (slot.send_rank != rank_)
        {
            peer_ranks.push_back(slot.send_rank);
        }
        if (slot.recv_rank != rank_)
        {
            peer_ranks.push_back(slot.recv_rank);
        }
    }
    std::sort(peer_ranks.begin(), peer_ranks.end());
    peer_ranks.erase(std::unique(peer_ranks.begin(), peer_ranks.end()), peer_ranks.end());

    std::vector<std::vector<ForceRecord> > send_buffers(peer_ranks.size());
    for (std::size_t iat = 0; iat < ghost_atoms.size(); ++iat)
    {
        LocalAtom& atom = ghost_atoms[iat];
        if (atom.owner_rank == rank_)
        {
            const std::map<std::pair<int, std::int64_t>, std::size_t>::const_iterator found
                = owned_lookup.find(std::make_pair(atom.type, atom.type_index));
            if (found == owned_lookup.end())
            {
                throw std::runtime_error("Cannot match a local ghost force to an owned atom.");
            }
            owned_atoms[found->second].force += atom.force;
        }
        else
        {
            ForceRecord record;
            record.type = atom.type;
            record.type_index = atom.type_index;
            record.force[0] = atom.force.x;
            record.force[1] = atom.force.y;
            record.force[2] = atom.force.z;
            const std::vector<int>::const_iterator peer = std::lower_bound(peer_ranks.begin(),
                                                                             peer_ranks.end(),
                                                                             atom.owner_rank);
            if (peer == peer_ranks.end() || *peer != atom.owner_rank)
            {
                throw std::runtime_error("Ghost force owner is outside the ghost communication stencil.");
            }
            send_buffers[static_cast<std::size_t>(peer - peer_ranks.begin())].push_back(record);
        }
        atom.force.set(0.0, 0.0, 0.0);
    }

    for (std::size_t ipeer = 0; ipeer < peer_ranks.size(); ++ipeer)
    {
        const std::vector<ForceRecord>& send_records = send_buffers[ipeer];
        const std::size_t bytes = send_records.size() * sizeof(ForceRecord);
        if (bytes > static_cast<std::size_t>(std::numeric_limits<int>::max()))
        {
            throw std::overflow_error("DomainDecomposition ghost force message exceeds MPI int range.");
        }
        const int send_bytes = static_cast<int>(bytes);
        int recv_bytes = 0;
        const int peer_rank = peer_ranks[ipeer];
        exchange_bytes_(&send_bytes, sizeof(send_bytes), &recv_bytes, sizeof(recv_bytes),
                        peer_rank, peer_rank, tag_ghost_force_size);
        if (recv_bytes < 0 || recv_bytes % static_cast<int>(sizeof(ForceRecord)) != 0)
        {
            throw std::runtime_error("Invalid ghost force message size.");
        }

        std::vector<ForceRecord> recv_records(
            static_cast<std::size_t>(recv_bytes / static_cast<int>(sizeof(ForceRecord))));
        exchange_bytes_(send_records.data(), send_bytes, recv_records.data(), recv_bytes,
                        peer_rank, peer_rank, tag_ghost_forces);

        for (std::size_t irecord = 0; irecord < recv_records.size(); ++irecord)
        {
            const ForceRecord& record = recv_records[irecord];
            const std::map<std::pair<int, std::int64_t>, std::size_t>::const_iterator found
                = owned_lookup.find(std::make_pair(record.type, record.type_index));
            if (found == owned_lookup.end())
            {
                throw std::runtime_error("Cannot match a received ghost force to an owned atom.");
            }
            LocalAtom& atom = owned_atoms[found->second];
            atom.force.x += record.force[0];
            atom.force.y += record.force[1];
            atom.force.z += record.force[2];
        }
    }
}

void DomainDecomposition::prepare_neighbors(MDCell& cell)
{
    if (cell.cutoff_ <= 0.0)
    {
        throw std::runtime_error("MDCell neighbors must be initialized before use.");
    }

    synchronize_geometry_(cell);
    ModuleBase::CommunicationDomain domain;
#ifdef __MPI
    domain.initialize(cell.communicator());
#endif
    bool rebuild = !cell.neighbor_layout_valid_ || cell.neighbor_reference_frac_.size() != cell.owned_atoms_.size();
    rebuild = domain.max(rebuild || !ghost_layout_valid_ ? 1 : 0) != 0;
    double local_max_displacement = 0.0;
    if (!rebuild)
    {
        for (std::size_t i = 0; i < cell.owned_atoms_.size(); ++i)
        {
            ModuleBase::Vector3<double> delta = cell.owned_atoms_[i].frac - cell.neighbor_reference_frac_[i];
            delta.x -= std::nearbyint(delta.x);
            delta.y -= std::nearbyint(delta.y);
            delta.z -= std::nearbyint(delta.z);
            local_max_displacement = std::max(local_max_displacement, (delta * cell.latvec_).norm() * cell.lat0_);
        }
        local_max_displacement = domain.max(local_max_displacement);
        rebuild = local_max_displacement >= cell.skin_ * 0.5;
    }

    if (rebuild)
    {
        migrate_owned_atoms(cell);
        cell.neighbor_search_.reset(new NeighborSearch);
        cell.neighbor_search_->init(cell, cell.cutoff_ + cell.skin_);
        cell.neighbor_search_->build_neighbors();
        cell.neighbor_search_->refresh_mdcell(cell, cell.cutoff_);
        cell.neighbor_reference_frac_.resize(cell.owned_atoms_.size());
        for (std::size_t i = 0; i < cell.owned_atoms_.size(); ++i)
        {
            cell.neighbor_reference_frac_[i] = cell.owned_atoms_[i].frac;
        }
        cell.neighbor_layout_valid_ = true;
        return;
    }

    update_ghost_atom_positions(cell);
    cell.neighbor_search_->refresh_mdcell(cell, cell.cutoff_);
}

const std::array<int, 3>& DomainDecomposition::dims() const
{
    return dims_;
}

const std::array<int, 3>& DomainDecomposition::coords() const
{
    return coords_;
}

int DomainDecomposition::owner_rank_from_frac(const ModuleBase::Vector3<double>& frac) const
{
    std::array<int, 3> owner_coords;
    const double values[3] = {
        wrap_fractional(frac.x),
        wrap_fractional(frac.y),
        wrap_fractional(frac.z)
    };
    for (int idim = 0; idim < 3; ++idim)
    {
        int index = static_cast<int>(std::floor(values[idim] * dims_[idim]));
        index = std::min(std::max(index, 0), dims_[idim] - 1);
        owner_coords[idim] = index;
    }
    return rank_from_coords(owner_coords);
}

void DomainDecomposition::update_geometry_(const ModuleBase::Matrix3& latvec,
                                            double lat0, double cutoff, double skin)
{
    if (lat0 <= 0.0 || std::abs(latvec.Det()) <= 0.0)
    {
        throw std::runtime_error("DomainDecomposition requires a positive lattice constant and nonzero volume.");
    }
    latvec_ = latvec;
    inv_latvec_ = latvec_.Inverse();
    lat0_ = lat0;
    cutoff_ = cutoff;
    skin_ = skin;
    ghost_slots_.clear();
    ghost_layout_valid_ = false;
    initialized_ = true;
    const ModuleBase::Vector3<double> a1(latvec_.e11, latvec_.e12, latvec_.e13);
    const ModuleBase::Vector3<double> a2(latvec_.e21, latvec_.e22, latvec_.e23);
    const ModuleBase::Vector3<double> a3(latvec_.e31, latvec_.e32, latvec_.e33);
    const ModuleBase::Vector3<double> a2xa3 = ModuleBase::cross(a2, a3);
    const ModuleBase::Vector3<double> a3xa1 = ModuleBase::cross(a3, a1);
    const ModuleBase::Vector3<double> a1xa2 = ModuleBase::cross(a1, a2);

    const double volume = std::abs(a1 * a2xa3);
    const double heights[3] = {
        volume / a2xa3.norm(),
        volume / a3xa1.norm(),
        volume / a1xa2.norm()
    };
    const double cutoff_lat0 = (cutoff_ + skin_) / lat0_;
    for (int idim = 0; idim < 3; ++idim)
    {
        margin_[idim] = cutoff_lat0 / heights[idim] + 1.0e-12;
    }
}

void DomainDecomposition::synchronize_geometry_(MDCell& cell)
{
    if (!initialized_)
    {
        throw std::runtime_error("DomainDecomposition must be initialized before use.");
    }
#ifdef __MPI
    if (comm_ != cell.communicator())
    {
        throw std::runtime_error("DomainDecomposition and MDCell must use the same communication domain.");
    }
#endif
    const ModuleBase::Matrix3& lattice = cell.latvec_;
    const bool changed = lattice.e11 != latvec_.e11
                         || lattice.e12 != latvec_.e12
                         || lattice.e13 != latvec_.e13
                         || lattice.e21 != latvec_.e21
                         || lattice.e22 != latvec_.e22
                         || lattice.e23 != latvec_.e23
                         || lattice.e31 != latvec_.e31
                         || lattice.e32 != latvec_.e32
                         || lattice.e33 != latvec_.e33
                         || cell.lat0_ != lat0_ || cell.cutoff_ != cutoff_ || cell.skin_ != skin_;
    if (changed)
    {
        update_geometry_(lattice, cell.lat0_, cell.cutoff_, cell.skin_);
        cell.neighbor_layout_valid_ = false;
    }
}

void DomainDecomposition::exchange_bytes_(const void* send, std::size_t send_bytes,
                                           void* recv, std::size_t recv_bytes,
                                           int send_rank, int recv_rank, int tag) const
{
    if (send_bytes > static_cast<std::size_t>(std::numeric_limits<int>::max())
        || recv_bytes > static_cast<std::size_t>(std::numeric_limits<int>::max()))
    {
        throw std::overflow_error("DomainDecomposition message exceeds MPI byte count range.");
    }
    if (send_rank == rank_ && recv_rank == rank_)
    {
        if (send_bytes != recv_bytes)
        {
            throw std::runtime_error("Local exchange requires matching buffer sizes.");
        }
        if (send_bytes != 0)
        {
            std::memcpy(recv, send, send_bytes);
        }
        return;
    }
#ifdef __MPI
    MPI_Sendrecv(send, static_cast<int>(send_bytes), MPI_BYTE, send_rank, tag,
                 recv, static_cast<int>(recv_bytes), MPI_BYTE, recv_rank, tag,
                 cart_comm_, MPI_STATUS_IGNORE);
#else
    static_cast<void>(tag);
    throw std::runtime_error("A serial communication domain cannot exchange with a remote rank.");
#endif
}

void DomainDecomposition::exchange_migration_(
    const std::array<int, 6>& neighbors,
    const std::array<std::vector<PackedAtom>, 6>& send_atoms,
    const std::array<int, 6>& send_counts,
    std::array<std::vector<PackedAtom>, 6>& recv_atoms) const
{
    const int direction_count = 6;
    if (size_ == 1)
    {
        for (int idir = 0; idir < direction_count; ++idir)
        {
            recv_atoms[idir] = send_atoms[idir ^ 1];
        }
        return;
    }
#ifdef __MPI
    std::array<int, direction_count> recv_counts;
    std::array<MPI_Request, 2 * direction_count> requests;
    for (int idir = 0; idir < direction_count; ++idir)
    {
        const int opposite = idir ^ 1;
        MPI_Irecv(&recv_counts[idir], 1, MPI_INT, neighbors[idir], tag_migration_count_base + opposite, comm_, &requests[idir]);
        MPI_Isend(&send_counts[idir], 1, MPI_INT, neighbors[idir], tag_migration_count_base + idir, comm_, &requests[direction_count + idir]);
    }
    MPI_Waitall(2 * direction_count, &requests[0], MPI_STATUSES_IGNORE);

    for (int idir = 0; idir < direction_count; ++idir)
    {
        if (recv_counts[idir] < 0 || recv_counts[idir] % static_cast<int>(sizeof(PackedAtom)) != 0)
        {
            throw std::runtime_error("Invalid DomainDecomposition migration receive count.");
        }
        recv_atoms[idir].resize(static_cast<std::size_t>(recv_counts[idir] / static_cast<int>(sizeof(PackedAtom))));
    }

    for (int idir = 0; idir < direction_count; ++idir)
    {
        const int opposite = idir ^ 1;
        MPI_Irecv(recv_atoms[idir].empty() ? NULL : reinterpret_cast<char*>(&recv_atoms[idir][0]),
                  recv_counts[idir], MPI_BYTE, neighbors[idir], tag_migration_atoms_base + opposite, comm_, &requests[idir]);
        MPI_Isend(send_atoms[idir].empty() ? NULL : reinterpret_cast<const char*>(&send_atoms[idir][0]),
                  send_counts[idir], MPI_BYTE, neighbors[idir], tag_migration_atoms_base + idir, comm_, &requests[direction_count + idir]);
    }
    MPI_Waitall(2 * direction_count, &requests[0], MPI_STATUSES_IGNORE);

#else
    static_cast<void>(neighbors);
    static_cast<void>(send_counts);
    throw std::runtime_error("A serial communication domain cannot migrate to a remote rank.");
#endif
}

double DomainDecomposition::wrap_fractional(double value)
{
    value -= std::floor(value);
    if (value >= 1.0 - 1.0e-12)
    {
        return 0.0;
    }
    if (value < 1.0e-12)
    {
        return 0.0;
    }
    return value;
}

int DomainDecomposition::floor_div(int value, int divisor)
{
    assert(divisor!=0);
    int quotient = value / divisor;
    const int remainder = value % divisor;
    if (remainder != 0 && ((remainder < 0) != (divisor < 0)))
    {
        --quotient;
    }
    return quotient;
}

int DomainDecomposition::positive_mod(int value, int divisor)
{
    int result = value % divisor;
    if (result < 0)
    {
        result += divisor;
    }
    return result;
}

ModuleBase::Vector3<double> DomainDecomposition::wrapped_frac_from_cart(
    const ModuleBase::Vector3<double>& cart) const
{
    const ModuleBase::Vector3<double> frac = cart * inv_latvec_;
    return ModuleBase::Vector3<double>(wrap_fractional(frac.x),
                                       wrap_fractional(frac.y),
                                       wrap_fractional(frac.z));
}

int DomainDecomposition::rank_from_coords(const std::array<int, 3>& coords) const
{
    // MPI_Cart_create uses reorder = 0 and the last coordinate varies fastest.
    // This also maps the single serial domain to rank zero.
    return (coords[0] * dims_[1] + coords[1]) * dims_[2] + coords[2];
}

std::vector<LocalAtom> DomainDecomposition::split_owned_atoms_from_ucell(const UnitCell& ucell) const
{
    std::vector<LocalAtom> owned_atoms;
    owned_atoms.clear();
    owned_atoms.reserve(static_cast<size_t>(ucell.nat / std::max(1, size_) + 1));

    for (int it = 0; it < ucell.ntype; ++it)
    {
        for (int ia = 0; ia < ucell.atoms[it].na; ++ia)
        {
            const ModuleBase::Vector3<double> original_cart = ucell.atoms[it].tau[ia];
            const ModuleBase::Vector3<double> frac = wrapped_frac_from_cart(original_cart);
            const int owner = owner_rank_from_frac(frac);
            if (owner == rank_)
            {
                const ModuleBase::Vector3<double> wrapped_cart = frac * latvec_;
                owned_atoms.push_back(LocalAtom(wrapped_cart,
                                                frac,
                                                ucell.atoms[it].vel[ia],
                                                ModuleBase::Vector3<double>(0.0, 0.0, 0.0),
                                                ucell.atoms[it].mbl[ia],
                                                ucell.atoms[it].mass / ModuleBase::AU_to_MASS,
                                                it,
                                                ia,
                                                owner));
            }
        }
    }
    return owned_atoms;
}

void DomainDecomposition::target_for_offset(const std::array<int, 3>& offset,
                                            std::array<int, 3>& target_coords,
                                            std::array<int, 3>& image_shift) const
{
    for (int idim = 0; idim < 3; ++idim)
    {
        const int unwrapped = coords_[idim] + offset[idim];
        const int period_shift = floor_div(unwrapped, dims_[idim]);
        target_coords[idim] = positive_mod(unwrapped, dims_[idim]);
        image_shift[idim] = -period_shift;
    }
}

int DomainDecomposition::neighbor_layer(int dim) const
{
    return std::max(1, static_cast<int>(std::ceil(margin_[dim] * dims_[dim])));
}

void DomainDecomposition::build_ghost_exchange_slots(std::vector<GhostExchangeSlot>& slots) const
{
    slots.clear();

    const int nlayer_x = neighbor_layer(0);
    const int nlayer_y = neighbor_layer(1);
    const int nlayer_z = neighbor_layer(2);

    slots.reserve(static_cast<std::size_t>((2 * nlayer_x + 1)
                                           * (2 * nlayer_y + 1)
                                           * (2 * nlayer_z + 1)
                                           - 1));
    for (int dx = -nlayer_x; dx <= nlayer_x; ++dx)
    {
        for (int dy = -nlayer_y; dy <= nlayer_y; ++dy)
        {
            for (int dz = -nlayer_z; dz <= nlayer_z; ++dz)
            {
                if (dx == 0 && dy == 0 && dz == 0)
                {
                    continue;
                }

                GhostExchangeSlot slot;
                slot.offset = {{dx, dy, dz}};
                const std::array<int, 3> recv_offset = {{-dx, -dy, -dz}};
                std::array<int, 3> recv_coords;
                target_for_offset(slot.offset, slot.target_coords, slot.image_shift);
                target_for_offset(recv_offset, recv_coords, slot.recv_image_shift);
                // target_for_offset gives our image in the receiving domain.
                // Incoming coordinates need the inverse shift back into ours.
                for (int dim = 0; dim < 3; ++dim)
                {
                    slot.recv_image_shift[dim] = -slot.recv_image_shift[dim];
                }
                slot.send_rank = rank_from_coords(slot.target_coords);
                slot.recv_rank = rank_from_coords(recv_coords);
                slots.push_back(slot);
            }
        }
    }
}

DomainDecomposition::PackedAtom DomainDecomposition::pack_atom(
    const LocalAtom& atom,
    const std::array<int, 3>& image_shift) const
{
    PackedAtom packed;
    packed.frac[0] = atom.frac.x;
    packed.frac[1] = atom.frac.y;
    packed.frac[2] = atom.frac.z;
    packed.vel[0] = atom.vel.x;
    packed.vel[1] = atom.vel.y;
    packed.vel[2] = atom.vel.z;
    packed.force[0] = atom.force.x;
    packed.force[1] = atom.force.y;
    packed.force[2] = atom.force.z;
    packed.mbl[0] = atom.mbl.x;
    packed.mbl[1] = atom.mbl.y;
    packed.mbl[2] = atom.mbl.z;
    packed.mass = atom.mass;
    packed.image_shift[0] = image_shift[0];
    packed.image_shift[1] = image_shift[1];
    packed.image_shift[2] = image_shift[2];
    packed.type = atom.type;
    packed.type_index = atom.type_index;
    packed.owner_rank = atom.owner_rank;
    return packed;
}

LocalAtom DomainDecomposition::unpack_ghost_atom(const PackedAtom& packed) const
{
    const ModuleBase::Vector3<double> frac(packed.frac[0], packed.frac[1], packed.frac[2]);
    const ModuleBase::Vector3<double> image_frac(packed.frac[0] + packed.image_shift[0],
                                                 packed.frac[1] + packed.image_shift[1],
                                                 packed.frac[2] + packed.image_shift[2]);
    const ModuleBase::Vector3<double> cart = image_frac * latvec_;
    const ModuleBase::Vector3<double> vel(packed.vel[0], packed.vel[1], packed.vel[2]);
    const ModuleBase::Vector3<double> force(packed.force[0], packed.force[1], packed.force[2]);
    const ModuleBase::Vector3<int> mbl(packed.mbl[0], packed.mbl[1], packed.mbl[2]);
    return LocalAtom(cart,
                     frac,
                     vel,
                     force,
                     mbl,
                     packed.mass,
                     packed.type,
                     packed.type_index,
                     packed.owner_rank);
}

LocalAtom DomainDecomposition::unpack_owned_atom(const PackedAtom& packed) const
{
    const ModuleBase::Vector3<double> frac(packed.frac[0], packed.frac[1], packed.frac[2]);
    const ModuleBase::Vector3<double> cart = frac * latvec_;
    const ModuleBase::Vector3<double> vel(packed.vel[0], packed.vel[1], packed.vel[2]);
    const ModuleBase::Vector3<double> force(packed.force[0], packed.force[1], packed.force[2]);
    const ModuleBase::Vector3<int> mbl(packed.mbl[0], packed.mbl[1], packed.mbl[2]);
    return LocalAtom(cart,
                     frac,
                     vel,
                     force,
                     mbl,
                     packed.mass,
                     packed.type,
                     packed.type_index,
                     packed.owner_rank);
}
