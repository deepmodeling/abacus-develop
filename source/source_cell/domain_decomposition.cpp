#include "source_cell/domain_decomposition.h"
#include "source_base/parallel_cell.h"
#include "source_cell/ghost_exchange_plan.h"
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
    bool local_periodic_crossing = false;
    if (!rebuild)
    {
        for (std::size_t i = 0; i < cell.owned_atoms_.size(); ++i)
        {
            ModuleBase::Vector3<double> delta = cell.owned_atoms_[i].frac - cell.neighbor_reference_frac_[i];
            if (std::nearbyint(delta.x) != 0.0 || std::nearbyint(delta.y) != 0.0
                || std::nearbyint(delta.z) != 0.0)
            {
                local_periodic_crossing = true;
            }
            delta.x -= std::nearbyint(delta.x);
            delta.y -= std::nearbyint(delta.y);
            delta.z -= std::nearbyint(delta.z);
            local_max_displacement = std::max(local_max_displacement, (delta * cell.latvec_).norm() * cell.lat0_);
        }
        local_max_displacement = domain.max(local_max_displacement);
        rebuild = local_max_displacement >= cell.skin_ * 0.5
                  || domain.max(local_periodic_crossing ? 1 : 0) != 0;
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
