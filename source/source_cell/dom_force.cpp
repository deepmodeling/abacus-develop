#include "source_cell/dom_decomp.h"
#include "source_cell/mdcell.h"

#include <algorithm>
#include <limits>
#include <map>
#include <stdexcept>
#include <utility>

namespace
{
const int tag_ghost_force_size = 9200;
const int tag_ghost_forces = 9201;
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
