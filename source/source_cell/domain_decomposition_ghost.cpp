#include "source_cell/domain_decomposition.h"
#include "source_cell/ghost_exchange_plan.h"
#include "source_cell/mdcell.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <stdexcept>
#include <utility>

namespace
{
const int tag_ghost_count = 9100;
const int tag_ghost_atoms = 9101;
const int tag_ghost_positions = 9110;
}

GhostExchangePlanView GhostExchangePlan::make_view()
{
    GhostExchangePlanView view;
    view.nswap = static_cast<int>(sendnum.size());
    view.sendnum = sendnum.empty() ? NULL : sendnum.data();
    view.recvnum = recvnum.empty() ? NULL : recvnum.data();
    view.firstrecv = firstrecv.empty() ? NULL : firstrecv.data();
    view.sendproc = sendproc.empty() ? NULL : sendproc.data();
    view.recvproc = recvproc.empty() ? NULL : recvproc.data();
    view.mapping = mapping.empty() ? NULL : mapping.data();
    view.sendlist_storage.resize(send_atom_indices.size());
    for (std::size_t iswap = 0; iswap < send_atom_indices.size(); ++iswap)
    {
        if (!send_atom_indices[iswap].empty())
        {
            view.sendlist_storage[iswap] = send_atom_indices[iswap].data();
        }
        else
        {
            view.sendlist_storage[iswap] = NULL;
        }
    }
    return view;
}

void DomainDecomposition::exchange_ghost_atoms(MDCell& cell)
{
    synchronize_geometry_(cell);
    const std::vector<LocalAtom>& owned_atoms = cell.owned_atoms_;
    std::vector<LocalAtom>& ghost_atoms = cell.ghost_atoms_;
    ghost_atoms.clear();
    ghost_layout_valid_ = false;
    cell.clear_ghost_exchange_plan_();
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
        slot.previous_fracs.clear();
        slot.dynamic_image_shifts.clear();
        slot.previous_fracs.reserve(slot.send_atom_indices.size());
        slot.dynamic_image_shifts.reserve(slot.send_atom_indices.size());
        for (const int index : slot.send_atom_indices)
        {
            const LocalAtom& atom = owned_atoms[static_cast<std::size_t>(index)];
            slot.previous_fracs.push_back({{atom.frac.x, atom.frac.y, atom.frac.z}});
            slot.dynamic_image_shifts.push_back(slot.image_shift);
        }
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
    cell.set_ghost_exchange_plan_(build_ghost_exchange_plan(cell));
    cell.clear_forces_(ghost_atoms);
    cell.neighbor_layout_valid_ = false;
}

GhostExchangePlan DomainDecomposition::build_ghost_exchange_plan(const MDCell& cell) const
{
    if (!initialized_ || !ghost_layout_valid_)
    {
        throw std::runtime_error("Ghost exchange plan is not available.");
    }

    GhostExchangePlan plan;
    plan.nprocs = size_;
    plan.sendnum.reserve(ghost_slots_.size());
    plan.recvnum.reserve(ghost_slots_.size());
    plan.firstrecv.reserve(ghost_slots_.size());
    plan.send_atom_indices.reserve(ghost_slots_.size());
    plan.sendproc.reserve(ghost_slots_.size());
    plan.recvproc.reserve(ghost_slots_.size());

    for (const GhostExchangeSlot& slot : ghost_slots_)
    {
        if (slot.send_atom_indices.size() > static_cast<std::size_t>(std::numeric_limits<int>::max()))
        {
            throw std::overflow_error("Ghost send list exceeds int range.");
        }
        if (slot.ghost_begin > static_cast<std::size_t>(std::numeric_limits<int>::max()))
        {
            throw std::overflow_error("Ghost receive offset exceeds int range.");
        }
        if (slot.ghost_count < 0)
        {
            throw std::runtime_error("Invalid ghost receive count.");
        }

        plan.sendnum.push_back(static_cast<int>(slot.send_atom_indices.size()));
        plan.recvnum.push_back(slot.ghost_count);
        const std::size_t firstrecv = cell.owned_atoms_.size() + slot.ghost_begin;
        if (firstrecv > static_cast<std::size_t>(std::numeric_limits<int>::max()))
        {
            throw std::overflow_error("Ghost receive index exceeds int range.");
        }
        plan.firstrecv.push_back(static_cast<int>(firstrecv));
        plan.send_atom_indices.push_back(slot.send_atom_indices);
        plan.sendproc.push_back(slot.send_rank);
        plan.recvproc.push_back(slot.recv_rank);
    }

    if (plan.nprocs == 1)
    {
        const std::size_t natom = cell.owned_atoms_.size() + cell.ghost_atoms_.size();
        plan.mapping.resize(natom, 0);
        std::map<std::pair<int, std::int64_t>, int> owned_lookup;
        for (std::size_t iat = 0; iat < cell.owned_atoms_.size(); ++iat)
        {
            const LocalAtom& atom = cell.owned_atoms_[iat];
            owned_lookup[std::make_pair(atom.type, atom.type_index)] = static_cast<int>(iat);
            plan.mapping[iat] = static_cast<int>(iat);
        }
        for (std::size_t iat = 0; iat < cell.ghost_atoms_.size(); ++iat)
        {
            const LocalAtom& atom = cell.ghost_atoms_[iat];
            const std::map<std::pair<int, std::int64_t>, int>::const_iterator owner
                = owned_lookup.find(std::make_pair(atom.type, atom.type_index));
            if (owner == owned_lookup.end())
            {
                throw std::runtime_error("Cannot match a single-rank ghost atom to an owned atom.");
            }
            plan.mapping[cell.owned_atoms_.size() + iat] = owner->second;
        }
    }

    return plan;
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
        GhostExchangeSlot& slot = ghost_slots_[islot];
        // Keep wrapped fractional coordinates for metadata and send the
        // dynamically selected image coordinates in the same message.
        std::vector<double> send_frac(6 * slot.send_atom_indices.size(), 0.0);
        for (std::size_t i = 0; i < slot.send_atom_indices.size(); ++i)
        {
            const LocalAtom& atom = owned_atoms[static_cast<std::size_t>(slot.send_atom_indices[i])];
            const std::array<int, 3> image_shift = image_shift_for_update(
                atom, slot.previous_fracs[i], slot.dynamic_image_shifts[i], slot);
            slot.previous_fracs[i] = {{atom.frac.x, atom.frac.y, atom.frac.z}};
            slot.dynamic_image_shifts[i] = image_shift;
            send_frac[3 * i] = atom.frac.x;
            send_frac[3 * i + 1] = atom.frac.y;
            send_frac[3 * i + 2] = atom.frac.z;
            send_frac[3 * slot.send_atom_indices.size() + 3 * i] = atom.frac.x + image_shift[0];
            send_frac[3 * slot.send_atom_indices.size() + 3 * i + 1] = atom.frac.y + image_shift[1];
            send_frac[3 * slot.send_atom_indices.size() + 3 * i + 2] = atom.frac.z + image_shift[2];
        }
        std::vector<double> recv_frac(6 * static_cast<std::size_t>(slot.ghost_count), 0.0);
        exchange_bytes_(send_frac.data(), send_frac.size() * sizeof(double),
                        recv_frac.data(), recv_frac.size() * sizeof(double),
                        slot.send_rank, slot.recv_rank, tag_ghost_positions);
        for (int i = 0; i < slot.ghost_count; ++i)
        {
            LocalAtom& ghost = ghost_atoms[slot.ghost_begin + static_cast<std::size_t>(i)];
            ghost.frac.set(recv_frac[3 * i], recv_frac[3 * i + 1], recv_frac[3 * i + 2]);
            const ModuleBase::Vector3<double> image_frac(
                recv_frac[3 * static_cast<std::size_t>(slot.ghost_count) + 3 * i],
                recv_frac[3 * static_cast<std::size_t>(slot.ghost_count) + 3 * i + 1],
                recv_frac[3 * static_cast<std::size_t>(slot.ghost_count) + 3 * i + 2]);
            ghost.cart = image_frac * latvec_;
            ghost.force.set(0.0, 0.0, 0.0);
        }
    }
    cell.clear_forces_(ghost_atoms);
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

std::array<int, 3> DomainDecomposition::image_shift_for_update(
    const LocalAtom& atom,
    const std::array<double, 3>& previous_frac,
    const std::array<int, 3>& previous_image_shift,
    const GhostExchangeSlot& slot) const
{
    std::array<int, 3> image_shift = previous_image_shift;
    const double frac[3] = {atom.frac.x, atom.frac.y, atom.frac.z};
    for (int idim = 0; idim < 3; ++idim)
    {
        if (slot.offset[idim] == 0)
        {
            continue;
        }

        const int crossing = static_cast<int>(std::nearbyint(frac[idim] - previous_frac[idim]));
        if (crossing == 0)
        {
            continue;
        }

        // A positive wrapped displacement changes the image index in the
        // opposite direction for a positive offset, and vice versa.  This
        // preserves the same physical image continuously across a boundary.
        image_shift[idim] -= slot.offset[idim] > 0 ? crossing : -crossing;
    }
    return image_shift;
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
