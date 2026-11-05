#include "source_cell/dom_decomp.h"
#include "source_cell/mdcell.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

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
    for (int axis = 0; axis < 3; ++axis)
    {
        packed.displacement[axis] = atom.displacement[axis];
        packed.boundary_shift[axis] = atom.boundary_shift[axis];
    }
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
    LocalAtom atom(cart,
                   frac,
                   vel,
                   force,
                   mbl,
                   packed.mass,
                   packed.type,
                   packed.type_index,
                   packed.owner_rank);
    for (int axis = 0; axis < 3; ++axis)
    {
        atom.displacement[axis] = packed.displacement[axis];
        atom.boundary_shift[axis] = packed.boundary_shift[axis];
    }
    return atom;
}

LocalAtom DomainDecomposition::unpack_owned_atom(const PackedAtom& packed) const
{
    const ModuleBase::Vector3<double> frac(packed.frac[0], packed.frac[1], packed.frac[2]);
    const ModuleBase::Vector3<double> cart = frac * latvec_;
    const ModuleBase::Vector3<double> vel(packed.vel[0], packed.vel[1], packed.vel[2]);
    const ModuleBase::Vector3<double> force(packed.force[0], packed.force[1], packed.force[2]);
    const ModuleBase::Vector3<int> mbl(packed.mbl[0], packed.mbl[1], packed.mbl[2]);
    LocalAtom atom(cart,
                   frac,
                   vel,
                   force,
                   mbl,
                   packed.mass,
                   packed.type,
                   packed.type_index,
                   packed.owner_rank);
    for (int axis = 0; axis < 3; ++axis)
    {
        atom.displacement[axis] = packed.displacement[axis];
        atom.boundary_shift[axis] = packed.boundary_shift[axis];
    }
    return atom;
}
