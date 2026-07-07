#include "source_cell/module_neighlist/domain_decomposition.h"

#ifdef __MPI

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstring>
#include <limits>
#include <stdexcept>

DomainDecomposition::DomainDecomposition()
    : comm_(MPI_COMM_NULL),
      cart_comm_(MPI_COMM_NULL),
      owns_cart_comm_(false),
      rank_(0),
      size_(1),
      dims_(),
      coords_(),
      margin_(),
      latvec_(),
      inv_latvec_(),
      lat0_(1.0),
      cutoff_(0.0),
      skin_(0.0)
{
    dims_[0] = dims_[1] = dims_[2] = 1;
    coords_[0] = coords_[1] = coords_[2] = 0;
    margin_[0] = margin_[1] = margin_[2] = 0.0;
}

DomainDecomposition::~DomainDecomposition()
{
    if (owns_cart_comm_ && cart_comm_ != MPI_COMM_NULL)
    {
        MPI_Comm_free(&cart_comm_);
    }
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

double DomainDecomposition::dot_product(const ModuleBase::Vector3<double>& a,
                                        const ModuleBase::Vector3<double>& b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

ModuleBase::Vector3<double> DomainDecomposition::cross_product(const ModuleBase::Vector3<double>& a,
                                                               const ModuleBase::Vector3<double>& b)
{
    return ModuleBase::Vector3<double>(a.y * b.z - a.z * b.y,
                                       a.z * b.x - a.x * b.z,
                                       a.x * b.y - a.y * b.x);
}

double DomainDecomposition::norm(const ModuleBase::Vector3<double>& value)
{
    return std::sqrt(dot_product(value, value));
}

void DomainDecomposition::init(MPI_Comm comm,
                               const ModuleBase::Matrix3& latvec,
                               double lat0,
                               double cutoff,
                               double skin)
{
    comm_ = comm;
    MPI_Comm_rank(comm_, &rank_);
    MPI_Comm_size(comm_, &size_);

    latvec_ = latvec;
    inv_latvec_ = latvec_.Inverse();
    lat0_ = lat0;
    cutoff_ = cutoff;
    skin_ = skin;

    int dims[3] = {0, 0, 0};
    MPI_Dims_create(size_, 3, dims);
    dims_[0] = std::max(1, dims[0]);
    dims_[1] = std::max(1, dims[1]);
    dims_[2] = std::max(1, dims[2]);

    int periods[3] = {1, 1, 1};
    if (owns_cart_comm_ && cart_comm_ != MPI_COMM_NULL)
    {
        MPI_Comm_free(&cart_comm_);
        cart_comm_ = MPI_COMM_NULL;
        owns_cart_comm_ = false;
    }
    MPI_Cart_create(comm_, 3, dims, periods, 0, &cart_comm_);
    owns_cart_comm_ = cart_comm_ != MPI_COMM_NULL;
    MPI_Comm_rank(cart_comm_, &rank_);
    int coords[3] = {0, 0, 0};
    MPI_Cart_coords(cart_comm_, rank_, 3, coords);
    coords_[0] = coords[0];
    coords_[1] = coords[1];
    coords_[2] = coords[2];

    const ModuleBase::Vector3<double> a1(latvec_.e11, latvec_.e12, latvec_.e13);
    const ModuleBase::Vector3<double> a2(latvec_.e21, latvec_.e22, latvec_.e23);
    const ModuleBase::Vector3<double> a3(latvec_.e31, latvec_.e32, latvec_.e33);
    const ModuleBase::Vector3<double> a2xa3 = cross_product(a2, a3);
    const ModuleBase::Vector3<double> a3xa1 = cross_product(a3, a1);
    const ModuleBase::Vector3<double> a1xa2 = cross_product(a1, a2);

    const double volume = std::abs(dot_product(a1, a2xa3));
    const double heights[3] = {
        volume / norm(a2xa3),
        volume / norm(a3xa1),
        volume / norm(a1xa2)
    };
    const double cutoff_lat0 = (cutoff_ + skin_) / lat0_;
    for (int idim = 0; idim < 3; ++idim)
    {
        margin_[idim] = cutoff_lat0 / heights[idim] + 1.0e-12;
    }
}

const std::array<int, 3>& DomainDecomposition::dims() const
{
    return dims_;
}

const std::array<int, 3>& DomainDecomposition::coords() const
{
    return coords_;
}

int DomainDecomposition::rank() const
{
    return rank_;
}

int DomainDecomposition::size() const
{
    return size_;
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
    int raw_coords[3] = {coords[0], coords[1], coords[2]};
    int rank = 0;
    MPI_Cart_rank(cart_comm_, raw_coords, &rank);
    return rank;
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

void DomainDecomposition::split_owned_atoms_from_ucell(const AtomProvider& ucell,
                                                       std::vector<LocalAtom>& owned_atoms) const
{
    owned_atoms.clear();
    owned_atoms.reserve(static_cast<size_t>(ucell.get_natom() / std::max(1, size_) + 1));

    ModuleNeighList::GlobalAtomId global_id = 0;
    for (int it = 0; it < ucell.get_ntype(); ++it)
    {
        for (int ia = 0; ia < ucell.get_na(it); ++ia)
        {
            const ModuleBase::Vector3<double> original_cart = ucell.get_tau(it, ia);
            const ModuleBase::Vector3<double> frac = wrapped_frac_from_cart(original_cart);
            const int owner = owner_rank_from_frac(frac);
            if (owner == rank_)
            {
                const ModuleBase::Vector3<double> wrapped_cart = frac * latvec_;
                owned_atoms.push_back(LocalAtom(wrapped_cart, frac, it, ia, global_id, owner, false));
            }
            ++global_id;
        }
    }
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

bool DomainDecomposition::atom_overlaps_target_halo(
    const LocalAtom& atom,
    const std::array<int, 3>& target_coords,
    const std::array<int, 3>& image_shift) const
{
    const double frac_values[3] = {
        atom.frac.x + image_shift[0],
        atom.frac.y + image_shift[1],
        atom.frac.z + image_shift[2]
    };
    for (int idim = 0; idim < 3; ++idim)
    {
        const double lo = static_cast<double>(target_coords[idim]) / dims_[idim];
        const double hi = static_cast<double>(target_coords[idim] + 1) / dims_[idim];
        if (frac_values[idim] < lo - margin_[idim] ||
            frac_values[idim] >= hi + margin_[idim])
        {
            return false;
        }
    }
    return true;
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
    packed.image_shift[0] = image_shift[0];
    packed.image_shift[1] = image_shift[1];
    packed.image_shift[2] = image_shift[2];
    packed.type = atom.type;
    packed.type_index = atom.type_index;
    packed.global_id = atom.global_id;
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
    return LocalAtom(cart,
                     frac,
                     packed.type,
                     packed.type_index,
                     packed.global_id,
                     packed.owner_rank,
                     true);
}

void DomainDecomposition::exchange_ghost_atoms(const std::vector<LocalAtom>& owned_atoms,
                                               std::vector<LocalAtom>& ghost_atoms) const
{
    ghost_atoms.clear();

    std::vector<GhostExchangeSlot> slots;
    build_ghost_exchange_slots(slots);

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
                }
            }
        }
    }

    for (std::size_t islot = 0; islot < slots.size(); ++islot)
    {
        const GhostExchangeSlot& slot = slots[islot];
        const std::vector<PackedAtom>& send_atoms = send_buffers[islot];

        if (send_atoms.size() > static_cast<std::size_t>(std::numeric_limits<int>::max()))
        {
            throw std::overflow_error("DomainDecomposition ghost send count exceeds int range.");
        }

        if (slot.send_rank == rank_ && slot.recv_rank == rank_)
        {
            for (size_t i = 0; i < send_atoms.size(); ++i)
            {
                ghost_atoms.push_back(unpack_ghost_atom(send_atoms[i]));
            }
            continue;
        }

        int send_count = static_cast<int>(send_atoms.size());
        int recv_count = 0;
        MPI_Sendrecv(&send_count,
                     1,
                     MPI_INT,
                     slot.send_rank,
                     9100,
                     &recv_count,
                     1,
                     MPI_INT,
                     slot.recv_rank,
                     9100,
                     cart_comm_,
                     MPI_STATUS_IGNORE);

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

        MPI_Sendrecv(send_atoms.empty() ? NULL : &send_atoms[0],
                     send_bytes,
                     MPI_BYTE,
                     slot.send_rank,
                     9101,
                     recv_atoms.empty() ? NULL : &recv_atoms[0],
                     recv_bytes,
                     MPI_BYTE,
                     slot.recv_rank,
                     9101,
                     cart_comm_,
                     MPI_STATUS_IGNORE);

        for (size_t i = 0; i < recv_atoms.size(); ++i)
        {
            ghost_atoms.push_back(unpack_ghost_atom(recv_atoms[i]));
        }
    }
}

#endif // __MPI
