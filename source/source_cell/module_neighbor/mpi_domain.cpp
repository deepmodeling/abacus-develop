#include "mpi_domain.h"

#include <algorithm>
#include <cmath>
#include <map>
#include <set>
#include <stdexcept>
#include <tuple>
#include <utility>

namespace ModuleNeighbor
{
namespace
{
const int kPackedAtomFields = 9;
const double kSingularTolerance = 1.0e-12;

bool in_half_open_range(const double value,
                        const double lower,
                        const double upper,
                        const bool include_upper)
{
    if (include_upper)
    {
        return value >= lower && value <= upper;
    }
    return value >= lower && value < upper;
}

double dot(const std::array<double, 3>& lhs, const std::array<double, 3>& rhs)
{
    return lhs[0] * rhs[0] + lhs[1] * rhs[1] + lhs[2] * rhs[2];
}

std::array<double, 3> cross(const std::array<double, 3>& lhs, const std::array<double, 3>& rhs)
{
    return std::array<double, 3>{{lhs[1] * rhs[2] - lhs[2] * rhs[1],
                                  lhs[2] * rhs[0] - lhs[0] * rhs[2],
                                  lhs[0] * rhs[1] - lhs[1] * rhs[0]}};
}

double norm(const std::array<double, 3>& vector)
{
    return std::sqrt(dot(vector, vector));
}

std::tuple<int, int, int, int, int, int> mpi_record_key(const MpiAtomRecord& atom)
{
    return std::tuple<int, int, int, int, int, int>(atom.global_index,
                                                    atom.type,
                                                    atom.natom,
                                                    atom.pbc_shift[0],
                                                    atom.pbc_shift[1],
                                                    atom.pbc_shift[2]);
}

void append_packed_atom(const MpiAtomRecord& atom, std::vector<double>& buffer)
{
    buffer.push_back(atom.x);
    buffer.push_back(atom.y);
    buffer.push_back(atom.z);
    buffer.push_back(static_cast<double>(atom.global_index));
    buffer.push_back(static_cast<double>(atom.type));
    buffer.push_back(static_cast<double>(atom.natom));
    buffer.push_back(static_cast<double>(atom.pbc_shift[0]));
    buffer.push_back(static_cast<double>(atom.pbc_shift[1]));
    buffer.push_back(static_cast<double>(atom.pbc_shift[2]));
}

MpiAtomRecord unpack_atom(const std::vector<double>& buffer, const int offset, const bool is_ghost)
{
    return MpiAtomRecord(buffer[offset + 0],
                         buffer[offset + 1],
                         buffer[offset + 2],
                         static_cast<int>(buffer[offset + 3]),
                         std::array<int, 3>{{static_cast<int>(buffer[offset + 6]),
                                             static_cast<int>(buffer[offset + 7]),
                                             static_cast<int>(buffer[offset + 8])}},
                         is_ghost,
                         static_cast<int>(buffer[offset + 4]),
                         static_cast<int>(buffer[offset + 5]));
}

NeighborMpiComm null_comm()
{
#ifdef __MPI
    return MPI_COMM_NULL;
#else
    return kSerialMpiCommNull;
#endif
}

NeighborMpiComm world_comm()
{
#ifdef __MPI
    return MPI_COMM_WORLD;
#else
    return kSerialMpiCommWorld;
#endif
}
} // namespace

MpiAtomRecord::MpiAtomRecord()
    : x(0.0), y(0.0), z(0.0), global_index(-1), type(-1), natom(-1), pbc_shift{{0, 0, 0}}, is_ghost(false)
{
}

MpiAtomRecord::MpiAtomRecord(const double x_in,
                             const double y_in,
                             const double z_in,
                             const int global_index_in,
                             const std::array<int, 3>& pbc_shift_in,
                             const bool is_ghost_in,
                             const int type_in,
                             const int natom_in)
    : x(x_in),
      y(y_in),
      z(z_in),
      global_index(global_index_in),
      type(type_in),
      natom(natom_in),
      pbc_shift(pbc_shift_in),
      is_ghost(is_ghost_in)
{
}

MpiDomain::MpiDomain()
    : cart_comm_(null_comm()),
      owns_comm_(false),
      initialized_(false),
      rank_(0),
      size_(1),
      dims_{{1, 1, 1}},
      coords_{{0, 0, 0}},
      periods_{{0, 0, 0}},
      origin_{{0.0, 0.0, 0.0}},
      lattice_vectors_{{std::array<double, 3>{{1.0, 0.0, 0.0}},
                        std::array<double, 3>{{0.0, 1.0, 0.0}},
                        std::array<double, 3>{{0.0, 0.0, 1.0}}}},
      reciprocal_rows_{{std::array<double, 3>{{1.0, 0.0, 0.0}},
                        std::array<double, 3>{{0.0, 1.0, 0.0}},
                        std::array<double, 3>{{0.0, 0.0, 1.0}}}},
      fractional_ghost_padding_{{0.0, 0.0, 0.0}},
      local_bounds_(),
      ghost_cutoff_(0.0)
{
    local_bounds_.lower = std::array<double, 3>{{0.0, 0.0, 0.0}};
    local_bounds_.upper = std::array<double, 3>{{0.0, 0.0, 0.0}};
}

MpiDomain::MpiDomain(MpiDomain&& other) : MpiDomain()
{
    *this = std::move(other);
}

MpiDomain& MpiDomain::operator=(MpiDomain&& other)
{
    if (this != &other)
    {
        reset();
        cart_comm_ = other.cart_comm_;
        owns_comm_ = other.owns_comm_;
        initialized_ = other.initialized_;
        rank_ = other.rank_;
        size_ = other.size_;
        dims_ = other.dims_;
        coords_ = other.coords_;
        periods_ = other.periods_;
        origin_ = other.origin_;
        lattice_vectors_ = other.lattice_vectors_;
        reciprocal_rows_ = other.reciprocal_rows_;
        fractional_ghost_padding_ = other.fractional_ghost_padding_;
        local_bounds_ = other.local_bounds_;
        ghost_cutoff_ = other.ghost_cutoff_;

        other.cart_comm_ = null_comm();
        other.owns_comm_ = false;
        other.initialized_ = false;
    }
    return *this;
}

MpiDomain::~MpiDomain()
{
    reset();
}

void MpiDomain::reset()
{
#ifdef __MPI
    if (owns_comm_ && cart_comm_ != MPI_COMM_NULL)
    {
        int mpi_finalized = 0;
        MPI_Finalized(&mpi_finalized);
        if (!mpi_finalized)
        {
            MPI_Comm_free(&cart_comm_);
        }
    }
#endif
    cart_comm_ = null_comm();
    owns_comm_ = false;
    initialized_ = false;
}

void MpiDomain::initialize(NeighborMpiComm parent_comm,
                           const std::array<double, 3>& global_lower,
                           const std::array<double, 3>& global_upper,
                           const double ghost_cutoff,
                           const bool pbc)
{
    std::array<std::array<double, 3>, 3> lattice_vectors;
    for (int axis = 0; axis < 3; ++axis)
    {
        if (!(global_upper[axis] > global_lower[axis]))
        {
            throw std::invalid_argument("MpiDomain global bounds must have positive length.");
        }
        lattice_vectors[axis] = std::array<double, 3>{{0.0, 0.0, 0.0}};
        lattice_vectors[axis][axis] = global_upper[axis] - global_lower[axis];
    }
    initialize_lattice(parent_comm, global_lower, lattice_vectors, ghost_cutoff, pbc);
}

void MpiDomain::initialize_lattice(NeighborMpiComm parent_comm,
                                   const std::array<double, 3>& origin,
                                   const std::array<std::array<double, 3>, 3>& lattice_vectors,
                                   const double ghost_cutoff,
                                   const bool pbc)
{
    reset();
    if (ghost_cutoff < 0.0)
    {
        throw std::invalid_argument("MpiDomain ghost cutoff must be non-negative.");
    }

    const std::array<double, 3> b_cross_c = cross(lattice_vectors[1], lattice_vectors[2]);
    const double determinant = dot(lattice_vectors[0], b_cross_c);
    if (std::abs(determinant) <= kSingularTolerance)
    {
        throw std::invalid_argument("MpiDomain lattice vectors must form a non-singular cell.");
    }

    origin_ = origin;
    lattice_vectors_ = lattice_vectors;
    reciprocal_rows_[0] = b_cross_c;
    reciprocal_rows_[1] = cross(lattice_vectors[2], lattice_vectors[0]);
    reciprocal_rows_[2] = cross(lattice_vectors[0], lattice_vectors[1]);
    for (int axis = 0; axis < 3; ++axis)
    {
        for (int component = 0; component < 3; ++component)
        {
            reciprocal_rows_[axis][component] /= determinant;
        }
        fractional_ghost_padding_[axis] = ghost_cutoff * norm(reciprocal_rows_[axis]);
    }

    ghost_cutoff_ = ghost_cutoff;
    periods_ = pbc ? std::array<int, 3>{{1, 1, 1}} : std::array<int, 3>{{0, 0, 0}};

#ifdef __MPI
    MPI_Comm_size(parent_comm, &size_);
    MPI_Comm_rank(parent_comm, &rank_);
    dims_ = std::array<int, 3>{{0, 0, 0}};
    MPI_Dims_create(size_, 3, dims_.data());
    MPI_Cart_create(parent_comm, 3, dims_.data(), periods_.data(), 0, &cart_comm_);
    if (cart_comm_ == MPI_COMM_NULL)
    {
        throw std::runtime_error("MPI_Cart_create returned MPI_COMM_NULL.");
    }
    owns_comm_ = true;
    MPI_Comm_rank(cart_comm_, &rank_);
    MPI_Cart_coords(cart_comm_, rank_, 3, coords_.data());
#else
    (void)parent_comm;
    size_ = 1;
    rank_ = 0;
    dims_ = std::array<int, 3>{{1, 1, 1}};
    coords_ = std::array<int, 3>{{0, 0, 0}};
    cart_comm_ = world_comm();
    owns_comm_ = false;
#endif

    compute_local_bounds();
    initialized_ = true;
}

void MpiDomain::compute_local_bounds()
{
    for (int axis = 0; axis < 3; ++axis)
    {
        const double width = 1.0 / static_cast<double>(dims_[axis]);
        local_bounds_.lower[axis] = width * static_cast<double>(coords_[axis]);
        local_bounds_.upper[axis] = coords_[axis] == dims_[axis] - 1
                                      ? 1.0
                                      : width * static_cast<double>(coords_[axis] + 1);
    }
}

bool MpiDomain::initialized() const
{
    return initialized_;
}

int MpiDomain::rank() const
{
    return rank_;
}

int MpiDomain::size() const
{
    return size_;
}

NeighborMpiComm MpiDomain::cart_comm() const
{
    return cart_comm_;
}

const std::array<int, 3>& MpiDomain::dims() const
{
    return dims_;
}

const std::array<int, 3>& MpiDomain::coords() const
{
    return coords_;
}

const std::array<int, 3>& MpiDomain::periods() const
{
    return periods_;
}

const DomainBounds& MpiDomain::local_bounds() const
{
    return local_bounds_;
}

const std::array<double, 3>& MpiDomain::fractional_ghost_padding() const
{
    return fractional_ghost_padding_;
}

double MpiDomain::ghost_cutoff() const
{
    return ghost_cutoff_;
}

std::array<double, 3> MpiDomain::cartesian_to_fractional(const double x,
                                                        const double y,
                                                        const double z) const
{
    const std::array<double, 3> displacement{{x - origin_[0], y - origin_[1], z - origin_[2]}};
    return std::array<double, 3>{{dot(reciprocal_rows_[0], displacement),
                                  dot(reciprocal_rows_[1], displacement),
                                  dot(reciprocal_rows_[2], displacement)}};
}

std::array<double, 3> MpiDomain::fractional_to_cartesian(const double fx,
                                                        const double fy,
                                                        const double fz) const
{
    return std::array<double, 3>{{origin_[0] + fx * lattice_vectors_[0][0] + fy * lattice_vectors_[1][0]
                                      + fz * lattice_vectors_[2][0],
                                  origin_[1] + fx * lattice_vectors_[0][1] + fy * lattice_vectors_[1][1]
                                      + fz * lattice_vectors_[2][1],
                                  origin_[2] + fx * lattice_vectors_[0][2] + fy * lattice_vectors_[1][2]
                                      + fz * lattice_vectors_[2][2]}};
}

std::array<double, 3> MpiDomain::normalize_fractional(const std::array<double, 3>& fractional) const
{
    std::array<double, 3> normalized = fractional;
    for (int axis = 0; axis < 3; ++axis)
    {
        if (periods_[axis])
        {
            normalized[axis] -= std::floor(normalized[axis]);
        }
    }
    return normalized;
}

bool MpiDomain::inside_local(const std::array<double, 3>& fractional) const
{
    for (int axis = 0; axis < 3; ++axis)
    {
        const bool include_upper = coords_[axis] == dims_[axis] - 1;
        if (!in_half_open_range(fractional[axis],
                                local_bounds_.lower[axis],
                                local_bounds_.upper[axis],
                                include_upper))
        {
            return false;
        }
    }
    return true;
}

DomainBounds MpiDomain::bounds_for_coords(const std::array<int, 3>& coords) const
{
    DomainBounds bounds;
    for (int axis = 0; axis < 3; ++axis)
    {
        const double width = 1.0 / static_cast<double>(dims_[axis]);
        bounds.lower[axis] = width * static_cast<double>(coords[axis]);
        bounds.upper[axis] = coords[axis] == dims_[axis] - 1
                                 ? 1.0
                                 : width * static_cast<double>(coords[axis] + 1);
    }
    return bounds;
}

bool MpiDomain::inside_local_for_bounds(const DomainBounds& bounds,
                                        const std::array<int, 3>& coords,
                                        const std::array<double, 3>& fractional) const
{
    for (int axis = 0; axis < 3; ++axis)
    {
        const bool include_upper = coords[axis] == dims_[axis] - 1;
        if (!in_half_open_range(fractional[axis], bounds.lower[axis], bounds.upper[axis], include_upper))
        {
            return false;
        }
    }
    return true;
}

bool MpiDomain::inside_expanded_for_bounds(const DomainBounds& bounds,
                                           const std::array<double, 3>& fractional) const
{
    for (int axis = 0; axis < 3; ++axis)
    {
        if (fractional[axis] < bounds.lower[axis] - fractional_ghost_padding_[axis]
            || fractional[axis] > bounds.upper[axis] + fractional_ghost_padding_[axis])
        {
            return false;
        }
    }
    return true;
}

bool MpiDomain::neighbor_from_direction(const std::array<int, 3>& direction,
                                        int& neighbor_rank,
                                        std::array<int, 3>& target_coords,
                                        std::array<int, 3>& image_shift) const
{
    target_coords = coords_;
    image_shift = std::array<int, 3>{{0, 0, 0}};
    for (int axis = 0; axis < 3; ++axis)
    {
        target_coords[axis] += direction[axis];
        if (target_coords[axis] < 0)
        {
            if (!periods_[axis])
            {
                return false;
            }
            target_coords[axis] += dims_[axis];
            image_shift[axis] = 1;
        }
        else if (target_coords[axis] >= dims_[axis])
        {
            if (!periods_[axis])
            {
                return false;
            }
            target_coords[axis] -= dims_[axis];
            image_shift[axis] = -1;
        }
    }

#ifdef __MPI
    MPI_Cart_rank(cart_comm_, target_coords.data(), &neighbor_rank);
#else
    neighbor_rank = 0;
#endif
    return neighbor_rank != rank_;
}

std::vector<int> MpiDomain::neighbor_ranks() const
{
    std::vector<int> ranks;
    const int dx_begin = dims_[0] == 1 ? 0 : -1;
    const int dx_end = dims_[0] == 1 ? 0 : 1;
    const int dy_begin = dims_[1] == 1 ? 0 : -1;
    const int dy_end = dims_[1] == 1 ? 0 : 1;
    const int dz_begin = dims_[2] == 1 ? 0 : -1;
    const int dz_end = dims_[2] == 1 ? 0 : 1;

    for (int dx = dx_begin; dx <= dx_end; ++dx)
    {
        for (int dy = dy_begin; dy <= dy_end; ++dy)
        {
            for (int dz = dz_begin; dz <= dz_end; ++dz)
            {
                if (dx == 0 && dy == 0 && dz == 0)
                {
                    continue;
                }
                int neighbor_rank = -1;
                std::array<int, 3> target_coords{{0, 0, 0}};
                std::array<int, 3> image_shift{{0, 0, 0}};
                if (neighbor_from_direction(std::array<int, 3>{{dx, dy, dz}},
                                            neighbor_rank,
                                            target_coords,
                                            image_shift))
                {
                    ranks.push_back(neighbor_rank);
                }
            }
        }
    }
    std::sort(ranks.begin(), ranks.end());
    ranks.erase(std::unique(ranks.begin(), ranks.end()), ranks.end());
    return ranks;
}

std::map<int, std::vector<MpiAtomRecord>>
MpiDomain::build_neighbor_send_records(const std::vector<MpiAtomRecord>& local_atoms) const
{
    std::map<int, std::vector<MpiAtomRecord>> send_records;
    std::map<int, std::set<std::tuple<int, int, int, int, int, int>>> seen;

    const int dx_begin = dims_[0] == 1 ? 0 : -1;
    const int dx_end = dims_[0] == 1 ? 0 : 1;
    const int dy_begin = dims_[1] == 1 ? 0 : -1;
    const int dy_end = dims_[1] == 1 ? 0 : 1;
    const int dz_begin = dims_[2] == 1 ? 0 : -1;
    const int dz_end = dims_[2] == 1 ? 0 : 1;

    for (int dx = dx_begin; dx <= dx_end; ++dx)
    {
        for (int dy = dy_begin; dy <= dy_end; ++dy)
        {
            for (int dz = dz_begin; dz <= dz_end; ++dz)
            {
                if (dx == 0 && dy == 0 && dz == 0)
                {
                    continue;
                }

                int neighbor_rank = -1;
                std::array<int, 3> target_coords{{0, 0, 0}};
                std::array<int, 3> image_shift{{0, 0, 0}};
                if (!neighbor_from_direction(std::array<int, 3>{{dx, dy, dz}},
                                             neighbor_rank,
                                             target_coords,
                                             image_shift))
                {
                    continue;
                }

                const DomainBounds target_bounds = bounds_for_coords(target_coords);
                for (const MpiAtomRecord& atom: local_atoms)
                {
                    // Build the image seen by the target rank, then keep only
                    // atoms that fall inside its fractional ghost shell.
                    MpiAtomRecord image = atom;
                    std::array<double, 3> fractional = cartesian_to_fractional(atom.x, atom.y, atom.z);
                    for (int axis = 0; axis < 3; ++axis)
                    {
                        fractional[axis] += static_cast<double>(image_shift[axis]);
                        image.x += static_cast<double>(image_shift[axis]) * lattice_vectors_[axis][0];
                        image.y += static_cast<double>(image_shift[axis]) * lattice_vectors_[axis][1];
                        image.z += static_cast<double>(image_shift[axis]) * lattice_vectors_[axis][2];
                        image.pbc_shift[axis] += image_shift[axis];
                    }
                    image.is_ghost = true;

                    if (!inside_expanded_for_bounds(target_bounds, fractional)
                        || inside_local_for_bounds(target_bounds, target_coords, fractional))
                    {
                        continue;
                    }

                    const auto key = mpi_record_key(image);
                    if (seen[neighbor_rank].insert(key).second)
                    {
                        send_records[neighbor_rank].push_back(image);
                    }
                }
            }
        }
    }

    return send_records;
}

bool MpiDomain::owns(const double x, const double y, const double z) const
{
    return inside_local(normalize_fractional(cartesian_to_fractional(x, y, z)));
}

std::vector<int> MpiDomain::select_local_atoms(const std::vector<MpiAtomRecord>& atoms) const
{
    std::vector<int> local_indices;
    for (int i = 0; i < static_cast<int>(atoms.size()); ++i)
    {
        if (owns(atoms[i].x, atoms[i].y, atoms[i].z))
        {
            local_indices.push_back(i);
        }
    }
    return local_indices;
}

std::vector<MpiAtomRecord> MpiDomain::exchange_ghost_atoms(const std::vector<MpiAtomRecord>& local_atoms) const
{
    return exchange_ghost_atoms(local_atoms, nullptr);
}

std::vector<MpiAtomRecord> MpiDomain::exchange_ghost_atoms(const std::vector<MpiAtomRecord>& local_atoms,
                                                          MpiGhostExchangeStats* stats) const
{
    if (!initialized_)
    {
        throw std::runtime_error("MpiDomain must be initialized before ghost exchange.");
    }
    if (stats != nullptr)
    {
        *stats = MpiGhostExchangeStats();
    }

#ifdef __MPI
    const std::map<int, std::vector<MpiAtomRecord>> send_records = build_neighbor_send_records(local_atoms);

    // Exchange with every adjacent rank, even when the payload is empty. This
    // keeps the nonblocking send/recv pairs symmetric across all ranks.
    const std::vector<int> neighbors = neighbor_ranks();
    if (stats != nullptr)
    {
        stats->neighbor_rank_count = static_cast<int>(neighbors.size());
    }

    std::vector<std::vector<double>> send_buffers(size_);
    for (const auto& entry: send_records)
    {
        std::vector<double>& buffer = send_buffers[entry.first];
        buffer.reserve(entry.second.size() * kPackedAtomFields);
        for (const MpiAtomRecord& atom: entry.second)
        {
            append_packed_atom(atom, buffer);
        }
        if (stats != nullptr && !entry.second.empty())
        {
            ++stats->nonempty_send_rank_count;
            stats->sent_atom_count += static_cast<int>(entry.second.size());
        }
    }

    std::vector<int> send_counts(size_, 0);
    std::vector<int> recv_counts(size_, 0);
    std::vector<MPI_Request> count_requests;
    count_requests.reserve(neighbors.size() * 2);
    for (const int neighbor_rank: neighbors)
    {
        send_counts[neighbor_rank] = static_cast<int>(send_buffers[neighbor_rank].size());
        if (stats != nullptr)
        {
            stats->sent_payload_count += send_counts[neighbor_rank];
        }
        MPI_Request recv_request;
        MPI_Request send_request;
        MPI_Irecv(&recv_counts[neighbor_rank], 1, MPI_INT, neighbor_rank, 4100, cart_comm_, &recv_request);
        MPI_Isend(&send_counts[neighbor_rank], 1, MPI_INT, neighbor_rank, 4100, cart_comm_, &send_request);
        count_requests.push_back(recv_request);
        count_requests.push_back(send_request);
    }
    if (!count_requests.empty())
    {
        MPI_Waitall(static_cast<int>(count_requests.size()), count_requests.data(), MPI_STATUSES_IGNORE);
    }

    std::vector<std::vector<double>> recv_buffers(size_);
    std::vector<MPI_Request> payload_requests;
    payload_requests.reserve(neighbors.size() * 2);
    for (const int neighbor_rank: neighbors)
    {
        recv_buffers[neighbor_rank].resize(recv_counts[neighbor_rank], 0.0);
        MPI_Request recv_request;
        MPI_Request send_request;
        MPI_Irecv(recv_buffers[neighbor_rank].empty() ? nullptr : recv_buffers[neighbor_rank].data(),
                  recv_counts[neighbor_rank],
                  MPI_DOUBLE,
                  neighbor_rank,
                  4101,
                  cart_comm_,
                  &recv_request);
        MPI_Isend(send_buffers[neighbor_rank].empty() ? nullptr : send_buffers[neighbor_rank].data(),
                  send_counts[neighbor_rank],
                  MPI_DOUBLE,
                  neighbor_rank,
                  4101,
                  cart_comm_,
                  &send_request);
        payload_requests.push_back(recv_request);
        payload_requests.push_back(send_request);
    }
    if (!payload_requests.empty())
    {
        MPI_Waitall(static_cast<int>(payload_requests.size()), payload_requests.data(), MPI_STATUSES_IGNORE);
    }

    std::vector<MpiAtomRecord> ghosts;
    for (const int neighbor_rank: neighbors)
    {
        const std::vector<double>& buffer = recv_buffers[neighbor_rank];
        if (stats != nullptr)
        {
            stats->received_payload_count += static_cast<int>(buffer.size());
        }
        for (int offset = 0; offset + kPackedAtomFields - 1 < static_cast<int>(buffer.size());
             offset += kPackedAtomFields)
        {
            ghosts.push_back(unpack_atom(buffer, offset, true));
        }
    }
    if (stats != nullptr)
    {
        stats->received_ghost_count = static_cast<int>(ghosts.size());
    }
    return ghosts;
#else
    (void)local_atoms;
    (void)stats;
    return std::vector<MpiAtomRecord>();
#endif
}

} // namespace ModuleNeighbor
