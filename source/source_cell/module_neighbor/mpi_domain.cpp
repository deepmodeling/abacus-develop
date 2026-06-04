#include "mpi_domain.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace ModuleNeighbor
{
namespace
{
const int kPackedAtomFields = 7;

std::array<double, 3> make_point(const double x, const double y, const double z)
{
    return std::array<double, 3>{{x, y, z}};
}

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
} // namespace

MpiAtomRecord::MpiAtomRecord()
    : x(0.0), y(0.0), z(0.0), global_index(-1), pbc_shift{{0, 0, 0}}, is_ghost(false)
{
}

MpiAtomRecord::MpiAtomRecord(const double x_in,
                             const double y_in,
                             const double z_in,
                             const int global_index_in,
                             const std::array<int, 3>& pbc_shift_in,
                             const bool is_ghost_in)
    : x(x_in), y(y_in), z(z_in), global_index(global_index_in), pbc_shift(pbc_shift_in), is_ghost(is_ghost_in)
{
}

MpiDomain::MpiDomain()
    : cart_comm_(MPI_COMM_NULL),
      owns_comm_(false),
      initialized_(false),
      rank_(0),
      size_(1),
      dims_{{1, 1, 1}},
      coords_{{0, 0, 0}},
      periods_{{0, 0, 0}},
      global_lower_{{0.0, 0.0, 0.0}},
      global_upper_{{0.0, 0.0, 0.0}},
      global_length_{{0.0, 0.0, 0.0}},
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
        global_lower_ = other.global_lower_;
        global_upper_ = other.global_upper_;
        global_length_ = other.global_length_;
        local_bounds_ = other.local_bounds_;
        ghost_cutoff_ = other.ghost_cutoff_;

        other.cart_comm_ = MPI_COMM_NULL;
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
    cart_comm_ = MPI_COMM_NULL;
    owns_comm_ = false;
    initialized_ = false;
}

void MpiDomain::initialize(MPI_Comm parent_comm,
                           const std::array<double, 3>& global_lower,
                           const std::array<double, 3>& global_upper,
                           const double ghost_cutoff,
                           const bool pbc)
{
    reset();
    if (ghost_cutoff < 0.0)
    {
        throw std::invalid_argument("MpiDomain ghost cutoff must be non-negative.");
    }
    for (int axis = 0; axis < 3; ++axis)
    {
        if (!(global_upper[axis] > global_lower[axis]))
        {
            throw std::invalid_argument("MpiDomain global bounds must have positive length.");
        }
    }

    global_lower_ = global_lower;
    global_upper_ = global_upper;
    ghost_cutoff_ = ghost_cutoff;
    periods_ = pbc ? std::array<int, 3>{{1, 1, 1}} : std::array<int, 3>{{0, 0, 0}};
    for (int axis = 0; axis < 3; ++axis)
    {
        global_length_[axis] = global_upper_[axis] - global_lower_[axis];
    }

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
    cart_comm_ = MPI_COMM_WORLD;
    owns_comm_ = false;
#endif

    compute_local_bounds();
    initialized_ = true;
}

void MpiDomain::compute_local_bounds()
{
    for (int axis = 0; axis < 3; ++axis)
    {
        const double width = global_length_[axis] / static_cast<double>(dims_[axis]);
        local_bounds_.lower[axis] = global_lower_[axis] + width * static_cast<double>(coords_[axis]);
        local_bounds_.upper[axis] = coords_[axis] == dims_[axis] - 1
                                      ? global_upper_[axis]
                                      : global_lower_[axis] + width * static_cast<double>(coords_[axis] + 1);
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

MPI_Comm MpiDomain::cart_comm() const
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

double MpiDomain::ghost_cutoff() const
{
    return ghost_cutoff_;
}

bool MpiDomain::inside_local(const std::array<double, 3>& point) const
{
    for (int axis = 0; axis < 3; ++axis)
    {
        const bool include_upper = coords_[axis] == dims_[axis] - 1;
        if (!in_half_open_range(point[axis], local_bounds_.lower[axis], local_bounds_.upper[axis], include_upper))
        {
            return false;
        }
    }
    return true;
}

bool MpiDomain::inside_expanded(const std::array<double, 3>& point) const
{
    for (int axis = 0; axis < 3; ++axis)
    {
        if (point[axis] < local_bounds_.lower[axis] - ghost_cutoff_
            || point[axis] > local_bounds_.upper[axis] + ghost_cutoff_)
        {
            return false;
        }
    }
    return true;
}

bool MpiDomain::owns(const double x, const double y, const double z) const
{
    return inside_local(make_point(x, y, z));
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

bool MpiDomain::append_ghost_images(const MpiAtomRecord& atom, std::vector<MpiAtomRecord>& ghosts) const
{
    bool appended = false;
    const int sx_begin = periods_[0] ? -1 : 0;
    const int sx_end = periods_[0] ? 1 : 0;
    const int sy_begin = periods_[1] ? -1 : 0;
    const int sy_end = periods_[1] ? 1 : 0;
    const int sz_begin = periods_[2] ? -1 : 0;
    const int sz_end = periods_[2] ? 1 : 0;

    for (int sx = sx_begin; sx <= sx_end; ++sx)
    {
        for (int sy = sy_begin; sy <= sy_end; ++sy)
        {
            for (int sz = sz_begin; sz <= sz_end; ++sz)
            {
                MpiAtomRecord image = atom;
                image.x += static_cast<double>(sx) * global_length_[0];
                image.y += static_cast<double>(sy) * global_length_[1];
                image.z += static_cast<double>(sz) * global_length_[2];
                image.pbc_shift[0] += sx;
                image.pbc_shift[1] += sy;
                image.pbc_shift[2] += sz;
                image.is_ghost = true;

                const std::array<double, 3> point = make_point(image.x, image.y, image.z);
                if (inside_expanded(point) && !inside_local(point))
                {
                    ghosts.push_back(image);
                    appended = true;
                }
            }
        }
    }
    return appended;
}

std::vector<MpiAtomRecord> MpiDomain::exchange_ghost_atoms(const std::vector<MpiAtomRecord>& local_atoms) const
{
    if (!initialized_)
    {
        throw std::runtime_error("MpiDomain must be initialized before ghost exchange.");
    }

#ifdef __MPI
    const int send_count = static_cast<int>(local_atoms.size()) * kPackedAtomFields;
    std::vector<double> send_buffer(send_count, 0.0);
    for (int i = 0; i < static_cast<int>(local_atoms.size()); ++i)
    {
        const int offset = i * kPackedAtomFields;
        send_buffer[offset + 0] = local_atoms[i].x;
        send_buffer[offset + 1] = local_atoms[i].y;
        send_buffer[offset + 2] = local_atoms[i].z;
        send_buffer[offset + 3] = static_cast<double>(local_atoms[i].global_index);
        send_buffer[offset + 4] = static_cast<double>(local_atoms[i].pbc_shift[0]);
        send_buffer[offset + 5] = static_cast<double>(local_atoms[i].pbc_shift[1]);
        send_buffer[offset + 6] = static_cast<double>(local_atoms[i].pbc_shift[2]);
    }

    std::vector<int> counts(size_, 0);
    MPI_Allgather(&send_count, 1, MPI_INT, counts.data(), 1, MPI_INT, cart_comm_);

    std::vector<int> displacements(size_, 0);
    int total_count = 0;
    for (int i = 0; i < size_; ++i)
    {
        displacements[i] = total_count;
        total_count += counts[i];
    }

    std::vector<double> recv_buffer(total_count, 0.0);
    MPI_Allgatherv(send_buffer.empty() ? 0 : send_buffer.data(),
                   send_count,
                   MPI_DOUBLE,
                   recv_buffer.empty() ? 0 : recv_buffer.data(),
                   counts.data(),
                   displacements.data(),
                   MPI_DOUBLE,
                   cart_comm_);

    std::vector<MpiAtomRecord> ghosts;
    for (int owner_rank = 0; owner_rank < size_; ++owner_rank)
    {
        if (owner_rank == rank_)
        {
            continue;
        }
        for (int offset = displacements[owner_rank]; offset < displacements[owner_rank] + counts[owner_rank];
             offset += kPackedAtomFields)
        {
            MpiAtomRecord atom(recv_buffer[offset + 0],
                               recv_buffer[offset + 1],
                               recv_buffer[offset + 2],
                               static_cast<int>(recv_buffer[offset + 3]),
                               std::array<int, 3>{{static_cast<int>(recv_buffer[offset + 4]),
                                                   static_cast<int>(recv_buffer[offset + 5]),
                                                   static_cast<int>(recv_buffer[offset + 6])}},
                               false);
            append_ghost_images(atom, ghosts);
        }
    }
    return ghosts;
#else
    (void)local_atoms;
    return std::vector<MpiAtomRecord>();
#endif
}

} // namespace ModuleNeighbor
