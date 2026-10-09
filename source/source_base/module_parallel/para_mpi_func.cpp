#include "para_mpi_func.h"

#include <algorithm>
#include <cstring>

namespace Parallel
{

#ifdef __MPI
namespace
{
inline MPI_Datatype mpi_type(int*)
{
    return MPI_INT;
}
inline MPI_Datatype mpi_type(double*)
{
    return MPI_DOUBLE;
}
inline MPI_Datatype mpi_type(std::complex<double>*)
{
    return MPI_DOUBLE;
} // 2 doubles
} // namespace
#endif

// ========== Broadcast ==========

void bcast_bool(bool& v, const ParaWorld& world, int root)
{
#ifdef __MPI
    if (!world.valid())
    {
        return;
    }
    int tmp = v ? 1 : 0;
    MPI_Bcast(&tmp, 1, MPI_INT, root, world.comm());
    v = (tmp != 0);
#endif
}

void bcast_int(int& v, const ParaWorld& world, int root)
{
#ifdef __MPI
    if (!world.valid())
    {
        return;
    }
    MPI_Bcast(&v, 1, MPI_INT, root, world.comm());
#endif
}

void bcast_double(double& v, const ParaWorld& world, int root)
{
#ifdef __MPI
    if (!world.valid())
    {
        return;
    }
    MPI_Bcast(&v, 1, MPI_DOUBLE, root, world.comm());
#endif
}

void bcast_complex(std::complex<double>& v, const ParaWorld& world, int root)
{
#ifdef __MPI
    if (!world.valid())
    {
        return;
    }
    MPI_Bcast(&v, 2, MPI_DOUBLE, root, world.comm());
#endif
}

void bcast_string(std::string& s, const ParaWorld& world, int root)
{
#ifdef __MPI
    if (!world.valid())
    {
        return;
    }
    int len = static_cast<int>(s.size());
    MPI_Bcast(&len, 1, MPI_INT, root, world.comm());
    if (world.rank() != root)
    {
        s.resize(len);
    }
    if (len > 0)
    {
        MPI_Bcast(&s[0], len, MPI_CHAR, root, world.comm());
    }
#endif
}

void bcast_int(int* v, int n, const ParaWorld& world, int root)
{
#ifdef __MPI
    if (!world.valid())
    {
        return;
    }
    MPI_Bcast(v, n, MPI_INT, root, world.comm());
#endif
}

void bcast_double(double* v, int n, const ParaWorld& world, int root)
{
#ifdef __MPI
    if (!world.valid())
    {
        return;
    }
    MPI_Bcast(v, n, MPI_DOUBLE, root, world.comm());
#endif
}

void bcast_complex(std::complex<double>* v, int n, const ParaWorld& world, int root)
{
#ifdef __MPI
    if (!world.valid())
    {
        return;
    }
    MPI_Bcast(v, 2 * n, MPI_DOUBLE, root, world.comm());
#endif
}

void bcast_char(char* v, int n, const ParaWorld& world, int root)
{
#ifdef __MPI
    if (!world.valid())
    {
        return;
    }
    MPI_Bcast(v, n, MPI_CHAR, root, world.comm());
#endif
}

void bcast_string(std::string* v, int n, const ParaWorld& world, int root)
{
#ifdef __MPI
    if (!world.valid())
    {
        return;
    }
    for (int i = 0; i < n; ++i)
    {
        bcast_string(v[i], world, root);
    }
#endif
}

// ========== Reduce ==========

void reduce_all(double& v, const ParaWorld& world)
{
#ifdef __MPI
    if (!world.valid())
    {
        return;
    }
    MPI_Allreduce(MPI_IN_PLACE, &v, 1, MPI_DOUBLE, MPI_SUM, world.comm());
#endif
}

void reduce_all(int& v, const ParaWorld& world)
{
#ifdef __MPI
    if (!world.valid())
    {
        return;
    }
    MPI_Allreduce(MPI_IN_PLACE, &v, 1, MPI_INT, MPI_SUM, world.comm());
#endif
}

void reduce_all(double* v, int n, const ParaWorld& world)
{
#ifdef __MPI
    if (!world.valid())
    {
        return;
    }
    MPI_Allreduce(MPI_IN_PLACE, v, n, MPI_DOUBLE, MPI_SUM, world.comm());
#endif
}

// ========== Min/Max ==========

void reduce_min(double& v, const ParaWorld& world)
{
#ifdef __MPI
    if (!world.valid())
    {
        return;
    }
    MPI_Allreduce(MPI_IN_PLACE, &v, 1, MPI_DOUBLE, MPI_MIN, world.comm());
#endif
}

void reduce_max(double& v, const ParaWorld& world)
{
#ifdef __MPI
    if (!world.valid())
    {
        return;
    }
    MPI_Allreduce(MPI_IN_PLACE, &v, 1, MPI_DOUBLE, MPI_MAX, world.comm());
#endif
}

void reduce_min(int& v, const ParaWorld& world)
{
#ifdef __MPI
    if (!world.valid())
    {
        return;
    }
    MPI_Allreduce(MPI_IN_PLACE, &v, 1, MPI_INT, MPI_MIN, world.comm());
#endif
}

void reduce_max(int& v, const ParaWorld& world)
{
#ifdef __MPI
    if (!world.valid())
    {
        return;
    }
    MPI_Allreduce(MPI_IN_PLACE, &v, 1, MPI_INT, MPI_MAX, world.comm());
#endif
}

// ========== Barrier ==========

void barrier(const ParaWorld& world)
{
#ifdef __MPI
    if (!world.valid())
    {
        return;
    }
    MPI_Barrier(world.comm());
#endif
}

// ========== Gather ==========

void allgather_int(int& v, int* all, const ParaWorld& world)
{
#ifdef __MPI
    if (!world.valid())
    {
        return;
    }
    MPI_Allgather(&v, 1, MPI_INT, all, 1, MPI_INT, world.comm());
#else
    all[0] = v;
#endif
}

void allgatherv_double(const double* send,
                       const int send_count,
                       double* recv,
                       const int* recv_counts,
                       const int* displs,
                       const ParaWorld& world)
{
    if (!world.valid())
    {
        return;
    }
#ifdef __MPI
    MPI_Allgatherv(send, send_count, MPI_DOUBLE, recv, recv_counts, displs, MPI_DOUBLE, world.comm());
#else
    static_cast<void>(recv_counts);
    if (send_count > 0)
    {
        std::copy(send, send + send_count, recv + displs[0]);
    }
#endif
}

} // namespace Parallel
