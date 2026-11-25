#include "parallel_common.h"

#include "source_base/parallel_reduce.h"

#ifdef __MPI
#include <mpi.h>
#endif

namespace Parallel_Common
{

#ifdef __MPI
/// Broadcast a trivially-copyable buffer of type T on MPI_COMM_WORLD from
/// rank 0. This is the single implementation behind all bcast_* wrappers.
template <typename T>
static void bcast_world_impl(T* object, const int n)
{
    MPI_Bcast(object, n, Parallel_Reduce::MPI_Type<T>::value, 0, MPI_COMM_WORLD);
}
#endif

void bcast_string(std::string& object) // Peize Lin fix bug 2019-03-18
{
#ifdef __MPI
    int size = object.size();
    MPI_Bcast(&size, 1, MPI_INT, 0, MPI_COMM_WORLD);

    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    if (0 != my_rank)
    {
        object.resize(size);
    }

    MPI_Bcast(&object[0], size, MPI_CHAR, 0, MPI_COMM_WORLD);
#endif
    return;
}

void bcast_string(std::string* object, const int n) // Peize Lin fix bug 2019-03-18
{
#ifdef __MPI
    for (int i = 0; i < n; i++)
        bcast_string(object[i]);
#endif
    return;
}

void bcast_complex_double(std::complex<double>& object)
{
#ifdef __MPI
    bcast_world_impl(&object, 1);
#endif
}

void bcast_complex_double(std::complex<double>* object, const int n)
{
#ifdef __MPI
    bcast_world_impl(object, n);
#endif
}

void bcast_double(double& object)
{
#ifdef __MPI
    bcast_world_impl(&object, 1);
#endif
}

void bcast_double(double* object, const int n)
{
#ifdef __MPI
    bcast_world_impl(object, n);
#endif
}

void bcast_int(int& object)
{
#ifdef __MPI
    bcast_world_impl(&object, 1);
#endif
}

void bcast_int(int* object, const int n)
{
#ifdef __MPI
    bcast_world_impl(object, n);
#endif
}

void bcast_bool(bool& object)
{
#ifdef __MPI
    bcast_bool(object, MPI_COMM_WORLD);
#endif
}

void bcast_bool(bool& object, const MPI_Comm& comm)
{
#ifdef __MPI
    int swap = object;
    int my_rank;
    MPI_Comm_rank(comm, &my_rank);
    MPI_Bcast(&swap, 1, MPI_INT, 0, comm);
    if (my_rank != 0)
    {
        object = static_cast<bool>(swap);
    }
#else
    // One rank, which is the root: the value is already correct.
    (void)object;
    (void)comm;
#endif
}

void bcast_char(char* object, const int n)
{
#ifdef __MPI
    MPI_Bcast(object, n, MPI_CHAR, 0, MPI_COMM_WORLD);
#endif
}

void bcast_int(int* object, const int n, const MPI_Comm mpi_comm, const int root)
{
#ifdef __MPI
    MPI_Bcast(object, n, MPI_INT, root, mpi_comm);
#else
    // A serial run has one rank, which is always the root, so the payload is
    // already in place everywhere. The other arguments are unused.
    (void)object;
    (void)n;
    (void)mpi_comm;
    (void)root;
#endif
}

void bcast_double(double* object, const int n, const MPI_Comm mpi_comm, const int root)
{
#ifdef __MPI
    MPI_Bcast(object, n, MPI_DOUBLE, root, mpi_comm);
#else
    (void)object;
    (void)n;
    (void)mpi_comm;
    (void)root;
#endif
}

int my_rank(const MPI_Comm mpi_comm)
{
#ifdef __MPI
    int rank = 0;
    MPI_Comm_rank(mpi_comm, &rank);
    return rank;
#else
    (void)mpi_comm;
    return 0;
#endif
}

int nproc(const MPI_Comm mpi_comm)
{
#ifdef __MPI
    int size = 1;
    MPI_Comm_size(mpi_comm, &size);
    return size;
#else
    (void)mpi_comm;
    return 1;
#endif
}

void barrier(const MPI_Comm mpi_comm)
{
#ifdef __MPI
    MPI_Barrier(mpi_comm);
#else
    (void)mpi_comm;
#endif
}

#ifdef __MPI
// bcast single on an explicit communicator, from the communicator's rank 0
void bcast_mpi_offset(MPI_Offset& object, const MPI_Comm& comm)
{
    MPI_Bcast(&object, 1, MPI_OFFSET, 0, comm);
}
#endif

} // namespace Parallel_Common
