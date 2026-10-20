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
    int swap = object;
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Bcast(&swap, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (my_rank != 0)
        object = static_cast<bool>(swap);
#endif
}

void bcast_char(char* object, const int n)
{
#ifdef __MPI
    MPI_Bcast(object, n, MPI_CHAR, 0, MPI_COMM_WORLD);
#endif
}

#ifdef __MPI
void bcast_int(int* object, const int n, const MPI_Comm mpi_comm, const int root)
{
    MPI_Bcast(object, n, MPI_INT, root, mpi_comm);
}

void bcast_double(double* object, const int n, const MPI_Comm mpi_comm, const int root)
{
    MPI_Bcast(object, n, MPI_DOUBLE, root, mpi_comm);
}

int my_rank(const MPI_Comm mpi_comm)
{
    int rank = 0;
    MPI_Comm_rank(mpi_comm, &rank);
    return rank;
}

int nproc(const MPI_Comm mpi_comm)
{
    int size = 1;
    MPI_Comm_size(mpi_comm, &size);
    return size;
}

void barrier(const MPI_Comm mpi_comm)
{
    MPI_Barrier(mpi_comm);
}
#endif

} // namespace Parallel_Common
