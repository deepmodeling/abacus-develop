#include "parallel_common.h"

#ifdef __MPI
#include <mpi.h>
#endif

#include <cstring>

#ifdef __MPI
void Parallel_Common::bcast_string(std::string& object, MPI_Comm comm) // Peize Lin fix bug 2019-03-18
{
    int size = object.size();
    MPI_Bcast(&size, 1, MPI_INT, 0, comm);

    int my_rank;
    MPI_Comm_rank(comm, &my_rank);

    if (0 != my_rank)
    {
        object.resize(size);
    }

    MPI_Bcast(&object[0], size, MPI_CHAR, 0, comm);
    return;
}

void Parallel_Common::bcast_string(std::string* object, const int n, MPI_Comm comm) // Peize Lin fix bug 2019-03-18
{
    for (int i = 0; i < n; i++)
        bcast_string(object[i], comm);
    return;
}

void Parallel_Common::bcast_complex_double(std::complex<double>& object, MPI_Comm comm)
{
    MPI_Bcast(&object, 1, MPI_DOUBLE_COMPLEX, 0, comm);
}

void Parallel_Common::bcast_complex_double(std::complex<double>* object, const int n, MPI_Comm comm)
{
    MPI_Bcast(object, n, MPI_DOUBLE_COMPLEX, 0, comm);
}

void Parallel_Common::bcast_double(double& object, MPI_Comm comm)
{
    MPI_Bcast(&object, 1, MPI_DOUBLE, 0, comm);
}

void Parallel_Common::bcast_double(double* object, const int n, MPI_Comm comm)
{
    MPI_Bcast(object, n, MPI_DOUBLE, 0, comm);
}

void Parallel_Common::bcast_int(int& object, MPI_Comm comm)
{
    MPI_Bcast(&object, 1, MPI_INT, 0, comm);
}

void Parallel_Common::bcast_int(int* object, const int n, MPI_Comm comm)
{
    MPI_Bcast(object, n, MPI_INT, 0, comm);
}

void Parallel_Common::bcast_bool(bool& object, MPI_Comm comm)
{
    int swap = object;
    int my_rank;
    MPI_Comm_rank(comm, &my_rank);
    if (my_rank == 0)
        swap = object;
    MPI_Bcast(&swap, 1, MPI_INT, 0, comm);
    if (my_rank != 0)
        object = static_cast<bool>(swap);
}

void Parallel_Common::bcast_char(char* object, const int n, MPI_Comm comm)
{
    MPI_Bcast(object, n, MPI_CHAR, 0, comm);
}

#endif
