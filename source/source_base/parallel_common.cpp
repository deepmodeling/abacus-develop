#include "parallel_common.h"

#ifdef __MPI
#include <mpi.h>
#endif

#include <cstring>

void Parallel_Common::bcast_string(std::string& object) // Peize Lin fix bug 2019-03-18
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

void Parallel_Common::bcast_string(std::string* object, const int n) // Peize Lin fix bug 2019-03-18
{
#ifdef __MPI
    for (int i = 0; i < n; i++)
        bcast_string(object[i]);
#endif
    return;
}

void Parallel_Common::bcast_complex_double(std::complex<double>& object)
{
#ifdef __MPI
    MPI_Bcast(&object, 1, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
#endif
}

void Parallel_Common::bcast_complex_double(std::complex<double>* object, const int n)
{
#ifdef __MPI
    MPI_Bcast(object, n, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
#endif
}

void Parallel_Common::bcast_double(double& object)
{
#ifdef __MPI
    MPI_Bcast(&object, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif
}

void Parallel_Common::bcast_double(double* object, const int n)
{
#ifdef __MPI
    MPI_Bcast(object, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif
}

void Parallel_Common::bcast_int(int& object)
{
#ifdef __MPI
    MPI_Bcast(&object, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif
}

void Parallel_Common::bcast_int(int* object, const int n)
{
#ifdef __MPI
    MPI_Bcast(object, n, MPI_INT, 0, MPI_COMM_WORLD);
#endif
}

void Parallel_Common::bcast_bool(bool& object)
{
#ifdef __MPI
    int swap = object;
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    if (my_rank == 0)
        swap = object;
    MPI_Bcast(&swap, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (my_rank != 0)
        object = static_cast<bool>(swap);
#endif
}

void Parallel_Common::bcast_char(char* object, const int n)
{
#ifdef __MPI
    MPI_Bcast(object, n, MPI_CHAR, 0, MPI_COMM_WORLD);
#endif
}
