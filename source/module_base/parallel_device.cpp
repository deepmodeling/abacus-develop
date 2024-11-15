#include "parallel_device.h"
#ifdef __MPI
namespace Parallel_Common
{
void bcast_data(std::complex<double>* object, const int& n, const MPI_Comm& comm)
{
    MPI_Bcast(object, n * 2, MPI_DOUBLE, 0, comm);
}
void bcast_data(std::complex<float>* object, const int& n, const MPI_Comm& comm)
{
    MPI_Bcast(object, n * 2, MPI_FLOAT, 0, comm);
}
void bcast_data(double* object, const int& n, const MPI_Comm& comm)
{
    MPI_Bcast(object, n, MPI_DOUBLE, 0, comm);
}
void bcast_data(float* object, const int& n, const MPI_Comm& comm)
{
    MPI_Bcast(object, n, MPI_FLOAT, 0, comm);
}
void reduce_data(std::complex<double>* object, const int& n, const MPI_Comm& comm)
{
    MPI_Allreduce(MPI_IN_PLACE, object, n * 2, MPI_DOUBLE, MPI_SUM, comm);
}
void reduce_data(std::complex<float>* object, const int& n, const MPI_Comm& comm)
{
    MPI_Allreduce(MPI_IN_PLACE, object, n * 2, MPI_FLOAT, MPI_SUM, comm);
}
void reduce_data(double* object, const int& n, const MPI_Comm& comm)
{
    MPI_Allreduce(MPI_IN_PLACE, object, n, MPI_DOUBLE, MPI_SUM, comm);
}
void reduce_data(float* object, const int& n, const MPI_Comm& comm)
{
    MPI_Allreduce(MPI_IN_PLACE, object, n, MPI_FLOAT, MPI_SUM, comm);
}
}
#endif