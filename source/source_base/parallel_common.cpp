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
void bcast_mpi_offset(MPI_Offset& object, const MPI_Comm& comm)
{
    MPI_Bcast(&object, 1, MPI_OFFSET, 0, comm);
}

void barrier(const MPI_Comm comm)
{
    MPI_Barrier(comm);
}

void gatherv_bytes(const void* send_buffer,
                   const int send_count,
                   void* receive_buffer,
                   const int* receive_counts,
                   const int* receive_displacements,
                   const int root,
                   const MPI_Comm comm)
{
    MPI_Gatherv(send_buffer,
                send_count,
                MPI_BYTE,
                receive_buffer,
                receive_counts,
                receive_displacements,
                MPI_BYTE,
                root,
                comm);
}

void allgatherv_bytes(const void* send_buffer,
                      const int send_count,
                      void* receive_buffer,
                      const int* receive_counts,
                      const int* receive_displacements,
                      const MPI_Comm comm)
{
    MPI_Allgatherv(send_buffer,
                   send_count,
                   MPI_BYTE,
                   receive_buffer,
                   receive_counts,
                   receive_displacements,
                   MPI_BYTE,
                   comm);
}

void scatterv_bytes(const void* send_buffer,
                    const int* send_counts,
                    const int* send_displacements,
                    void* receive_buffer,
                    const int receive_count,
                    const int root,
                    const MPI_Comm comm)
{
    MPI_Scatterv(send_buffer,
                 send_counts,
                 send_displacements,
                 MPI_BYTE,
                 receive_buffer,
                 receive_count,
                 MPI_BYTE,
                 root,
                 comm);
}

void allreduce_sum(int& value, const MPI_Comm comm)
{
    MPI_Allreduce(MPI_IN_PLACE, &value, 1, MPI_INT, MPI_SUM, comm);
}

void allreduce_sum(float& value, const MPI_Comm comm)
{
    MPI_Allreduce(MPI_IN_PLACE, &value, 1, MPI_FLOAT, MPI_SUM, comm);
}

void allreduce_sum(double& value, const MPI_Comm comm)
{
    MPI_Allreduce(MPI_IN_PLACE, &value, 1, MPI_DOUBLE, MPI_SUM, comm);
}

void allreduce_sum(std::complex<float>& value, const MPI_Comm comm)
{
    float values[2] = {value.real(), value.imag()};
    MPI_Allreduce(MPI_IN_PLACE, values, 2, MPI_FLOAT, MPI_SUM, comm);
    value = std::complex<float>(values[0], values[1]);
}

void allreduce_sum(std::complex<double>& value, const MPI_Comm comm)
{
    double values[2] = {value.real(), value.imag()};
    MPI_Allreduce(MPI_IN_PLACE, values, 2, MPI_DOUBLE, MPI_SUM, comm);
    value = std::complex<double>(values[0], values[1]);
}

void allreduce_sum(std::complex<float>* values, const int count, const MPI_Comm comm)
{
    MPI_Allreduce(
        MPI_IN_PLACE, values, count, MPI_C_FLOAT_COMPLEX, MPI_SUM, comm);
}

void allreduce_sum(std::complex<double>* values, const int count, const MPI_Comm comm)
{
    MPI_Allreduce(
        MPI_IN_PLACE, values, count, MPI_C_DOUBLE_COMPLEX, MPI_SUM, comm);
}

void allreduce_max(double& value, const MPI_Comm comm)
{
    MPI_Allreduce(MPI_IN_PLACE, &value, 1, MPI_DOUBLE, MPI_MAX, comm);
}
#endif

} // namespace Parallel_Common
