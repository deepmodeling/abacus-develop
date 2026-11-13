#ifndef PARALLEL_COMMON_H
#define PARALLEL_COMMON_H

#ifdef __MPI
#include "mpi.h"
#endif

#include <algorithm>
#include <complex>
#include <string>

namespace Parallel_Common
{

#ifdef __MPI
inline int comm_size(const MPI_Comm comm)
{
    int size = 0;
    MPI_Comm_size(comm, &size);
    return size;
}

inline int comm_rank(const MPI_Comm comm)
{
    int rank = 0;
    MPI_Comm_rank(comm, &rank);
    return rank;
}
#else
inline int comm_size()
{
    return 1;
}

inline int comm_rank()
{
    return 0;
}
#endif

//(1) bcast array
void bcast_complex_double(std::complex<double>* object, const int n);
void bcast_string(std::string* object, const int n);
void bcast_double(double* object, const int n);
void bcast_int(int* object, const int n);
void bcast_char(char* object, const int n);

//(2) bcast single
void bcast_complex_double(std::complex<double>& object);
void bcast_string(std::string& object);
void bcast_double(double& object);
void bcast_int(int& object);
void bcast_bool(bool& object);

#ifdef __MPI
//(3) bcast single on an explicit communicator, from the communicator's rank 0
void bcast_mpi_offset(MPI_Offset& object, const MPI_Comm& comm);

void barrier(const MPI_Comm comm);

void gatherv_bytes(const void* send_buffer,
                   int send_count,
                   void* receive_buffer,
                   const int* receive_counts,
                   const int* receive_displacements,
                   int root,
                   MPI_Comm comm);

void allgatherv_bytes(const void* send_buffer,
                      int send_count,
                      void* receive_buffer,
                      const int* receive_counts,
                      const int* receive_displacements,
                      MPI_Comm comm);

void scatterv_bytes(const void* send_buffer,
                    const int* send_counts,
                    const int* send_displacements,
                    void* receive_buffer,
                    int receive_count,
                    int root,
                    MPI_Comm comm);

void allreduce_sum(int& value, MPI_Comm comm);
void allreduce_sum(float& value, MPI_Comm comm);
void allreduce_sum(double& value, MPI_Comm comm);
void allreduce_sum(std::complex<float>& value, MPI_Comm comm);
void allreduce_sum(std::complex<double>& value, MPI_Comm comm);
void allreduce_sum(std::complex<float>* values, int count, MPI_Comm comm);
void allreduce_sum(std::complex<double>* values, int count, MPI_Comm comm);
void allreduce_max(double& value, MPI_Comm comm);
#else
inline void barrier()
{
}

inline void gatherv_bytes(const void* send_buffer,
                          int send_count,
                          void* receive_buffer,
                          const int*,
                          const int*,
                          int,
                          int)
{
    if (receive_buffer != nullptr && send_count > 0)
    {
        std::copy(static_cast<const unsigned char*>(send_buffer),
                  static_cast<const unsigned char*>(send_buffer) + send_count,
                  static_cast<unsigned char*>(receive_buffer));
    }
}

inline void allgatherv_bytes(const void* send_buffer,
                             int send_count,
                             void* receive_buffer,
                             const int*,
                             const int*,
                             int)
{
    if (receive_buffer != nullptr && send_count > 0)
    {
        std::copy(static_cast<const unsigned char*>(send_buffer),
                  static_cast<const unsigned char*>(send_buffer) + send_count,
                  static_cast<unsigned char*>(receive_buffer));
    }
}

inline void scatterv_bytes(const void* send_buffer,
                           const int*,
                           const int*,
                           void* receive_buffer,
                           int receive_count,
                           int,
                           int)
{
    if (receive_buffer != nullptr && receive_count > 0)
    {
        std::copy(static_cast<const unsigned char*>(send_buffer),
                  static_cast<const unsigned char*>(send_buffer) + receive_count,
                  static_cast<unsigned char*>(receive_buffer));
    }
}
#endif

} // namespace Parallel_Common

#endif
