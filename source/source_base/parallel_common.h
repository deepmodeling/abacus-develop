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

void gatherv_bytes(const void* sendbuf,
                   int send,
                   void* recvbuf,
                   const int* recvcounts,
                   const int* rdispls,
                   int root,
                   MPI_Comm comm);

void allgatherv_bytes(const void* sendbuf,
                      int send,
                      void* recvbuf,
                      const int* recvcounts,
                      const int* rdispls,
                      MPI_Comm comm);

void scatterv_bytes(const void* sendbuf,
                    const int* sendcounts,
                    const int* sdispls,
                    void* recvbuf,
                    int recv,
                    int root,
                    MPI_Comm comm);

void allreduce_sum(int& value, MPI_Comm comm);
void allreduce_sum(float& value, MPI_Comm comm);
void allreduce_sum(double& value, MPI_Comm comm);
void allreduce_sum(std::complex<float>& value, MPI_Comm comm);
void allreduce_sum(std::complex<double>& value, MPI_Comm comm);
void allreduce_max(double& value, MPI_Comm comm);
#else
inline void barrier()
{
}

inline void gatherv_bytes(const void* sendbuf,
                          int send,
                          void* recvbuf,
                          const int*,
                          const int*,
                          int,
                          int)
{
    if (recvbuf != nullptr && send > 0)
    {
        std::copy(static_cast<const unsigned char*>(sendbuf),
                  static_cast<const unsigned char*>(sendbuf) + send,
                  static_cast<unsigned char*>(recvbuf));
    }
}

inline void allgatherv_bytes(const void* sendbuf,
                             int send,
                             void* recvbuf,
                             const int*,
                             const int*,
                             int)
{
    if (recvbuf != nullptr && send > 0)
    {
        std::copy(static_cast<const unsigned char*>(sendbuf),
                  static_cast<const unsigned char*>(sendbuf) + send,
                  static_cast<unsigned char*>(recvbuf));
    }
}

inline void scatterv_bytes(const void* sendbuf,
                           const int*,
                           const int*,
                           void* recvbuf,
                           int recv,
                           int,
                           int)
{
    if (recvbuf != nullptr && recv > 0)
    {
        std::copy(static_cast<const unsigned char*>(sendbuf),
                  static_cast<const unsigned char*>(sendbuf) + recv,
                  static_cast<unsigned char*>(recvbuf));
    }
}
#endif

} // namespace Parallel_Common

#endif
