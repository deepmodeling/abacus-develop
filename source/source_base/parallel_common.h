#ifndef W_ABACUS_DEVELOP_ABACUS_DEVELOP_SOURCE_MODULE_BASE_PARALLEL_COMMON_H
#define W_ABACUS_DEVELOP_ABACUS_DEVELOP_SOURCE_MODULE_BASE_PARALLEL_COMMON_H

#ifdef __MPI
#include "mpi.h"
#endif
#include <complex>
#include <string>

namespace Parallel_Common
{
//(1) bcast array
void bcast_complex_double(std::complex<double>* object, const int n
#ifdef __MPI
    , MPI_Comm comm = MPI_COMM_WORLD
#endif
    );
void bcast_string(std::string* object, const int n
#ifdef __MPI
    , MPI_Comm comm = MPI_COMM_WORLD
#endif
    );
void bcast_double(double* object, const int n
#ifdef __MPI
    , MPI_Comm comm = MPI_COMM_WORLD
#endif
    );
void bcast_int(int* object, const int n
#ifdef __MPI
    , MPI_Comm comm = MPI_COMM_WORLD
#endif
    );
void bcast_char(char* object, const int n
#ifdef __MPI
    , MPI_Comm comm = MPI_COMM_WORLD
#endif
    );

//(2) bcast single
void bcast_complex_double(std::complex<double>& object
#ifdef __MPI
    , MPI_Comm comm = MPI_COMM_WORLD
#endif
    );
void bcast_string(std::string& object
#ifdef __MPI
    , MPI_Comm comm = MPI_COMM_WORLD
#endif
    );
void bcast_double(double& object
#ifdef __MPI
    , MPI_Comm comm = MPI_COMM_WORLD
#endif
    );
void bcast_int(int& object
#ifdef __MPI
    , MPI_Comm comm = MPI_COMM_WORLD
#endif
    );
void bcast_bool(bool& object
#ifdef __MPI
    , MPI_Comm comm = MPI_COMM_WORLD
#endif
    );

} // namespace Parallel_Common

#endif
