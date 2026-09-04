#ifndef PARA_MPI_FUNC_H
#define PARA_MPI_FUNC_H

#include <complex>
#include <string>

#include "para_world.h"

namespace Parallel
{

// Domain-aware MPI communication functions.
// Each function takes the target communication domain (const ParaWorld&)
// explicitly instead of hardcoding MPI_COMM_WORLD/POOL_WORLD.
// In serial builds all functions are no-ops (gather_int copies locally);
// invalid/empty domains are safely skipped.

// ========== Broadcast ==========

void bcast_bool(bool& v, const ParaWorld& world, int root = 0);
void bcast_int(int& v, const ParaWorld& world, int root = 0);
void bcast_double(double& v, const ParaWorld& world, int root = 0);
void bcast_complex(std::complex<double>& v, const ParaWorld& world, int root = 0);
void bcast_string(std::string& s, const ParaWorld& world, int root = 0);

void bcast_int(int* v, int n, const ParaWorld& world, int root = 0);
void bcast_double(double* v, int n, const ParaWorld& world, int root = 0);
void bcast_complex(std::complex<double>* v, int n, const ParaWorld& world, int root = 0);
void bcast_char(char* v, int n, const ParaWorld& world, int root = 0);
void bcast_string(std::string* v, int n, const ParaWorld& world, int root = 0);

// ========== Reduce (Allreduce, result on all ranks) ==========

void reduce_all(double& v, const ParaWorld& world);
void reduce_all(int& v, const ParaWorld& world);
void reduce_all(double* v, int n, const ParaWorld& world);

// ========== Reduce min/max ==========

void reduce_min(double& v, const ParaWorld& world);
void reduce_max(double& v, const ParaWorld& world);
void reduce_min(int& v, const ParaWorld& world);
void reduce_max(int& v, const ParaWorld& world);

// ========== Barrier ==========

void barrier(const ParaWorld& world);

// ========== Gather ==========

void gather_int(int& v, int* all, const ParaWorld& world);

} // namespace Parallel

#endif // PARA_MPI_FUNC_H
