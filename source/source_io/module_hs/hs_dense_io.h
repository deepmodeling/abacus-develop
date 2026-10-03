#ifndef HS_DENSE_IO_H
#define HS_DENSE_IO_H

#include <string>

#include "source_base/parallel_2d.h"

// mohan note: this file holds the dense square-matrix writer
// (formerly save_mat in write_hs.h / write_hs.hpp).

namespace ModuleIO
{

/// @brief save a square matrix, such as H(k) and S(k)
/// @param[in] istep : the step of the calculation
/// @param[in] mat : the local matrix
/// @param[in] dim : the dimension of the square matrix
/// @param[in] bit : true for binary, false for decimal
/// @param[in] precision : the precision of the decimal output
/// @param[in] tri : true for upper triangle, false for full matrix
/// @param[in] app : true for append, false for overwrite
/// @param[in] file_name : the name of the output file
/// @param[in] pv : the 2d-block parallelization information
/// @param[in] drank : the rank of the current process in the diagonalization world
/// @param[in] ks_solver : the name of the eigensolver (decides row/column major)
/// @param[in] reduce : whether to reduce across ranks before writing
template <typename T>
void save_mat(const int istep,
              const T* mat,
              const int dim,
              const bool bit,
              const int precision,
              const bool tri,
              const bool app,
              const std::string& file_name,
              const Parallel_2D& pv,
              const int drank,
              const std::string& ks_solver,
              const bool reduce = true);

} // namespace ModuleIO

#endif // HS_DENSE_IO_H
