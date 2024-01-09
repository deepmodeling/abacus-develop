#ifndef WRITE_HS_H
#define WRITE_HS_H

#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_hamilt_lcao/hamilt_lcaodft/LCAO_matrix.h"

#include <string>

// mohan add this file 2010-09-10
namespace ModuleIO
{
    /// @brief save a square matrix
    /// @param[in] istep : the step of the calculation
    /// @param[in] mat : the local matrix
    /// @param[in] bit : true for binary, false for decimal
    /// @param[in] tri : true for upper triangle, false for full matrix
    /// @param[in] app : true for append, false for overwrite
    /// @param[in] label : the symbol of the matrix, like "H", "S"
    /// @param[in] file_name : the name of the output file
    /// @param[in] pv : the 2d-block parallelization information
    /// @param[in] drank : the rank of the current process
    template<typename T>
    void save_mat(const int istep,
        const T* mat,
        const int dim,
        const bool bit,
        const bool tri,
        const bool app,
        const std::string label,
        const std::string& file_name,
        const Parallel_2D& pv,
        const int drank,
        const int ndigits);

    // comment out this function for not used
    // void save_HSR_tr(const int current_spin, LCAO_Matrix& lm); // LiuXh add 2019-07-15

// mohan comment out 2021-02-10
// void save_HS_ccf(const int &iter, const int &Hnnz, const int *colptr_H, const int *rowind_H, 
// const double *nzval_H, const double *nzval_S, bool bit);
}
#include "write_HS.hpp"
#endif
