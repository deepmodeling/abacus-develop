#ifndef HS_MATRIX_H
#define HS_MATRIX_H

#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"
#include "src_lcao/LCAO_matrix.h"

#include <string>

// mohan add this file 2010-09-10
namespace HS_Matrix
{
    template <typename T>
    void save_HS(
        const int ik,
        const T* H,
        const T* S,
        const bool write_binary,
        const std::string &file_name,
        const Parallel_Orbitals &pv);

    void save_HSR_tr(const int current_spin, LCAO_Matrix &lm); //LiuXh add 2019-07-15

    // jingan add 2021-6-4, modify 2021-12-2
    void save_HSR_sparse(
        LCAO_Matrix &lm,
        const double& sparse_threshold,
        const bool &binary
    );
    void save_SR_sparse(
        LCAO_Matrix &lm,
        const double& sparse_threshold,
        const bool &binary,  
        const std::string &SR_filename
    );
    void output_single_R(std::ofstream &ofs, const std::map<size_t, std::map<size_t, double>> &XR, const double &sparse_threshold, const bool &binary, const Parallel_Orbitals &pv);
    void output_soc_single_R(std::ofstream &ofs, const std::map<size_t, std::map<size_t, std::complex<double>>> &XR, const double &sparse_threshold, const bool &binary, const Parallel_Orbitals &pv);

}

#endif
