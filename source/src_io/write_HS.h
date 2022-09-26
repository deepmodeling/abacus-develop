#ifndef HS_MATRIX_H
#define HS_MATRIX_H

#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"
#include "src_lcao/LCAO_matrix.h"

#include <string>

// mohan add this file 2010-09-10
namespace HS_Matrix
{
    void saving_HS(const double *Hloc, const double* Sloc, const bool bit, const int &out_hs, const std::string &file_name, const Parallel_Orbitals &pv);

    void save_HS(const double *H, const double *S, const bool bit, const std::string &file_name, const Parallel_Orbitals &pv);

    void save_HS_complex(const std::complex<double> *H, const std::complex<double> *S, const bool bit, const std::string &file__name, const Parallel_Orbitals &pv);

    void save_HSR_tr(const int current_spin, LCAO_Matrix &lm); //LiuXh add 2019-07-15

    // jingan add 2021-6-4, modify 2021-12-2
    void save_HSR_sparse(
        const int &istep,
        LCAO_Matrix &lm,
        const double& sparse_threshold,
        const bool &binary,  
        const std::string &SR_filename, 
        const std::string &HR_filename_up, 
        const std::string &HR_filename_down
    );
    void save_SR_sparse(
        LCAO_Matrix &lm,
        const double& sparse_threshold,
        const bool &binary,  
        const std::string &SR_filename
    );
    void output_single_R(std::ofstream &ofs, const std::map<size_t, std::map<size_t, double>> &XR, const double &sparse_threshold, const bool &binary, const Parallel_Orbitals &pv);
    void output_soc_single_R(std::ofstream &ofs, const std::map<size_t, std::map<size_t, std::complex<double>>> &XR, const double &sparse_threshold, const bool &binary, const Parallel_Orbitals &pv);

// mohan comment out 2021-02-10
// void save_HS_ccf(const int &iter, const int &Hnnz, const int *colptr_H, const int *rowind_H, 
// const double *nzval_H, const double *nzval_S, bool bit);

    void saving_HS(std::complex<double> *Hloc, std::complex<double>* Sloc, bool bit, const int &out_hs, const std::string &file__name, const Parallel_Orbitals &pv); //LiuXh, 2017-03-21

    void save_HS_complex(std::complex<double> *H, std::complex<double> *S, const bool bit, const std::string &file__name, const Parallel_Orbitals &pv); //LiuXh, 2017-03-21
}

#endif
