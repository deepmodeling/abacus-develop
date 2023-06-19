#ifndef WRITE_HS_H
#define WRITE_HS_H

#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_hamilt_lcao/hamilt_lcaodft/LCAO_matrix.h"

#include <string>

// mohan add this file 2010-09-10
namespace ModuleIO
{
void saving_HS(const int istep,
               const double* Hloc,
               const double* Sloc,
               const bool bit,
               const int& out_hs,
               const std::string& file_name,
               const Parallel_Orbitals& pv,
               bool tri);
// tri : if set 1 , output triangle matrix ; if set 0 , output complete matrix .

void save_HS_triangle(const int istep,
                      const double* H,
                      const double* S,
                      const bool bit,
                      const std::string& file_name,
                      const Parallel_Orbitals& pv);
void save_HS_complete(const int istep,
                      const double* H,
                      const double* S,
                      const bool bit,
                      const std::string& file_name,
                      const Parallel_Orbitals& pv);

void save_HS_complex(const int istep,
                     const std::complex<double>* H,
                     const std::complex<double>* S,
                     const bool bit,
                     const std::string& file__name,
                     const Parallel_Orbitals& pv,
                     bool tri);

void save_HSR_tr(const int current_spin, LCAO_Matrix& lm); // LiuXh add 2019-07-15

// mohan comment out 2021-02-10
// void save_HS_ccf(const int &iter, const int &Hnnz, const int *colptr_H, const int *rowind_H, 
// const double *nzval_H, const double *nzval_S, bool bit);

void saving_HS(const int istep,
               std::complex<double>* Hloc,
               std::complex<double>* Sloc,
               bool bit,
               const int& out_hs,
               const std::string& file__name,
               const Parallel_Orbitals& pv,
               bool tri); // LiuXh, 2017-03-21

void save_HS_complex_triangle(const int istep,
                              std::complex<double>* H,
                              std::complex<double>* S,
                              const bool bit,
                              const std::string& file__name,
                              const Parallel_Orbitals& pv); // LiuXh, 2017-03-21

void save_HS_complex_complete(const int istep,
                              std::complex<double>* H,
                              std::complex<double>* S,
                              const bool bit,
                              const std::string& file__name,
                              const Parallel_Orbitals& pv); // LiuXh, 2017-03-21
}

#endif
