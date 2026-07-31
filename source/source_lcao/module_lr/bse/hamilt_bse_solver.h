#pragma once

#include "source_base/parallel_2d.h"
#include "source_base/parallel_reduce.h"
#include "source_base/module_external/blacs_connector.h"
#include "source_base/module_external/scalapack_connector.h"
#include "source_lcao/module_lr/utils/lr_util.h"
#define HAVE_SKEWSYMMETRIC // for elpa_skew_eigenvectors
#include <elpa/elpa.h>
#ifdef I // to avoid conflict with macro I defined in elpa
#undef I
#endif
#include <omp.h>

namespace BSE
{
template <typename T>
void printM(const std::vector<T>& A, int m, int n, std::string file, std::string name);

/// @brief result={{A,B},{-A*,-B*}}
void arrayFlatten1(const std::vector<std::complex<double>>& A,
                   const std::vector<std::complex<double>>& B,
                   std::vector<std::complex<double>>& result,
                   const Parallel_2D& pA,
                   const Parallel_2D& pM);

/// @brief result={{Re(A+B), Im(A-B)},{-Im(A+B), Re(A-B)}}
void arrayFlatten2(const std::vector<std::complex<double>>& A,
                   const std::vector<std::complex<double>>& B,
                   std::vector<double>& result,
                   const Parallel_2D& pA,
                   const Parallel_2D& pM);

/// @brief solve full BSE through matrix M = {{Re(A+B), Im(A-B)},{-Im(A+B), Re(A-B)}}
/// @attention A_part and B_part will be destroyed (swapped with empty vectors)
///          after constructing M. Pass copies if the caller needs them afterwards.
void solve_full(const int my_rank,
                std::vector<std::complex<double>>& A_part,
                std::vector<std::complex<double>>& B_part,
                const Parallel_2D& pA,
                const Parallel_2D& pM,
                std::vector<double>& ev,
                std::vector<std::complex<double>>& global_v);
/// @brief template solve TDA BSE, implemented in hamilt_bse_solver.hpp
template <typename T>
void solve_tda(const int my_rank,
               const std::vector<T>& A,
               const Parallel_2D& pA,
               std::vector<double>& ev,
               std::vector<T>& v);
}

#include "hamilt_bse_solver.hpp"