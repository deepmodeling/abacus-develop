#ifndef BASIC_FUNCS_H
#define BASIC_FUNCS_H

#include <cmath>
#include <vector>

#include "module_base/vector3.h"

/**
 * @brief Find the maximum absolute value in a 2D array.
 */
double maxval_abs_2d(const std::vector<std::vector<double>>& array);

/**
 * @brief Find the maximum absolute value in a 2D array and its index.
 */
void maxloc_abs_2d(const std::vector<std::vector<double>>& array, std::vector<int>& result);

/**
 * @brief sum of all elements in a 2D array.
 */
double sum_2d(const std::vector<std::vector<double>>& array);

/**
 * @brief scalar multiply a 2D array.
 */
void scalar_multiply_2d(const std::vector<std::vector<double>>& array,
                        double scalar,
                        std::vector<std::vector<double>>& result);

/**
 * @brief array_1 + scalar * array_2.
 */
void add_scalar_multiply_2d(const std::vector<ModuleBase::Vector3<double>>& array_1,
                            const std::vector<ModuleBase::Vector3<double>>& array_2,
                            double scalar,
                            std::vector<ModuleBase::Vector3<double>>& result);

/**
 * @brief array_1 - array_2.
 */
void subtract_2d(const std::vector<std::vector<double>>& array_1,
                 const std::vector<std::vector<double>>& array_2,
                 std::vector<std::vector<double>>& result);

/**
 * @brief fill a 2D array with a scalar.
 */
void fill_scalar_2d(double scalar, std::vector<std::vector<double>>& result);

/**
 * @brief fill a 2D array with a scalar if the corresponding element is equal to mask.
 */
void where_fill_scalar_2d(const std::vector<std::vector<int>>& array_mask,
                          int mask,
                          double scalar,
                          std::vector<std::vector<double>>& result);

void where_fill_scalar_else_2d(const std::vector<ModuleBase::Vector3<int>>& array_mask,
                               int mask,
                               double scalar,
                               const std::vector<ModuleBase::Vector3<double>>& rest,
                               std::vector<ModuleBase::Vector3<double>>& result);

void print_2d(const std::vector<std::vector<double>>& array);

void optimize_lambda(const std::vector<std::vector<double>>& M_CONSTR,
                     const std::vector<std::vector<int>>& CONSTRL,
                     const int NIONS,
                     const int NTYP,
                     const std::vector<int>& NITYP,
                     const double INISC,
                     const double SCDIFF,
                     const std::vector<double>& SCCONV_GRAD,
                     const int NSC,
                     const int NSCMIN,
                     const double SCCUT,
                     const int N,
                     std::vector<std::vector<double>>& MW,
                     std::vector<std::vector<double>>& OUT_LAMBDA);

#endif // BASIC_FUNCS_H