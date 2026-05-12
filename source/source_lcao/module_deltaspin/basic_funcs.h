#ifndef BASIC_FUNCS_H
#define BASIC_FUNCS_H

#include <cmath>
#include <complex>
#include <vector>
#include <ostream>

#include "source_base/vector3.h"

/**
 * @file basic_funcs.h
 * @brief Utility vector/array operations for per-atom 3D vector arrays.
 *
 * @par Data structure
 * All functions operate on std::vector<ModuleBase::Vector3<double>>, which
 * represents a 2D array of shape [nat][3] where:
 * - First index: atom index (iat = 0 to nat-1)
 * - Second index: component (x=0, y=1, z=2)
 *
 * These are NumPy-style element-wise operations used throughout the lambda
 * optimization loop for manipulating magnetic moments, lambda values,
 * search directions, and constraint masks.
 */

/**
 * @brief Find the maximum absolute value across all atoms and components.
 * @return max(|array[iat][ic]|) for all iat, ic
 */
double maxval_abs_2d(const std::vector<ModuleBase::Vector3<double>>& array);

/**
 * @brief Find the maximum absolute value and its (atom, component) index.
 * @return pair<iat, ic> of the element with maximum absolute value
 */
std::pair<int,int> maxloc_abs_2d(const std::vector<ModuleBase::Vector3<double>>& array);

/**
 * @brief Sum of all elements across all atoms and components.
 * @tparam T Numeric type (int or double)
 * @return sum(array[iat][ic]) for all iat, ic
 */
template <typename T>
T sum_2d(const std::vector<ModuleBase::Vector3<T>>& array);

/**
 * @brief Element-wise scalar multiplication: result = scalar * array.
 */
void scalar_multiply_2d(const std::vector<ModuleBase::Vector3<double>>& array,
                        double scalar,
                        std::vector<ModuleBase::Vector3<double>>& result);

/**
 * @brief Element-wise fused multiply-add: result = array_1 + scalar * array_2.
 */
void add_scalar_multiply_2d(const std::vector<ModuleBase::Vector3<double>>& array_1,
                            const std::vector<ModuleBase::Vector3<double>>& array_2,
                            double scalar,
                            std::vector<ModuleBase::Vector3<double>>& result);

/**
 * @brief Element-wise subtraction: result = array_1 - array_2.
 */
void subtract_2d(const std::vector<ModuleBase::Vector3<double>>& array_1,
                 const std::vector<ModuleBase::Vector3<double>>& array_2,
                 std::vector<ModuleBase::Vector3<double>>& result);

/**
 * @brief Fill all elements with a scalar value.
 */
void fill_scalar_2d(double scalar, std::vector<ModuleBase::Vector3<double>>& result);

/**
 * @brief Conditional fill: if mask[iat][ic] == value, set result[iat][ic] = scalar.
 *
 * Used to mask unconstrained components to zero:
 *   where_fill_scalar_2d(constrain_, 0, 0.0, delta_spin)
 * sets delta_spin[ia][ic] = 0 where constrain[ia][ic] == 0.
 */
void where_fill_scalar_2d(const std::vector<ModuleBase::Vector3<int>>& array_mask,
                          int mask,
                          double scalar,
                          std::vector<ModuleBase::Vector3<double>>& result);

/**
 * @brief Conditional fill with else branch: if mask == value, set scalar; otherwise copy from rest.
 *
 * Used to create masked copies:
 *   where_fill_scalar_else_2d(constrain_, 0, 0.0, lambda_, initial_lambda)
 * sets initial_lambda[ia][ic] = 0 if unconstrained, else lambda_[ia][ic].
 */
void where_fill_scalar_else_2d(const std::vector<ModuleBase::Vector3<int>>& array_mask,
                               int mask,
                               double scalar,
                               const std::vector<ModuleBase::Vector3<double>>& rest,
                               std::vector<ModuleBase::Vector3<double>>& result);

/**
 * @brief Formatted print of a 2D array.
 * @param info Header string
 * @param array Data to print
 * @param nspin Spin type: 2=z-only, 4=xyz
 * @param unit_convert Multiplicative factor (e.g., Ry_to_eV for unit conversion)
 * @param ofs Output stream (default: stdout)
 */
void print_2d(const std::string info, const std::vector<ModuleBase::Vector3<double>> &array, const int nspin, const double unit_convert = 1.0, std::ostream& ofs = std::cout);

#endif // BASIC_FUNCS_H
