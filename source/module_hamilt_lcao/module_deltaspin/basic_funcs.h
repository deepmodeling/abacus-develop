#ifndef BASIC_FUNCS_H
#define BASIC_FUNCS_H

#include <cmath>
#include <vector>
#include <ostream>

#include "module_base/vector3.h"

/**
 * @brief Find the maximum absolute value in a 2D array.
 */
double maxval_abs_2d(const std::vector<ModuleBase::Vector3<double>>& array);

/**
 * @brief Find the maximum absolute value in a 2D array and its index.
 */
std::pair<int,int> maxloc_abs_2d(const std::vector<ModuleBase::Vector3<double>>& array);

/**
 * @brief sum of all elements in a 2D array.
 */
template <typename T>
T sum_2d(const std::vector<ModuleBase::Vector3<T>>& array);

/**
 * @brief scalar multiply a 2D array.
 */
void scalar_multiply_2d(const std::vector<ModuleBase::Vector3<double>>& array,
                        double scalar,
                        std::vector<ModuleBase::Vector3<double>>& result);

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
void subtract_2d(const std::vector<ModuleBase::Vector3<double>>& array_1,
                 const std::vector<ModuleBase::Vector3<double>>& array_2,
                 std::vector<ModuleBase::Vector3<double>>& result);

/**
 * @brief fill a 2D array with a scalar.
 */
void fill_scalar_2d(double scalar, std::vector<ModuleBase::Vector3<double>>& result);

/**
 * @brief fill a 2D array with a scalar if the corresponding element is equal to mask.
 */
void where_fill_scalar_2d(const std::vector<ModuleBase::Vector3<int>>& array_mask,
                          int mask,
                          double scalar,
                          std::vector<ModuleBase::Vector3<double>>& result);

void where_fill_scalar_else_2d(const std::vector<ModuleBase::Vector3<int>>& array_mask,
                               int mask,
                               double scalar,
                               const std::vector<ModuleBase::Vector3<double>>& rest,
                               std::vector<ModuleBase::Vector3<double>>& result);

void print_2d(std::string info, const std::vector<ModuleBase::Vector3<double>> &array, std::ostream& ofs = std::cout);

#endif // BASIC_FUNCS_H