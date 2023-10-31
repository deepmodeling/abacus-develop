#include "basic_funcs.h"

#include <iostream>

double maxval_abs_2d(const std::vector<ModuleBase::Vector3<double>>& array)
{
    double max = 0;
    for (const auto& value: array)
    {
            max = std::max(max, std::abs(value.x));
            max = std::max(max, std::abs(value.y));
            max = std::max(max, std::abs(value.z));
    }
    return max;
}

std::pair<int,int> maxloc_abs_2d(const std::vector<ModuleBase::Vector3<double>>& array)
{
    double max = 0;
    int i_max = 0;
    int j_max = 0;
    for (int i = 0; i < array.size(); i++)
    {
        for (int j = 0; j < 3; j++)
        {
            if ((max < std::abs(array[i][j])))
            {
                max = std::abs(array[i][j]);
                i_max = i;
                j_max = j;
            }
        }
    }
    return std::make_pair(i_max, j_max);
}

template <typename T>
T sum_2d(const std::vector<ModuleBase::Vector3<T>>& array)
{
    ModuleBase::Vector3<T> sum;
    for (const auto& element: array)
    {
            sum += element;
    }
    T final_sum = sum.x + sum.y + sum.z;
    return final_sum;
}

// Explicit instantiation
template int sum_2d<int>(const std::vector<ModuleBase::Vector3<int>>& array);
template double sum_2d<double>(const std::vector<ModuleBase::Vector3<double>>& array);

void scalar_multiply_2d(const std::vector<ModuleBase::Vector3<double>>& array,
                        double scalar,
                        std::vector<ModuleBase::Vector3<double>>& result)
{
    int size = array.size();
    result.reserve(size);
    for (int i = 0; i < size; i++)
    {
        result[i] = scalar * array[i];
    }
}

void add_scalar_multiply_2d(const std::vector<ModuleBase::Vector3<double>>& array_1,
                            const std::vector<ModuleBase::Vector3<double>>& array_2,
                            double scalar,
                            std::vector<ModuleBase::Vector3<double>>& result)
{
    int size = array_1.size();
    result.reserve(size);
    for (int i = 0; i < size; i++)
    {
        result[i] = array_1[i] + scalar * array_2[i];
    }
}

void subtract_2d(const std::vector<ModuleBase::Vector3<double>>& array_1,
                 const std::vector<ModuleBase::Vector3<double>>& array_2,
                 std::vector<ModuleBase::Vector3<double>>& result)
{
    int size = array_1.size();
    result.reserve(size);
    for (int i = 0; i < size; i++)
    {
            result[i] = array_1[i] - array_2[i];
    }
}

void fill_scalar_2d(double scalar, std::vector<ModuleBase::Vector3<double>>& result)
{
    for (auto& row: result)
    {
        row.x = scalar;
        row.y = scalar;
        row.z = scalar;
    }
}

void where_fill_scalar_2d(const std::vector<ModuleBase::Vector3<int>>& array_mask,
                          int mask,
                          double scalar,
                          std::vector<ModuleBase::Vector3<double>>& result)
{
    int size = array_mask.size();
    result.resize(size);
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            if (array_mask[i][j] == mask)
            {
                result[i][j] = scalar;
            }
        }
    }
}

void where_fill_scalar_else_2d(const std::vector<ModuleBase::Vector3<int>>& array_mask,
                               int mask,
                               double scalar,
                               const std::vector<ModuleBase::Vector3<double>>& rest,
                               std::vector<ModuleBase::Vector3<double>>& result)
{
    int size = array_mask.size();
    result.resize(size);
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            result[i][j] = (array_mask[i][j] == mask) ? scalar : rest[i][j];
        }
    }
}

void print_2d(std::string info, const std::vector<ModuleBase::Vector3<double>> &array, std::ostream& ofs)
{
    ofs << info << std::endl;
    for (const auto &row : array)
    {
        ofs << std::setprecision(16) << row.x << " " << row.y << " " << row.z << std::endl;
    }
}