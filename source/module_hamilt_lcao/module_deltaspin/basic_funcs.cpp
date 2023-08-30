#include "basic_funcs.h"

#include <iostream>

double maxval_abs_2d(const std::vector<std::vector<double>>& array)
{
    double max = 0;
    for (const auto& row: array)
    {
        for (double value: row)
        {
            max = std::max(max, std::abs(value));
        }
    }
    return max;
}

void maxloc_abs_2d(const std::vector<std::vector<double>>& array, std::vector<int>& result)
{
    double max = 0;
    int size_1 = array.size();
    int size_2 = array[0].size();
    for (int i = 0; i < size_1; i++)
    {
        for (int j = 0; j < size_2; j++)
        {
            if ((max < abs(array[i][j])))
            {
                max = abs(array[i][j]);
                result[0] = i;
                result[1] = j;
            }
        }
    }
}

double sum_2d(const std::vector<ModuleBase::Vector3<double>>& array)
{
    double sum = 0;
    for (const auto& element: array)
    {
            sum += element.x;
            sum += element.y;
            sum += element.z;
    }
    return sum;
}

void scalar_multiply_2d(const std::vector<std::vector<double>>& array,
                        double scalar,
                        std::vector<std::vector<double>>& result)
{
    int size_1 = array.size();
    int size_2 = array[0].size();
    result.resize(size_1);
    for (int i = 0; i < size_1; i++)
    {
        result[i].resize(size_2);
        for (int j = 0; j < size_2; j++)
        {
            result[i][j] = scalar * array[i][j];
        }
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

void subtract_2d(const std::vector<std::vector<double>>& array_1,
                 const std::vector<std::vector<double>>& array_2,
                 std::vector<std::vector<double>>& result)
{
    int size_1 = array_1.size();
    int size_2 = array_1[0].size();
    result.resize(size_1);
    for (int i = 0; i < size_1; i++)
    {
        result[i].resize(size_2);
        for (int j = 0; j < size_2; j++)
        {
            result[i][j] = array_1[i][j] - array_2[i][j];
        }
    }
}

void fill_scalar_2d(double scalar, std::vector<std::vector<double>>& result)
{
    for (auto& row: result)
    {
        std::fill(row.begin(), row.end(), scalar);
    }
}

void where_fill_scalar_2d(const std::vector<std::vector<int>>& array_mask,
                          int mask,
                          double scalar,
                          std::vector<std::vector<double>>& result)
{
    int size_1 = array_mask.size();
    int size_2 = array_mask[0].size();
    result.resize(size_1);
    for (int i = 0; i < size_1; i++)
    {
        result[i].resize(size_2);
        for (int j = 0; j < size_2; j++)
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
    result.reserve(size);
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            result[i][j] = (array_mask[i][j] == mask) ? scalar : rest[i][j];
        }
    }
}

void print_2d(const std::vector<std::vector<double>>& array)
{
    for (const auto& row: array)
    {
        for (const auto& element: row)
        {
            std::cout << element << " ";
        }
        std::cout << std::endl;
    }
}