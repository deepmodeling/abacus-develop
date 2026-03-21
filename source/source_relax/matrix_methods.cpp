#include "matrix_methods.h"




std::vector<double> ReshapeMToV(std::vector<ModuleBase::Vector3<double>>& matrix) 
{
    assert(!matrix.empty());
    int size = matrix.size();
    std::vector<double> result;
    result.reserve(3*size);
    for (const auto& v : matrix) 
    {
        result.push_back(v.x);
        result.push_back(v.y);
        result.push_back(v.z);
    }
    return result;
}

std::vector<std::vector<double>> MAddM(std::vector<std::vector<double>>& a, 
                                             std::vector<std::vector<double>>& b) 
{
    assert(!a.empty() && !b.empty());
    assert(a.size() == b.size() && a[0].size() == b[0].size());
    std::vector<std::vector<double>> result = std::vector<std::vector<double>>(a.size(), std::vector<double>(a[0].size(), 0.0));
    for(int i = 0; i < a.size(); i++)
    {
        for(int j = 0; j < a[0].size(); j++)
        {
            result[i][j] = a[i][j] + b[i][j];
        }
    }
    return result;
}

std::vector<double> VSubV(std::vector<double>& a, std::vector<double>& b) 
{
    assert(a.size() == b.size());
    std::vector<double> result = std::vector<double>(a.size(), 0.0);
    for(int i = 0; i < a.size(); i++)
    {
        result[i] = a[i] - b[i];
    }
    return result;
}

std::vector<ModuleBase::Vector3<double>> ReshapeVToM(std::vector<double>& matrix) 
{
    assert(matrix.size() % 3 == 0);
    std::vector<ModuleBase::Vector3<double>> result = std::vector<ModuleBase::Vector3<double>>(matrix.size() / 3, ModuleBase::Vector3<double>(0.0, 0.0, 0.0));
    for(int i = 0; i < result.size(); i++)
    {
        result[i].x = matrix[i*3 ];
        result[i].y = matrix[i*3 + 1];
        result[i].z = matrix[i*3 + 2];
    }
    return result;
}

std::vector<double> DotInMAndV1(std::vector<std::vector<double>>& matrix, std::vector<double>& vec) 
{
    assert(!matrix.empty());
    assert(matrix[0].size() == vec.size());
    std::vector<double> result(matrix.size(), 0.0);
    for(int i = 0; i < result.size(); i++)
    {
        for(int j = 0; j < vec.size(); j++)
        {
            result[i] += matrix[i][j] * vec[j];
        }
    }
    return result;
}
std::vector<double> DotInMAndV2(std::vector<std::vector<double>>& matrix, std::vector<double>& vec) 
{
    assert(!matrix.empty());
    assert(matrix.size() == vec.size());
    std::vector<double> result(matrix.size(), 0.0);
    for(int i = 0; i < result.size(); i++)
    {
        for(int j = 0; j < vec.size(); j++)
        {
            result[i] += matrix[j][i] * vec[j];
        }
    }
    return result;
}

double DotInVAndV(std::vector<double>& vec1, std::vector<double>& vec2) 
{
    assert(vec1.size() == vec2.size());
    double result = 0.0;
    for(int i = 0; i < vec1.size(); i++)
    {
        result += vec1[i] * vec2[i];
    }
    return result;
}

std::vector<std::vector<double>> OuterVAndV(std::vector<double>& a, std::vector<double>& b) 
{
    assert(a.size() == b.size());
    std::vector<std::vector<double>> result = std::vector<std::vector<double>>(a.size(), std::vector<double>(b.size(), 0.0));
    for(int i = 0; i < a.size(); i++)
    {
        for(int j = 0; j < b.size(); j++)
        {
            result[i][j] = a[i] * b[j];
        }
    }
    return result;
}

std::vector<std::vector<double>> MPlus(std::vector<std::vector<double>>& a, double b)
{
    assert(!a.empty());
    assert(b != 0);
    std::vector<std::vector<double>> result = std::vector<std::vector<double>>(a.size(), std::vector<double>(a[0].size(), 0.0));
    for(int i = 0; i < a.size(); i++)
    {
        for(int j = 0; j < a[0].size(); j++)
        {
            result[i][j] = a[i][j] / b;
        }
    }
    return result;
}

std::vector<std::vector<double>> MSubM(std::vector<std::vector<double>>& a, std::vector<std::vector<double>>& b)
{
    assert(!a.empty() && !b.empty());
    assert(a.size() == b.size() && a[0].size() == b[0].size());
    std::vector<std::vector<double>> result = std::vector<std::vector<double>>(a.size(), std::vector<double>(a[0].size(), 0.0));
    for(int i = 0; i < a.size(); i++)
    {
        for(int j = 0; j < a[0].size(); j++)
        {
            result[i][j] = a[i][j] - b[i][j];
        }
    }
    return result;
}

std::vector<double> DotInVAndFloat(std::vector<double>& vec, double b) 
{
    assert(b != 0);
    std::vector<double> result(vec.size(), 0.0);
    for(int i = 0; i < vec.size(); i++)
    {
        result[i] = vec[i] * b;
    }
    return result;
}

std::vector<double> VAddV(std::vector<double>& a, std::vector<double>& b) 
{
    assert(a.size() == b.size());
    std::vector<double> result = std::vector<double>(a.size(), 0.0);
    for(int i = 0; i < a.size(); i++)
    {
        result[i] = a[i] + b[i];
    }
    return result;
}
