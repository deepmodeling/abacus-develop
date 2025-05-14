#pragma once

#include <memory>
#include <cmath>
#include "gint_type.h"
#include "module_base/timer.h"

template <typename T>
std::shared_ptr<const T> toConstSharedPtr(std::shared_ptr<T> ptr) {
    return std::static_pointer_cast<const T>(ptr);
}


inline int index3Dto1D(const int id_x, const int id_y, const int id_z,
                        const int dim_x, const int dim_y, const int dim_z)
{
    return id_z + id_y * dim_z + id_x * dim_y * dim_z;
};

inline Vec3i index1Dto3D(const int index_1d,
                            const int dim_x, const int dim_y, const int dim_z)
{
    int id_x = index_1d / (dim_y * dim_z);
    int id_y = (index_1d - id_x * dim_y * dim_z) / dim_z;
    int id_z = index_1d % dim_z;
    return Vec3i(id_x, id_y, id_z);
};

// if exponent is an integer between 0 and 5 (the most common cases in gint) and
// and exp is a variable that cannot be determined at compile time (which means the compiler cannot optimize the code),
// pow_int is much faster than std::pow
inline double pow_int(const double base, const int exp)
{
    switch (exp)
    {
    case 0:
        return 1.0;
    case 1:
        return base;
    case 2:
        return base * base;
    case 3:
        return base * base * base;
    case 4:
        return base * base * base * base;
    case 5:
        return base * base * base * base * base;
    default:
        double result = std::pow(base, exp);
        return result;
    }
};

inline int floor_div(const int a, const int b)
{
    // a ^ b < 0 means a and b have different signs
    return a / b - (a % b != 0 && (a ^ b) < 0);
};

inline int ceil_div(const int a, const int b)
{
    return a / b + (a % b != 0 && (a ^ b) > 0); 
};