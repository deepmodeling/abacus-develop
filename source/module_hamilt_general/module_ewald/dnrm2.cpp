#include "dnrm2.h"

#include<math.h>
#include <iostream>

double dnrm2(const int n, const double *x, const int incx)
{
    // compute Euclidean length (12 norm) of std::vector x,
    if (n < 0 || incx <= 0)
    {
        std::cerr << "\n error in dnrm2, n < 0 or incx <= 0, ";
        return 0;
    }
    if (n == 0)
    {
        return 0;
    }

    double norm2=0.0;
    for (int ix=0; ix<n; ix+=incx)
    {
        norm2 += x[ix] * x[ix];
    }

    return sqrt(norm2);
}