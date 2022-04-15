/* myfunc.cpp */
// from LPACK
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <complex>

using namespace std;
#include "myfunc.h"
#include "../module_base/blas_connector.h"
#include "global.h"

// daxpy compute y := alpha * x + y where alpha is a scalar and x and y
void daxpy(const int n, const double &alpha, const double *x, const int incx, double *y, const int incy)
{
    if (n < 1 || incy <= 0 || incx <= 0)
    {
        std::cerr << "\n error in daxpy, n < 1 or incx <= 0 or incy <= 0, ";
        return;
    }
    for (int ix = 0, iy = 0;ix < n && iy < n;ix += incx, iy += incy)
    {
        y[iy] += alpha * x[ix];
    }
    return;
}

//-----------------------------------------------------------------
void dcopy(const ModuleBase::ComplexMatrix &a,
           int i,
           std::complex < double> *y)
{
    // dcopy Copy a(i,:) to y where a is std::complex matrix, and y are n-vectors.
    const int nr = a.nr;
    const int nc = a.nc;

    if (nr < 1 || nc < 1 || i < 0 || i >= nr)
    {
        std::cerr << "\n error in dcopy(ModuleBase::ComplexMatrix a, int i, std::complex < double> *),"
             << " nr or nc < 1 or i out of range ";
        return;
    }

    for (int j = 0;j < nc; j++)
    {
        y[j] = a(i, j);
    }
} // end dcopy

// ------------------------------------
void dcopy(double *x, ModuleBase::matrix &b, int i)
{
    // copy x to ith row of b where b is a matrix and x is a std::vector
    int nr, nc;
    nr = b.nr;
    nc = b.nc;

    if (nr < 1 || nc < 1 || i < 0 || i >= nr)
    {
        std::cerr << "\n error in dcopy(double *x, matrix &b, int i), "
             << "nr or nc < 1 or i out of range ";
        return;
    }

    for (int j = 0;j < nc; j++)
    {
        b(i, j) = x[j];
    }
}

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