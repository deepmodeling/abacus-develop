#ifndef MYFUNC_H
#define MYFUNC_H
using namespace std;

#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"
#include "../module_base/matrix.h"
#include "../module_base/vector3.h"
#include "../module_base/complexmatrix.h"
//#include "../src_algorithms/mymath.h"     //  only wgauss(),wlgauss()
#include <complex>

// compute y = alpha * x + y where alpha is a scalar and x and y are n-vectors
void daxpy(const int n, const double &alpha, const double *x, const int incx, double *y, const int incy);
// copy ith row of a to y where a is a matrix and y is a std::vector
void dcopy(const ModuleBase::ComplexMatrix &a,int i,std::complex<double> *y);
// copy x to ith row of b where b is a matrix and x is a std::vector
void dcopy(double *x,ModuleBase::matrix &b,int i);

// compute the Euclidean length (12 norm) of std::vector x, with scaling of
// input to avoid destructive underflow and overflow
double dnrm2(const int n, const double *x, const int incx) ;
#endif // NYFUNC



