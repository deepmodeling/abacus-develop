#include "./opt_test_tools.h"
#include <math.h>

LinearEqu::LinearEqu()
{
    // initial A, x, b, p, Ap

    // A = [[2,1,0], [1,2,1], [0,1,2]]
    A = new double*[this->nx];
    for (int i = 0; i < this->nx; ++i) 
    {
        A[i] = new double[this->nx];
    }
    A[0][0] = 2; A[0][1] = 1; A[0][2] = 0;
    A[1][0] = 1; A[1][1] = 2; A[1][2] = 1;
    A[2][0] = 0; A[2][1] = 1; A[2][2] = 2;
    // A[0][0] = 1; A[0][1] = 1; A[0][2] = 0;
    // A[1][0] = 1; A[1][1] = 1; A[1][2] = 1;
    // A[2][0] = 0; A[2][1] = 1; A[2][2] = 1;

    b = new double[this->nx];
    b[0] = 1; b[1] = 2; b[2] = this->nx;
}

LinearEqu::~LinearEqu()
{
    delete[] b;
    for (int i = 0; i < this->nx; ++i)
    {
        delete[] A[i];
    }
    delete[] A;
}

void LinearEqu::get_Ap(double **A, double *p, double *Ap, int nx, int ny)
{
    for (int i = 0; i < this->nx; ++i)
    {
        Ap[i] = 0;
        for (int j = 0; j < ny; ++j)
        {
            Ap[i] += A[i][j] * p[j];
        }
    }
}

// f(x) = xAx/2 - bx
// A must be symmetrical positive definite matrix
double LinearEqu::func(double *x)
{
    double *Ax = new double[this->nx];
    this->get_Ap(A, x, Ax, this->nx, this->nx);
    double result = 0;
    for (int i = 0; i < this->nx; ++i)
    {
        result += x[i] * Ax[i] / 2 - b[i] * x[i];
    }
    delete[] Ax;
    return result;
}

// df(x)/dx = Ax - b
void LinearEqu::dfuncdx(double *x, double *gradient)
{
    double *Ax = new double[this->nx];
    this->get_Ap(A, x, Ax, this->nx, this->nx);
    for (int i = 0; i < this->nx; ++i)
    {
        gradient[i] = Ax[i] - b[i];
    }
    delete[] Ax;
}

// x = x + ap
// df(x)/da = df(x)/dx * dx/da = (Ax - b)*p
double LinearEqu::dfuncdstp(double *x, double *p)
{
    double *Ax = new double[this->nx];
    get_Ap(A, x, Ax, this->nx, this->nx);
    double result = 0;
    for (int i = 0; i < this->nx; ++i)
    {
        result += (Ax[i] - b[i]) * p[i];
    }
    delete[] Ax;
    return result;
}

namespace ModuleESolver
{
// f(x) = xAx/2 - bx
// A must be symmetrical positive definite matrix
double ESolver_OF::func(double *x)
{
    double result = 0.;
    result += pow(x[0] - x[1] - 2, 4.);
    result += pow((x[0] * x[1] - x[2] + 1), 2.);
    result += pow(x[0] - 4, 2.);
    // result += pow(x[0] - 2, 4.);
    // result += pow(x[1] - 2, 2.);
    // result += pow(x[2] - 4, 2.);
    return result;
}

// df(x)/dx = Ax - b
void ESolver_OF::dfuncdx(double *x, double *gradient)
{
    gradient[0] = 4 * pow(x[0] - x[1] - 2, 3) + 2 * (x[0] * x[1] - x[2] + 1) * x[1] + 2 * (x[0] - 4);
    gradient[1] = -4 * pow(x[0] - x[1] - 2, 3) + 2 * (x[0] * x[1] - x[2] + 1) * x[0];
    gradient[2] = -2 * (x[0] * x[1] - x[2] + 1);
    // gradient[0] = 4 * pow(x[0] - 2, 3.);
    // gradient[1] = 2 * (x[1] - 2);
    // gradient[2] = 2 * (x[2] - 4);
}

// x = x + ap
// df(x)/da = df(x)/dx * dx/da = gradient*p
double ESolver_OF::dfuncdstp(double *x, double *p)
{
    double *grad = new double[3];
    dfuncdx(x, grad);
    double result = 0;
    for (int i = 0; i < 3; ++i) result += grad[i] * p[i];
    delete[] grad;
    return result;
}
}