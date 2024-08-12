#include "sto_che.h"
#include "module_base/blas_connector.h"

template <typename REAL>
StoChe<REAL>::~StoChe()
{
    delete p_che;
    delete[] spolyv;
}

template <typename REAL>
StoChe<REAL>::StoChe(const int& nche, const int& method, const REAL& emax_sto, const REAL& emin_sto)
{
    this->nche = nche;
    this->method_sto = method;
    p_che = new ModuleBase::Chebyshev<REAL>(nche);
    if (method == 1)
    {
        spolyv = new REAL[nche];
    }
    else
    {
        spolyv = new REAL[nche * nche];
    }

    this->emax_sto = emax_sto;
    this->emin_sto = emin_sto;
}

template class StoChe<double>;
// template class StoChe<float>;

double vTMv(const double* v, const double* M, const int n)
{
    const char normal = 'N';
    const double one = 1;
    const int inc = 1;
    const double zero = 0;
    double* y = new double[n];
    dgemv_(&normal, &n, &n, &one, M, &n, v, &inc, &zero, y, &inc);
    double result = BlasConnector::dot(n, y, 1, v, 1);
    delete[] y;
    return result;
}