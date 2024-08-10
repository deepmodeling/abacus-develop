#include "sto_che.h"

template <typename REAL>
StoChe<REAL>::~StoChe()
{
    delete p_che;
    delete[] spolyv;
}

template <typename REAL>
StoChe<REAL>::StoChe(const int& nche, const int& method)
{
    this->nche = nche;
    p_che = new ModuleBase::Chebyshev<REAL>(nche);
    if (method == 1)
    {
        spolyv = new REAL[nche];
    }
    else
    {
        spolyv = new REAL[nche * nche];
    }
}

template class StoChe<double>;
// template class StoChe<float>;