#ifndef STO_CHE_H
#define STO_CHE_H
#include "module_base/math_chebyshev.h"

template <typename REAL>
class StoChe
{
  public:
    StoChe(const int& nche, const int& method, const REAL& emax_sto, const REAL& emin_sto);
    ~StoChe();

  public:
    int nche = 0;           ///< order of Chebyshev expansion
    REAL* spolyv = nullptr; ///< coefficients of Chebyshev expansion
    int method_sto = 0;     ///< method for the stochastic calculation

    // Chebyshev expansion
    // It stores the plan of FFTW and should be initialized at the beginning of the calculation
    ModuleBase::Chebyshev<REAL>* p_che = nullptr;

    REAL emax_sto = 0.0; ///< maximum energy for normalization
    REAL emin_sto = 0.0; ///< minimum energy for normalization
};

/**
 * @brief calculate v^T*M*v
 * 
 * @param v v
 * @param M M
 * @param n the dimension of v
 * @return double 
 */
double vTMv(const double* v, const double* M, const int n);

#endif