#ifndef STO_CHE_H
#define STO_CHE_H
#include "module_base/math_chebyshev.h"

template <typename REAL>
class StoChe
{
  public:
    StoChe(const int& nche, const int& method);
    ~StoChe();

  public:
    int nche = 0;                                     ///< order of Chebyshev expansion
    REAL* spolyv = nullptr;                         ///< coefficients of Chebyshev expansion
    
    // Chebyshev expansion
    // It stores the plan of FFTW and should be initialized at the beginning of the calculation
    ModuleBase::Chebyshev<REAL>* p_che = nullptr;     
    
};

#endif