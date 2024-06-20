#ifndef LCAO_HS_ARRAYS_H
#define LCAO_HS_ARRAYS_H

#include <vector>
#include <complex>

class LCAO_HS_Arrays
{
    public:

    LCAO_HS_Arrays(){};
    ~LCAO_HS_Arrays(){};

    //------------------------------
    // Store H(mu,nu')
    // nu' : nu in near unitcell R.
    // used in kpoint algorithm.
    // these matrixed are used
    // for 'folding_matrix' in lcao_nnr,
    // HlocR -> Hloc2,
    // SlocR -> Sloc2,
    //------------------------------
    std::vector<double> Hloc_fixedR;
};

#endif
