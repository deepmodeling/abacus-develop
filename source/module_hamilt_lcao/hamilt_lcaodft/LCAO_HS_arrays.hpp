#ifndef LCAO_HS_ARRAYS_H
#define LCAO_HS_ARRAYS_H

#include "module_base/abfs-vector3_order.h"

#include <complex>
#include <vector>

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

    std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, std::complex<double>>>> dHRx_soc_sparse;
    std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, std::complex<double>>>> dHRy_soc_sparse;
    std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, std::complex<double>>>> dHRz_soc_sparse;
};

#endif
