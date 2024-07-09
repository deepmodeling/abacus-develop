#ifndef LCAO_HS_ARRAYS_H
#define LCAO_HS_ARRAYS_H

#include "module_base/abfs-vector3_order.h"

#include <complex>
#include <vector>

class LCAO_HS_Arrays {
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

    // For HR_sparse[2], when nspin=1, only 0 is valid, when nspin=2, 0 means
    // spin up, 1 means spin down
    std::map<Abfs::Vector3_Order<int>,
             std::map<size_t, std::map<size_t, double>>>
        HR_sparse[2];
    std::map<Abfs::Vector3_Order<int>,
             std::map<size_t, std::map<size_t, double>>>
        SR_sparse;
    std::map<Abfs::Vector3_Order<int>,
             std::map<size_t, std::map<size_t, double>>>
        TR_sparse;

    std::map<Abfs::Vector3_Order<int>,
             std::map<size_t, std::map<size_t, double>>>
        dHRx_sparse[2];
    std::map<Abfs::Vector3_Order<int>,
             std::map<size_t, std::map<size_t, double>>>
        dHRy_sparse[2];
    std::map<Abfs::Vector3_Order<int>,
             std::map<size_t, std::map<size_t, double>>>
        dHRz_sparse[2];

    // For nspin = 4
    std::map<Abfs::Vector3_Order<int>,
             std::map<size_t, std::map<size_t, std::complex<double>>>>
        HR_soc_sparse;
    std::map<Abfs::Vector3_Order<int>,
             std::map<size_t, std::map<size_t, std::complex<double>>>>
        SR_soc_sparse;

    std::map<Abfs::Vector3_Order<int>,
             std::map<size_t, std::map<size_t, std::complex<double>>>>
        dHRx_soc_sparse;
    std::map<Abfs::Vector3_Order<int>,
             std::map<size_t, std::map<size_t, std::complex<double>>>>
        dHRy_soc_sparse;
    std::map<Abfs::Vector3_Order<int>,
             std::map<size_t, std::map<size_t, std::complex<double>>>>
        dHRz_soc_sparse;

    // Records the R direct coordinates of HR and SR output, This variable will
    // be filled with data when HR and SR files are output.
    std::set<Abfs::Vector3_Order<int>> output_R_coor;

    // Record all R direct coordinate information, even if HR or SR is a zero
    // matrix
    std::set<Abfs::Vector3_Order<int>> all_R_coor;
};

#endif
