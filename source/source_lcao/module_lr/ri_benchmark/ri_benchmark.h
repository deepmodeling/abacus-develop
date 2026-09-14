#ifdef __EXX
#pragma once
#include "source_cell/unitcell.h"
#include "source_psi/psi.h"
#include "source_lcao/module_lr/utils/lr_io.h"
#include "source_basis/module_ao/parallel_orbitals.h"
#include <RI/global/Tensor.h>
namespace RI_Benchmark
{
    using TA = int;
    using TC = std::array<int, 3>;
    using TAC = std::pair<int, TC>;

    template <typename T>
    using TLRI = std::map<int, std::map<TAC, RI::Tensor<T>>>;
    template<typename T>
    using TLRIX = std::map<int, std::map<TAC, std::vector<T>>>;

    template <typename TK, typename TR>
    void benchmark_driver_A(std::string& file_Cs, std::string& file_Vs, std::string& file_kswfc, const int nocc, const int nvirt);
    template <typename TK, typename TR>
    void benchmark_driver_AX(std::string& file_Cs, std::string& file_Vs, std::string& file_kswfc, const int nocc, const int nvirt);

    // 1. full-matrix version

    /// calculate $Cs_{Ua\alpha,Vi}^{mo} = \sum_{\mu\nu} C_{U\mu\alpha,V\nu}^{ao} c^*_{\mu a} c_{\nu i}$
    /// if occ_first, calculate $Cs_{Ui\alpha,Va}^{mo} = \sum_{\mu\nu} C_{U\mu\alpha,V\nu}^{ao} c^*_{\mu i} c_{\nu a}$
    template <typename TK, typename TR>
    TLRI<TK> cal_Cs_mo(const UnitCell& ucell,
        const TLRI<TR>& Cs_ao,
        const psi::Psi<TK>& wfc_ks,
        const int& nocc,
        const int& nvirt,
        const int& occ_first=false);

    /// A=CVC, sum over atom quads
    template <typename TK, typename TR>
    std::vector<TK> cal_Amat_full(const TLRI<TK>& Cs_a,
        const TLRI<TK>& Cs_b,
        const TLRI<TR>& Vs);

    // 2. AX version
    template <typename TK>
    TLRIX<TK> cal_CsX(const TLRI<TK>& Cs_mo, TK* X);

    template <typename TK, typename TR>
    TLRI<TK> cal_CV(const TLRI<TK>& Cs_a,
        const TLRI<TR>& Vs);

    /// AX=CV(CX), sum over atom quads
    template <typename TK, typename TR>
    void cal_AX(const TLRI<TK>& Cs_a,
        const TLRIX<TK>& Cs_bX,
        const TLRI<TR>& Vs,
        TK* AX,
        const double& scale = 2.0);
    /// AX=（CV)(CX), sum over atom quads
    template <typename TK>
    void cal_AX(const TLRI<TK>& CV,
        const TLRIX<TK>& Cs_bX,
        TK* AX,
        const double& scale = 2.0);

    // 3. read/write tools    
    template<typename FPTYPE>
    std::vector<FPTYPE> read_aims_ebands(const std::string& file, const int nocc, const int nvirt, int& ncore);

    template <typename TK>
    void read_aims_eigenvectors(psi::Psi<TK>& wfc_ks, const std::string& file, const int ncore, const int nbands, const int nbasis);

    template <typename TR>
    bool compare_Vs(const TLRI<TR>& Vs1, const TLRI<TR>& Vs2, const double thr = 1e-4);
    template <typename TR>
    std::vector<TLRI<TR>> split_Ds(const std::vector<std::vector<TR>>& Ds, const std::vector<int>& aims_nbasis, const UnitCell& ucell);
}
#include "ri_benchmark.hpp"
#endif