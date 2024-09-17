#ifdef __EXX
#pragma once
#include "module_cell/unitcell.h"
#include "module_psi/psi.h"

#include <RI/global/Tensor.h>
namespace RI_Benchmark
{
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
        const int& occ_first=false,
        const bool& read_from_aims=false,
        const std::vector<int>& aims_nbasis={});

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

    template<typename FPTYPE>
    std::vector<FPTYPE> read_bands(const std::string& file, const int nocc, const int nvirt, int& ncore);
    template <typename TK>
    void read_aims_eigenvectors(psi::Psi<TK>& wfc_ks, const std::string& file, const int ncore, const int nbands, const int nbasis);
    /// only for blocking by atom pairs
    template <typename TR>
    TLRI<TR> read_coulomb_mat(const std::string& file, const TLRI<TR>& Cs);
    /// for any way of blocking
    template <typename TR>
    TLRI<TR> read_coulomb_mat_general(const std::string& file, const TLRI<TR>& Cs);
    template <typename TR>
    bool compare_Vs(const TLRI<TR>& Vs1, const TLRI<TR>& Vs2, const double thr = 1e-4);
    template <typename TR>
    std::vector<TLRI<TR>> split_Ds(const std::vector<std::vector<TR>>& Ds, const std::vector<int>& aims_nbasis, const UnitCell& ucell);
}
#include "ri_benchmark.hpp"
#endif