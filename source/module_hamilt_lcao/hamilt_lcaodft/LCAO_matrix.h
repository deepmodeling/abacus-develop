#ifndef LCAO_MATRIX_H
#define LCAO_MATRIX_H

#include "module_base/complexmatrix.h"
#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_base/vector3.h"
#include "module_basis/module_ao/parallel_orbitals.h"
#include "module_hamilt_lcao/hamilt_lcaodft/LCAO_HS_arrays.hpp"

// add by jingan for map<> in 2021-12-2, will be deleted in the future
#include "module_base/abfs-vector3_order.h"
#ifdef __EXX
#include <RI/global/Tensor.h>
#endif

class LCAO_Matrix {
  public:
    LCAO_Matrix();
    ~LCAO_Matrix();

    // folding the fixed Hamiltonian (T+Vnl) if
    // k-point algorithm is used.
    void folding_fixedH(const int& ik,
                        const std::vector<ModuleBase::Vector3<double>>& kvec_d,
                        bool cal_syns = false);

    Parallel_Orbitals* ParaV;

#ifdef __EXX
    using TAC = std::pair<int, std::array<int, 3>>;
    std::vector<std::map<int, std::map<TAC, RI::Tensor<double>>>>* Hexxd;
    std::vector<std::map<int, std::map<TAC, RI::Tensor<std::complex<double>>>>>*
        Hexxc;
    /// @brief Hexxk for all k-points, only for the 1st scf loop ofrestart load
    std::vector<std::vector<double>> Hexxd_k_load;
    std::vector<std::vector<std::complex<double>>> Hexxc_k_load;
#endif

  public:
    //------------------------------
    // H, S, Hfixed
    // used in gamma only algorithm.
    // thse matrix are used to
    // diagonalize.
    //------------------------------
    std::vector<double> Hloc;
    std::vector<double> Sloc;
    std::vector<double> Hloc_fixed;

    //------------------------------
    // 1. Hamiltonian(vl),
    // 2. overlap matrix Sloc2
    // 3. fixed (vna+T+Vnl) matrix.
    // used in kpoint algorithm.
    // these matrix are used to
    // diagonalize.
    //------------------------------
    std::vector<std::complex<double>> Hloc2;
    std::vector<std::complex<double>> Sloc2;
    std::vector<std::complex<double>> Hloc_fixed2;

    // Record all R direct coordinate information, even if HR or SR is a zero
    // matrix
    std::set<Abfs::Vector3_Order<int>> all_R_coor;

    // Records the R direct coordinates of HR and SR output, This variable will
    // be filled with data when HR and SR files are output.
    std::set<Abfs::Vector3_Order<int>> output_R_coor;

    template <typename T>
    static void set_mat2d(const int& global_ir,
                          const int& global_ic,
                          const T& v,
                          const Parallel_Orbitals& pv,
                          T* mat);

    void set_HSgamma(const int& iw1_all,
                     const int& iw2_all,
                     const double& v,
                     double* HSloc);

    void zeros_HSgamma(const char& mtype);

    void zeros_HSk(const char& mtype);
};

#include "LCAO_matrix.hpp"

#endif
