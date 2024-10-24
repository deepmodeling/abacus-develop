#pragma once
#include <cstddef>
#include <vector>
#include <utility>
#include "module_base/matrix.h"
#include "module_base/complexmatrix.h"
#include "module_basis/module_ao/parallel_2d.h"
#include "module_psi/psi.h"
#include <ATen/core/tensor.h>
#include "module_basis/module_pw/pw_basis.h"

using DAT = container::DataType;
using DEV = container::DeviceType;

namespace LR_Util
{
    /// =====================PHYSICS====================

    /// @brief calculate the number of electrons
    /// @tparam TCell 
    /// @param ucell 
    template <typename TCell>
    int cal_nelec(const TCell& ucell);
    
    /// @brief calculate the number of occupied orbitals
    /// @param nelec 
    int cal_nocc(int nelec);
    
    /// @brief  set the index map: ix to (ic, iv) and vice versa
    /// by diagonal traverse the c-v pairs
    /// leftdown -> rightup for mode 0, rightup -> leftdown for mode 1
    /// @param mode  0: homo-1 -> lumo first; 1: homo -> lumo+1 first
    /// @param nc number of occupied bands
    /// @param nv number of virtual bands
    /// @return [iciv2ix, ix2iciv]
    std::pair<ModuleBase::matrix, std::vector<std::pair<int, int>>>
        set_ix_map_diagonal(bool mode, int nc, int nv);

#ifdef  USE_LIBXC
    /// operators to calculate XC kernels
    void grad(const double* rhor,
        ModuleBase::Vector3<double>* gdr,
        const ModulePW::PW_Basis& rho_basis,
        const double& tpiba);
    void laplace(const double* rhor,
        double* lapn,
        const ModulePW::PW_Basis& rho_basis,
        const double& tpiba2);
#endif
    /// =================ALGORITHM====================

    //====== newers and deleters========
    /// @brief  delete 2d pointer  
    template <typename T>
    void delete_p2(T** p2, size_t size);
    /// @brief  new 2d pointer  
    template <typename T>
    void new_p2(T**& p2, size_t size1, size_t size2);

    template<typename T> ct::Tensor newTensor(const ct::TensorShape& shape)
    {
        return ct::Tensor(ct::DataTypeToEnum<T>::value, ct::DeviceTypeToEnum<base_device::DEVICE_CPU>::value, shape);
    }

    ///================ BLAS ======================
    /// calculate (A+A^T)/2
    template<typename T>
    void matsym(const T* in, const int n, T* out);
    /// calculate (A+A^T)/2 (in-place version)
    template<typename T>
    void matsym(T* inout, const int n);
#ifdef __MPI
    template<typename T>
    void matsym(const T* in, const int n, const Parallel_2D& pmat, T* out);
    template<typename T>
    void matsym(T* inout, const int n, const Parallel_2D& pmat);
#endif
    ///======== Tensor - Matrix transformer==========
    container::Tensor mat2ten_double(ModuleBase::matrix& m);
    std::vector<container::Tensor> mat2ten_double(std::vector<ModuleBase::matrix>& m);
    ModuleBase::matrix ten2mat_double(container::Tensor& t);
    std::vector<ModuleBase::matrix> ten2mat_double(std::vector<container::Tensor>& t);
    container::Tensor mat2ten_complex(ModuleBase::ComplexMatrix& m);
    std::vector<container::Tensor> mat2ten_complex(std::vector<ModuleBase::ComplexMatrix>& m);
    ModuleBase::ComplexMatrix ten2mat_complex(container::Tensor& t);
    std::vector<ModuleBase::ComplexMatrix> ten2mat_complex(std::vector<container::Tensor>& t);

    ModuleBase::matrix vec2mat(const std::vector<double>& v, const int nr, const int nc);
    ModuleBase::ComplexMatrix vec2mat(const std::vector<std::complex<double>>& v, const int nr, const int nc);
    std::vector<ModuleBase::matrix> vec2mat(const std::vector<std::vector<double>>& v, const int nr, const int nc);
    std::vector<ModuleBase::ComplexMatrix> vec2mat(const std::vector<std::vector<std::complex<double>>>& v, const int nr, const int nc);

    ///===================Psi wrapper=================
    /// psi(nk=1, nbands=nb, nk * nbasis) -> psi(nb, nk, nbasis) without memory copy
    template<typename T, typename Device>
    psi::Psi<T, Device> k1_to_bfirst_wrapper(const psi::Psi<T, Device>& psi_kfirst, int nk_in, int nbasis_in);
    ///  psi(nb, nk, nbasis) -> psi(nk=1, nbands=nb, nk * nbasis)  without memory copy
    template<typename T, typename Device>
    psi::Psi<T, Device> bfirst_to_k1_wrapper(const psi::Psi<T, Device>& psi_bfirst);

    ///=================2D-block Parallel===============
    // pack the process to setup 2d divion reusing blacs_ctxt of a new 2d-matrix
    void setup_2d_division(Parallel_2D& pv, int nb, int gr, int gc);

#ifdef __MPI
    // pack the process to setup 2d divion reusing blacs_ctxt of an existing 2d-matrix
    void setup_2d_division(Parallel_2D& pv, int nb, int gr, int gc, const int& blacs_ctxt_in);
    /// @brief  gather 2d matrix to full matrix
    /// the defination of row and col is consistent with setup_2d_division
    template <typename T>
    void gather_2d_to_full(const Parallel_2D& pv, const T* submat, T* fullmat, bool col_first, int global_nrow, int global_ncol);
#endif

    ///=================diago-lapack====================
    /// @brief  diagonalize a hermitian matrix
    void diag_lapack(const int& n, double* mat, double* eig);
    void diag_lapack(const int& n, std::complex<double>* mat, double* eig);
    /// @brief  diagonalize a general matrix
    void diag_lapack_nh(const int& n, double* mat, std::complex<double>* eig);
    void diag_lapack_nh(const int& n, std::complex<double>* mat, std::complex<double>* eig);
}
#include "lr_util.hpp"
