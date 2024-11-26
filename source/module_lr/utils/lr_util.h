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

#ifndef TO_COMPLEX_H
#define TO_COMPLEX_H
template <typename T> struct ToComplex;
template <> struct ToComplex<double> { using type = std::complex<double>; };
template <> struct ToComplex<std::complex<double>> { using type = std::complex<double>; };
template <> struct ToComplex<float> { using type = std::complex<float>; };
template <> struct ToComplex<std::complex<float>> { using type = std::complex<float>; };
#endif

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

    // Operators to calculate xc kernel have been moved into lr_util_xc.hpp.
    /// =================ALGORITHM====================

    //====== newers and deleters========
    /// @brief  delete 2d pointer  
    template <typename T>
    void _deallocate_2order_nested_ptr(T** p2, size_t size);
    /// @brief  new 2d pointer  
    template <typename T>
    void _allocate_2order_nested_ptr(T**& p2, size_t size1, size_t size2);

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
