#pragma once
#include <cstddef>
#include <complex>
#include <vector>
#include <utility>
#include "source_base/matrix.h"
#include "source_base/complexmatrix.h"
#include "source_base/parallel_2d.h"
#include "source_psi/psi.h"
#include <ATen/core/tensor.h>

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

    inline double inner_product(const double* vec1, const double* vec2, const int& size)
    {
        double result = 0.0;
        for (int i = 0; i < size; ++i)
        {
            result += vec1[i] * vec2[i];
        }
        return result;
    }

    inline std::complex<double> inner_product(const std::complex<double>* vec1,
                                              const std::complex<double>* vec2,
                                              const int& size)
    {
        std::complex<double> result = 0.0;
        for (int i = 0; i < size; ++i)
        {
            result += std::conj(vec1[i]) * vec2[i];
        }
        return result;
    }

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
    template<typename T>
    bool is_hermitian(const T* mat, const Parallel_2D& pmat, const double threshold, const int my_rank);

    template<typename T>
    bool is_symmetric(const T* mat, const Parallel_2D& pmat, const double threshold, const int my_rank);

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
    
    /// @brief assign global X to 2d-matrix, its col is band(excition state), and row is { spin, k-point, occ, virt }
    /// @attention pX is 2d-blocked as {occ, virt}, this assignment is used to calculate transition density matrix c_b X_{bj} c_j
    template <typename T>
    void global2local_X(T* local_X, const T* global_X, const int nband, const int nk, 
        const std::vector<int>& nocc, const std::vector<int>& nvirt, const std::vector<Parallel_2D>& pX,
        const bool openshell);
#ifdef __MPI
    // pack the process to setup 2d divion reusing blacs_ctxt of an existing 2d-matrix
    void setup_2d_division(Parallel_2D& pv, int nb, int gr, int gc, const int& blacs_ctxt_in);

    /// @brief Struct to get MPI_traits for different data types
    template <typename T>
    struct MPIType {
        static MPI_Datatype value() { return MPI_DATATYPE_NULL; }
    };
    template <>
    struct MPIType<int> {
        static MPI_Datatype value() { return MPI_INT; }
    };
    template <>
    struct MPIType<double> {
        static MPI_Datatype value() { return MPI_DOUBLE; }
    };
    template <>
    struct MPIType<std::complex<double>> {
        static MPI_Datatype value() { return MPI_DOUBLE_COMPLEX; }
    };
    /// @brief assign X in pA to X in pX
    template <typename T>
    void pA2pX(T* X_pX, const T* X_pA, const int nband, const int nk, 
        const std::vector<int>& nocc, const std::vector<int>& nvirt,
        const std::vector<Parallel_2D>& pX, const Parallel_2D& pA,
        const int row_offset, const int col_offset, const bool openshell,
        const int my_rank, const int nproc);

    /// @brief  gather 2d matrix to full matrix
    /// the defination of row and col is consistent with setup_2d_division
    template <typename T>
    void gather_2d_to_full(const Parallel_2D& pv, const T* submat, T* fullmat,
        const bool row_major, const std::size_t global_nrow, const std::size_t global_ncol);
#endif

    ///=================diago-lapack====================
    /// @brief  diagonalize a hermitian matrix
    void diag_lapack(const int& n, double* mat, double* eig);
    void diag_lapack(const int& n, std::complex<double>* mat, double* eig);
    /// @brief  diagonalize a general matrix
    void diag_lapack_nh(const int& n, double* mat, std::complex<double>* eig);
    void diag_lapack_nh(const int& n, std::complex<double>* mat, std::complex<double>* eig);

    ///=================string option====================
    std::string tolower(const std::string& str);
    std::string toupper(const std::string& str);
}
#include "lr_util.hpp"
