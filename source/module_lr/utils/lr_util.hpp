#pragma once
#include <cstddef>
#include "lr_util.h"
#include <algorithm>
#include "module_cell/unitcell.h"
#include "module_base/constants.h"
#include "module_hamilt_general/module_xc/xc_functional.h"
namespace LR_Util
{
    /// =================PHYSICS====================

    template <typename TCell>
    int cal_nelec(const TCell& ucell) {
        int nelec = 0;
        for (int it = 0; it < ucell.ntype; ++it) {
            nelec += ucell.atoms[it].ncpp.zv * ucell.atoms[it].na;
}
        return nelec;
    }

    /// =================ALGORITHM====================
    ///for lack of make_unique in c++11 
    template<typename T, typename... Args>
    std::unique_ptr<T> make_unique(Args &&... args)
    {
        return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
    }

    //====== newers and deleters========
    /// @brief  new 2d pointer
    /// @tparam T
    /// @param size1
    /// @param size2
    template <typename T>
    void _allocate_2order_nested_ptr(T**& p2, size_t size1, size_t size2)
    {
        p2 = new T * [size1];
        for (size_t i = 0; i < size1; ++i)
        {
            p2[i] = new T[size2];
        }
    };

    /// @brief  delete 2d pointer 
    /// @tparam T 
    /// @param p2 
    /// @param size 
    template <typename T>
    void _deallocate_2order_nested_ptr(T** p2, size_t size)
    {
        if (p2 != nullptr)
        {
            for (size_t i = 0; i < size; ++i)
            {
                if (p2[i] != nullptr) { delete[] p2[i]; }
            }
            delete[] p2;
        }
    };

    inline double get_conj(const double& x)
    {
        return x;
    }
    inline std::complex<double> get_conj(const std::complex<double>& x)
    {
        return std::conj(x);
    }
    template <typename T>
    void matsym(const T* in, const int n, T* out)
    {
        for (int i = 0; i < n; ++i) {
            out[i * n + i] = 0.5 * in[i * n + i] + 0.5 * get_conj(in[i * n + i]);
}
        for (int i = 0;i < n;++i) {
            for (int j = i + 1;j < n;++j)
            {
                out[i * n + j] = 0.5 * (in[i * n + j] + get_conj(in[j * n + i]));
                out[j * n + i] = get_conj(out[i * n + j]);
            }
}
    }
    template <typename T>
    void matsym(T* inout, const int n)
    {
        for (int i = 0; i < n; ++i) {
            inout[i * n + i] = 0.5 * (inout[i * n + i] + get_conj(inout[i * n + i]));
}
        for (int i = 0;i < n;++i) {
            for (int j = i + 1;j < n;++j)
            {
                inout[i * n + j] = 0.5 * (inout[i * n + j] + get_conj(inout[j * n + i]));
                inout[j * n + i] = get_conj(inout[i * n + j]);
            }
}
    }

    /// get the Psi wrapper of the selected spin from the Psi object
    template<typename T>
    psi::Psi<T> get_psi_spin(const psi::Psi<T>& psi_in, const int& is, const int& nk)
    {
        return psi::Psi<T>(&psi_in(is * nk, 0, 0), psi_in, nk, psi_in.get_nbands());
    }

    /// psi(nk=1, nbands=nb, nk * nbasis) -> psi(nb, nk, nbasis) without memory copy
    template<typename T, typename Device>
    psi::Psi<T, Device> k1_to_bfirst_wrapper(const psi::Psi<T, Device>& psi_kfirst, int nk_in, int nbasis_in)
    {
        assert(psi_kfirst.get_nk() == 1);
        assert(nk_in * nbasis_in == psi_kfirst.get_nbasis());
        int ib_now = psi_kfirst.get_current_b();
        psi_kfirst.fix_b(0);    // for get_pointer() to get the head pointer
        psi::Psi<T, Device> psi_bfirst(psi_kfirst.get_pointer(), nk_in, psi_kfirst.get_nbands(), nbasis_in, psi_kfirst.get_ngk_pointer(), false);
        psi_kfirst.fix_b(ib_now);
        return psi_bfirst;
    }

    ///  psi(nb, nk, nbasis) -> psi(nk=1, nbands=nb, nk * nbasis)  without memory copy
    template<typename T, typename Device>
    psi::Psi<T, Device> bfirst_to_k1_wrapper(const psi::Psi<T, Device>& psi_bfirst)
    {
        int ib_now = psi_bfirst.get_current_b();
        int ik_now = psi_bfirst.get_current_k();
        psi_bfirst.fix_kb(0, 0);    // for get_pointer() to get the head pointer
        psi::Psi<T, Device> psi_kfirst(psi_bfirst.get_pointer(), 1, psi_bfirst.get_nbands(), psi_bfirst.get_nk() * psi_bfirst.get_nbasis(), psi_bfirst.get_ngk_pointer(), true);
        psi_bfirst.fix_kb(ik_now, ib_now);
        return psi_kfirst;
    }

#ifdef __MPI
    template <typename T>
    void gather_2d_to_full(const Parallel_2D& pv, const T* submat, T* fullmat, bool col_first, int global_nrow, int global_ncol)
    {
        ModuleBase::TITLE("LR_Util", "gather_2d_to_full");
        auto get_mpi_datatype = []() -> MPI_Datatype {
            if (std::is_same<T, int>::value) { return MPI_INT; }
            if (std::is_same<T, float>::value) { return MPI_FLOAT; }
            else if (std::is_same<T, double>::value) { return MPI_DOUBLE; }
            if (std::is_same<T, std::complex<float>>::value) { return MPI_COMPLEX; }
            else if (std::is_same<T, std::complex<double>>::value) { return MPI_DOUBLE_COMPLEX; }
            else { throw std::runtime_error("gather_2d_to_full: unsupported type"); }
            };

        // zeros
        for (int i = 0;i < global_nrow * global_ncol;++i) { fullmat[i] = 0.0;
}
        //copy
        for (int i = 0;i < pv.get_row_size();++i) {
            for (int j = 0;j < pv.get_col_size();++j) {
                if (col_first) {
                    fullmat[pv.local2global_row(i) * global_ncol + pv.local2global_col(j)] = submat[i * pv.get_col_size() + j];
                } else {
                    fullmat[pv.local2global_col(j) * global_nrow + pv.local2global_row(i)] = submat[j * pv.get_row_size() + i];
}
}
}

        //reduce to root
        MPI_Allreduce(MPI_IN_PLACE, fullmat, global_nrow * global_ncol, get_mpi_datatype(), MPI_SUM, pv.comm());
    };
#endif

}
