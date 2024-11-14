#pragma once
#include <vector>
#include <functional>
#include "module_base/module_device/types.h"
#include "module_base/module_device/memory_op.h"
#include "module_hsolver/kernels/math_kernel_op.h"
namespace hsolver
{
    template <typename T>
    using Real = typename GetTypeReal<T>::type;

    /// @brief Transform a single value
    namespace fval
    {
        template <typename T> T none(const T& x) { return x; }
        template <typename T> T qe_pw(const T& x) { return 0.5 * (1.0 + x + sqrt(1 + (x - 1.0) * (x - 1.0))); }
    }

    /// @brief Transform vectors
    namespace fvec
    {
        /// @brief To be called in the iterative eigensolver.
        /// Users can add other types of operation than the following ones at one's need.
        /// fixed parameters: object vector, eigenvalue, leading dimension, number of vectors

        ///---------------------------------------------------------------------------------------------
        /// type 1: directly divide each vector by the precondition vector
        ///---------------------------------------------------------------------------------------------
        template <typename T, typename Device = base_device::DEVICE_CPU>
        void div_prevec(T* ptr, const size_t& dim, const size_t& nvec,
            const Real<T>* const pre)
        {
            for (int m = 0; m < nvec; m++)
            {
                T* const ptr_m = ptr + m * dim;
                vector_div_vector_op<T, Device>()({}, dim, ptr_m, ptr_m, pre);
            }
        }
        /// calling intereface in the eigensolver
        template <typename T>
        using Div = std::function<void(T*, const size_t&, const size_t&)>;
        // Kernel function full of dependence
        template <typename T, typename Device = base_device::DEVICE_CPU>
        using DivKernel = std::function<decltype(div_prevec<T, Device>)>;

        ///---------------------------------------------------------------------------------------------
        ///type2: Divide transval(precon_vec - eigen_subspace[m]) for each vector[m]
        ///$X \to (A-\lambda I)^{-1} X$
        ///---------------------------------------------------------------------------------------------
        template <typename T, typename Device = base_device::DEVICE_CPU>
        void div_trans_prevec_minus_eigen(T* ptr, const Real<T>* eig, const size_t& dim, const size_t& nvec,
            const Real<T>* const pre, Real<T>* const d_pre = nullptr, const std::function<Real<T>(const Real<T>&)>& transval = fval::none<Real<T>>)
        {
            using syncmem_var_h2d_op = base_device::memory::synchronize_memory_op<Real<T>, Device, base_device::DEVICE_CPU>;
            std::vector<Real<T>> pre_trans(dim, 0.0);
            const auto device = base_device::get_device_type<Device>({});

            for (int m = 0; m < nvec; m++)
            {
                T* const ptr_m = ptr + m * dim;
                for (size_t i = 0; i < dim; i++) { pre_trans[i] = transval(pre[i] - eig[m]); }
#if defined(__CUDA) || defined(__ROCM)
                if (device == base_device::GpuDevice)
                {
                    assert(d_pre);
                    syncmem_var_h2d_op()({}, {}, d_pre, pre_trans.data(), dim);
                    vector_div_vector_op<T, Device>()({}, dim, ptr_m, ptr_m, d_pre);
                }
                else
#endif
                {
                    vector_div_vector_op<T, Device>()({}, dim, ptr_m, ptr_m, pre_trans.data());
                }
            }
        }
        /// calling intereface in the eigensolver
        template <typename T>
        using DivTransMinusEig = std::function<void(T*, const Real<T>*, const size_t&, const size_t&)>;
        // Kernel function full of dependence
        template <typename T, typename Device = base_device::DEVICE_CPU>
        using DivTransMinusEigKernel = std::function<decltype(div_trans_prevec_minus_eigen<T, Device>)>;
    }

    /// @brief A operator-like class of precondition function
    ///     to encapsulate the pre-allocation of memory on different devices before starting the iterative eigensolver.
    ///     One can pass the operatr() function of this class, or other custom lambdas/functions to eigensolvers.
    template <typename T, typename Device = base_device::DEVICE_CPU, typename Kernel_t = fvec::DivKernel<T, Device>>
    struct PreOP
    {
        using FVal_t = std::function<Real<T>(const Real<T>&)>;  //single-value transformer
        using resmem_real_op = base_device::memory::resize_memory_op<Real<T>, Device>;
        using delmem_real_op = base_device::memory::delete_memory_op<Real<T>, Device>;
        PreOP(const std::vector<Real<T>>& prevec,
            const Kernel_t& transvec = fvec::div_prevec<T, Device>,
            const FVal_t& transval = fval::none<Real<T>>)
            : PreOP<T, Device, Kernel_t>(prevec.data(), prevec.size(), transvec, transval) {}
        PreOP(const Real<T>* const prevec, const int& dim,
            const Kernel_t& transvec = fvec::div_prevec<T, Device>,
            const FVal_t& transval = fval::none<Real<T>>)
            : prevec_(prevec), dim_(dim), transvec_(transvec), transval_(transval),
            dev_(base_device::get_device_type<Device>({}))
        {
#if defined(__CUDA) || defined(__ROCM)
            if (this->dev_ == base_device::GpuDevice) { resmem_real_op()({}, this->d_prevec_, dim_); }
#endif
        }
        PreOP(const PreOP&) = delete;
        ~PreOP() {
#if defined(__CUDA) || defined(__ROCM)
            if (this->dev_ == base_device::GpuDevice) { delmem_real_op()({}, this->d_prevec_); }
#endif
        }

        /// @brief Bind a PreOP object to a function
        template<typename U = Kernel_t, typename std::enable_if<std::is_same<U, fvec::DivKernel<T, Device>>::value, bool>::type = 0>
        fvec::Div<T> get() const
        {
            return std::bind(PreOP<T, Device, Kernel_t>::transvec_,
                std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, this->prevec_);
        }
        template<typename U = Kernel_t, typename std::enable_if<std::is_same<U, fvec::DivTransMinusEigKernel<T, Device>>::value, bool>::type = 0>
        fvec::DivTransMinusEig<T> get() const
        {
            return std::bind(PreOP<T, Device, Kernel_t>::transvec_,
                std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4,
                this->prevec_, this->d_prevec_, this->transval_);
        }

    private:
        const Real<T>* const prevec_ = nullptr;
        const int dim_;
        Real<T>* d_prevec_ = nullptr;
        const Kernel_t transvec_;
        const FVal_t transval_;
        const base_device::AbacusDevice_t dev_;
    };
}
