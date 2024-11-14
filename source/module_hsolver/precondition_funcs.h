#include <functional>
#include "module_base/module_device/types.h"
#include "module_hsolver/kernels/math_kernel_op.h"
#include <iostream> // for debugging
#include <vector>
namespace hsolver
{
    /// @brief Transforming a single value, 
    namespace transfunc
    {
        template <typename T> T none(const T& x) { return x; }
        template <typename T> T qe_pw(const T& x) { return 0.5 * (1.0 + x + sqrt(1 + (x - 1.0) * (x - 1.0))); }
    }

    template <typename T>
    using Real = typename GetTypeReal<T>::type;

    /// @brief to be called in the iterative eigensolver.
    /// fixed parameters: object vector, eigenvalue, leading dimension, number of vectors
    template <typename T>
    using PreFunc = const std::function<void(T*, const Real<T>*, const size_t&, const size_t&)>;
    // using PreFunc = std::function<void(T*, const Real<T>*, const int&, const int&)>;

    /// type1: Divide transfunc(precon_vec - eigen_subspace[m]) for each vector[m]
    ///$X \to (A-\lambda I)^{-1} X$
    // There may be other types of operation than this one.
    template <typename T, typename Device = base_device::DEVICE_CPU>
    void div_trans_prevec_minus_eigen(T* ptr, const Real<T>* eig, const size_t& dim, const size_t& nvec,
        const Real<T>* const pre, Real<T>* const d_pre = nullptr, const std::function<Real<T>(const Real<T>&)>& transfunc = transfunc::none<Real<T>>)
    {
        using syncmem_var_h2d_op = base_device::memory::synchronize_memory_op<Real<T>, Device, base_device::DEVICE_CPU>;
        std::vector<Real<T>> pre_trans(dim, 0.0);
        const auto device = base_device::get_device_type<Device>({});

        for (int m = 0; m < nvec; m++)
        {
            T* const ptr_m = ptr + m * dim;
            for (size_t i = 0; i < dim; i++) { pre_trans[i] = transfunc(pre[i] - eig[m]); }
            std::cout << std::endl;
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

    /// @brief A operator-like class of precondition function
    ///     to encapsulate the pre-allocation of memory on different devices before starting the iterative eigensolver.
    ///     One can pass the operatr() function of this class, or other custom lambdas/functions to eigensolvers.
    template <typename T, typename Device = base_device::DEVICE_CPU>
    struct PreOP
    {
        PreOP(const std::vector<Real<T>>& prevec, const std::function<Real<T>(const Real<T>&)>& transfunc = transfunc::none)
            : PreOP<T, Device>(prevec.data(), prevec.size(), transfunc) {}
        PreOP(const Real<T>* const prevec, const int& dim, const std::function<Real<T>(const Real<T>&)>& transfunc = transfunc::none)
            : prevec_(prevec), dim_(dim), transfunc_(transfunc),
            dev_(base_device::get_device_type<Device>({}))
        {
#if defined(__CUDA) || defined(__ROCM)
            if (this->dev_ == base_device::GpuDevice) { resmem_real_op<T, Device>()({}, this->d_prevec_, dim_); }
#endif
        }
        PreOP(const PreOP& other) = delete;
        ~PreOP() {
#if defined(__CUDA) || defined(__ROCM)
            if (this->dev_ == base_device::GpuDevice) { delmem_real_op<T, Device>()({}, this->d_precondition); }
#endif
        }
        void operator()(T* ptr, const Real<T>* eig, const size_t& dim, const size_t& nvec) const
        {
            assert(dim <= dim_);
            div_trans_prevec_minus_eigen<T, Device>(ptr, eig, dim, nvec, prevec_, d_prevec_, transfunc_);
        }
    private:
        const Real<T>* const prevec_;
        const int dim_;
        Real<T>* d_prevec_;
        const std::function<Real<T>(const Real<T>&)> transfunc_;
        const base_device::AbacusDevice_t dev_;
    };

    /// @brief Bind a PreOP object to a function
    template <typename T, typename Device>
    PreFunc<T> bind_pre_op(const PreOP<T, Device>& pre_op)
    {
        return std::bind(&PreOP<T, Device>::operator(), &pre_op, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
    }
}
