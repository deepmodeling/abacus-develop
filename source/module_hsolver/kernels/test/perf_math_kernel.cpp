#include "module_base/blas_connector.h"
#include "module_base/constants.h"
#include "module_base/module_device/memory_op.h"
#include "module_hsolver/kernels/math_kernel_op.h"

#include <complex>
#include <benchmark/benchmark.h>
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <chrono>

/************************************************
*  performace test of class math_kernel_op
***********************************************/

/**
 * Tested function: 
 *      - zdot_real_cpu_op
 *      - vector_div_constant_op_cpu
 *      - vector_mul_vector_op_cpu
 *      - vector_div_vector_op_cpu
 *      - constantvector_addORsub_constantVector_op_cpu
 *      - axpy_cpu
 *      - scal_cpu

 *      - zdot_real_gpu_op
 *      - vector_div_constant_op_gpu
 *      - vector_mul_vector_op_gpu
 *      - vector_div_vector_op_gpu
 *      - constantvector_addORsub_constantVector_op_gpu
 *      - axpy_gpu
 *      - scal_gpu
 */

class PerfModuleHsolverMathKernel : public benchmark::Fixture {
    public:

    // DEVICE SYMBOL
    const base_device::DEVICE_CPU* cpu_ctx = {};

    int dim_vector = 1;

    std::complex<double>* test_zvector_a = nullptr;
    std::complex<double>* test_zvector_b = nullptr;
    std::complex<double>* result_zvector = nullptr;

    double* test_dvector_a = nullptr;

    double dconstant_a = 1.0;
    double dconstant_b = 1.0;

    std::complex<double> zconstant_a = {1.0,1.0};

 #if __CUDA || __UT_USE_CUDA || __ROCM || __UT_USE_ROCM
    const base_device::DEVICE_GPU * gpu_ctx = {};

    // from CPU to GPU
    using synchronize_memory_op
        = base_device::memory::synchronize_memory_op<std::complex<double>, base_device::DEVICE_GPU, base_device::DEVICE_CPU>;

    // form GPU to CPU
    using synchronize_memory_op_gpu
        = base_device::memory::synchronize_memory_op<std::complex<double>, base_device::DEVICE_CPU, base_device::DEVICE_GPU>;

    using resize_memory_op = base_device::memory::resize_memory_op<std::complex<double>, base_device::DEVICE_GPU>;
    using delete_memory_op = base_device::memory::delete_memory_op<std::complex<double>, base_device::DEVICE_GPU>;
    using resize_memory_op_double = base_device::memory::resize_memory_op<double, base_device::DEVICE_GPU>;
    using delete_memory_op_double = base_device::memory::delete_memory_op<double, base_device::DEVICE_GPU>;
    using synchronize_memory_op_double = base_device::memory::synchronize_memory_op<double, base_device::DEVICE_GPU, base_device::DEVICE_CPU>;

    using set_memory_op = base_device::memory::set_memory_op<std::complex<double>, base_device::DEVICE_GPU>;
    using set_memory_op_double = base_device::memory::set_memory_op<double, base_device::DEVICE_GPU>;

    std::complex<double>* test_zvector_a_gpu = nullptr;
    std::complex<double>* test_zvector_b_gpu = nullptr;
    std::complex<double>* result_zvector_gpu = nullptr;

    double* test_dvector_a_gpu = nullptr;

 #endif // __CUDA || __UT_USE_CUDA || __ROCM || __UT_USE_ROCM

    void SetUp(const benchmark::State& state){
        dim_vector = state.range(0); // Generate vectors with different diminsions (1,10,100,...,1e6)

        // This should be complex vectors
        test_zvector_a = new std::complex<double>[dim_vector + 1];
        test_zvector_b = new std::complex<double>[dim_vector + 1];
        result_zvector = new std::complex<double>[dim_vector + 1];

        // The following is double vectors
        test_dvector_a = new double[dim_vector + 1];


        for (int i=0;i<dim_vector;i++){
            // Generate vectors using random
            test_zvector_a[i] = std::complex<double>{(double)rand()+(double)rand()/(RAND_MAX+1.0),(double)rand()+(double)rand()/(RAND_MAX+1.0)};
            test_zvector_b[i] = std::complex<double>{(double)rand()+(double)rand()/(RAND_MAX+1.0),(double)rand()+(double)rand()/(RAND_MAX+1.0)};
            test_dvector_a[i] = (double)rand()+(double)rand()/(RAND_MAX+1.0);
        }

        dconstant_a = (double)rand()+(double)rand()/(RAND_MAX+1.0);
        dconstant_b = (double)rand()+(double)rand()/(RAND_MAX+1.0);

        zconstant_a = std::complex<double>{(double)rand()+(double)rand()/(RAND_MAX+1.0),(double)rand()+(double)rand()/(RAND_MAX+1.0)};
#if __CUDA || __UT_USE_CUDA || __ROCM || __UT_USE_ROCM

        resize_memory_op()(test_zvector_a_gpu, dim_vector);
        resize_memory_op()(test_zvector_b_gpu, dim_vector);
        synchronize_memory_op()(test_zvector_a_gpu, test_zvector_a, dim_vector);
        synchronize_memory_op()(test_zvector_b_gpu, test_zvector_b, dim_vector);

        resize_memory_op()(result_zvector_gpu, dim_vector);
        resize_memory_op_double()(test_dvector_a_gpu, dim_vector);
        synchronize_memory_op_double()(test_dvector_a_gpu, test_dvector_a, dim_vector);

        hsolver::createGpuBlasHandle();


#endif // __CUDA || __UT_USE_CUDA || __ROCM || __UT_USE_ROCM
    }
    void TearDown(const benchmark::State& state){
        delete[] test_zvector_a;
        delete[] test_zvector_b;
        delete[] result_zvector;
        delete[] test_dvector_a;
#if __CUDA || __UT_USE_CUDA || __ROCM || __UT_USE_ROCM
        hsolver::destoryBLAShandle();
#endif // __CUDA || __UT_USE_CUDA || __ROCM || __UT_USE_ROCM
    }


    // OPs need benchmark
    // CPU operator
    using zdot_real_cpu_op = hsolver::dot_real_op<std::complex<double>, base_device::DEVICE_CPU>;
    
    using vector_div_constant_op_cpu = hsolver::vector_div_constant_op<std::complex<double>, base_device::DEVICE_CPU>;
    using vector_mul_vector_op_cpu = hsolver::vector_mul_vector_op<std::complex<double>, base_device::DEVICE_CPU>;
    using vector_div_vector_op_cpu = hsolver::vector_div_vector_op<std::complex<double>, base_device::DEVICE_CPU>;
    using constantvector_addORsub_constantVector_op_cpu
        = hsolver::constantvector_addORsub_constantVector_op<std::complex<double>, base_device::DEVICE_CPU>;
    using axpy_op_cpu = hsolver::axpy_op<std::complex<double>, base_device::DEVICE_CPU>;
    using scal_op_cpu = hsolver::scal_op<double, base_device::DEVICE_CPU>;
    using gemv_op_cpu = hsolver::gemv_op<std::complex<double>, base_device::DEVICE_CPU>;

#if __CUDA || __UT_USE_CUDA || __ROCM || __UT_USE_ROCM

    // GPU operator
    using zdot_real_gpu_op = hsolver::dot_real_op<std::complex<double>, base_device::DEVICE_GPU>;

    using vector_div_constant_op_gpu = hsolver::vector_div_constant_op<std::complex<double>, base_device::DEVICE_GPU>;
    using vector_mul_vector_op_gpu = hsolver::vector_mul_vector_op<std::complex<double>, base_device::DEVICE_GPU>;
    using vector_div_vector_op_gpu = hsolver::vector_div_vector_op<std::complex<double>, base_device::DEVICE_GPU>;
    using constantvector_addORsub_constantVector_op_gpu
        = hsolver::constantvector_addORsub_constantVector_op<std::complex<double>, base_device::DEVICE_GPU>;
    using axpy_op_gpu = hsolver::axpy_op<std::complex<double>, base_device::DEVICE_GPU>;
    using scal_op_gpu = hsolver::scal_op<double, base_device::DEVICE_GPU>;

#endif // __CUDA || __UT_USE_CUDA || __ROCM || __UT_USE_ROCM
};


BENCHMARK_DEFINE_F(PerfModuleHsolverMathKernel, BM_zdot_real_cpu_op)(benchmark::State& state) {
    for (auto _ : state) {
        double result = zdot_real_cpu_op()(cpu_ctx, dim_vector, test_zvector_a, test_zvector_b, false);
    }
}


BENCHMARK_DEFINE_F(PerfModuleHsolverMathKernel, BM_vector_div_constant_op_cpu)(benchmark::State& state) {
    for (auto _ : state) {
        vector_div_constant_op_cpu()(cpu_ctx, dim_vector, result_zvector, test_zvector_a, dconstant_a);
    }
}

BENCHMARK_DEFINE_F(PerfModuleHsolverMathKernel, BM_vector_mul_vector_op_cpu)(benchmark::State& state) {
    for (auto _ : state) {
        vector_mul_vector_op_cpu()(cpu_ctx, dim_vector, result_zvector, test_zvector_a, test_dvector_a);
    }
}

BENCHMARK_DEFINE_F(PerfModuleHsolverMathKernel, BM_vector_div_vector_op_cpu)(benchmark::State& state) {
    for (auto _ : state) {
        vector_div_vector_op_cpu()(cpu_ctx, dim_vector, result_zvector, test_zvector_a, test_dvector_a);
    }
}

BENCHMARK_DEFINE_F(PerfModuleHsolverMathKernel, BM_constantvector_addORsub_constantVector_op_cpu)(benchmark::State& state) {
    for (auto _ : state) {
        constantvector_addORsub_constantVector_op_cpu()(cpu_ctx, dim_vector, result_zvector, test_zvector_a, dconstant_a ,test_zvector_b, dconstant_b);
    }
}

BENCHMARK_DEFINE_F(PerfModuleHsolverMathKernel, BM_axpy_op_cpu)(benchmark::State& state) {
    for (auto _ : state) {
        axpy_op_cpu()(cpu_ctx, dim_vector, &zconstant_a, test_zvector_a, 1 ,test_zvector_b, 1);
    }
}

BENCHMARK_DEFINE_F(PerfModuleHsolverMathKernel, BM_scal_op_cpu)(benchmark::State& state) {
    for (auto _ : state) {
        scal_op_cpu()(cpu_ctx, dim_vector, &zconstant_a, test_zvector_a, 1);
    }
}


BENCHMARK_REGISTER_F(PerfModuleHsolverMathKernel, BM_zdot_real_cpu_op)->RangeMultiplier(10)->Range(1,10e6)->Unit(benchmark::kMicrosecond);
BENCHMARK_REGISTER_F(PerfModuleHsolverMathKernel, BM_vector_div_constant_op_cpu)->RangeMultiplier(10)->Range(1,10e6)->Unit(benchmark::kMicrosecond);
BENCHMARK_REGISTER_F(PerfModuleHsolverMathKernel, BM_vector_mul_vector_op_cpu)->RangeMultiplier(10)->Range(1,10e6)->Unit(benchmark::kMicrosecond);
BENCHMARK_REGISTER_F(PerfModuleHsolverMathKernel, BM_vector_div_vector_op_cpu)->RangeMultiplier(10)->Range(1,10e6)->Unit(benchmark::kMicrosecond);
BENCHMARK_REGISTER_F(PerfModuleHsolverMathKernel, BM_constantvector_addORsub_constantVector_op_cpu)->RangeMultiplier(10)->Range(1,10e6)->Unit(benchmark::kMicrosecond);
BENCHMARK_REGISTER_F(PerfModuleHsolverMathKernel, BM_axpy_op_cpu)->RangeMultiplier(10)->Range(1,10e6)->Unit(benchmark::kMicrosecond);
BENCHMARK_REGISTER_F(PerfModuleHsolverMathKernel, BM_scal_op_cpu)->RangeMultiplier(10)->Range(1,10e6)->Unit(benchmark::kMicrosecond);


#if __CUDA || __UT_USE_CUDA || __ROCM || __UT_USE_ROCM

// If you want to use manual timer, you can refer to this.
/*
BENCHMARK_DEFINE_F(PerfModuleHsolverMathKernel, BM_zdot_real_gpu_op)(benchmark::State& state) {
    for (auto _ : state) {
        auto start = std::chrono::high_resolution_clock::now();

        double result = zdot_real_gpu_op()(gpu_ctx, dim_vector, test_zvector_a_gpu, test_zvector_b_gpu, false);

        auto end = std::chrono::high_resolution_clock::now();
        auto elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
        state.SetIterationTime(elapsed_seconds.count());
    }
}
*/

BENCHMARK_DEFINE_F(PerfModuleHsolverMathKernel, BM_zdot_real_gpu_op)(benchmark::State& state) {
    for (auto _ : state) {
        double result = zdot_real_gpu_op()(gpu_ctx, dim_vector, test_zvector_a_gpu, test_zvector_b_gpu, false);
    }
}

BENCHMARK_DEFINE_F(PerfModuleHsolverMathKernel, BM_vector_div_constant_op_gpu)(benchmark::State& state) {
    for (auto _ : state) {
        vector_div_constant_op_gpu()(gpu_ctx, dim_vector, result_zvector_gpu, test_zvector_a_gpu, dconstant_a);
    }
}

BENCHMARK_DEFINE_F(PerfModuleHsolverMathKernel, BM_vector_mul_vector_op_gpu)(benchmark::State& state) {
    for (auto _ : state) {
        vector_mul_vector_op_gpu()(gpu_ctx, dim_vector, result_zvector_gpu, test_zvector_a_gpu, test_dvector_a_gpu);
    }
}

BENCHMARK_DEFINE_F(PerfModuleHsolverMathKernel, BM_vector_div_vector_op_gpu)(benchmark::State& state) {
    for (auto _ : state) {
        vector_div_vector_op_gpu()(gpu_ctx, dim_vector, result_zvector_gpu, test_zvector_a_gpu, test_dvector_a_gpu);
    }
}

BENCHMARK_DEFINE_F(PerfModuleHsolverMathKernel, BM_constantvector_addORsub_constantVector_op_gpu)(benchmark::State& state) {
    for (auto _ : state) {
        constantvector_addORsub_constantVector_op_gpu()(gpu_ctx, dim_vector, result_zvector_gpu, test_zvector_a_gpu, dconstant_a ,test_zvector_b_gpu, dconstant_b);
    }
}

BENCHMARK_DEFINE_F(PerfModuleHsolverMathKernel, BM_axpy_op_gpu)(benchmark::State& state) {
    for (auto _ : state) {
        axpy_op_gpu()(gpu_ctx, dim_vector, &zconstant_a, test_zvector_a_gpu, 1 ,test_zvector_b_gpu, 1);
    }
}

BENCHMARK_DEFINE_F(PerfModuleHsolverMathKernel, BM_scal_op_gpu)(benchmark::State& state) {
    for (auto _ : state) {
        scal_op_gpu()(gpu_ctx, dim_vector, &zconstant_a, test_zvector_a_gpu, 1);
    }
}

// If you want to use manual timer, you can refer to this.
// BENCHMARK_REGISTER_F(PerfModuleHsolverMathKernel, BM_zdot_real_gpu_op)->RangeMultiplier(10)->Range(1,10e6)->UseManualTime()->Unit(benchmark::kMicrosecond);

BENCHMARK_REGISTER_F(PerfModuleHsolverMathKernel, BM_zdot_real_gpu_op)->RangeMultiplier(10)->Range(1,10e6)->Unit(benchmark::kMicrosecond);
BENCHMARK_REGISTER_F(PerfModuleHsolverMathKernel, BM_vector_div_constant_op_gpu)->RangeMultiplier(10)->Range(1,10e6)->Unit(benchmark::kMicrosecond);
BENCHMARK_REGISTER_F(PerfModuleHsolverMathKernel, BM_vector_mul_vector_op_gpu)->RangeMultiplier(10)->Range(1,10e6)->Unit(benchmark::kMicrosecond);
BENCHMARK_REGISTER_F(PerfModuleHsolverMathKernel, BM_vector_div_vector_op_gpu)->RangeMultiplier(10)->Range(1,10e6)->Unit(benchmark::kMicrosecond);
BENCHMARK_REGISTER_F(PerfModuleHsolverMathKernel, BM_constantvector_addORsub_constantVector_op_gpu)->RangeMultiplier(10)->Range(1,10e6)->Unit(benchmark::kMicrosecond);
BENCHMARK_REGISTER_F(PerfModuleHsolverMathKernel, BM_axpy_op_gpu)->RangeMultiplier(10)->Range(1,10e6)->Unit(benchmark::kMicrosecond);
BENCHMARK_REGISTER_F(PerfModuleHsolverMathKernel, BM_scal_op_gpu)->RangeMultiplier(10)->Range(1,10e6)->Unit(benchmark::kMicrosecond);

#endif // __CUDA || __UT_USE_CUDA || __ROCM || __UT_USE_ROCM


BENCHMARK_MAIN(); 