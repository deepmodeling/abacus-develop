#ifndef __PWTEST
#define __PWTEST
#include "gtest/gtest.h"
#include <iostream>
#include "module_base/module_device/memory_op.h"
#include "module_basis/module_pw/kernels/pw_op.h"
using namespace std;
extern int nproc_in_pool, rank_in_pool;
extern string precision_flag, device_flag;

class PWTEST: public testing::Test
{
public:
    static void SetUpTestCase()
    {
        if(rank_in_pool == 0)
        {
            cout<<"\033[32m"<<"============================"<<"\033[0m"<<endl;
            cout<<"\033[32m"<<"=    PW GPU MODULE TEST START  ="<<"\033[0m"<<endl;
            cout<<"\033[32m"<<"============================"<<"\033[0m"<<endl;
        }
    }
    static void TearDownTestCase()
    {
        if(rank_in_pool == 0)
        {
            cout<<"\033[32m"<<"============================"<<"\033[0m"<<endl;
            cout<<"\033[32m"<<"=     PW GPU MODULE TEST END   ="<<"\033[0m"<<endl;
            cout<<"\033[32m"<<"============================"<<"\033[0m"<<endl;
        }
    }
    void SetUp()
    {
        cout<<"\033[32m"<<"[    CASE  ]"<<"\033[0m"<<" ";
    }
    void TearDown(){}
};
    using set_3d_fft_box_cpu_op = ModulePW::set_3d_fft_box_op<double, base_device::DEVICE_CPU>;
    using set_3d_fft_box_gpu_op = ModulePW::set_3d_fft_box_op<double, base_device::DEVICE_GPU>;
    using set_recip_to_real_output_cpu_op = ModulePW::set_recip_to_real_output_op<double, base_device::DEVICE_CPU>;
    using set_recip_to_real_output_gpu_op = ModulePW::set_recip_to_real_output_op<double, base_device::DEVICE_GPU>;
    using set_real_to_recip_output_cpu_op = ModulePW::set_real_to_recip_output_op<double, base_device::DEVICE_CPU>;
    using set_real_to_recip_output_gpu_op = ModulePW::set_real_to_recip_output_op<double, base_device::DEVICE_GPU>;

    using resize_memory_complex_gpu_op
        = base_device::memory::resize_memory_op<std::complex<double>, base_device::DEVICE_GPU>;
    using delete_memory_complex_gpu_op
        = base_device::memory::delete_memory_op<std::complex<double>, base_device::DEVICE_GPU>;
    using synchronize_memory_complex_h2d_op = base_device::memory::
        synchronize_memory_op<std::complex<double>, base_device::DEVICE_GPU, base_device::DEVICE_CPU>;
    using synchronize_memory_complex_d2h_op = base_device::memory::
        synchronize_memory_op<std::complex<double>, base_device::DEVICE_CPU, base_device::DEVICE_GPU>;

    using resize_memory_double_gpu_op = base_device::memory::resize_memory_op<double, base_device::DEVICE_GPU>;
    using delete_memory_double_gpu_op = base_device::memory::delete_memory_op<double, base_device::DEVICE_GPU>;
    using synchronize_memory_double_h2d_op
        = base_device::memory::synchronize_memory_op<double, base_device::DEVICE_GPU, base_device::DEVICE_CPU>;
    using synchronize_memory_double_d2h_op
        = base_device::memory::synchronize_memory_op<double, base_device::DEVICE_CPU, base_device::DEVICE_GPU>;

    using delete_memory_int_gpu_op = base_device::memory::delete_memory_op<int, base_device::DEVICE_GPU>;
    using resize_memory_int_gpu_op = base_device::memory::resize_memory_op<int, base_device::DEVICE_GPU>;
    using synchronize_memory_int_h2d_op
        = base_device::memory::synchronize_memory_op<int, base_device::DEVICE_GPU, base_device::DEVICE_CPU>;
#endif