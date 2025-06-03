#include "test.cuh"
#include <cstdio>
#include <cufft.h>
#include <cuda_runtime.h>
#include <iostream>
#include <complex>

namespace ModulePW
{
    void test()
    {
    ModulePW::FFT_Bundle fft_bundle;
    fft_bundle.setfft("gpu", "double");
    fft_bundle.initfft(16, 16, 16, 64, 64, 1, 1, 1, false, false, false);
    fft_bundle.setupFFT();
        // Test the FFT functionality here
        // This is a placeholder for the actual test implementation
    const int BATCH1 = 4;       // 批处理大小
    const int BATCH = 1;       // 批处理大小
    const int NX = 16;          // 第一维大小
    const int NY = 16;          // 第二维大小
    const int NZ = 16;          // 第三维大小
    const int R2C_OUPUT_SIZE = NX * NY * NZ; // 复数输出大小（利用对称性）
    // const int size = BATCH * NX * NY * NZ;
    // 分配主机内存
    std::complex<double>* h_input = new std::complex<double>[BATCH * NX * NY * NZ];
    cufftDoubleComplex*   h_output = new cufftDoubleComplex[BATCH * R2C_OUPUT_SIZE];
    std::complex<double>* h_output_complex = new std::complex<double>[BATCH * NX * NY * NZ];
    // 初始化输入数据（示例）
    for (int b = 0; b < BATCH; ++b) {
        for (int i = 0; i < NX*NY*NZ; ++i) {
            h_input[b * NX*NY*NZ + i] = std::complex<double>(2,0) ;
            h_output_complex[b * NX*NY*NZ + i] = std::complex<double>(0,0) ;
        }
    }

    // 分配设备内存
    // std::complex<double>* d_input;
    // cufftDoubleComplex* d_output;
    // std::complex<double>* d_output_complex;
    // cudaMalloc(&d_input, BATCH * NX * NY * NZ * sizeof(std::complex<double>));
    // cudaMalloc(&d_output, BATCH * R2C_OUPUT_SIZE * sizeof(cufftDoubleComplex));
    // cudaMalloc(&d_output_complex, BATCH * R2C_OUPUT_SIZE * sizeof(std::complex<double>));

    // 拷贝数据到设备
    cudaMemcpy(fft_bundle.get_auxr_3d_data<double>(), h_input, BATCH * NX * NY * NZ * sizeof(std::complex<double>), cudaMemcpyHostToDevice);

    // 创建 cuFFT 计划
    cufftHandle plan;
    int rank = 3;                  // 3D 变换
    int n[3] = {NX, NY, NZ};       // 各维度大小
    int inembed[3] = {NX, NY, NZ}; // 输入步长（紧密排列）
    int onembed[3] = {NX, NY, NZ}; // 输出步长
    int istride = 1;                // 输入连续元素间距（1 = 无间隔）
    int ostride = 1;                // 输出连续元素间距
    int idist = NX * NY * NZ;       // 输入批处理间距
    int odist = R2C_OUPUT_SIZE;     // 输出批处理间距
    cufftPlanMany(&plan, rank, n,
                  inembed, istride, idist, // 输入描述
                  onembed, ostride, odist, // 输出描述
                  CUFFT_Z2Z, BATCH1);       // 类型和批处理量

    // 执行 R2C 变换
    // cufftExecZ2Z(plan, 
    //             reinterpret_cast<cufftDoubleComplex*>(d_input), 
    //             d_output,CUFFT_FORWARD);
    // cufftExecZ2Z(plan, reinterpret_cast<cufftDoubleComplex*>(fft_bundle.get_auxr_3d_data<double>()), 
    //                   reinterpret_cast<cufftDoubleComplex*>(fft_bundle.get_auxr_3d_data<double>()),
    //                   CUFFT_FORWARD);
    fft_bundle.fft3D_forward(fft_bundle.get_auxr_3d_data<double>(), fft_bundle.get_auxr_3d_data<double>());
    // 同步等待完成
    // cudaDeviceSynchronize();

    // 将结果拷贝回主机
    // cudaMemcpy(h_output, d_output, BATCH * R2C_OUPUT_SIZE * sizeof(cufftDoubleComplex), cudaMemcpyDeviceToHost);
    cudaMemcpy(h_output_complex, fft_bundle.get_auxr_3d_data<double>(), BATCH * R2C_OUPUT_SIZE * sizeof(std::complex<double>), cudaMemcpyDeviceToHost);
    // for (int b = 0; b < BATCH; ++b) {
    //     for (int i = 0; i < R2C_OUPUT_SIZE; ++i) {
    //         if (h_output[i].x != 0)
    //         {
    //             std::cout << "BATCH is " << b<< "  " << i << ": ";
    //             std::cout << "  " << h_output[b * R2C_OUPUT_SIZE + i].x << " + " 
    //                       << h_output[b * R2C_OUPUT_SIZE + i].y << "i" << std::endl;
    //         }
    //     }
    // }

    for (int b = 0; b < BATCH; ++b) {
        for (int i = 0; i < R2C_OUPUT_SIZE; ++i) {
            if (h_output_complex[i].real() > 0)
            {
                std::cout << "BATCH is " << b<< "  " << i << ": ";
                std::cout << "  " << h_output_complex[b * R2C_OUPUT_SIZE + i].real() << " + " 
                          << h_output_complex[b * R2C_OUPUT_SIZE + i].imag() << "i" << std::endl;
            }
        }
    }
    // 清理资源
    cufftDestroy(plan);
    // cudaFree(d_input);
    // cudaFree(d_output);
    delete[] h_input;
    // delete[] h_output;
    delete[] h_output_complex;
    }
} // namespace ModulePW