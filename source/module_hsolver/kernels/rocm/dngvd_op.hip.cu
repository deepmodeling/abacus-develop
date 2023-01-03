#include "module_hsolver/include/dngvd_op.h"

#include <hip/hip_runtime.h>

namespace hsolver {

void createCUSOLVERhandle() {
    return;
}

void destoryCUSOLVERhandle() {
    return;
}

template <>
void dngvx_op<double, psi::DEVICE_GPU>::operator()(const psi::DEVICE_GPU* ctx,
                                                   const int nstart,
                                                   const int ldh,
                                                   const std::complex<double>* _hcc, // hcc
                                                   const std::complex<double>* _scc, // scc
                                                   const int nbands, // nbands
                                                   double* _eigenvalue,  // eigenvalue
                                                   std::complex<double>* _vcc) // vcc
{
    std::vector<std::complex<double>> hcc(nstart * nstart, {0, 0});
    std::vector<std::complex<double>> scc(nstart * nstart, {0, 0});
    std::vector<std::complex<double>> vcc(nstart * nstart, {0, 0});
    std::vector<double> eigenvalue(nbands, 0);
    hipMemcpy(hcc.data(), _hcc, sizeof(std::complex<double>) * hcc.size(), hipMemcpyDeviceToHost);
    hipMemcpy(scc.data(), _scc, sizeof(std::complex<double>) * scc.size(), hipMemcpyDeviceToHost);
    psi::DEVICE_CPU * cpu_ctx = {};
    dngvx_op<double, psi::DEVICE_CPU>()(cpu_ctx, nstart, ldh, hcc.data(), scc.data(), nbands, eigenvalue.data(), vcc.data());
    hipMemcpy(_vcc, vcc.data(), sizeof(std::complex<double>) * vcc.size(), hipMemcpyHostToDevice);
    hipMemcpy(_eigenvalue, eigenvalue.data(), sizeof(double) * eigenvalue.size(), hipMemcpyHostToDevice);
};

template <>
void dngv_op<double, psi::DEVICE_GPU>::operator()(const psi::DEVICE_GPU* ctx,
                                                  const int nstart,
                                                  const int ldh,
                                                  const std::complex<double>* _hcc,
                                                  const std::complex<double>* _scc,
                                                  double* _eigenvalue,
                                                  std::complex<double>* _vcc)
                                                  
{
    std::vector<std::complex<double>> hcc(nstart * nstart, {0, 0});
    std::vector<std::complex<double>> scc(nstart * nstart, {0, 0});
    std::vector<std::complex<double>> vcc(nstart * nstart, {0, 0});
    std::vector<double> eigenvalue(nstart, 0);
    hipMemcpy(hcc.data(), _hcc, sizeof(std::complex<double>) * hcc.size(), hipMemcpyDeviceToHost);
    hipMemcpy(scc.data(), _scc, sizeof(std::complex<double>) * scc.size(), hipMemcpyDeviceToHost);
    psi::DEVICE_CPU * cpu_ctx = {};
    dngv_op<double, psi::DEVICE_CPU>()(cpu_ctx, nstart, ldh, hcc.data(), scc.data(), eigenvalue.data(), vcc.data());
    hipMemcpy(_vcc, vcc.data(), sizeof(std::complex<double>) * vcc.size(), hipMemcpyHostToDevice);
    hipMemcpy(_eigenvalue, eigenvalue.data(), sizeof(double) * eigenvalue.size(), hipMemcpyHostToDevice);
}

template <>
void dngvd_op<double, psi::DEVICE_GPU>::operator()(const psi::DEVICE_GPU* ctx,
                                                   const int nstart,
                                                   const int ldh,
                                                   const std::complex<double>* _hcc,
                                                   const std::complex<double>* _scc,
                                                   double* _eigenvalue,
                                                   std::complex<double>* _vcc)
{
    std::vector<std::complex<double>> hcc(nstart * nstart, {0, 0});
    std::vector<std::complex<double>> scc(nstart * nstart, {0, 0});
    std::vector<std::complex<double>> vcc(nstart * nstart, {0, 0});
    std::vector<double> eigenvalue(nstart, 0);
    hipMemcpy(hcc.data(), _hcc, sizeof(std::complex<double>) * hcc.size(), hipMemcpyDeviceToHost);
    hipMemcpy(scc.data(), _scc, sizeof(std::complex<double>) * scc.size(), hipMemcpyDeviceToHost);
    psi::DEVICE_CPU * cpu_ctx = {};
    dngvd_op<double, psi::DEVICE_CPU>()(cpu_ctx, nstart, ldh, hcc.data(), scc.data(), eigenvalue.data(), vcc.data());
    hipMemcpy(_vcc, vcc.data(), sizeof(std::complex<double>) * vcc.size(), hipMemcpyHostToDevice);
    hipMemcpy(_eigenvalue, eigenvalue.data(), sizeof(double) * eigenvalue.size(), hipMemcpyHostToDevice);
}

template <>
void dnevx_op<double, psi::DEVICE_GPU>::operator()(const psi::DEVICE_GPU* ctx,
                                                   const int nstart,
                                                   const int ldh,
                                                   const std::complex<double>* _hcc,
                                                   const int m,
                                                   double* _eigenvalue,
                                                   std::complex<double>* _vcc)
{
    std::vector<std::complex<double>> hcc(ldh * ldh, {0, 0});
    std::vector<std::complex<double>> vcc(ldh * ldh, {0, 0});
    std::vector<double> eigenvalue(ldh, 0);
    hipMemcpy(hcc.data(), _hcc, sizeof(std::complex<double>) * hcc.size(), hipMemcpyDeviceToHost);
    psi::DEVICE_CPU * cpu_ctx = {};
    dnevx_op<double, psi::DEVICE_CPU>()(cpu_ctx, nstart, ldh, hcc.data(), m, eigenvalue.data(), vcc.data());
    hipMemcpy(_vcc, vcc.data(), sizeof(std::complex<double>) * vcc.size(), hipMemcpyHostToDevice);
    hipMemcpy(_eigenvalue, eigenvalue.data(), sizeof(double) * eigenvalue.size(), hipMemcpyHostToDevice);
}

} // namespace hsolver