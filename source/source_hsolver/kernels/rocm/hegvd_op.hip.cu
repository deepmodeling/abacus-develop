#include "source_hsolver/kernels/hegvd_op.h"

#include <hip/hip_runtime.h>
#include <base/macros/macros.h>

namespace hsolver {

// NOTE: mimicked from ../cuda/hegvd_op.cu for three hegvd_op

static hipsolverHandle_t hipsolver_H = nullptr;
// Test on DCU platform. When nstart is greater than 234, code on DCU performs better.
const int N_DCU = 234;

void createGpuSolverHandle() {
    if (hipsolver_H == nullptr)
    {
        hipsolverErrcheck(hipsolverCreate(&hipsolver_H));
    }
}

void destroyGpuSolverHandle() {
    if (hipsolver_H != nullptr)
    {
        hipsolverErrcheck(hipsolverDestroy(hipsolver_H));
        hipsolver_H = nullptr;
    }
}

#ifdef __LCAO
template <>
void hegvd_op<double, base_device::DEVICE_GPU>::operator()(const base_device::DEVICE_GPU* ctx,
                                                           const int nstart,
                                                           const int ldh,
                                                           const double* _hcc,
                                                           double* _scc,
                                                           double* _eigenvalue,
                                                           double* _vcc)
{
    // copied from ../cuda/hegvd_op.cu, "hegvd_op"
    assert(nstart == ldh);

    if (nstart > N_DCU){
        hipErrcheck(hipMemcpy(_vcc, _hcc, sizeof(double) * ldh * nstart, hipMemcpyDeviceToDevice));
        // now vcc contains hcc

        // prepare some values for hipsolverDnZhegvd_bufferSize
        int * devInfo = nullptr;
        int lwork = 0, info_gpu = 0;
        double * work = nullptr;
        hipErrcheck(hipMalloc((void**)&devInfo, sizeof(int)));
        hipsolverFillMode_t uplo = HIPSOLVER_FILL_MODE_UPPER;

        // calculate the sizes needed for pre-allocated buffer.
        hipsolverErrcheck(hipsolverDnDsygvd_bufferSize(
            hipsolver_H, HIPSOLVER_EIG_TYPE_1, HIPSOLVER_EIG_MODE_VECTOR, uplo,
            nstart,
            _vcc, ldh,
            _scc, ldh,
            _eigenvalue,
            &lwork));

        // allocate memery
        hipErrcheck(hipMalloc((void**)&work, sizeof(double) * lwork));

        // compute eigenvalues and eigenvectors.
        hipsolverErrcheck(hipsolverDnDsygvd(
            hipsolver_H, HIPSOLVER_EIG_TYPE_1, HIPSOLVER_EIG_MODE_VECTOR, uplo,
            nstart,
            _vcc, ldh,
            const_cast<double *>(_scc), ldh,
            _eigenvalue,
            work, lwork, devInfo));

        hipErrcheck(hipMemcpy(&info_gpu, devInfo, sizeof(int), hipMemcpyDeviceToHost));

        // free the buffer
        hipErrcheck(hipFree(work));
        hipErrcheck(hipFree(devInfo));
    }
    // if(fail_info != nullptr) *fail_info = info_gpu;
    else{
        std::vector<double> hcc(nstart * nstart, 0.0);
        std::vector<double> scc(nstart * nstart, 0.0);
        std::vector<double> vcc(nstart * nstart, 0.0);
        std::vector<double> eigenvalue(nstart, 0);
        hipErrcheck(hipMemcpy(hcc.data(), _hcc, sizeof(double) * hcc.size(), hipMemcpyDeviceToHost));
        hipErrcheck(hipMemcpy(scc.data(), _scc, sizeof(double) * scc.size(), hipMemcpyDeviceToHost));
        base_device::DEVICE_CPU* cpu_ctx = {};
        hegvd_op<double, base_device::DEVICE_CPU>()(cpu_ctx,
                                                   nstart,
                                                   ldh,
                                                   hcc.data(),
                                                   scc.data(),
                                                   eigenvalue.data(),
                                                   vcc.data());
        hipErrcheck(hipMemcpy(_vcc, vcc.data(), sizeof(double) * vcc.size(), hipMemcpyHostToDevice));
        hipErrcheck(hipMemcpy(_eigenvalue, eigenvalue.data(), sizeof(double) * eigenvalue.size(), hipMemcpyHostToDevice));
    }


}
#endif // __LCAO

template <>
void hegvd_op<std::complex<float>, base_device::DEVICE_GPU>::operator()(const base_device::DEVICE_GPU* ctx,
                                                                        const int nstart,
                                                                        const int ldh,
                                                                        const std::complex<float>* _hcc,
                                                                        const std::complex<float>* _scc,
                                                                        float* _eigenvalue,
                                                                        std::complex<float>* _vcc)
{
    // copied from ../cuda/hegvd_op.cu, "hegvd_op"
    assert(nstart == ldh);

    if (nstart > N_DCU){
        hipErrcheck(hipMemcpy(_vcc, _hcc, sizeof(std::complex<float>) * ldh * nstart, hipMemcpyDeviceToDevice));
        // now vcc contains hcc

        // prepare some values for hipsolverDnZhegvd_bufferSize
        int * devInfo = nullptr;
        int lwork = 0, info_gpu = 0;
        float2 * work = nullptr;
        hipErrcheck(hipMalloc((void**)&devInfo, sizeof(int)));
        hipsolverFillMode_t uplo = HIPSOLVER_FILL_MODE_UPPER;

        // calculate the sizes needed for pre-allocated buffer.
        hipsolverErrcheck(hipsolverDnChegvd_bufferSize(
            hipsolver_H, HIPSOLVER_EIG_TYPE_1, HIPSOLVER_EIG_MODE_VECTOR, uplo,
            nstart,
            reinterpret_cast<const float2 *>(_vcc), ldh,
            reinterpret_cast<const float2 *>(_scc), ldh,
            _eigenvalue,
            &lwork));

        // allocate memery
        hipErrcheck(hipMalloc((void**)&work, sizeof(float2) * lwork));

        // compute eigenvalues and eigenvectors.
        hipsolverErrcheck(hipsolverDnChegvd(
            hipsolver_H, HIPSOLVER_EIG_TYPE_1, HIPSOLVER_EIG_MODE_VECTOR, uplo,
            nstart,
            reinterpret_cast<float2 *>(_vcc), ldh,
            const_cast<float2 *>(reinterpret_cast<const float2 *>(_scc)), ldh,
            _eigenvalue,
            work, lwork, devInfo));

        hipErrcheck(hipMemcpy(&info_gpu, devInfo, sizeof(int), hipMemcpyDeviceToHost));
        // free the buffer
        hipErrcheck(hipFree(work));
        hipErrcheck(hipFree(devInfo));
    }
    // if(fail_info != nullptr) *fail_info = info_gpu;
    else{
        std::vector<std::complex<float>> hcc(nstart * nstart, {0, 0});
        std::vector<std::complex<float>> scc(nstart * nstart, {0, 0});
        std::vector<std::complex<float>> vcc(nstart * nstart, {0, 0});
        std::vector<float> eigenvalue(nstart, 0);
        hipErrcheck(hipMemcpy(hcc.data(), _hcc, sizeof(std::complex<float>) * hcc.size(), hipMemcpyDeviceToHost));
        hipErrcheck(hipMemcpy(scc.data(), _scc, sizeof(std::complex<float>) * scc.size(), hipMemcpyDeviceToHost));
        base_device::DEVICE_CPU* cpu_ctx = {};
        hegvd_op<std::complex<float>, base_device::DEVICE_CPU>()(cpu_ctx,
                                                                nstart,
                                                                ldh,
                                                                hcc.data(),
                                                                scc.data(),
                                                                eigenvalue.data(),
                                                                vcc.data());
        hipErrcheck(hipMemcpy(_vcc, vcc.data(), sizeof(std::complex<float>) * vcc.size(), hipMemcpyHostToDevice));
        hipErrcheck(hipMemcpy(_eigenvalue, eigenvalue.data(), sizeof(float) * eigenvalue.size(), hipMemcpyHostToDevice));
    }


}

template <>
void hegvd_op<std::complex<double>, base_device::DEVICE_GPU>::operator()(const base_device::DEVICE_GPU* ctx,
                                                                         const int nstart,
                                                                         const int ldh,
                                                                         const std::complex<double>* _hcc,
                                                                         const std::complex<double>* _scc,
                                                                         double* _eigenvalue,
                                                                         std::complex<double>* _vcc
                                                                        )
{
    // copied from ../cuda/hegvd_op.cu, "hegvd_op"
    // assert(nstart == ldh);

    // save a copy of scc in case the diagonalization fails
    if (nstart > N_DCU){
        std::vector<std::complex<double>> scc(nstart * nstart, {0, 0});
        hipErrcheck(hipMemcpy(scc.data(), _scc, sizeof(std::complex<double>) * scc.size(), hipMemcpyDeviceToHost));

        hipErrcheck(hipMemcpy(_vcc, _hcc, sizeof(std::complex<double>) * ldh * nstart, hipMemcpyDeviceToDevice));

        // now vcc contains hcc

        // prepare some values for hipsolverDnZhegvd_bufferSize
        int * devInfo = nullptr;
        int lwork = 0, info_gpu = 0;
        double2 * work = nullptr;
        hipErrcheck(hipMalloc((void**)&devInfo, sizeof(int)));
        hipsolverFillMode_t uplo = HIPSOLVER_FILL_MODE_UPPER;

        // calculate the sizes needed for pre-allocated buffer.
        hipsolverErrcheck(hipsolverDnZhegvd_bufferSize(
            hipsolver_H, HIPSOLVER_EIG_TYPE_1, HIPSOLVER_EIG_MODE_VECTOR, uplo,
            nstart,
            reinterpret_cast<const double2 *>(_vcc), ldh,
            reinterpret_cast<const double2 *>(_scc), ldh,
            _eigenvalue,
            &lwork));

        // allocate memery
        hipErrcheck(hipMalloc((void**)&work, sizeof(double2) * lwork));

        // compute eigenvalues and eigenvectors.
        hipsolverErrcheck(hipsolverDnZhegvd(
            hipsolver_H, HIPSOLVER_EIG_TYPE_1, HIPSOLVER_EIG_MODE_VECTOR, uplo,
            nstart,
            reinterpret_cast<double2 *>(_vcc), ldh,
            const_cast<double2 *>(reinterpret_cast<const double2 *>(_scc)), ldh,
            _eigenvalue,
            work, lwork, devInfo));

        hipErrcheck(hipMemcpy(&info_gpu, devInfo, sizeof(int), hipMemcpyDeviceToHost));
        // free the buffer
        hipErrcheck(hipFree(work));
        hipErrcheck(hipFree(devInfo));
    }
    // if(fail_info != nullptr) *fail_info = info_gpu;
    else{
        std::vector<std::complex<double>> hcc(nstart * nstart, {0, 0});
        std::vector<std::complex<double>> scc(nstart * nstart, {0, 0});
        std::vector<std::complex<double>> vcc(nstart * nstart, {0, 0});
        std::vector<double> eigenvalue(nstart, 0);
        hipErrcheck(hipMemcpy(hcc.data(), _hcc, sizeof(std::complex<double>) * hcc.size(), hipMemcpyDeviceToHost));
        hipErrcheck(hipMemcpy(scc.data(), _scc, sizeof(std::complex<double>) * scc.size(), hipMemcpyDeviceToHost));
        base_device::DEVICE_CPU* cpu_ctx = {};
        hegvd_op<std::complex<double>, base_device::DEVICE_CPU>()(cpu_ctx,
                                                                nstart,
                                                                ldh,
                                                                hcc.data(),
                                                                scc.data(),
                                                                eigenvalue.data(),
                                                                vcc.data());
        hipErrcheck(hipMemcpy(_vcc, vcc.data(), sizeof(std::complex<double>) * vcc.size(), hipMemcpyHostToDevice));
        hipErrcheck(hipMemcpy(_eigenvalue, eigenvalue.data(), sizeof(double) * eigenvalue.size(), hipMemcpyHostToDevice));
    }







}

#ifdef __LCAO
template <>
void heevx_op<double, base_device::DEVICE_GPU>::operator()(const base_device::DEVICE_GPU* ctx,
                                                           const int nstart,
                                                           const int ldh,
                                                           const double* _hcc,
                                                           const int m,
                                                           double* _eigenvalue,
                                                           double* _vcc)
{
    std::vector<double> hcc(ldh * ldh, 0.0);
    std::vector<double> vcc(ldh * ldh, 0.0);
    std::vector<double> eigenvalue(ldh, 0);
    hipErrcheck(hipMemcpy(hcc.data(), _hcc, sizeof(double) * hcc.size(), hipMemcpyDeviceToHost));
    base_device::DEVICE_CPU* cpu_ctx = {};
    heevx_op<double, base_device::DEVICE_CPU>()(cpu_ctx, nstart, ldh, hcc.data(), m, eigenvalue.data(), vcc.data());
    hipErrcheck(hipMemcpy(_vcc, vcc.data(), sizeof(double) * vcc.size(), hipMemcpyHostToDevice));
    hipErrcheck(hipMemcpy(_eigenvalue, eigenvalue.data(), sizeof(double) * eigenvalue.size(), hipMemcpyHostToDevice));
}
#endif // __LCAO

template <>
void heevx_op<std::complex<float>, base_device::DEVICE_GPU>::operator()(const base_device::DEVICE_GPU* ctx,
                                                                        const int nstart,
                                                                        const int ldh,
                                                                        const std::complex<float>* _hcc,
                                                                        const int m,
                                                                        float* _eigenvalue,
                                                                        std::complex<float>* _vcc)
{
    std::vector<std::complex<float>> hcc(ldh * ldh, {0, 0});
    std::vector<std::complex<float>> vcc(ldh * ldh, {0, 0});
    std::vector<float> eigenvalue(ldh, 0);
    hipErrcheck(hipMemcpy(hcc.data(), _hcc, sizeof(std::complex<float>) * hcc.size(), hipMemcpyDeviceToHost));
    base_device::DEVICE_CPU* cpu_ctx = {};
    heevx_op<std::complex<float>, base_device::DEVICE_CPU>()(cpu_ctx,
                                                             nstart,
                                                             ldh,
                                                             hcc.data(),
                                                             m,
                                                             eigenvalue.data(),
                                                             vcc.data());
    hipErrcheck(hipMemcpy(_vcc, vcc.data(), sizeof(std::complex<float>) * vcc.size(), hipMemcpyHostToDevice));
    hipErrcheck(hipMemcpy(_eigenvalue, eigenvalue.data(), sizeof(float) * eigenvalue.size(), hipMemcpyHostToDevice));
}

template <>
void heevx_op<std::complex<double>, base_device::DEVICE_GPU>::operator()(const base_device::DEVICE_GPU* ctx,
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
    hipErrcheck(hipMemcpy(hcc.data(), _hcc, sizeof(std::complex<double>) * hcc.size(), hipMemcpyDeviceToHost));
    base_device::DEVICE_CPU* cpu_ctx = {};
    heevx_op<std::complex<double>, base_device::DEVICE_CPU>()(cpu_ctx,
                                                              nstart,
                                                              ldh,
                                                              hcc.data(),
                                                              m,
                                                              eigenvalue.data(),
                                                              vcc.data());
    hipErrcheck(hipMemcpy(_vcc, vcc.data(), sizeof(std::complex<double>) * vcc.size(), hipMemcpyHostToDevice));
    hipErrcheck(hipMemcpy(_eigenvalue, eigenvalue.data(), sizeof(double) * eigenvalue.size(), hipMemcpyHostToDevice));
}

template <>
void hegvx_op<std::complex<float>, base_device::DEVICE_GPU>::operator()(const base_device::DEVICE_GPU* d,
                                                                        const int nbase,
                                                                        const int ldh,
                                                                        std::complex<float>* hcc,
                                                                        std::complex<float>* scc,
                                                                        const int m,
                                                                        float* eigenvalue,
                                                                        std::complex<float>* vcc)
{
}

template <>
void hegvx_op<std::complex<double>, base_device::DEVICE_GPU>::operator()(const base_device::DEVICE_GPU* d,
                                                                         const int nbase,
                                                                         const int ldh,
                                                                         std::complex<double>* hcc,
                                                                         std::complex<double>* scc,
                                                                         const int m,
                                                                         double* eigenvalue,
                                                                         std::complex<double>* vcc)
{
}

#ifdef __LCAO
template <>
void hegvx_op<double, base_device::DEVICE_GPU>::operator()(const base_device::DEVICE_GPU* d,
                                                           const int nbase,
                                                           const int ldh,
                                                           double* hcc,
                                                           double* scc,
                                                           const int m,
                                                           double* eigenvalue,
                                                           double* vcc)
{
}
#endif // __LCAO

} // namespace hsolver