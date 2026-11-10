// Auxiliary real/reciprocal-space utilities shared by the gradient
// correction and meta-GGA stress:
// 1. grad_rho      : gradient of the density
// 2. grad_dot      : divergence of a vector field
// 3. laplacian_rho : Laplacian of the density
// 4. noncolin_rho  : diagonalizes the spin density matrix and gives the
//                    spin up and spin down components of the charge.

#include "xc_functional.h"
#include "source_base/constants.h"
#include "source_base/timer.h"

#if defined(__CUDA)
#include "source_basis/module_pw/kernels/pw_op.h"
#include "source_hamilt/module_xc/kernels/xc_functional_op.h"

#include <ATen/core/tensor.h>
#include <ATen/core/tensor_map.h>
#include <ATen/core/tensor_types.h>
#endif

#if defined(__CUDA)
namespace
{
bool xc_gpu_transform_eligible(const ModulePW::PW_Basis& basis)
{
    return basis.get_device() == "gpu" && basis.poolnproc == 1 && !basis.gamma_only
           && basis.fft3d_bundle.template get_auxr_3d_data<double>() != nullptr;
}

template <typename T>
ct::Tensor xc_copy_to_gpu(const T* source, int size)
{
    ct::Tensor result(ct::DataTypeToEnum<T>::value, ct::DeviceType::GpuDevice, {size});
    base_device::memory::synchronize_memory_op<T, base_device::DEVICE_GPU, base_device::DEVICE_CPU>()(
        result.data<T>(), source, size);
    return result;
}

ct::Tensor xc_density_grid_indices(const ModulePW::PW_Basis& basis)
{
    ct::Tensor indices(ct::DataType::DT_INT, ct::DeviceType::CpuDevice, {basis.npw});
    auto* index = indices.data<int>();
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
    for (int ig = 0; ig < basis.npw; ++ig)
    {
        const int iz = basis.ig2isz[ig] % basis.nz;
        const int is = basis.ig2isz[ig] / basis.nz;
        index[ig] = basis.is2fftixy[is] * basis.nz + iz;
    }
    return indices.to_device<ct::DEVICE_GPU>();
}

void xc_density_gradient_gpu(const std::complex<double>* rhog, ModuleBase::Vector3<double>* gradient,
                             const ModulePW::PW_Basis& basis, double scale)
{
    ModuleBase::timer::start("XC_Functional", "grad_rho_gpu");
    static_assert(sizeof(ModuleBase::Vector3<double>) == 3 * sizeof(double), "Vector3 must be contiguous");
    auto indices = xc_density_grid_indices(basis);
    auto reciprocal = xc_copy_to_gpu(rhog, basis.npw);
    auto gcar = xc_copy_to_gpu(reinterpret_cast<const double*>(basis.gcar), 3 * basis.npw);
    ct::Tensor result(ct::DataType::DT_DOUBLE, ct::DeviceType::GpuDevice, {basis.nrxx, 3});
    auto* grid = basis.fft3d_bundle.template get_auxr_3d_data<double>();
    for (int component = 0; component < 3; ++component)
    {
        base_device::memory::set_memory_op<std::complex<double>, base_device::DEVICE_GPU>()(grid, 0, basis.nxyz);
        hamilt::xc_gradient_coefficients_gpu(basis.npw, component, indices.data<int>(), gcar.data<double>(),
                                             reciprocal.data<std::complex<double>>(), grid);
        basis.fft3d_bundle.fft3D_backward(grid, grid);
        hamilt::xc_extract_real_gpu(basis.nrxx, 3, component, scale, grid, result.data<double>());
    }
    base_device::memory::synchronize_memory_op<double, base_device::DEVICE_CPU, base_device::DEVICE_GPU>()(
        reinterpret_cast<double*>(gradient), result.data<double>(), 3 * basis.nrxx);
    ModuleBase::timer::end("XC_Functional", "grad_rho_gpu");
}

void xc_density_divergence_gpu(const ModuleBase::Vector3<double>* vectors, double* divergence,
                               const ModulePW::PW_Basis& basis, double scale)
{
    ModuleBase::timer::start("XC_Functional", "grad_dot_gpu");
    static_assert(sizeof(ModuleBase::Vector3<double>) == 3 * sizeof(double), "Vector3 must be contiguous");
    auto indices = xc_density_grid_indices(basis);
    auto input = xc_copy_to_gpu(reinterpret_cast<const double*>(vectors), 3 * basis.nrxx);
    auto gcar = xc_copy_to_gpu(reinterpret_cast<const double*>(basis.gcar), 3 * basis.npw);
    ct::Tensor reciprocal(ct::DataType::DT_COMPLEX_DOUBLE, ct::DeviceType::GpuDevice, {basis.npw});
    ct::Tensor result(ct::DataType::DT_DOUBLE, ct::DeviceType::GpuDevice, {basis.nrxx});
    auto* grid = basis.fft3d_bundle.template get_auxr_3d_data<double>();
    for (int component = 0; component < 3; ++component)
    {
        hamilt::xc_pack_component_gpu(basis.nrxx, component, input.data<double>(), grid);
        basis.fft3d_bundle.fft3D_forward(grid, grid);
        hamilt::xc_accumulate_divergence_gpu(basis.npw, basis.nxyz, component, indices.data<int>(),
                                            gcar.data<double>(), grid, reciprocal.data<std::complex<double>>());
    }
    base_device::memory::set_memory_op<std::complex<double>, base_device::DEVICE_GPU>()(grid, 0, basis.nxyz);
    ModulePW::set_3d_fft_box_op<double, base_device::DEVICE_GPU>()(
        basis.npw, indices.data<int>(), reciprocal.data<std::complex<double>>(), grid);
    basis.fft3d_bundle.fft3D_backward(grid, grid);
    hamilt::xc_extract_real_gpu(basis.nrxx, 1, 0, scale, grid, result.data<double>());
    base_device::memory::synchronize_memory_op<double, base_device::DEVICE_CPU, base_device::DEVICE_GPU>()(
        divergence, result.data<double>(), basis.nrxx);
    ModuleBase::timer::end("XC_Functional", "grad_dot_gpu");
}
} // namespace
#endif

void XC_Functional::grad_rho(
    const std::complex<double>* rhog,
    ModuleBase::Vector3<double>* gdr,
    const ModulePW::PW_Basis* rho_basis,
    const double tpiba)
{
#if defined(__CUDA)
    if (xc_gpu_transform_eligible(*rho_basis))
    {
        xc_density_gradient_gpu(rhog, gdr, *rho_basis, tpiba);
        return;
    }
#endif
    std::vector<std::complex<double>> gdrtmp(rho_basis->nmaxgr);

    // the formula is : rho(r)^prime = int iG * rho(G)e^{iGr} dG
    for(int i = 0 ; i < 3 ; ++i)
    {
        // calculate the charge density gradient in reciprocal space.
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
        for(int ig=0; ig<rho_basis->npw; ig++)
        {
            gdrtmp[ig] = ModuleBase::IMAG_UNIT * rhog[ig] * rho_basis->gcar[ig][i];
        }

        // bring the gdr from G --> R
        rho_basis->recip2real(gdrtmp.data(), gdrtmp.data());

        // remember to multily 2pi/a0, which belongs to G vectors.
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
        for(int ir=0; ir<rho_basis->nrxx; ir++)
        {
            gdr[ir][i] = gdrtmp[ir].real() * tpiba;
        }
    }

    return;
}


void XC_Functional::grad_dot(
    const ModuleBase::Vector3<double>* h,
    double* dh,
    const ModulePW::PW_Basis* rho_basis,
    const double tpiba)
{
#if defined(__CUDA)
    if (xc_gpu_transform_eligible(*rho_basis))
    {
        xc_density_divergence_gpu(h, dh, *rho_basis, tpiba);
        return;
    }
#endif
    std::vector<std::complex<double>> aux(rho_basis->nmaxgr);
    std::vector<std::complex<double>> gaux(rho_basis->npw);

    for(int i = 0 ; i < 3 ; ++i)
    {
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
        for(int ir = 0; ir < rho_basis->nrxx; ++ir)
        {
            aux[ir] = std::complex<double>( h[ir][i], 0.0);
        }

        // bring to G space.
        rho_basis->real2recip(aux.data(), aux.data());
        if (i == 0)
        {
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
            for(int ig = 0; ig < rho_basis->npw; ++ig)
            {
                gaux[ig] =  ModuleBase::IMAG_UNIT * aux[ig] * rho_basis->gcar[ig][i];
            }
        }
        else
        {
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
            for(int ig = 0; ig < rho_basis->npw; ++ig)
            {
                gaux[ig] +=  ModuleBase::IMAG_UNIT * aux[ig] * rho_basis->gcar[ig][i];
            }
        }
    }

    // bring back to R space
    rho_basis->recip2real(gaux.data(), aux.data());

#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
    for(int ir=0; ir<rho_basis->nrxx; ir++)
    {
        dh[ir] = aux[ir].real() * tpiba;
    }

    return;
}

void XC_Functional::laplacian_rho(
    const std::complex<double>* rhog,
    double* lapl,
    const ModulePW::PW_Basis* rho_basis,
    const double tpiba)
{
    std::vector<std::complex<double>> lapl_tmp(rho_basis->nmaxgr);

    for(int ig=0; ig<rho_basis->npw; ig++)
    {
        double g2 = 0.0;
        for(int i=0; i<3; i++)
        {
            g2 += rho_basis->gcar[ig][i] * rho_basis->gcar[ig][i];
        }
        lapl_tmp[ig] = -rhog[ig] * g2;
    }
    rho_basis->recip2real(lapl_tmp.data(), lapl_tmp.data());
    for(int ir=0; ir<rho_basis->nrxx; ir++)
    {
        lapl[ir] = lapl_tmp[ir].real() * tpiba * tpiba;
    }
}

void XC_Functional::noncolin_rho(
    double *rhoout1,
    double *rhoout2,
    double *neg,
    const double*const*const rho,
    const int nrxx,
    const double* ux_,
    const bool lsign_)
{
    //this function diagonalizes the spin density matrix and gives as output the
    //spin up and spin down components of the charge.
    //If lsign is true up and dw are with respect to the fixed quantization axis
    //ux, otherwise rho + |m| is always rhoup and rho-|m| is always rhodw.
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
    for(int ir = 0;ir<nrxx;ir++)
    {
        neg[ir] = 1.0;
    }
    if(lsign_)
    {
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
        for(int ir = 0;ir<nrxx;ir++)
        {
            if(rho[1][ir]*ux_[0] + rho[2][ir]*ux_[1] + rho[3][ir]*ux_[2]>0)
            {
                neg[ir] = 1.0;
            }
            else
            {
                neg[ir] = -1.0;
            }
        }
    }
    // Note: with lsign_ = true the up/down channels are defined w.r.t. the
    // global axis ux_ through sign(m . ux_), evaluated point by point above;
    // with lsign_ = false, "up" is always the local |m| (neg stays +1
    // everywhere). amag = |m(r)| is non-negative; neg[ir] carries that sign.
    // The global axis ux_ itself is derived from the STRU initial moments in
    // cal_ux(). Because GGA also consumes the *gradients* of rhoout1/rhoout2,
    // a flipping sign field changes the gradient terms too, so the two
    // branches can give noticeably different v_xc even for predominantly
    // single-sign mz. (Locked by GgaGradTools.NoncolinRhoGlobalAxis in
    // module_xc/test/test_xc5.cpp.)
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(int ir = 0;ir<nrxx;ir++)
    {
        double amag = sqrt(pow(rho[1][ir],2)+pow(rho[2][ir],2)+pow(rho[3][ir],2));
        rhoout1[ir] = 0.5 * (rho[0][ir] + neg[ir] * amag);
        rhoout2[ir] = 0.5 * (rho[0][ir] - neg[ir] * amag);
    }
    return;
}
