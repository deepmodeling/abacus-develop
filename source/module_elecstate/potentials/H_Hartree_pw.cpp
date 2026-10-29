#include "H_Hartree_pw.h"

#include "module_parameter/parameter.h"
#include "module_base/constants.h"
#include "module_base/timer.h"
#include "module_base/parallel_reduce.h"
#include <algorithm>
#if defined(__CUDA)
#include "module_basis/module_pw/kernels/pw_op.h"
#include "module_elecstate/kernels/hartree_op.h"
#include <ATen/core/tensor.h>
#endif

namespace elecstate
{

#if defined(__CUDA)
namespace
{

ModuleBase::matrix hartree_gpu(const UnitCell& cell, const ModulePW::PW_Basis& basis,
                               int nspin, const double* const* rho, double& energy)
{
    using Complex = std::complex<double>;
    const base_device::DEVICE_GPU* gpu = nullptr;
    const base_device::DEVICE_CPU* cpu = nullptr;
    const int channels = nspin == 2 ? 2 : 1;
    const int blocks = hartree_reduction_size_gpu(basis.npw);
    ct::Tensor density(ct::DataType::DT_DOUBLE, ct::DeviceType::GpuDevice, {channels, basis.nrxx});
    ct::Tensor gg(ct::DataType::DT_DOUBLE, ct::DeviceType::GpuDevice, {basis.npw});
    ct::Tensor coefficients(ct::DataType::DT_COMPLEX_DOUBLE, ct::DeviceType::GpuDevice, {basis.npw});
    ct::Tensor sums(ct::DataType::DT_DOUBLE, ct::DeviceType::GpuDevice, {blocks});
    for (int spin = 0; spin < channels; ++spin)
    {
        base_device::memory::synchronize_memory_op<double, base_device::DEVICE_GPU, base_device::DEVICE_CPU>()(
            gpu, cpu, density.data<double>() + spin * basis.nrxx, rho[spin], basis.nrxx);
    }
    base_device::memory::synchronize_memory_op<double, base_device::DEVICE_GPU, base_device::DEVICE_CPU>()(
        gpu, cpu, gg.data<double>(), basis.gg, basis.npw);
    auto* grid = basis.ft.get_auxr_3d_data<double>();
    hartree_pack_density_gpu(basis.nrxx, channels, density.data<double>(), grid);
    basis.ft.fft3D_forward(gpu, grid, grid);
    hartree_coefficients_gpu(basis.npw, basis.d_ig2fftbox, gg.data<double>(), 1.0 / basis.nxyz,
                             ModuleBase::e2 * ModuleBase::FOUR_PI, cell.tpiba2, grid,
                             coefficients.data<Complex>(), sums.data<double>());
    base_device::memory::set_memory_op<Complex, base_device::DEVICE_GPU>()(gpu, grid, 0, basis.nxyz);
    ModulePW::set_3d_fft_box_op<double, base_device::DEVICE_GPU>()(
        gpu, basis.npw, basis.d_ig2fftbox, coefficients.data<Complex>(), grid);
    basis.ft.fft3D_backward(gpu, grid, grid);
    ModulePW::extract_real_grid_gpu(basis.nrxx, false, 1.0, grid, density.data<double>());
    ModuleBase::matrix potential(nspin, basis.nrxx);
    base_device::memory::synchronize_memory_op<double, base_device::DEVICE_CPU, base_device::DEVICE_GPU>()(
        cpu, gpu, potential.c, density.data<double>(), basis.nrxx);
    if (nspin != 4)
    {
        for (int spin = 1; spin < nspin; ++spin)
        {
            std::copy(potential.c, potential.c + basis.nrxx, potential.c + spin * basis.nrxx);
        }
    }
    auto host_sums = sums.to_device<ct::DEVICE_CPU>();
    energy = 0;
    for (int block = 0; block < blocks; ++block) { energy += host_sums.data<double>()[block]; }
    energy *= 0.5 * cell.omega;
    return potential;
}

} // namespace
#endif

double H_Hartree_pw::hartree_energy = 0.0;

//--------------------------------------------------------------------
// Transform charge density to hartree potential.
//--------------------------------------------------------------------
ModuleBase::matrix H_Hartree_pw::v_hartree(const UnitCell &cell,
                                           ModulePW::PW_Basis *rho_basis,
                                           const int &nspin,
                                           const double *const *const rho)
{
    ModuleBase::TITLE("H_Hartree_pw", "v_hartree");
    ModuleBase::timer::tick("H_Hartree_pw", "v_hartree");

#if defined(__CUDA)
    if (rho_basis->d_ig2fftbox != nullptr)
    {
        auto potential = hartree_gpu(cell, *rho_basis, nspin, rho, H_Hartree_pw::hartree_energy);
        ModuleBase::timer::tick("H_Hartree_pw", "v_hartree");
        return potential;
    }
#endif

    //  Hartree potential VH(r) from n(r)
    std::vector<std::complex<double>> Porter(rho_basis->nmaxgr);
    const int nspin0 = (nspin == 2) ? 2 : 1;
    for (int is = 0; is < nspin0; is++)
    {
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 256)
#endif
        for (int ir = 0; ir < rho_basis->nrxx; ir++)
            Porter[ir] += std::complex<double>(rho[is][ir], 0.0);
    }
    //=============================
    //  bring rho (aux) to G space
    //=============================
    rho_basis->real2recip(Porter.data(), Porter.data());

    //=======================================================
    // calculate hartree potential in G-space (NB: V(G=0)=0 )
    //=======================================================

    double ehart = 0.0;

    std::vector<std::complex<double>> vh_g(rho_basis->npw);
#ifdef _OPENMP
#pragma omp parallel for reduction(+:ehart)
#endif
    for (int ig = 0; ig < rho_basis->npw; ig++)
    {
        if (rho_basis->gg[ig] >= 1.0e-8) // LiuXh 20180410
        {
            const double fac = ModuleBase::e2 * ModuleBase::FOUR_PI / (cell.tpiba2 * rho_basis->gg[ig]);

            ehart += (conj(Porter[ig]) * Porter[ig]).real() * fac;
            vh_g[ig] = fac * Porter[ig];
        }
    }

    Parallel_Reduce::reduce_pool(ehart);
    ehart *= 0.5 * cell.omega;
    // std::cout << " ehart=" << ehart << std::endl;
    H_Hartree_pw::hartree_energy = ehart;

    //==========================================
    // transform hartree potential to real space
    //==========================================
    rho_basis->recip2real(vh_g.data(), Porter.data());

    //==========================================
    // Add hartree potential to the xc potential
    //==========================================
    ModuleBase::matrix v(nspin, rho_basis->nrxx);
    if (nspin == 4)
    {
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 512)
#endif
        for (int ir = 0; ir < rho_basis->nrxx; ir++)
            v(0, ir) = Porter[ir].real();
    }
    else
    {
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static, 512)
#endif
        for (int is = 0; is < nspin; is++)
            for (int ir = 0; ir < rho_basis->nrxx; ir++)
                v(is, ir) = Porter[ir].real();
    }

    ModuleBase::timer::tick("H_Hartree_pw", "v_hartree");
    return v;
} // end subroutine v_h

PotHartree::PotHartree(const ModulePW::PW_Basis* rho_basis_in)
{
    this->rho_basis_ = rho_basis_in;
    this->dynamic_mode = true;
    this->fixed_mode = false;
}

void PotHartree::cal_v_eff(const Charge*const chg, const UnitCell*const ucell, ModuleBase::matrix& v_eff)
{
    if(PARAM.inp.use_paw)
    {
        double ** rho_tmp;
        rho_tmp = new double*[chg->nspin];
        for(int is = 0; is < chg->nspin; is++)
        {
            rho_tmp[is] = new double[rho_basis_->nrxx];
            for(int ir = 0; ir < rho_basis_->nrxx; ir++)
            {
                rho_tmp[is][ir] = chg->rho[is][ir] + chg->nhat[is][ir];
            }
        }
        v_eff += H_Hartree_pw::v_hartree(*ucell, const_cast<ModulePW::PW_Basis*>(this->rho_basis_), v_eff.nr, rho_tmp);

        for(int is = 0; is < chg->nspin; is++)
        {
            delete[] rho_tmp[is];
        }
        delete[] rho_tmp;
    }
    else
    {
        v_eff += H_Hartree_pw::v_hartree(*ucell, const_cast<ModulePW::PW_Basis*>(this->rho_basis_), v_eff.nr, chg->rho);
    }
    return;
}

} // namespace elecstate
