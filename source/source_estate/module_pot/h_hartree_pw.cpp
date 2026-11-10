#include "h_hartree_pw.h"

#include "source_io/module_parameter/parameter.h"
#include "source_base/constants.h"
#include "source_base/timer.h"
#include "source_base/parallel_reduce.h"

#include <algorithm>
#if defined(__CUDA)
#include "source_basis/module_pw/kernels/pw_op.h"
#include "source_estate/kernels/hartree_op.h"
#endif

namespace elecstate
{

double H_Hartree_pw::hartree_energy = 0.0;

#if defined(__CUDA)
namespace
{

bool hartree_gpu_eligible(const ModulePW::PW_Basis& basis)
{
    return basis.get_device() == "gpu" && basis.poolnproc == 1 && !basis.gamma_only
           && basis.ig2ixyz_gpu != nullptr;
}

ModuleBase::matrix hartree_gpu(const UnitCell& cell, const ModulePW::PW_Basis& basis,
                               int nspin, const double* const* rho, double& energy)
{
    using Complex = std::complex<double>;
    const int channels = nspin == 2 ? 2 : 1;
    const int blocks = hartree_reduction_size_gpu(basis.npw);
    std::vector<double> density(channels * basis.nrxx);
    for (int spin = 0; spin < channels; ++spin)
    {
        std::copy(rho[spin], rho[spin] + basis.nrxx, density.data() + spin * basis.nrxx);
    }
    std::vector<double> host_gg(basis.npw);
    std::copy(basis.gg, basis.gg + basis.npw, host_gg.begin());

    double* d_density = nullptr;
    double* d_gg = nullptr;
    Complex* d_coefficients = nullptr;
    double* d_sums = nullptr;
    base_device::memory::resize_memory_op<double, base_device::DEVICE_GPU>()(d_density, channels * basis.nrxx);
    base_device::memory::resize_memory_op<double, base_device::DEVICE_GPU>()(d_gg, basis.npw);
    base_device::memory::resize_memory_op<Complex, base_device::DEVICE_GPU>()(d_coefficients, basis.npw);
    base_device::memory::resize_memory_op<double, base_device::DEVICE_GPU>()(d_sums, blocks);
    base_device::memory::synchronize_memory_op<double, base_device::DEVICE_GPU, base_device::DEVICE_CPU>()(
        d_density, density.data(), channels * basis.nrxx);
    base_device::memory::synchronize_memory_op<double, base_device::DEVICE_GPU, base_device::DEVICE_CPU>()(
        d_gg, host_gg.data(), basis.npw);

    auto* grid = basis.fft_bundle.template get_auxr_3d_data<double>();
    hartree_pack_density_gpu(basis.nrxx, channels, d_density, grid);
    basis.fft_bundle.fft3D_forward(grid, grid);
    hartree_coefficients_gpu(basis.npw, basis.ig2ixyz_gpu, d_gg, 1.0 / basis.nxyz,
                             ModuleBase::e2 * ModuleBase::FOUR_PI, cell.tpiba2, grid,
                             d_coefficients, d_sums);
    base_device::memory::set_memory_op<Complex, base_device::DEVICE_GPU>()(grid, 0, basis.nxyz);
    ModulePW::set_3d_fft_box_op<double, base_device::DEVICE_GPU>()(
        basis.npw, basis.ig2ixyz_gpu, d_coefficients, grid);
    basis.fft_bundle.fft3D_backward(grid, grid);
    ModulePW::extract_real_grid_gpu(basis.nrxx, false, 1.0, grid, d_density);
    ModuleBase::matrix potential(nspin, basis.nrxx);
    base_device::memory::synchronize_memory_op<double, base_device::DEVICE_CPU, base_device::DEVICE_GPU>()(
        potential.c, d_density, basis.nrxx);
    if (nspin != 4)
    {
        for (int spin = 1; spin < nspin; ++spin)
        {
            std::copy(potential.c, potential.c + basis.nrxx, potential.c + spin * basis.nrxx);
        }
    }
    std::vector<double> host_sums(blocks);
    base_device::memory::synchronize_memory_op<double, base_device::DEVICE_CPU, base_device::DEVICE_GPU>()(
        host_sums.data(), d_sums, blocks);
    energy = 0;
    for (int block = 0; block < blocks; ++block)
    {
        energy += host_sums[block];
    }
    energy *= 0.5 * cell.omega;

    base_device::memory::delete_memory_op<double, base_device::DEVICE_GPU>()(d_density);
    base_device::memory::delete_memory_op<double, base_device::DEVICE_GPU>()(d_gg);
    base_device::memory::delete_memory_op<Complex, base_device::DEVICE_GPU>()(d_coefficients);
    base_device::memory::delete_memory_op<double, base_device::DEVICE_GPU>()(d_sums);
    return potential;
}

} // namespace
#endif

//--------------------------------------------------------------------
// Transform charge density to hartree potential.
//--------------------------------------------------------------------
ModuleBase::matrix H_Hartree_pw::v_hartree(const UnitCell &cell,
                                           ModulePW::PW_Basis *rho_basis,
                                           const int &nspin,
                                           const double *const *const rho)
{
    ModuleBase::TITLE("H_Hartree_pw", "v_hartree");
    ModuleBase::timer::start("H_Hartree_pw", "v_hartree");

#if defined(__CUDA)
    if (hartree_gpu_eligible(*rho_basis))
    {
        auto potential = hartree_gpu(cell, *rho_basis, nspin, rho, H_Hartree_pw::hartree_energy);
        ModuleBase::timer::end("H_Hartree_pw", "v_hartree");
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
    const int ig0 = rho_basis->ig_gge0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+:ehart)
#endif
    for (int ig = 0; ig < rho_basis->npw; ig++)
    {
        if (ig == ig0) 
        {
            continue; // skip G=0
        }
        const double fac = ModuleBase::e2 * ModuleBase::FOUR_PI / (cell.tpiba2 * rho_basis->gg[ig]);
        ehart += (conj(Porter[ig]) * Porter[ig]).real() * fac;
        vh_g[ig] = fac * Porter[ig];
        
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

    ModuleBase::timer::end("H_Hartree_pw", "v_hartree");
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
    v_eff += H_Hartree_pw::v_hartree(*ucell, const_cast<ModulePW::PW_Basis*>(this->rho_basis_), v_eff.nr, chg->rho);
    return;
}

} // namespace elecstate
