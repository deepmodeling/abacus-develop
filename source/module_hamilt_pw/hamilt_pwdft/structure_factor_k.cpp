#include "module_base/memory.h"
#include "module_base/timer.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_hamilt_pw/hamilt_pwdft/kernels/wf_op.h"
#include "module_psi/kernels/device.h"
#include "structure_factor.h"
std::complex<double>* Structure_Factor::get_sk(const int ik,
                                               const int it,
                                               const int ia,
                                               const ModulePW::PW_Basis_K* wfc_basis) const
{
    ModuleBase::timer::tick("Structure_Factor", "get_sk");
    const double arg = (wfc_basis->kvec_c[ik] * GlobalC::ucell.atoms[it].tau[ia]) * ModuleBase::TWO_PI;
    const std::complex<double> kphase = std::complex<double>(cos(arg), -sin(arg));
    const int npw = wfc_basis->npwk[ik];
    std::complex<double> *sk = new std::complex<double>[npw];
    const int nx = wfc_basis->nx, ny = wfc_basis->ny, nz = wfc_basis->nz;
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int igl = 0; igl < npw; ++igl)
    {
        const int isz = wfc_basis->getigl2isz(ik, igl);
        int iz = isz % nz;
        const int is = isz / nz;
        const int ixy = wfc_basis->is2fftixy[is];
        int ix = ixy / wfc_basis->fftny;
        int iy = ixy % wfc_basis->fftny;
        if (ix >= int(nx / 2) + 1)
            ix -= nx;
        if (iy >= int(ny / 2) + 1)
            iy -= ny;
        if (iz >= int(nz / 2) + 1)
            iz -= nz;
        ix += this->rho_basis->nx;
        iy += this->rho_basis->ny;
        iz += this->rho_basis->nz;
        const int iat = GlobalC::ucell.itia2iat(it, ia);
        sk[igl] = kphase * this->eigts1(iat, ix) * this->eigts2(iat, iy) * this->eigts3(iat, iz);
    }
    ModuleBase::timer::tick("Structure_Factor", "get_sk");
    return sk;
}

template <typename FPTYPE, typename Device>
void Structure_Factor::get_sk(Device* ctx,
                              const int ik,
                              const ModulePW::PW_Basis_K* wfc_basis,
                              std::complex<FPTYPE>* sk) const
{
    ModuleBase::timer::tick("Structure_Factor", "get_sk");

    psi::DEVICE_CPU *cpu_ctx = {};
    psi::AbacusDevice_t device = psi::device::get_device_type<Device>(ctx);
    using cal_sk_op = hamilt::cal_sk_op<FPTYPE, Device>;
    using resmem_int_op = psi::memory::resize_memory_op<int, Device>;
    using delmem_int_op = psi::memory::delete_memory_op<int, Device>;
    using syncmem_int_op = psi::memory::synchronize_memory_op<int, Device, psi::DEVICE_CPU>;

    using resmem_var_op = psi::memory::resize_memory_op<FPTYPE, Device>;
    using delmem_var_op = psi::memory::delete_memory_op<FPTYPE, Device>;
    using syncmem_var_op = psi::memory::synchronize_memory_op<FPTYPE, Device, psi::DEVICE_CPU>;

    int iat = 0, _npw = wfc_basis->npwk[ik], eigts1_nc = this->eigts1.nc, eigts2_nc = this->eigts2.nc,
            eigts3_nc = this->eigts3.nc;
    int *igl2isz = nullptr, *is2fftixy = nullptr, *atom_na = nullptr, *h_atom_na = new int[GlobalC::ucell.ntype];
    FPTYPE *atom_tau = nullptr, *h_atom_tau = new FPTYPE[GlobalC::ucell.nat * 3], *kvec = wfc_basis->get_kvec_c_data<FPTYPE>();
    std::complex<FPTYPE> *eigts1 = this->get_eigts1_data<FPTYPE>(), *eigts2 = this->get_eigts2_data<FPTYPE>(),
            *eigts3 = this->get_eigts3_data<FPTYPE>();
    for (int it = 0; it < GlobalC::ucell.ntype; it++)
    {
        h_atom_na[it] = GlobalC::ucell.atoms[it].na;
    }
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int iat = 0; iat < GlobalC::ucell.nat; iat++)
    {
        int it = GlobalC::ucell.iat2it[iat];
        int ia = GlobalC::ucell.iat2ia[iat];
        auto *tau = reinterpret_cast<double *>(GlobalC::ucell.atoms[it].tau);
        h_atom_tau[iat * 3 + 0] = static_cast<FPTYPE>(tau[ia * 3 + 0]);
        h_atom_tau[iat * 3 + 1] = static_cast<FPTYPE>(tau[ia * 3 + 1]);
        h_atom_tau[iat * 3 + 2] = static_cast<FPTYPE>(tau[ia * 3 + 2]);
    }
    if (device == psi::GpuDevice)
    {
        resmem_int_op()(ctx, atom_na, GlobalC::ucell.ntype);
        syncmem_int_op()(ctx, cpu_ctx, atom_na, h_atom_na, GlobalC::ucell.ntype);

        resmem_var_op()(ctx, atom_tau, GlobalC::ucell.nat * 3);
        syncmem_var_op()(ctx, cpu_ctx, atom_tau, h_atom_tau, GlobalC::ucell.nat * 3);

        igl2isz = wfc_basis->d_igl2isz_k;
        is2fftixy = wfc_basis->d_is2fftixy;
    }
    else
    {
        atom_na = h_atom_na;
        atom_tau = h_atom_tau;
        igl2isz = wfc_basis->igl2isz_k;
        is2fftixy = wfc_basis->is2fftixy;
    }

    cal_sk_op()(ctx,
                ik,
                GlobalC::ucell.ntype,
                wfc_basis->nx,
                wfc_basis->ny,
                wfc_basis->nz,
                this->rho_basis->nx,
                this->rho_basis->ny,
                this->rho_basis->nz,
                _npw,
                wfc_basis->npwk_max,
                wfc_basis->fftny,
                eigts1_nc,
                eigts2_nc,
                eigts3_nc,
                atom_na,
                igl2isz,
                is2fftixy,
                ModuleBase::TWO_PI,
                kvec,
                atom_tau,
                eigts1,
                eigts2,
                eigts3,
                sk);
    if (device == psi::GpuDevice)
    {
        delmem_int_op()(ctx, atom_na);
        delmem_var_op()(ctx, atom_tau);
    }
    delete[] h_atom_na;
    delete[] h_atom_tau;
    ModuleBase::timer::tick("Structure_Factor", "get_sk");
}

std::complex<double>* Structure_Factor::get_skq(int ik,
                                                const int it,
                                                const int ia,
                                                const ModulePW::PW_Basis_K* wfc_basis,
                                                ModuleBase::Vector3<double> q) // pengfei 2016-11-23
{
    const int npw = wfc_basis->npwk[ik];
    std::complex<double> *skq = new std::complex<double>[npw];

    for (int ig = 0; ig < npw; ig++)
    {
        ModuleBase::Vector3<double> qkq = wfc_basis->getgpluskcar(ik, ig) + q;
        double arg = (qkq * GlobalC::ucell.atoms[it].tau[ia]) * ModuleBase::TWO_PI;
        skq[ig] = std::complex<double>(cos(arg), -sin(arg));
    }

    return skq;
}

template void Structure_Factor::get_sk<float, psi::DEVICE_CPU>(psi::DEVICE_CPU*,
                                                               int,
                                                               const ModulePW::PW_Basis_K*,
                                                               std::complex<float>*) const;
template void Structure_Factor::get_sk<double, psi::DEVICE_CPU>(psi::DEVICE_CPU*,
                                                                int,
                                                                const ModulePW::PW_Basis_K*,
                                                                std::complex<double>*) const;
#if defined(__CUDA) || defined(__ROCM)
template void Structure_Factor::get_sk<float, psi::DEVICE_GPU>(psi::DEVICE_GPU*,
                                                               int,
                                                               const ModulePW::PW_Basis_K*,
                                                               std::complex<float>*) const;
template void Structure_Factor::get_sk<double, psi::DEVICE_GPU>(psi::DEVICE_GPU*,
                                                                int,
                                                                const ModulePW::PW_Basis_K*,
                                                                std::complex<double>*) const;
#endif