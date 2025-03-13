#include "velocity_pw.h"

#include "module_base/kernels/math_kernel_op.h"
#include "module_base/parallel_reduce.h"
#include "module_base/timer.h"
namespace hamilt
{

template <typename FPTYPE, typename Device>
Velocity<FPTYPE, Device>::Velocity(const ModulePW::PW_Basis_K* wfcpw_in,
                                   const int* isk_in,
                                   pseudopot_cell_vnl* ppcell_in,
                                   const UnitCell* ucell_in,
                                   const bool nonlocal_in)
{
    if (wfcpw_in == nullptr || isk_in == nullptr || ppcell_in == nullptr || ucell_in == nullptr)
    {
        ModuleBase::WARNING_QUIT("Velocity", "Constuctor of Operator::Velocity is failed, please check your code!");
    }
    this->wfcpw = wfcpw_in;
    this->isk = isk_in;
    this->ppcell = ppcell_in;
    this->ucell = ucell_in;
    this->nonlocal = nonlocal_in;
    this->tpiba = ucell_in->tpiba;
    if (this->nonlocal)
    {
        this->ppcell->initgradq_vnl(*this->ucell);
    }
}

template <typename FPTYPE, typename Device>
Velocity<FPTYPE, Device>::~Velocity()
{
    delmem_var_op()(this->gx_);
    delmem_var_op()(this->gy_);
    delmem_var_op()(this->gz_);
    delmem_complex_op()(vkb_);
    delmem_complex_op()(gradvkb_);
    delmem_var_op()(deeq_);
}

template <typename FPTYPE, typename Device>
void Velocity<FPTYPE, Device>::init(const int ik_in)
{
    this->ik = ik_in;
    // init G+K
    const int npw = this->wfcpw->npwk[ik_in];
    const int npwk_max = this->wfcpw->npwk_max;
    std::vector<FPTYPE> gtmp(npw);
    std::vector<FPTYPE*> gtmp_ptr = {this->gx_, this->gy_, this->gz_};
    for(int i=0; i<3; ++i)
    {
        resmem_var_op()(gtmp_ptr[i], npw);
        for (int ig = 0; ig < npw; ++ig)
        {
            const ModuleBase::Vector3<double> tmpg = wfcpw->getgpluskcar(this->ik, ig);
            gtmp[ig] = static_cast<FPTYPE>(tmpg[i] * tpiba);
        }
        syncmem_var_h2d_op()(gtmp_ptr[i], gtmp.data(), npw);
    }

    // Calculate nonlocal pseudopotential vkb
    if (this->ppcell->nkb > 0 && this->nonlocal)
    {
        this->ppcell->getgradq_vnl(*this->ucell, ik_in);
    }

    // sync to device
    if(std::is_same<Device, base_device::DEVICE_GPU>::value || std::is_same<FPTYPE, float>::value)
    {
        const int nkb = this->ppcell->nkb;
        // vkb
        resmem_complex_op()(vkb_, nkb * npwk_max);
        for(int ib = 0; ib < nkb; ++ib)
        {
            castmem_complex_h2d_op()(vkb_ + ib * npw, this->ppcell->vkb.c + ib * npwk_max, npw);
        }

        // gradvkb
        resmem_complex_op()(gradvkb_, 3 * nkb * npwk_max);
        for(int ib = 0; ib < 3 * nkb; ++ib)
        {
            castmem_complex_h2d_op()(gradvkb_ + ib * npw, this->ppcell->gradvkb.ptr + ib * npwk_max, npw);
        }

        // deep
        int dim_deeq = 0;
        for (int it = 0; it < this->ucell->ntype; ++it)
        {
            const int Nprojs = this->ucell->atoms[it].ncpp.nh;
            dim_deeq += Nprojs * Nprojs * this->ucell->atoms[it].na;
        }
        resmem_var_op()(deeq_, dim_deeq);
        const int current_spin = this->isk[ik];
        int ipp = 0;
        int iat = 0;
        for (int it = 0; it < this->ucell->ntype; it++)
        {
            const int Nprojs = this->ucell->atoms[it].ncpp.nh;
            const int Nprojs2 = Nprojs * Nprojs;
            castmem_var_h2d_op()(deeq_ + ipp, &(this->ppcell->deeq(current_spin, iat, 0, 0)), Nprojs2);
            iat += this->ucell->atoms[it].na;
            ipp += Nprojs2;
        }
    }

}

template <typename FPTYPE, typename Device>
void Velocity<FPTYPE, Device>::act(const psi::Psi<std::complex<FPTYPE>, Device>* psi_in,
                                   const int n_npwx,
                                   const std::complex<FPTYPE>* psi0,
                                   std::complex<FPTYPE>* vpsi,
                                   const bool add) const
{
    ModuleBase::timer::tick("Operator", "Velocity");

    const int npw = psi_in->get_current_nbas();

    const int max_npw = psi_in->get_nbasis() / psi_in->get_npol();
    const int npol = psi_in->get_npol();
    
    std::vector<FPTYPE*> gtmp_ptr = {this->gx_, this->gy_, this->gz_};
    // -------------
    //       p
    // -------------
    for (int id = 0; id < 3; ++id)
    {
        const Complex* tmpsi_in = psi0;
        Complex* tmpvpsi = vpsi + id * n_npwx * max_npw;
        for (int ib = 0; ib < n_npwx; ++ib)
        {
            ModuleBase::vector_mul_vector_op<Complex, Device>()(npw, tmpvpsi, tmpsi_in, gtmp_ptr[id], add);
            tmpvpsi += max_npw;
            tmpsi_in += max_npw;
        }
    }

    // ---------------------------------------------
    // i[V_NL, r] = (\nabla_q+\nabla_q')V_{NL}(q,q')
    // |\beta><\beta|\psi>
    // ---------------------------------------------
    if (this->ppcell->nkb <= 0 || !this->nonlocal)
    {
        ModuleBase::timer::tick("Operator", "Velocity");
        return;
    }

    // 1. <\beta|\psi>
    Complex* becp1_ = nullptr; ///<[Device, n_npwx * nkb] <\beta|\psi>
    Complex* becp2_ = nullptr; ///<[Device, n_npwx * 3*nkb] <\nabla\beta|\psi>
    Complex* ps1_ = nullptr; ///<[Device, nkb * n_npwx] sum of becp1
    Complex* ps2_ = nullptr; ///<[Device, 3*nkb * n_npwx] sum of becp2
    resmem_complex_op()(ps1_, this->ppcell->nkb * n_npwx);
    resmem_complex_op()(ps2_, 3 * this->ppcell->nkb * n_npwx);
    resmem_complex_op()(becp1_, this->ppcell->nkb * n_npwx);
    resmem_complex_op()(becp2_, 3 * this->ppcell->nkb * n_npwx);

    const int nkb = this->ppcell->nkb;
    const int nkb3 = 3 * nkb;
    Complex one = 1.0;
    Complex zero = 0.0;
    if (n_npwx == 1)
    {
        int inc = 1;
        ModuleBase::gemv_op<Complex, Device>()('C', npw, nkb, &one, vkb_, max_npw, psi0, inc, &zero, becp1_, inc);
        ModuleBase::gemv_op<Complex, Device>()('C', npw, nkb3, &one, gradvkb_, max_npw, psi0, inc, &zero, becp2_, inc);
    }
    else
    {
        ModuleBase::gemm_op<Complex, Device>()('C',
                                               'N',
                                               nkb,
                                               n_npwx,
                                               npw,
                                               &one,
                                               vkb_,
                                               max_npw,
                                               psi0,
                                               max_npw,
                                               &zero,
                                               becp1_,
                                               nkb);
        ModuleBase::gemm_op<Complex, Device>()('C',
                                               'N',
                                               nkb3,
                                               n_npwx,
                                               npw,
                                               &one,
                                               gradvkb_,
                                               max_npw,
                                               psi0,
                                               max_npw,
                                               &zero,
                                               becp2_,
                                               nkb3);
    }

    Complex* becp1_cpu = nullptr;
    Complex* becp2_cpu = nullptr;
    Complex* ps1_cpu = nullptr;
    Complex* ps2_cpu = nullptr;
    std::vector<Complex> tmp_space1, tmp_space2;
    if(std::is_same<Device, base_device::DEVICE_GPU>::value)
    {
        tmp_space1.resize(nkb * n_npwx + nkb3 * n_npwx);
        becp1_cpu = tmp_space1.data();
        becp2_cpu = becp1_cpu + nkb * n_npwx;
        syncmem_complex_d2h_op()(becp1_cpu, becp1_, nkb * n_npwx);
        syncmem_complex_d2h_op()(becp2_cpu, becp2_, nkb3 * n_npwx);
        
        tmp_space2.resize(nkb * n_npwx + nkb3 * n_npwx, 0.0);
        ps1_cpu = tmp_space2.data();
        ps2_cpu = ps1_cpu + nkb * n_npwx;
    }
    else
    {
        Parallel_Reduce::reduce_pool(becp1_, nkb * n_npwx);
        Parallel_Reduce::reduce_pool(becp2_, nkb3 * n_npwx);
        becp1_cpu = becp1_;
        becp2_cpu = becp2_;
        ps1_cpu = ps1_;
        ps2_cpu = ps2_;
    }

    // 2. <\beta \psi><psi|
    int sum = 0;
    int iat = 0;
    if (npol == 1)
    {
        const int current_spin = 0;
        for (int it = 0; it < this->ucell->ntype; it++)
        {
            const int nproj = this->ucell->atoms[it].ncpp.nh;
            for (int ia = 0; ia < this->ucell->atoms[it].na; ia++)
            {
                for (int ip = 0; ip < nproj; ip++)
                {
                    for (int ip2 = 0; ip2 < nproj; ip2++)
                    {
                        for (int ib = 0; ib < n_npwx; ++ib)
                        {
                            double dij = this->ppcell->deeq(current_spin, iat, ip, ip2);
                            int sumip2 = sum + ip2;
                            int sumip = sum + ip;
                            ps1_cpu[sumip2 * n_npwx + ib] += dij * becp1_cpu[ib * nkb + sumip];
                            ps2_cpu[sumip2 * n_npwx + ib] += dij * becp2_cpu[ib * nkb3 + sumip];
                            ps2_cpu[(sumip2 + nkb) * n_npwx + ib] += dij * becp2_cpu[ib * nkb3 + sumip + nkb];
                            ps2_cpu[(sumip2 + 2 * nkb) * n_npwx + ib] += dij * becp2_cpu[ib * nkb3 + sumip + 2 * nkb];
                        }
                    }
                }
                sum += nproj;
                ++iat;
            }
        }
    }
    else
    {
        ModuleBase::WARNING_QUIT("Velocity", "Velocity operator does not support the non-collinear case yet!");
    }

    if(std::is_same<Device, base_device::DEVICE_GPU>::value)
    {
        syncmem_complex_h2d_op()(ps1_, ps1_cpu, nkb * n_npwx);
        syncmem_complex_h2d_op()(ps2_, ps2_cpu, nkb3 * n_npwx);
    }

    if (n_npwx == 1)
    {
        int inc = 1;
        for (int id = 0; id < 3; ++id)
        {
            int vkbshift = id * max_npw * nkb;
            int ps2shift = id * nkb;
            int npwshift = id * max_npw;
            ModuleBase::gemv_op<Complex, Device>()('N',
                                                   npw,
                                                   nkb,
                                                   &one,
                                                   gradvkb_ + vkbshift,
                                                   max_npw,
                                                   ps1_,
                                                   inc,
                                                   &one,
                                                   vpsi + npwshift,
                                                   inc);
            ModuleBase::gemv_op<Complex, Device>()('N',
                                                   npw,
                                                   nkb,
                                                   &one,
                                                   vkb_,
                                                   max_npw,
                                                   ps2_ + ps2shift,
                                                   inc,
                                                   &one,
                                                   vpsi + npwshift,
                                                   inc);
        }
    }
    else
    {
        for (int id = 0; id < 3; ++id)
        {
            int vkbshift = id * max_npw * nkb;
            int ps2shift = id * n_npwx * nkb;
            int npwshift = id * max_npw * n_npwx;
            ModuleBase::gemm_op<Complex, Device>()('N',
                                                   'T',
                                                   npw,
                                                   n_npwx,
                                                   nkb,
                                                   &one,
                                                   gradvkb_ + vkbshift,
                                                   max_npw,
                                                   ps1_,
                                                   n_npwx,
                                                   &one,
                                                   vpsi + npwshift,
                                                   max_npw);
            ModuleBase::gemm_op<Complex, Device>()('N',
                                                   'T',
                                                   npw,
                                                   n_npwx,
                                                   nkb,
                                                   &one,
                                                   vkb_,
                                                   max_npw,
                                                   ps2_ + ps2shift,
                                                   n_npwx,
                                                   &one,
                                                   vpsi + npwshift,
                                                   max_npw);
        }
    }
    delmem_complex_op()(ps1_);
    delmem_complex_op()(ps2_);
    delmem_complex_op()(becp1_);
    delmem_complex_op()(becp2_);
    ModuleBase::timer::tick("Operator", "Velocity");
    return;
}

template class Velocity<double, base_device::DEVICE_CPU>;
template class Velocity<float, base_device::DEVICE_CPU>;
#if ((defined __CUDA) || (defined __ROCM))
template class Velocity<double, base_device::DEVICE_GPU>;
template class Velocity<float, base_device::DEVICE_GPU>;
#endif

} // namespace hamilt