#include "op_pw_vel.h"

#include "source_base/kernels/math_kernel_op.h"
#include "source_base/parallel_reduce.h"
#include "source_base/timer.h"
#include "source_hamilt/module_xc/xc_functional.h"
#include "source_pw/module_pwdft/kernels/meta_op.h"
namespace hamilt
{

template <typename FPTYPE, typename Device>
Velocity<FPTYPE, Device>::Velocity(const ModulePW::PW_Basis_K* wfcpw_in,
                                   const int* isk_in,
                                   pseudopot_cell_vnl* ppcell_in,
                                   const UnitCell* ucell_in,
                                   const bool nonlocal_in,
                                   const typename GetTypeReal<FPTYPE>::type* vtau_in,
                                   const int vtau_col_in,
                                   const int vtau_row_in)
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
    this->vtau_ = vtau_in;
    this->vtau_col_ = vtau_col_in;
    this->vtau_row_ = vtau_row_in;
    if (this->nonlocal)
    {
        this->ppcell->ensure_grad_table(*this->ucell);
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
    delmem_complex_op()(porter1_);
    delmem_complex_op()(porter2_);
}

template <typename FPTYPE, typename Device>
void Velocity<FPTYPE, Device>::init(const int ik_in)
{
    this->init(ik_in, ModuleBase::Vector3<double>(0.0, 0.0, 0.0));
}

template <typename FPTYPE, typename Device>
void Velocity<FPTYPE, Device>::init(const int ik_in, const ModuleBase::Vector3<double>& vector_potential)
{
    ModuleBase::timer::start("Velocity", "init");
    this->ik = ik_in;
    this->tpiba = this->ucell->tpiba;
    // init G+K
    const int npw = this->wfcpw->npwk[ik_in];
    const int npwk_max = this->wfcpw->npwk_max;
    std::vector<FPTYPE> gtmp(npw);
    if (npw > momentum_capacity_)
    {
        resmem_var_op()(gx_, npw);
        resmem_var_op()(gy_, npw);
        resmem_var_op()(gz_, npw);
        momentum_capacity_ = npw;
    }
    std::vector<FPTYPE*> gtmp_ptr = {this->gx_, this->gy_, this->gz_};
    for (int i = 0; i < 3; ++i)
    {
        for (int ig = 0; ig < npw; ++ig)
        {
            const ModuleBase::Vector3<double> tmpg = wfcpw->getgpluskcar(this->ik, ig);
            gtmp[ig] = static_cast<FPTYPE>(tmpg[i] * tpiba + vector_potential[i]);
        }
        if (npw > 0)
        {
            syncmem_var_h2d_op()(gtmp_ptr[i], gtmp.data(), npw);
        }
    }

    // Calculate nonlocal pseudopotential vkb
    if (npw > 0 && this->ppcell->nkb > 0 && this->nonlocal)
    {

        const int nkb = this->ppcell->nkb;
        if (nkb * npwk_max > projector_capacity_)
        {
            resmem_complex_op()(vkb_, nkb * npwk_max);
            resmem_complex_op()(gradvkb_, 3 * nkb * npwk_max);
            projector_capacity_ = nkb * npwk_max;
        }
        this->ppcell->getvnl(this->ctx, *this->ucell, ik_in, vector_potential, vkb_);
        this->gradient_.calculate(this->ppcell, *this->ucell, *this->wfcpw, ik_in, vector_potential, gradvkb_);
    }
    ModuleBase::timer::end("Velocity", "init");
}

template <typename FPTYPE, typename Device>
void Velocity<FPTYPE, Device>::act(const psi::Psi<std::complex<FPTYPE>, Device>* psi_in,
                                   const int n_npwx,
                                   const std::complex<FPTYPE>* psi0,
                                   std::complex<FPTYPE>* vpsi,
                                   const bool add) const
{
    ModuleBase::timer::start("Operator", "Velocity");

    const int npw = this->wfcpw->npwk[this->ik];
    const int max_npw = this->wfcpw->npwk_max;
    const int npol = psi_in->get_npol();
    using Real = typename GetTypeReal<FPTYPE>::type;

    std::vector<FPTYPE*> gtmp_ptr = {this->gx_, this->gy_, this->gz_};
    // -------------
    //       p
    // -------------
    for (int id = 0; id < 3; ++id)
    {
        const Complex* tmpsi_in = psi0;
        Complex* tmpvpsi = vpsi + id * n_npwx * max_npw;
        for (int ib = 0; npw > 0 && ib < n_npwx; ++ib)
        {
            ModuleBase::vector_mul_vector_op<Complex, Device, FPTYPE>()(npw, tmpvpsi, tmpsi_in, gtmp_ptr[id], add);
            tmpvpsi += max_npw;
            tmpsi_in += max_npw;
        }
    }

    // ---------------------------------------------
    // meta-GGA velocity correction
    // V_tau = -(1/2) div(v_tau grad), whose plane-wave matrix element is
    // <k+G|V_tau|k+G'> = (1/2) v_tau(G-G') (k+G) dot (k+G').
    // Therefore
    // i[V_tau, r_\alpha]_{G,G'} =
    // (1/2) v_tau(G-G') [2k_alpha + G_alpha + G'_alpha].
    // In real space this is implemented as
    // -i/2 [\partial_\alpha(v_tau psi) + v_tau \partial_\alpha psi].
    // ---------------------------------------------
    // The row count also enables the correction on empty real-space partitions.
    if (this->vtau_row_ > 0 && XC_Functional::get_ked_flag())
    {
        if (this->vtau_col_ != this->wfcpw->nrxx || (this->vtau_col_ > 0 && this->vtau_ == nullptr))
        {
            ModuleBase::WARNING_QUIT("Velocity", "Invalid local potential for meta-GGA velocity correction.");
        }
        if (this->wfcpw->nmaxgr > porter_capacity_)
        {
            resmem_complex_op()(this->porter1_, this->wfcpw->nmaxgr);
            resmem_complex_op()(this->porter2_, this->wfcpw->nmaxgr);
            porter_capacity_ = this->wfcpw->nmaxgr;
        }
        int current_spin = 0;
        if (this->vtau_row_ > 1)
        {
            current_spin = this->isk[this->ik];
            if (current_spin < 0 || current_spin >= this->vtau_row_)
            {
                ModuleBase::WARNING_QUIT("Velocity", "invalid spin index for meta-GGA velocity correction");
            }
        }
        const Real* vtau_spin = this->vtau_col_ > 0 ? this->vtau_ + current_spin * this->vtau_col_ : nullptr;
        Complex minus_half_i(0.0, -0.5);
        for (int ib = 0; ib < n_npwx; ++ib)
        {
            const Complex* bandpsi = psi0 + ib * max_npw;
            this->wfcpw->recip_to_real(this->ctx, bandpsi, this->porter1_, this->ik);
            if (this->vtau_col_ > 0)
            {
                ModuleBase::vector_mul_vector_op<Complex, Device, FPTYPE>()(this->vtau_col_,
                                                                            this->porter1_,
                                                                            this->porter1_,
                                                                            vtau_spin,
                                                                            false);
            }
            this->wfcpw->real_to_recip(this->ctx, this->porter1_, this->porter1_, this->ik);
            for (int id = 0; id < 3; ++id)
            {
                Complex* vpsi_slice = vpsi + id * n_npwx * max_npw + ib * max_npw;
                Complex one = 1.0;
                // term1: partial_id (v_tau * psi)
                if (npw > 0)
                {
                    meta_pw_op<Real, Device>()(this->ctx,
                                               this->ik,
                                               id,
                                               npw,
                                               max_npw,
                                               this->tpiba,
                                               this->wfcpw->template get_gcar_data<Real>(),
                                               this->wfcpw->template get_kvec_c_data<Real>(),
                                               this->porter1_,
                                               this->porter2_,
                                               false);
                    ModuleBase::scal_op<Real, Device>()(npw, &minus_half_i, this->porter2_, 1);
                    ModuleBase::axpy_op<Complex, Device>()(npw, &one, this->porter2_, 1, vpsi_slice, 1);

                    // term2: v_tau * partial_id psi
                    meta_pw_op<Real, Device>()(this->ctx,
                                               this->ik,
                                               id,
                                               npw,
                                               max_npw,
                                               this->tpiba,
                                               this->wfcpw->template get_gcar_data<Real>(),
                                               this->wfcpw->template get_kvec_c_data<Real>(),
                                               bandpsi,
                                               this->porter2_,
                                               false);
                }
                this->wfcpw->recip_to_real(this->ctx, this->porter2_, this->porter2_, this->ik);
                if (this->vtau_col_ > 0)
                {
                    ModuleBase::vector_mul_vector_op<Complex, Device, FPTYPE>()(this->vtau_col_,
                                                                                this->porter2_,
                                                                                this->porter2_,
                                                                                vtau_spin,
                                                                                false);
                }
                this->wfcpw->real_to_recip(this->ctx, this->porter2_, this->porter2_, this->ik);
                if (npw > 0)
                {
                    ModuleBase::scal_op<Real, Device>()(npw, &minus_half_i, this->porter2_, 1);
                    ModuleBase::axpy_op<Complex, Device>()(npw, &one, this->porter2_, 1, vpsi_slice, 1);
                }
            }
        }
    }

    // ---------------------------------------------
    // i[V_NL, r] = (\nabla_q+\nabla_q')V_{NL}(q,q')
    // |\beta><\beta|\psi>
    // ---------------------------------------------
    if (this->ppcell->nkb <= 0 || !this->nonlocal)
    {
        ModuleBase::timer::end("Operator", "Velocity");
        return;
    }

    // 1. <\beta|\psi>
    const int block = this->ppcell->nkb * n_npwx;
    Complex* becp1_ = this->contraction_.prepare(block);
    Complex* becp2_ = becp1_ + block;
    Complex* ps1_ = becp1_ + 4 * block;
    Complex* ps2_ = ps1_ + block;

    const int nkb = this->ppcell->nkb;
    const int nkb3 = 3 * nkb;
    Complex one = 1.0;
    Complex zero = 0.0;

    const Complex* vkb_d = vkb_;
    const Complex* gradvkb_d = gradvkb_;

    if (npw == 0)
    {
        // Keep the reduction collective even when this rank has no local plane waves.
        base_device::memory::set_memory_op<Complex, Device>()(becp1_, 0, 4 * block);
    }
    else if (n_npwx == 1)
    {
        int inc = 1;
        ModuleBase::gemv_op<Complex, Device>()('C', npw, nkb, &one, vkb_d, max_npw, psi0, inc, &zero, becp1_, inc);
        ModuleBase::gemv_op<Complex, Device>()('C', npw, nkb3, &one, gradvkb_d, max_npw, psi0, inc, &zero, becp2_, inc);
    }
    else
    {
        ModuleBase::gemm_op<Complex, Device>()('C', 'N', nkb, n_npwx, npw, &one, vkb_d, max_npw, psi0, max_npw, &zero, becp1_, nkb);
        ModuleBase::gemm_op<Complex, Device>()('C', 'N', nkb3, n_npwx, npw, &one, gradvkb_d, max_npw, psi0, max_npw, &zero, becp2_, nkb3);
    }

    if (npol != 1)
    {
        ModuleBase::WARNING_QUIT("Velocity", "Non-collinear velocity is not supported.");
    }
    // deeq stores energies in Ry; velocity uses Hartree atomic units like p.
    this->contraction_.contract(*this->ucell, *this->ppcell, this->isk[this->ik], n_npwx, FPTYPE(0.5), *this->wfcpw, becp1_);

    if (npw == 0)
    {
        ModuleBase::timer::end("Operator", "Velocity");
        return;
    }

    if (n_npwx == 1)
    {
        int inc = 1;
        for (int id = 0; id < 3; ++id)
        {
            int vkbshift = id * max_npw * nkb;
            int ps2shift = id * nkb;
            int npwshift = id * max_npw;
            ModuleBase::gemv_op<Complex,
                                Device>()('N', npw, nkb, &one, gradvkb_d + vkbshift, max_npw, ps1_, inc, &one, vpsi + npwshift, inc);
            ModuleBase::gemv_op<Complex, Device>()('N', npw, nkb, &one, vkb_d, max_npw, ps2_ + ps2shift, inc, &one, vpsi + npwshift, inc);
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
                                                   gradvkb_d + vkbshift,
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
                                                   vkb_d,
                                                   max_npw,
                                                   ps2_ + ps2shift,
                                                   n_npwx,
                                                   &one,
                                                   vpsi + npwshift,
                                                   max_npw);
        }
    }

    ModuleBase::timer::end("Operator", "Velocity");
    return;
}

template class Velocity<double, base_device::DEVICE_CPU>;
template class Velocity<float, base_device::DEVICE_CPU>;
#if ((defined __CUDA) || (defined __ROCM))
template class Velocity<double, base_device::DEVICE_GPU>;
template class Velocity<float, base_device::DEVICE_GPU>;
#endif

} // namespace hamilt
