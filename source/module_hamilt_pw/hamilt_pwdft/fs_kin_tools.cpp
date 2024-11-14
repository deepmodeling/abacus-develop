#include "fs_kin_tools.h"

#include "module_base/parallel_reduce.h"
namespace hamilt
{
template <typename FPTYPE, typename Device>
FS_Kin_tools<FPTYPE, Device>::FS_Kin_tools(const UnitCell& ucell_in,
                                           const K_Vectors* p_kv,
                                           const ModulePW::PW_Basis_K* wfc_basis_in,
                                           const ModuleBase::matrix& wg)
    : ucell_(ucell_in), nksbands_(wg.nc), wg(wg.c), wk(p_kv->wk.data())
{
    this->device = base_device::get_device_type<Device>(this->ctx);
    this->wfc_basis_ = wfc_basis_in;
    const int npwk_max = this->wfc_basis_->npwk_max;
    const int nks = this->wfc_basis_->nks;
    const int npol = ucell_in.get_npol();

    this->gk3_.resize(3 * npwk_max);
    this->gk.resize(3);
    for (int i = 0; i < 3; ++i)
    {
        this->gk[i] = &this->gk3_[i * npwk_max];
    }
    this->kfac.resize(npwk_max);
    this->s_kin.resize(9, 0.0);

    if (this->device == base_device::GpuDevice)
    {
        resmem_var_op()(this->ctx, d_gk, 3 * npwk_max);
        resmem_var_op()(this->ctx, d_kfac, npwk_max);
    }
    else
    {
        d_gk = gk3_.data();
        d_kfac = kfac.data();
    }
}

template <typename FPTYPE, typename Device>
FS_Kin_tools<FPTYPE, Device>::~FS_Kin_tools()
{
    if (this->device == base_device::GpuDevice)
    {
        delmem_var_op()(this->ctx, d_gk);
        delmem_var_op()(this->ctx, d_kfac);
    }
}

template <typename FPTYPE, typename Device>
void FS_Kin_tools<FPTYPE, Device>::cal_gk(const int& ik)
{
    const double tpiba = ModuleBase::TWO_PI / this->ucell_.lat0;
    const double twobysqrtpi = 2.0 / std::sqrt(ModuleBase::PI);
    const int npw = wfc_basis_->npwk[ik];
    const int npwk_max = wfc_basis_->npwk_max;
    for (int i = 0; i < npw; ++i)
    {
        gk[0][i] = wfc_basis_->getgpluskcar(ik, i)[0] * tpiba;
        gk[1][i] = wfc_basis_->getgpluskcar(ik, i)[1] * tpiba;
        gk[2][i] = wfc_basis_->getgpluskcar(ik, i)[2] * tpiba;
        if (wfc_basis_->erf_height > 0)
        {
            double gk2 = gk[0][i] * gk[0][i] + gk[1][i] * gk[1][i] + gk[2][i] * gk[2][i];
            double arg = (gk2 - wfc_basis_->erf_ecut) / wfc_basis_->erf_sigma;
            kfac[i] = 1.0 + wfc_basis_->erf_height / wfc_basis_->erf_sigma * twobysqrtpi * std::exp(-arg * arg);
        }
        else
        {
            kfac[i] = 1.0;
        }
    }
    if (this->device == base_device::GpuDevice)
    {
        syncmem_var_h2d_op()(this->ctx, this->cpu_ctx, d_gk, gk[0], 3 * npwk_max);
        syncmem_var_h2d_op()(this->ctx, this->cpu_ctx, d_kfac, kfac.data(), npwk_max);
    }
}

template <typename FPTYPE, typename Device>
void FS_Kin_tools<FPTYPE, Device>::cal_stress_kin(const int& ik,
                                                  const int& npm,
                                                  const bool& occ,
                                                  const std::complex<FPTYPE>* psi)
{
    if (npm == 0)
    {
        return;
    }
    const int npw = wfc_basis_->npwk[ik];
    const int npwk_max = wfc_basis_->npwk_max;
    const int npol = this->ucell_.get_npol();
    for (int ib = 0; ib < npm; ib++)
    {
        const std::complex<FPTYPE>* ppsi = psi + ib * npwk_max * npol;
        const std::complex<FPTYPE>* ppsi2 = ppsi + npwk_max;
        FPTYPE fac = 0.0;
        if (occ)
        {
            fac = wg[ik * this->nksbands_ + ib];
            if (fac == 0.0)
            {
                continue;
            }
        }
        else
        {
            fac = wk[ik];
        }
        for (int l = 0; l < 3; l++)
        {
            const FPTYPE* d_gkl = d_gk + l * npwk_max;
            for (int m = 0; m < l + 1; m++)
            {
                const FPTYPE* d_gkm = d_gk + m * npwk_max;
                FPTYPE sum = 0;
                sum += cal_multi_dot_op()(npw, fac, d_gkl, d_gkm, d_kfac, ppsi);
                if (npol == 2)
                {
                    sum += cal_multi_dot_op()(npw, fac, d_gkl, d_gkm, d_kfac, ppsi2);
                }
                s_kin[l * 3 + m] += sum;
            }
        }
    }
}

template <typename FPTYPE, typename Device>
void FS_Kin_tools<FPTYPE, Device>::symmetrize_stress(ModuleSymmetry::Symmetry* p_symm, ModuleBase::matrix& sigma)
{
    for (int l = 0; l < 3; ++l)
    {
        for (int m = 0; m < l; ++m)
        {
            s_kin[m * 3 + l] = s_kin[l * 3 + m];
        }
    }

    Parallel_Reduce::reduce_all(s_kin.data(), 9);

    for (int l = 0; l < 3; ++l)
    {
        for (int m = 0; m < 3; ++m)
        {
            sigma(l, m) = s_kin[l * 3 + m] * ModuleBase::e2 / this->ucell_.omega;
        }
    }
    // do symmetry
    if (ModuleSymmetry::Symmetry::symm_flag == 1)
    {
        p_symm->symmetrize_mat3(sigma, this->ucell_.lat);
    }
}

template class FS_Kin_tools<double, base_device::DEVICE_CPU>;
#if ((defined __CUDA) || (defined __ROCM))
template class FS_Kin_tools<double, base_device::DEVICE_GPU>;
#endif

} // namespace hamilt