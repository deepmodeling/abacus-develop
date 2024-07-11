#ifndef HAMILTPW_NONLOCAL_MATHS_H
#define HAMILTPW_NONLOCAL_MATHS_H

#include "module_base/module_device/device.h"
#include "module_basis/module_pw/pw_basis_k.h"
#include "module_cell/klist.h"
#include "module_cell/unitcell.h"
#include "module_hamilt_pw/hamilt_pwdft/VNL_in_pw.h"
#include "module_hamilt_pw/hamilt_pwdft/kernels/stress_op.h"
#include "module_hsolver/kernels/math_kernel_op.h"

namespace hamilt
{

template <typename FPTYPE, typename Device>
class Nonlocal_maths
{
  public:
    Nonlocal_maths(const pseudopot_cell_vnl* nlpp_in, const UnitCell* ucell_in)
    {
        this->device = base_device::get_device_type<Device>(this->ctx);
        this->nlpp_ = nlpp_in;
        this->ucell_ = ucell_in;
    }

  private:
    const pseudopot_cell_vnl* nlpp_;
    const UnitCell* ucell_;

    Device* ctx = {};
    base_device::DEVICE_CPU* cpu_ctx = {};
    base_device::AbacusDevice_t device = {};

  public:
    // functions
    /// calculate the G+K vectors
    std::vector<FPTYPE> cal_gk(int ik, const ModulePW::PW_Basis_K* wfc_basis);
    /// calculate the sperical bessel function for projections
    void cal_ylm(int lmax, int npw, const FPTYPE* gk_in, FPTYPE* ylm);
    /// calculate the derivate of the sperical bessel function for projections
    void cal_ylm_deri(int lmax, int npw, const FPTYPE* gk_in, FPTYPE* ylm_deri);
    /// calculate the (-i)^l factors
    std::vector<complex<FPTYPE>> cal_pref(int it);
    /// calculate the vkb matrix for this atom
    /// vkb = sum_lm (-i)^l * ylm(g^) * vq(g^) * sk(g^)
    void cal_vkb(int it,
                 int ia,
                 int npw,
                 const FPTYPE* vq_in,
                 const FPTYPE* ylm_in,
                 const complex<FPTYPE>* sk_in,
                 const complex<FPTYPE>* pref_in,
                 complex<FPTYPE>* vkb_out);
    /// calculate the dvkb matrix for this atom
    void cal_vkb_deri(int it,
                      int ia,
                      int npw,
                      int ipol,
                      int jpol,
                      const FPTYPE* vq_in,
                      const FPTYPE* vq_deri_in,
                      const FPTYPE* ylm_in,
                      const FPTYPE* ylm_deri_in,
                      const complex<FPTYPE>* sk_in,
                      const complex<FPTYPE>* pref_in,
                      const FPTYPE* gk_in,
                      complex<FPTYPE>* vkb_out);

    /// calculate the ptr used in vkb_op
    void prepare_vkb_ptr(int nbeta,
                         double* nhtol,
                         int nhtol_nc,
                         int npw,
                         int it,
                         std::complex<FPTYPE>* vkb_out,
                         std::complex<FPTYPE>** vkb_ptrs,
                         FPTYPE* ylm_in,
                         FPTYPE** ylm_ptrs,
                         FPTYPE* vq_in,
                         FPTYPE** vq_ptrs);

    /// calculate the indexes used in vkb_deri_op
    /// indexes save (lm, nb, dylm_lm_ipol, dylm_lm_jpol) for nh
    void cal_dvkb_index(const int nbeta,
                        const double* nhtol,
                        const int nhtol_nc,
                        const int npw,
                        const int it,
                        const int ipol,
                        const int jpol,
                        int* indexes);

    static void dylmr2(const int nylm, const int ngy, const FPTYPE* gk, FPTYPE* dylm, const int ipol);
    /// polynomial interpolation tool for calculate derivate of vq
    static FPTYPE Polynomial_Interpolation_nl(const ModuleBase::realArray& table,
                                              const int& dim1,
                                              const int& dim2,
                                              const FPTYPE& table_interval,
                                              const FPTYPE& x);
};

// cal_gk
template <typename FPTYPE, typename Device>
std::vector<FPTYPE> Nonlocal_maths<FPTYPE, Device>::cal_gk(int ik, const ModulePW::PW_Basis_K* wfc_basis)
{
    int npw = wfc_basis->npwk[ik];
    std::vector<FPTYPE> gk(npw * 5);
    ModuleBase::Vector3<FPTYPE> tmp;
    for (int ig = 0; ig < npw; ++ig)
    {
        tmp = wfc_basis->getgpluskcar(ik, ig);
        gk[ig * 3] = tmp.x;
        gk[ig * 3 + 1] = tmp.y;
        gk[ig * 3 + 2] = tmp.z;
        FPTYPE norm = sqrt(tmp.norm2());
        gk[3 * npw + ig] = norm * this->ucell_->tpiba;
        gk[4 * npw + ig] = norm < 1e-8 ? 0.0 : 1.0 / norm * this->ucell_->tpiba;
    }
    return gk;
}

// cal_ylm
template <typename FPTYPE, typename Device>
void Nonlocal_maths<FPTYPE, Device>::cal_ylm(int lmax, int npw, const FPTYPE* gk_in, FPTYPE* ylm)
{

    const int x1 = (lmax + 1) * (lmax + 1);

    if (this->device == base_device::GpuDevice)
    {
        using syncmem_var_h2d_op = base_device::memory::synchronize_memory_op<FPTYPE, Device, base_device::DEVICE_CPU>;
        std::vector<FPTYPE> ylm_cpu(x1 * npw);
        ModuleBase::YlmReal::Ylm_Real(cpu_ctx, x1, npw, gk_in, ylm_cpu.data());
        syncmem_var_h2d_op()(this->ctx, this->cpu_ctx, ylm, ylm_cpu.data(), ylm_cpu.size());
    }
    else
    {
        ModuleBase::YlmReal::Ylm_Real(cpu_ctx, x1, npw, gk_in, ylm);
    }

    return;
}
// cal_ylm_deri
template <typename FPTYPE, typename Device>
void Nonlocal_maths<FPTYPE, Device>::cal_ylm_deri(int lmax, int npw, const FPTYPE* gk_in, FPTYPE* ylm_deri)
{
    const int x1 = (lmax + 1) * (lmax + 1);

    if (this->device == base_device::GpuDevice)
    {
        std::vector<FPTYPE> dylm(3 * x1 * npw);
        for (int ipol = 0; ipol < 3; ipol++)
        {
            Nonlocal_maths<FPTYPE, Device>::dylmr2(x1, npw, gk_in, &dylm[ipol * x1 * npw], ipol);
        }
        using syncmem_var_h2d_op = base_device::memory::synchronize_memory_op<FPTYPE, Device, base_device::DEVICE_CPU>;
        syncmem_var_h2d_op()(this->ctx, this->cpu_ctx, ylm_deri, dylm.data(), dylm.size());
    }
    else
    {
        for (int ipol = 0; ipol < 3; ipol++)
        {
            Nonlocal_maths<FPTYPE, Device>::dylmr2(x1, npw, gk_in, &ylm_deri[ipol * x1 * npw], ipol);
        }
    }

    return;
}
// cal_pref
template <typename FPTYPE, typename Device>
std::vector<std::complex<FPTYPE>> Nonlocal_maths<FPTYPE, Device>::cal_pref(int it)
{
    const int nh = this->ucell_->atoms[it].ncpp.nh;
    std::vector<std::complex<FPTYPE>> pref(nh);
    for (int ih = 0; ih < nh; ih++)
    {
        pref[ih] = std::pow(std::complex<FPTYPE>(0.0, -1.0), this->nlpp_->nhtol(it, ih));
    }
    return pref;
}

// cal_vkb
// cpu version first, gpu version later
template <typename FPTYPE, typename Device>
void Nonlocal_maths<FPTYPE, Device>::cal_vkb(int it,
                                             int ia,
                                             int npw,
                                             const FPTYPE* vq_in,
                                             const FPTYPE* ylm_in,
                                             const std::complex<FPTYPE>* sk_in,
                                             const std::complex<FPTYPE>* pref_in,
                                             std::complex<FPTYPE>* vkb_out)
{
    int ih = 0;
    // loop over all beta functions
    for (int nb = 0; nb < this->ucell_->atoms[it].ncpp.nbeta; nb++)
    {
        int l = this->nlpp_->nhtol(it, ih);
        // loop over all m angular momentum
        for (int m = 0; m < 2 * l + 1; m++)
        {
            int lm = l * l + m;
            std::complex<FPTYPE>* vkb_ptr = &vkb_out[ih * npw];
            const FPTYPE* ylm_ptr = &ylm_in[lm * npw];
            const FPTYPE* vq_ptr = &vq_in[nb * npw];
            // loop over all G-vectors
            for (int ig = 0; ig < npw; ig++)
            {
                vkb_ptr[ig] = ylm_ptr[ig] * vq_ptr[ig] * sk_in[ig] * pref_in[ih];
            }
            ih++;
        }
    }
}

// cal_vkb
// cpu version first, gpu version later
template <typename FPTYPE, typename Device>
void Nonlocal_maths<FPTYPE, Device>::cal_vkb_deri(int it,
                                                  int ia,
                                                  int npw,
                                                  int ipol,
                                                  int jpol,
                                                  const FPTYPE* vq_in,
                                                  const FPTYPE* vq_deri_in,
                                                  const FPTYPE* ylm_in,
                                                  const FPTYPE* ylm_deri_in,
                                                  const std::complex<FPTYPE>* sk_in,
                                                  const std::complex<FPTYPE>* pref_in,
                                                  const FPTYPE* gk_in,
                                                  std::complex<FPTYPE>* vkb_out)
{
    const int x1 = (this->nlpp_->lmaxkb + 1) * (this->nlpp_->lmaxkb + 1);
    int ih = 0;
    // loop over all beta functions
    for (int nb = 0; nb < this->ucell_->atoms[it].ncpp.nbeta; nb++)
    {
        const int l = this->nlpp_->nhtol(it, ih);
        // loop over all m angular momentum
        for (int m = 0; m < 2 * l + 1; m++)
        {
            const int lm = l * l + m;
            std::complex<FPTYPE>* vkb_ptr = &vkb_out[ih * npw];
            const FPTYPE* ylm_ptr = &ylm_in[lm * npw];
            const FPTYPE* vq_ptr = &vq_in[nb * npw];
            // set vkb to zero
            for (int ig = 0; ig < npw; ig++)
            {
                vkb_ptr[ig] = std::complex<FPTYPE>(0.0, 0.0);
            }
            // first term: ylm * vq * sk * pref
            // loop over all G-vectors
            if (ipol == jpol)
            {
                for (int ig = 0; ig < npw; ig++)
                {
                    vkb_ptr[ig] -= ylm_ptr[ig] * vq_ptr[ig] * sk_in[ig] * pref_in[ih];
                }
            }
            // second term: ylm_deri * vq_deri * sk * pref
            //  loop over all G-vectors
            const FPTYPE* ylm_deri_ptr1 = &ylm_deri_in[(ipol * x1 + lm) * npw];
            const FPTYPE* ylm_deri_ptr2 = &ylm_deri_in[(jpol * x1 + lm) * npw];
            const FPTYPE* vq_deri_ptr = &vq_deri_in[nb * npw];
            const FPTYPE* gkn = &gk_in[4 * npw];
            for (int ig = 0; ig < npw; ig++)
            {
                vkb_ptr[ig] -= (gk_in[ig * 3 + ipol] * ylm_deri_ptr2[ig] + gk_in[ig * 3 + jpol] * ylm_deri_ptr1[ig])
                               * vq_ptr[ig] * sk_in[ig] * pref_in[ih];
            }
            // third term: ylm * vq_deri * sk * pref
            //  loop over all G-vectors
            for (int ig = 0; ig < npw; ig++)
            {
                vkb_ptr[ig] -= 2.0 * ylm_ptr[ig] * vq_deri_ptr[ig] * sk_in[ig] * pref_in[ih] * gk_in[ig * 3 + ipol]
                               * gk_in[ig * 3 + jpol] * gkn[ig];
            }
            ih++;
        }
    }
}

template <typename FPTYPE, typename Device>
void Nonlocal_maths<FPTYPE, Device>::prepare_vkb_ptr(int nbeta,
                                                     double* nhtol,
                                                     int nhtol_nc,
                                                     int npw,
                                                     int it,
                                                     std::complex<FPTYPE>* vkb_out,
                                                     std::complex<FPTYPE>** vkb_ptrs,
                                                     FPTYPE* ylm_in,
                                                     FPTYPE** ylm_ptrs,
                                                     FPTYPE* vq_in,
                                                     FPTYPE** vq_ptrs)
{
    // std::complex<FPTYPE>** vkb_ptrs[nh];
    // const FPTYPE** ylm_ptrs[nh];
    // const FPTYPE** vq_ptrs[nh];
    int ih = 0;
    for (int nb = 0; nb < nbeta; nb++)
    {
        int l = nhtol[it * nhtol_nc + ih];
        for (int m = 0; m < 2 * l + 1; m++)
        {
            int lm = l * l + m;
            vkb_ptrs[ih] = &vkb_out[ih * npw];
            ylm_ptrs[ih] = &ylm_in[lm * npw];
            vq_ptrs[ih] = &vq_in[nb * npw];
            ih++;
        }
    }
}

template <typename FPTYPE, typename Device>
void Nonlocal_maths<FPTYPE, Device>::cal_dvkb_index(const int nbeta,
                                                    const double* nhtol,
                                                    const int nhtol_nc,
                                                    const int npw,
                                                    const int it,
                                                    const int ipol,
                                                    const int jpol,
                                                    int* indexes)
{
    int ih = 0;
    const int x1 = (this->nlpp_->lmaxkb + 1) * (this->nlpp_->lmaxkb + 1);
    for (int nb = 0; nb < nbeta; nb++)
    {
        int l = nhtol[it * nhtol_nc + ih];
        for (int m = 0; m < 2 * l + 1; m++)
        {
            int lm = l * l + m;
            indexes[ih * 4] = lm;
            indexes[ih * 4 + 1] = nb;
            indexes[ih * 4 + 2] = (ipol * x1 + lm);
            indexes[ih * 4 + 3] = (jpol * x1 + lm);

            ih++;
        }
    }
}

template <typename FPTYPE, typename Device>
void Nonlocal_maths<FPTYPE, Device>::dylmr2(const int nylm,
                                            const int ngy,
                                            const FPTYPE* gk,
                                            FPTYPE* dylm,
                                            const int ipol)
{
    //-----------------------------------------------------------------------
    //
    //     compute \partial Y_lm(G) \over \partial (G)_ipol
    //     using simple numerical derivation (SdG)
    //     The spherical harmonics are calculated in ylmr2
    //
    // int nylm, ngy, ipol;
    // number of spherical harmonics
    // the number of g vectors to compute
    // desired polarization
    // FPTYPE g (3, ngy), gg (ngy), dylm (ngy, nylm)
    // the coordinates of g vectors
    // the moduli of g vectors
    // the spherical harmonics derivatives
    //
    const FPTYPE delta = 1e-6;
    const FPTYPE small = 1e-15;

    ModuleBase::matrix ylmaux;
    // dg is the finite increment for numerical derivation:
    // dg = delta |G| = delta * sqrt(gg)
    // dgi= 1 /(delta * sqrt(gg))
    // gx = g +/- dg

    std::vector<FPTYPE> gx(ngy * 3);

    std::vector<FPTYPE> dg(ngy);
    std::vector<FPTYPE> dgi(ngy);

    ylmaux.create(nylm, ngy);

    ModuleBase::GlobalFunc::ZEROS(dylm, nylm * ngy);
    ylmaux.zero_out();

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int ig = 0; ig < 3 * ngy; ig++)
    {
        gx[ig] = gk[ig];
    }
    //$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ig)
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int ig = 0; ig < ngy; ig++)
    {
        const int igx = ig * 3, igy = ig * 3 + 1, igz = ig * 3 + 2;
        FPTYPE norm2 = gx[igx] * gx[igx] + gx[igy] * gx[igy] + gx[igz] * gx[igz];
        dg[ig] = delta * sqrt(norm2);
        if (dg[ig] > small)
        {
            dgi[ig] = 1.0 / dg[ig];
        }
        else
        {
            dgi[ig] = 0.0;
        }
    }
    //$OMP END PARALLEL DO

    //$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ig)
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int ig = 0; ig < ngy; ig++)
    {
        const int index = ig * 3 + ipol;
        gx[index] = gk[index] + dg[ig];
    }
    //$OMP END PARALLEL DO

    base_device::DEVICE_CPU* cpu = {};
    ModuleBase::YlmReal::Ylm_Real(cpu, nylm, ngy, gx.data(), dylm);
    //$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ig)
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int ig = 0; ig < ngy; ig++)
    {
        const int index = ig * 3 + ipol;
        gx[index] = gk[index] - dg[ig];
    }
    //$OMP END PARALLEL DO

    ModuleBase::YlmReal::Ylm_Real(cpu, nylm, ngy, gx.data(), ylmaux.c);

    //  zaxpy ( - 1.0, ylmaux, 1, dylm, 1);
#ifdef _OPENMP
#pragma omp parallel for collapse(2)
#endif
    for (int lm = 0; lm < nylm; lm++)
    {
        for (int ig = 0; ig < ngy; ig++)
        {
            dylm[lm * ngy + ig] -= ylmaux(lm, ig);
            dylm[lm * ngy + ig] *= 0.5 * dgi[ig];
        }
    }

    return;
}

template <typename FPTYPE, typename Device>
FPTYPE Nonlocal_maths<FPTYPE, Device>::Polynomial_Interpolation_nl(const ModuleBase::realArray& table,
                                                                   const int& dim1,
                                                                   const int& dim2,
                                                                   const FPTYPE& table_interval,
                                                                   const FPTYPE& x // input value
)
{

    assert(table_interval > 0.0);
    const FPTYPE position = x / table_interval;
    const int iq = static_cast<int>(position);

    const FPTYPE x0 = position - static_cast<FPTYPE>(iq);
    const FPTYPE x1 = 1.0 - x0;
    const FPTYPE x2 = 2.0 - x0;
    const FPTYPE x3 = 3.0 - x0;
    const FPTYPE y = (table(dim1, dim2, iq) * (-x2 * x3 - x1 * x3 - x1 * x2) / 6.0
                      + table(dim1, dim2, iq + 1) * (+x2 * x3 - x0 * x3 - x0 * x2) / 2.0
                      - table(dim1, dim2, iq + 2) * (+x1 * x3 - x0 * x3 - x0 * x1) / 2.0
                      + table(dim1, dim2, iq + 3) * (+x1 * x2 - x0 * x2 - x0 * x1) / 6.0)
                     / table_interval;

    return y;
}
} // namespace hamilt

#endif