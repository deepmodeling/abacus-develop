#include "module_base/math_polyint.h"
#include "module_parameter/parameter.h"
#include "module_base/math_ylmreal.h"
#include "module_base/memory.h"
#include "module_base/module_device/device.h"
#include "module_base/timer.h"
#include "module_hamilt_pw/hamilt_pwdft/fs_nonlocal_tools.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_hamilt_pw/hamilt_pwdft/nonlocal_maths.hpp"
#include "stress_func.h"
// calculate the nonlocal pseudopotential stress in PW
template <typename FPTYPE, typename Device>
void Stress_Func<FPTYPE, Device>::stress_nl(ModuleBase::matrix& sigma,
                                            const ModuleBase::matrix& wg,
                                            const ModuleBase::matrix& ekb,
                                            Structure_Factor* p_sf,
                                            K_Vectors* p_kv,
                                            ModuleSymmetry::Symmetry* p_symm,
                                            ModulePW::PW_Basis_K* wfc_basis,
                                            const psi::Psi<complex<FPTYPE>, Device>* psi_in,
                                            pseudopot_cell_vnl* nlpp_in,
                                            const UnitCell& ucell_in)
{
    ModuleBase::TITLE("Stress_Func", "stress_nl");
    // skip nkb==0
    if (nlpp_in->nkb == 0 || psi_in == nullptr || wfc_basis == nullptr)
    {
        return;
    }
    ModuleBase::timer::tick("Stress_Func", "stress_nl");

    FPTYPE* stress_device = nullptr;
    resmem_var_op()(this->ctx, stress_device, 9);
    setmem_var_op()(this->ctx, stress_device, 0, 9);
    std::vector<FPTYPE> sigmanlc(9, 0.0);

    hamilt::FS_Nonlocal_tools<FPTYPE, Device> nl_tools(nlpp_in, &ucell_in, p_kv, wfc_basis, p_sf, wg, &ekb);

    const int nks = p_kv->get_nks();
    const int max_nbands = wg.nc;
    for (int ik = 0; ik < nks; ik++) // loop k points
    {
        // skip zero weights to speed up
        int nbands_occ = wg.nc;
        while (wg(ik, nbands_occ - 1) == 0.0)
        {
            nbands_occ--;
            if (nbands_occ == 0)
            {
                break;
            }
        }
        const int npm = nbands_occ;

        nl_tools.cal_vkb(ik, max_nbands);
        // calculate becp = <psi|beta> for all beta functions
        nl_tools.cal_becp(ik, npm, &psi_in[0](ik,0,0));
        nl_tools.reduce_pool_becp(max_nbands);
        // calculate dbecp = <psi|d(beta)/dR> for all beta functions
        // calculate stress = \sum <psi|d(beta_j)/dR> * <psi|beta_i> * D_{ij}
        for (int ipol = 0; ipol < 3; ipol++)
        {
            for (int jpol = 0; jpol <= ipol; jpol++)
            {
                nl_tools.cal_vkb_deri_s(ik, max_nbands, ipol, jpol);
                nl_tools.cal_dbecp_s(ik, npm, &psi_in[0](ik,0,0));
                nl_tools.cal_stress(ik, npm, true, ipol, jpol, stress_device);
            }
        }
    }
    // transfer stress from device to host
    syncmem_var_d2h_op()(this->cpu_ctx, this->ctx, sigmanlc.data(), stress_device, 9);
    delmem_var_op()(this->ctx, stress_device);
    // sum up forcenl from all processors
    for (int l = 0; l < 3; l++)
    {
        for (int m = 0; m < 3; m++)
        {
            if (m > l)
            {
                sigmanlc[l * 3 + m] = sigmanlc[m * 3 + l];
            }
        }
    }
    // sum up forcenl from all processors
    Parallel_Reduce::reduce_all(sigmanlc.data(), 9);

    for (int ipol = 0; ipol < 3; ipol++)
    {
        for (int jpol = 0; jpol < 3; jpol++)
        {
            sigma(ipol, jpol) = sigmanlc[ipol * 3 + jpol] / ucell_in.omega;
        }
    }
    // do symmetry
    if (ModuleSymmetry::Symmetry::symm_flag == 1)
    {
        p_symm->symmetrize_mat3(sigma, ucell_in.lat);
    } // end symmetry

    ModuleBase::timer::tick("Stress_Func", "stress_nl");
}

template <typename FPTYPE, typename Device>
void Stress_Func<FPTYPE, Device>::get_dvnl1(ModuleBase::ComplexMatrix& vkb,
                                            const int ik,
                                            const int ipol,
                                            Structure_Factor* p_sf,
                                            ModulePW::PW_Basis_K* wfc_basis)
{
    if (PARAM.inp.test_pp) {
        ModuleBase::TITLE("Stress_Func", "get_dvnl1");
}

    const int npw = wfc_basis->npwk[ik];
    const int lmaxkb = nlpp->lmaxkb;
    if (lmaxkb < 0)
    {
        return;
    }

    const int nhm = nlpp->nhm;
    ModuleBase::matrix vkb1(nhm, npw);
    vkb1.zero_out();
    FPTYPE* vq = new FPTYPE[npw];
    const int x1 = (lmaxkb + 1) * (lmaxkb + 1);

    ModuleBase::matrix dylm(x1, npw);
    ModuleBase::Vector3<FPTYPE>* gk = new ModuleBase::Vector3<FPTYPE>[npw];
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int ig = 0; ig < npw; ig++)
    {
        gk[ig] = wfc_basis->getgpluskcar(ik, ig);
    }

    hamilt::Nonlocal_maths<FPTYPE, Device>::dylmr2(x1, npw, reinterpret_cast<double*>(gk), dylm.c, ipol);

    const int imag_pow_period = 4;
    // result table of pow(0-1i, int)
    static const std::complex<FPTYPE> pref_tab[imag_pow_period] = {{1, 0}, {0, -1}, {-1, 0}, {0, 1}};
    int jkb = 0;
    for (int it = 0; it < this->ucell->ntype; it++)
    {
        if (PARAM.inp.test_pp > 1) {
            ModuleBase::GlobalFunc::OUT("it", it);
}
        // calculate beta in G-space using an interpolation table
        const int nbeta = this->ucell->atoms[it].ncpp.nbeta;
        const int nh = this->ucell->atoms[it].ncpp.nh;

        if (PARAM.inp.test_pp > 1) {
            ModuleBase::GlobalFunc::OUT("nbeta", nbeta);
}

        for (int nb = 0; nb < nbeta; nb++)
        {
            if (PARAM.inp.test_pp > 1) {
                ModuleBase::GlobalFunc::OUT("ib", nb);
}
#ifdef _OPENMP
#pragma omp parallel for
#endif
            for (int ig = 0; ig < npw; ig++)
            {
                const FPTYPE gnorm = gk[ig].norm() * this->ucell->tpiba;

                // cout << "\n gk[ig] = " << gk[ig].x << " " << gk[ig].y << " " << gk[ig].z;
                // cout << "\n gk.norm = " << gnorm;

                vq[ig] = ModuleBase::PolyInt::Polynomial_Interpolation(nlpp->tab,
                                                                       it,
                                                                       nb,
                                                                       PARAM.globalv.nqx,
                                                                       PARAM.globalv.dq,
                                                                       gnorm);

            } // enddo

            // add spherical harmonic part
            for (int ih = 0; ih < nh; ih++)
            {
                if (nb == nlpp->indv(it, ih))
                {
                    const int lm = static_cast<int>(nlpp->nhtolm(it, ih));
#ifdef _OPENMP
#pragma omp parallel for
#endif
                    for (int ig = 0; ig < npw; ig++)
                    {
                        vkb1(ih, ig) = dylm(lm, ig) * vq[ig];
                    }
                }

            } // end ih

        } // end nbeta

        // vkb1 contains all betas including angular part for type nt
        // now add the structure factor and factor (-i)^l
        for (int ia = 0; ia < this->ucell->atoms[it].na; ia++)
        {
            std::complex<FPTYPE>* sk = p_sf->get_sk(ik, it, ia, wfc_basis);
#ifdef _OPENMP
#pragma omp parallel for collapse(2)
#endif
            for (int ih = 0; ih < nh; ih++)
            {
                for (int ig = 0; ig < npw; ig++)
                {
                    std::complex<FPTYPE> pref = pref_tab[int(nlpp->nhtol(it, ih)) % imag_pow_period]; //?
                    vkb(jkb + ih, ig) = vkb1(ih, ig) * sk[ig] * pref;
                }

            } // end ih
            jkb += nh;
            delete[] sk;
        } // end ia
    }     // enddo
    delete[] gk;
    delete[] vq;
    return;
} // end get_dvnl1

template <typename FPTYPE, typename Device>
void Stress_Func<FPTYPE, Device>::get_dvnl2(ModuleBase::ComplexMatrix& vkb,
                                            const int ik,
                                            Structure_Factor* p_sf,
                                            ModulePW::PW_Basis_K* wfc_basis)
{
    if (PARAM.inp.test_pp) {
        ModuleBase::TITLE("Stress", "get_dvnl2");
}
    //	ModuleBase::timer::tick("Stress","get_dvnl2");
    const int npw = wfc_basis->npwk[ik];
    const int lmaxkb = nlpp->lmaxkb;
    if (lmaxkb < 0)
    {
        return;
    }

    const int nhm = nlpp->nhm;
    ModuleBase::matrix vkb1(nhm, npw);
    FPTYPE* vq = new FPTYPE[npw];
    const int x1 = (lmaxkb + 1) * (lmaxkb + 1);

    ModuleBase::matrix ylm(x1, npw);
    ModuleBase::Vector3<FPTYPE>* gk = new ModuleBase::Vector3<FPTYPE>[npw];
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int ig = 0; ig < npw; ig++)
    {
        gk[ig] = wfc_basis->getgpluskcar(ik, ig);
    }
    ModuleBase::YlmReal::Ylm_Real(x1, npw, gk, ylm);

    const int imag_pow_period = 4;
    // result table of pow(0-1i, int)
    static const std::complex<FPTYPE> pref_tab[imag_pow_period] = {{1, 0}, {0, -1}, {-1, 0}, {0, 1}};
    int jkb = 0;
    for (int it = 0; it < this->ucell->ntype; it++)
    {
        if (PARAM.inp.test_pp > 1) {
            ModuleBase::GlobalFunc::OUT("it", it);
}
        // calculate beta in G-space using an interpolation table
        const int nbeta = this->ucell->atoms[it].ncpp.nbeta;
        const int nh = this->ucell->atoms[it].ncpp.nh;

        if (PARAM.inp.test_pp > 1) {
            ModuleBase::GlobalFunc::OUT("nbeta", nbeta);
}

        for (int nb = 0; nb < nbeta; nb++)
        {
            if (PARAM.inp.test_pp > 1) {
                ModuleBase::GlobalFunc::OUT("ib", nb);
}
#ifdef _OPENMP
#pragma omp parallel for
#endif
            for (int ig = 0; ig < npw; ig++)
            {
                const FPTYPE gnorm = gk[ig].norm() * this->ucell->tpiba;
                // cout << "\n gk[ig] = " << gk[ig].x << " " << gk[ig].y << " " << gk[ig].z;
                // cout << "\n gk.norm = " << gnorm;
                vq[ig] = hamilt::Nonlocal_maths<FPTYPE, Device>::Polynomial_Interpolation_nl(nlpp->tab,
                                                                                             it,
                                                                                             nb,
                                                                                             PARAM.globalv.dq,
                                                                                             gnorm);

            } // enddo

            // add spherical harmonic part
            for (int ih = 0; ih < nh; ih++)
            {
                if (nb == nlpp->indv(it, ih))
                {
                    const int lm = static_cast<int>(nlpp->nhtolm(it, ih));
#ifdef _OPENMP
#pragma omp parallel for
#endif
                    for (int ig = 0; ig < npw; ig++)
                    {
                        vkb1(ih, ig) = ylm(lm, ig) * vq[ig];
                    }
                }
            } // end ih
        }     // end nbeta

        // vkb1 contains all betas including angular part for type nt
        // now add the structure factor and factor (-i)^l
        for (int ia = 0; ia < this->ucell->atoms[it].na; ia++)
        {
            std::complex<FPTYPE>* sk = p_sf->get_sk(ik, it, ia, wfc_basis);
#ifdef _OPENMP
#pragma omp parallel for collapse(2)
#endif
            for (int ih = 0; ih < nh; ih++)
            {
                for (int ig = 0; ig < npw; ig++)
                {
                    std::complex<FPTYPE> pref = pref_tab[int(nlpp->nhtol(it, ih)) % imag_pow_period]; //?
                    vkb(jkb + ih, ig) = vkb1(ih, ig) * sk[ig] * pref;
                }

            } // end ih
            jkb += nh;
            delete[] sk;
        } // end ia
    }     // enddo

    delete[] gk;
    delete[] vq;
    //	ModuleBase::timer::tick("Stress","get_dvnl2");

    return;
}

template <typename FPTYPE, typename Device>
FPTYPE Stress_Func<FPTYPE, Device>::Polynomial_Interpolation_nl(const ModuleBase::realArray& table,
                                                                const int& dim1,
                                                                const int& dim2,
                                                                const int& dim3,
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
    const FPTYPE y = (table(dim1, dim2, dim3, iq) * (-x2 * x3 - x1 * x3 - x1 * x2) / 6.0
                      + table(dim1, dim2, dim3, iq + 1) * (+x2 * x3 - x0 * x3 - x0 * x2) / 2.0
                      - table(dim1, dim2, dim3, iq + 2) * (+x1 * x3 - x0 * x3 - x0 * x1) / 2.0
                      + table(dim1, dim2, dim3, iq + 3) * (+x1 * x2 - x0 * x2 - x0 * x1) / 6.0)
                     / table_interval;

    return y;
}

template class Stress_Func<double, base_device::DEVICE_CPU>;
#if ((defined __CUDA) || (defined __ROCM))
template class Stress_Func<double, base_device::DEVICE_GPU>;
#endif