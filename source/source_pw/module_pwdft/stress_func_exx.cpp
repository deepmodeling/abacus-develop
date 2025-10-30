#include "global.h"
#include "operator_pw/op_exx_pw.h"
#include "source_base/parallel_common.h"
#include "stress_pw.h"

template <typename FPTYPE, typename Device>
void Stress_PW<FPTYPE, Device>::stress_exx(ModuleBase::matrix& sigma,
                                           const ModuleBase::matrix& wg,
                                           ModulePW::PW_Basis* rhopw,
                                           ModulePW::PW_Basis_K* wfcpw,
                                           const K_Vectors *p_kv,
                                           const psi::Psi <std::complex<FPTYPE>, Device>* d_psi_in, const UnitCell& ucell)
{
    bool gamma_extrapolation = PARAM.inp.exx_gamma_extrapolation;
    bool is_mp = p_kv->get_is_mp();
#ifdef __MPI
    Parallel_Common::bcast_bool(is_mp);
#endif
    if (!is_mp)
    {
        gamma_extrapolation = false;
    }

    // T is complex of FPTYPE, if FPTYPE is double, T is std::complex<double>
    // but if FPTYPE is std::complex<double>, T is still std::complex<double>
    using T = std::complex<FPTYPE>;
    using Real = FPTYPE;
    using setmem_complex_op = base_device::memory::set_memory_op<T, Device>;
    using resmem_complex_op = base_device::memory::resize_memory_op<T, Device>;
    using resmem_real_op = base_device::memory::resize_memory_op<Real, Device>;
    using delmem_complex_op = base_device::memory::delete_memory_op<T, Device>;
    using delmem_real_op = base_device::memory::delete_memory_op<Real, Device>;
    using syncmem_complex_op = base_device::memory::synchronize_memory_op<T, Device, Device>;

    int nks = wfcpw->nks;
    int nqs = wfcpw->nks; // currently q-points downsampling is not supported
    double omega = ucell.omega;
    double tpiba = ucell.tpiba;
    double tpiba2 = ucell.tpiba2;
    double omega_inv = 1.0 / omega;

    // allocate space
    T* psi_nk_real = nullptr;
    T* psi_mq_real = nullptr;
    T* density_real = nullptr;
    T* density_recip = nullptr;
    Real* pot = nullptr; // This factor is 2x of the potential in 10.1103/PhysRevB.73.125120
    Real* pot_stress = nullptr;

    resmem_complex_op()(psi_nk_real, wfcpw->nrxx);
    resmem_complex_op()(psi_mq_real, wfcpw->nrxx);
    resmem_complex_op()(density_real, rhopw->nrxx);
    resmem_complex_op()(density_recip, rhopw->npw);
    resmem_real_op()(pot, rhopw->npw);
    resmem_real_op()(pot_stress, rhopw->npw);

    // hamilt::get_exx_potential<Real, Device>(p_kv, wfcpw, rhopw, pot, tpiba, gamma_extrapolation, omega);
    // hamilt::get_exx_stress_potential<Real, Device>(p_kv, wfcpw, rhopw, pot_stress, tpiba, gamma_extrapolation, omega);

    // calculate the stress

    // for nk, mq
    for (int ik = 0; ik < nks; ik++)
    {
        for (int nband = 0; nband < d_psi_in->get_nbands(); nband++)
        {
            if (wg(ik, nband) < 1e-12) continue;
            // psi_nk in real space
            d_psi_in->fix_kb(ik, nband);
            T* psi_nk = d_psi_in->get_pointer();
            wfcpw->recip2real(psi_nk, psi_nk_real, ik);

            for (int iq = 0; iq < nqs; iq++)
            {
                hamilt::get_exx_potential<Real, Device>(p_kv, wfcpw, rhopw, pot, tpiba, gamma_extrapolation, omega, ik, iq, true);
                hamilt::get_exx_stress_potential<Real, Device>(p_kv, wfcpw, rhopw, pot_stress, tpiba, gamma_extrapolation, omega, ik, iq);
                for (int mband = 0; mband < d_psi_in->get_nbands(); mband++)
                {
                    // psi_mq in real space
                    d_psi_in->fix_kb(iq, mband);
                    T* psi_mq = d_psi_in->get_pointer();
                    wfcpw->recip2real(psi_mq, psi_mq_real, iq);

                    // overlap density in real space
                    setmem_complex_op()(density_real, 0.0, rhopw->nrxx);
                    for (int ig = 0; ig < rhopw->nrxx; ig++)
                    {
                        density_real[ig] = psi_nk_real[ig] * std::conj(psi_mq_real[ig]) * omega_inv;
                    }

                    // density in reciprocal space
                    rhopw->real2recip(density_real, density_recip);

                    // really calculate the stress

                    // for alpha beta
                    for (int alpha = 0; alpha < 3; alpha++)
                    {
                        for (int beta = alpha; beta < 3; beta++)
                        {
                            int delta_ab = (alpha == beta) ? 1 : 0;
                            double sigma_ab_loc = 0.0;
                            #ifdef _OPENMP
                            #pragma omp parallel for schedule(static) reduction(+:sigma_ab_loc)
                            #endif
                            for (int ig = 0; ig < rhopw->npw; ig++)
                            {
                                const ModuleBase::Vector3<double> kqg = wfcpw->kvec_c[ik] - wfcpw->kvec_c[iq] + rhopw->gcar[ig];
                                double kqg_alpha = kqg[alpha] * tpiba;
                                double kqg_beta = kqg[beta] * tpiba;
                                // equation 10 of 10.1103/PhysRevB.73.125120
                                double density_recip2 = std::real(density_recip[ig] * std::conj(density_recip[ig]));
                                const int idx = ig;
                                double pot_local = pot[idx];
                                double pot_stress_local = pot_stress[idx];
                                sigma_ab_loc += density_recip2 * pot_local * (kqg_alpha * kqg_beta * pot_stress_local - delta_ab) ;

                            }

                            // 0.5 in the following line is caused by 2x in the pot
                            sigma(alpha, beta) -= GlobalC::exx_info.info_global.hybrid_alpha
                                                  * 0.25 * sigma_ab_loc
                                                  * wg(ik, nband) * wg(iq, mband) / nqs / p_kv->wk[ik];
                        }
                    }
                }
            }
        }
    }

    for (int l = 0; l < 3; l++)
    {
        for (int m = l + 1; m < 3; m++)
        {
            sigma(m, l) = sigma(l, m);
        }
    }

    Parallel_Reduce::reduce_all(sigma.c, sigma.nr * sigma.nc);


    delmem_complex_op()(psi_nk_real);
    delmem_complex_op()(psi_mq_real);
    delmem_complex_op()(density_real);
    delmem_complex_op()(density_recip);
    delmem_real_op()(pot);
}

template class Stress_PW<double, base_device::DEVICE_CPU>;
#if ((defined __CUDA) || (defined __ROCM))
template class Stress_PW<double, base_device::DEVICE_GPU>;
#endif
