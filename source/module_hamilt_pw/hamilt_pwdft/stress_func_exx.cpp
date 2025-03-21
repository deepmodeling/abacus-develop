#include "stress_pw.h"
#include "global.h"

template <typename FPTYPE, typename Device>
void Stress_PW<FPTYPE, Device>::stress_exx(ModuleBase::matrix& sigma,
                                           const ModuleBase::matrix& wg,
                                           ModulePW::PW_Basis* rhopw,
                                           ModulePW::PW_Basis_K* wfcpw,
                                           const K_Vectors *p_kv,
                                           const psi::Psi<complex<FPTYPE>, Device>* d_psi_in, const UnitCell& ucell)
{
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

    resmem_complex_op()(psi_nk_real, wfcpw->nrxx);
    resmem_complex_op()(psi_mq_real, wfcpw->nrxx);
    resmem_complex_op()(density_real, rhopw->nrxx);
    resmem_complex_op()(density_recip, rhopw->npw);
    resmem_real_op()(pot, rhopw->npw * nks * nks);

    // prepare the coefficients
    double exx_div = 0;

    // pasted from op_exx_pw.cpp
    {
        if (GlobalC::exx_info.info_lip.lambda == 0.0)
        {
            return;
        }

        // here we follow the exx_divergence subroutine in q-e (PW/src/exx_base.f90)
        double alpha = 10.0 / wfcpw->gk_ecut;
        double div = 0;

        // this is the \sum_q F(q) part
        // temporarily for all k points, should be replaced to q points later
        for (int ik = 0; ik < wfcpw->nks; ik++)
        {
            auto k = wfcpw->kvec_c[ik];
#ifdef _OPENMP
#pragma omp parallel for reduction(+:div)
#endif
            for (int ig = 0; ig < rhopw->npw; ig++)
            {
                auto q = k + rhopw->gcar[ig];
                double qq = q.norm2();
                if (qq <= 1e-8) continue;
                else if (GlobalC::exx_info.info_global.ccp_type == Conv_Coulomb_Pot_K::Ccp_Type::Erfc)
                {
                    double hse_omega = GlobalC::exx_info.info_global.hse_omega;
                    double omega2 = hse_omega * hse_omega;
                    div += std::exp(-alpha * qq) / qq * (1.0 - std::exp(-qq*tpiba2 / 4.0 / omega2));
                }
                else
                {
                    div += std::exp(-alpha * qq) / qq;
                }
            }
        }

        Parallel_Reduce::reduce_pool(div);
        // std::cout << "EXX div: " << div << std::endl;

        // if (PARAM.inp.dft_functional == "hse")
        if (GlobalC::exx_info.info_global.ccp_type == Conv_Coulomb_Pot_K::Ccp_Type::Erfc)
        {
            double hse_omega = GlobalC::exx_info.info_global.hse_omega;
            div += tpiba2 / 4.0 / hse_omega / hse_omega; // compensate for the finite value when qq = 0
        }
        else
        {
            div -= alpha;
        }

        div *= ModuleBase::e2 * ModuleBase::FOUR_PI / tpiba2 / wfcpw->nks;

        // numerically value the mean value of F(q) in the reciprocal space
        // This means we need to calculate the average of F(q) in the first brillouin zone
        alpha /= tpiba2;
        int nqq = 100000;
        double dq = 5.0 / std::sqrt(alpha) / nqq;
        double aa = 0.0;
        // if (PARAM.inp.dft_functional == "hse")
        if (GlobalC::exx_info.info_global.ccp_type == Conv_Coulomb_Pot_K::Ccp_Type::Erfc)
        {
            double hse_omega = GlobalC::exx_info.info_global.hse_omega;
            double omega2 = hse_omega * hse_omega;
#ifdef _OPENMP
#pragma omp parallel for reduction(+:aa)
#endif
            for (int i = 0; i < nqq; i++)
            {
                double q = dq * (i+0.5);
                aa -= exp(-alpha * q * q) * exp(-q*q / 4.0 / omega2) * dq;
            }
        }
        aa *= 8 / ModuleBase::FOUR_PI;
        aa += 1.0 / std::sqrt(alpha * ModuleBase::PI);
        div -= ModuleBase::e2 * omega * aa;
        exx_div = div * wfcpw->nks;
        // std::cout << "EXX divergence: " << exx_div << std::endl;

    }

    // prepare for the potential
    for (int ik = 0; ik < nks; ik++)
    {
        for (int iq = 0; iq < nqs; iq++)
        {
            auto k = wfcpw->kvec_c[ik];
            auto q = wfcpw->kvec_c[iq];
            #ifdef _OPENMP
            #pragma omp parallel for schedule(static)
            #endif
            for (int ig = 0; ig < rhopw->npw; ig++)
            {
                FPTYPE qq = (k - q + rhopw->gcar[ig]).norm2() * tpiba2;
                FPTYPE fac = -ModuleBase::FOUR_PI * ModuleBase::e2 / qq;
                if (qq < 1e-8)
                {
                    pot[ig + iq * rhopw->npw + ik * rhopw->npw * nqs] = exx_div;
                }
                else
                {
                    pot[ig + iq * rhopw->npw + ik * rhopw->npw * nqs] = fac;
                }
            }
        }
    }

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
                                auto kqg = wfcpw->kvec_c[ik] - wfcpw->kvec_c[iq] + rhopw->gcar[ig];
                                double kqg_alpha = kqg[alpha] * tpiba;
                                double kqg_beta = kqg[beta] * tpiba;
                                // equation 10 of 10.1103/PhysRevB.73.125120
                                double density_recip2 = std::real(density_recip[ig] * std::conj(density_recip[ig]));
                                double pot_local = pot[ig + iq * rhopw->npw + ik * rhopw->npw * nqs];
                                double _4pi_e2 = ModuleBase::FOUR_PI * ModuleBase::e2;
                                sigma_ab_loc += density_recip2 * pot_local * (kqg_alpha * kqg_beta * (-pot_local) / _4pi_e2 - delta_ab) ;
//                                if (std::abs(pot_local + 22.235163511253440) < 1e-2)
//                                {
//                                    std::cout << "delta_ab: " << delta_ab << std::endl;
//                                    std::cout << "density_recip2: " << density_recip2 << std::endl;
//                                    std::cout << "pot_local: " << pot_local << std::endl;
//                                    std::cout << "kqg_alpha: " << kqg_alpha << std::endl;
//                                    std::cout << "kqg_beta: " << kqg_beta << std::endl;
//
//                                }
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

////    print sigma
//    for (int i = 0; i < 3; i++)
//    {
//        for (int j = 0; j < 3; j++)
//        {
//            std::cout << sigma(i, j) * ModuleBase::RYDBERG_SI / pow(ModuleBase::BOHR_RADIUS_SI, 3) * 1.0e-8 << " ";
//            sigma(i, j) = 0;
//        }
//        std::cout << std::endl;
//    }


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
