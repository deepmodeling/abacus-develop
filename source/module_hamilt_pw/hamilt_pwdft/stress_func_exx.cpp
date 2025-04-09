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
    double nqs_half1 = 0.5 * p_kv->nmp[0];
    double nqs_half2 = 0.5 * p_kv->nmp[1];
    double nqs_half3 = 0.5 * p_kv->nmp[2];
    bool gamma_extrapolation = PARAM.inp.exx_gamma_extrapolation;
    if (!p_kv->get_is_mp())
    {
        gamma_extrapolation = false;
    }
    auto isint = [](double x)
    {
        double epsilon = 1e-6; // this follows the isint judgement in q-e
        return std::abs(x - std::round(x)) < epsilon;
    };

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
    resmem_real_op()(pot, rhopw->npw * nks * nks);
    resmem_real_op()(pot_stress, rhopw->npw * nks * nks);

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
            const ModuleBase::Vector3<double> k_c = wfcpw->kvec_c[ik];
            const ModuleBase::Vector3<double> k_d = wfcpw->kvec_d[ik];
#ifdef _OPENMP
#pragma omp parallel for reduction(+:div)
#endif
            for (int ig = 0; ig < rhopw->npw; ig++)
            {
                const ModuleBase::Vector3<double> q_c = k_c + rhopw->gcar[ig];
                const ModuleBase::Vector3<double> q_d = k_d + rhopw->gdirect[ig];
                double qq = q_c.norm2();
                // For gamma_extrapolation (https://doi.org/10.1103/PhysRevB.79.205114)
                // 7/8 of the points in the grid are "activated" and 1/8 are disabled.
                // grid_factor is designed for the 7/8 of the grid to function like all of the points
                double grid_factor = 1;
                double extrapolate_grid = 8.0/7.0;
                if (gamma_extrapolation)
                {
                    if (isint(q_d[0] * nqs_half1) &&
                        isint(q_d[1] * nqs_half2) &&
                        isint(q_d[2] * nqs_half3))
                    {
                        grid_factor = 0;
                    }
                    else
                    {
                        grid_factor = extrapolate_grid;
                    }
                }


                if (qq <= 1e-8) continue;
                else if (GlobalC::exx_info.info_global.ccp_type == Conv_Coulomb_Pot_K::Ccp_Type::Erfc)
                {
                    double hse_omega = GlobalC::exx_info.info_global.hse_omega;
                    double omega2 = hse_omega * hse_omega;
                    div += std::exp(-alpha * qq) / qq * (1.0 - std::exp(-qq*tpiba2 / 4.0 / omega2)) * grid_factor;
                }
                else
                {
                    div += std::exp(-alpha * qq) / qq * grid_factor;
                }
            }
        }

        Parallel_Reduce::reduce_pool(div);
        // std::cout << "EXX div: " << div << std::endl;

        // if (PARAM.inp.dft_functional == "hse")
        if (!gamma_extrapolation)
        {
            if (GlobalC::exx_info.info_global.ccp_type == Conv_Coulomb_Pot_K::Ccp_Type::Erfc)
            {
                double omega = GlobalC::exx_info.info_global.hse_omega;
                div += tpiba2 / 4.0 / omega / omega; // compensate for the finite value when qq = 0
            }
            else
            {
                div -= alpha;
            }

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
        for (int iq = 0; iq < nks; iq++)
        {
            const ModuleBase::Vector3<double> k_c = wfcpw->kvec_c[ik];
            const ModuleBase::Vector3<double> k_d = wfcpw->kvec_d[ik];
            const ModuleBase::Vector3<double> q_c = wfcpw->kvec_c[iq];
            const ModuleBase::Vector3<double> q_d = wfcpw->kvec_d[iq];

            #ifdef _OPENMP
            #pragma omp parallel for schedule(static)
            #endif
            for (int ig = 0; ig < rhopw->npw; ig++)
            {
                const ModuleBase::Vector3<double> g_d = rhopw->gdirect[ig];
                const ModuleBase::Vector3<double> kqg_d = k_d - q_d + g_d;

                // For gamma_extrapolation (https://doi.org/10.1103/PhysRevB.79.205114)
                // 7/8 of the points in the grid are "activated" and 1/8 are disabled.
                // grid_factor is designed for the 7/8 of the grid to function like all of the points
                Real grid_factor = 1;
                if (gamma_extrapolation)
                {
                    double extrapolate_grid = 8.0/7.0;
                    if (isint(kqg_d[0] * nqs_half1) &&
                        isint(kqg_d[1] * nqs_half2) &&
                        isint(kqg_d[2] * nqs_half3))
                    {
                        grid_factor = 0;
                    }
                    else
                    {
                        grid_factor = extrapolate_grid;
                    }
                }

                const int ig_kq = ik * nks * rhopw->npw + iq * rhopw->npw + ig;

                Real gg = (k_c - q_c + rhopw->gcar[ig]).norm2() * tpiba2;
                Real hse_omega2 = GlobalC::exx_info.info_global.hse_omega * GlobalC::exx_info.info_global.hse_omega;
                // if (kqgcar2 > 1e-12) // vasp uses 1/40 of the smallest (k spacing)**2
                if (gg >= 1e-8)
                {
                    Real fac = -ModuleBase::FOUR_PI * ModuleBase::e2 / gg;
                    // if (PARAM.inp.dft_functional == "hse")
                    if (GlobalC::exx_info.info_global.ccp_type == Conv_Coulomb_Pot_K::Ccp_Type::Erfc)
                    {
                        pot[ig_kq] = fac * (1.0 - std::exp(-gg / 4.0 / hse_omega2)) * grid_factor;
                        pot_stress[ig_kq] = (1.0 - (1.0 + gg / 4.0 / hse_omega2) * std::exp(-gg / 4.0 / hse_omega2)) / (1.0 - std::exp(-gg / 4.0 / hse_omega2)) / gg;
                    }
                    else if (GlobalC::exx_info.info_global.ccp_type == Conv_Coulomb_Pot_K::Ccp_Type::Erf)
                    {
                        ModuleBase::WARNING("Stress_PW", "Stress for Erf is not implemented yet");
                        pot[ig_kq] = fac * grid_factor;
                        pot_stress[ig_kq] = 1.0 / gg;
                    }
                    else if (GlobalC::exx_info.info_global.ccp_type == Conv_Coulomb_Pot_K::Ccp_Type::Hf)
                    {
                        pot[ig_kq] = fac * grid_factor;
                        pot_stress[ig_kq] = 1.0 / gg;
                    }
                }
                // }
                else
                {
                    // if (PARAM.inp.dft_functional == "hse")
                    if (GlobalC::exx_info.info_global.ccp_type == Conv_Coulomb_Pot_K::Ccp_Type::Erfc && !gamma_extrapolation)
                    {
                        pot[ig_kq] = - ModuleBase::PI * ModuleBase::e2 / hse_omega2; // maybe we should add a exx_div here, but q-e does not do that
                        pot_stress[ig_kq] = 1 / 4.0 / hse_omega2;
                    }
                    else
                    {
                        pot[ig_kq] = exx_div;
                        pot_stress[ig_kq] = 0;
                    }
                }
                // assert(is_finite(density_recip[ig]));
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
                                const ModuleBase::Vector3<double> kqg = wfcpw->kvec_c[ik] - wfcpw->kvec_c[iq] + rhopw->gcar[ig];
                                double kqg_alpha = kqg[alpha] * tpiba;
                                double kqg_beta = kqg[beta] * tpiba;
                                // equation 10 of 10.1103/PhysRevB.73.125120
                                double density_recip2 = std::real(density_recip[ig] * std::conj(density_recip[ig]));
                                const int idx = ig + iq * rhopw->npw + ik * rhopw->npw * nqs;
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
