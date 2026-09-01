#include "to_w90_pw.h"

#include "source_io/module_parameter/parameter.h"
#include "source_base/math_integral.h"
#include "source_base/math_polyint.h"
#include "source_base/math_sphbes.h"
#include "source_base/math_ylmreal.h"
#include "source_base/parallel_reduce.h"

void toW90_PW::gen_radial_function_in_q(std::vector<ModuleBase::matrix> &radial_in_q)
{
    // The radial function is Fourier transformed into the q-space.
    radial_in_q.resize(num_wannier);

    std::vector<double> r(mesh_r);
    std::vector<double> dr(mesh_r);
    std::vector<double> psi(mesh_r);
    std::vector<double> psir(mesh_r);

    for (int wannier_index = 0; wannier_index < num_wannier; wannier_index++)
    {
        double x = 0;
        for (int ir = 0; ir < mesh_r; ir++)
        {
            x = x_min + ir * dx;
            r[ir] = exp(x) / alfa[wannier_index];
            dr[ir] = dx * r[ir];
        }

        double alfa32 = pow(alfa[wannier_index], 3.0 / 2.0);
        double alfa_new = alfa[wannier_index];

        if (rvalue[wannier_index] == 1)
        {
            for (int ir = 0; ir < mesh_r; ir++)
            {
                psi[ir] = 2.0 * alfa32 * exp(-alfa_new * r[ir]);
            }
        }

        if (rvalue[wannier_index] == 2)
        {
            for (int ir = 0; ir < mesh_r; ir++)
            {
                psi[ir] = 1.0 / sqrt(8.0) * alfa32 * (2.0 - alfa_new * r[ir]) * exp(-alfa_new * r[ir] * 0.5);
            }
        }

        if (rvalue[wannier_index] == 3)
        {
            for (int ir = 0; ir < mesh_r; ir++)
            {
                psi[ir] = sqrt(4.0 / 27.0) * alfa32
                        * (1.0 - 2.0 / 3.0 * alfa_new * r[ir] + 2.0 / 27.0 * pow(alfa_new, 2.0) * r[ir] * r[ir])
                        * exp(-alfa_new * r[ir] * 1.0 / 3.0);
            }
        }

        for (int ir = 0; ir < mesh_r; ir++)
        {
            psir[ir] = psi[ir] * r[ir];
        }

        auto &tmp_radial = radial_in_q[wannier_index];
        if (L[wannier_index] >= 0)
        {
            tmp_radial.create(1, nqx_);
            integral(mesh_r, psir.data(), r.data(), dr.data(), L[wannier_index], tmp_radial.c);
        }
        else
        {
            int tmp_size = 0;

            if (L[wannier_index] == -1 || L[wannier_index] == -2 || L[wannier_index] == -3) tmp_size = 2;

            if (L[wannier_index] == -4 || L[wannier_index] == -5) tmp_size = 3;

            tmp_radial.create(tmp_size, nqx_);

            for (int tmp_L = 0; tmp_L < tmp_size; tmp_L++)
            {
                integral(mesh_r, psir.data(), r.data(), dr.data(), tmp_L, tmp_radial.c+tmp_L*nqx_);
            }
        }

    }

}

void toW90_PW::produce_trial_in_pw(
    const psi::Psi<std::complex<double>>& psi_pw,
    const int& ik,
    const ModulePW::PW_Basis_K* wfcpw,
    const std::vector<ModuleBase::matrix> &radial_in_q,
    ModuleBase::ComplexMatrix& trial_orbitals_k
)
{
    const int npw = wfcpw->npwk[ik];
    const int npwx = wfcpw->npwk_max;

    trial_orbitals_k.create(num_wannier, npwx);

    const int total_lm = 16;
    ModuleBase::matrix ylm(total_lm, npw);

    std::vector<ModuleBase::Vector3<double>> gk(npw);
    for (int ig = 0; ig < npw; ig++)
    {
        gk[ig] = wfcpw->getgpluskcar(ik, ig);
    }

    ModuleBase::YlmReal::Ylm_Real(total_lm, npw, gk.data(), ylm);

    // Keep it consistent with the definition of spherical harmonic functions in Wannier90 
    std::vector<int> need_inv = {2, 3, 5, 6, 14, 15};
    for (auto index : need_inv)
    {
        for (int ig = 0; ig < npw; ig++)
        {
            ylm(index, ig) = -1.0 * ylm(index, ig);
        }
    }

    const double bs2 = 1.0 / sqrt(2.0);
    const double bs3 = 1.0 / sqrt(3.0);
    const double bs6 = 1.0 / sqrt(6.0);
    const double bs12 = 1.0 / sqrt(12.0);

    std::vector<std::complex<double>> sf(npw);
    for (int wannier_index = 0; wannier_index < num_wannier; wannier_index++)
    {
        double* radial_c = radial_in_q[wannier_index].c;
        const int nqx = nqx_;

        if (L[wannier_index] >= 0)
        {
            get_trial_orbitals_lm_k(L[wannier_index],
                                    m[wannier_index],
                                    ylm,
                                    gk.data(),
                                    npw,
                                    radial_c,
                                    trial_orbitals_k.c + wannier_index*npwx);

            for (int ig = 0; ig < npw; ig++)
            {
                const double arg = (gk[ig] * R_centre[wannier_index]) * ModuleBase::TWO_PI;
                sf[ig] = std::complex<double>(cos(arg), -sin(arg));

                trial_orbitals_k(wannier_index, ig) = sf[ig] * trial_orbitals_k(wannier_index, ig);
            }

        }
        else
        {
            if (L[wannier_index] == -1)
            {
                double tmp_bs2 = 0;
                if (m[wannier_index] == 0) tmp_bs2 = bs2;
                if (m[wannier_index] == 1) tmp_bs2 = -bs2;

                std::vector<std::complex<double>> orb_s(npw);
                std::vector<std::complex<double>> orb_px(npw);

                get_trial_orbitals_lm_k(0, 0, ylm, gk.data(), npw, radial_c, orb_s.data());
                get_trial_orbitals_lm_k(1, 1, ylm, gk.data(), npw, radial_c + nqx, orb_px.data());

                for (int ig = 0; ig < npw; ig++)
                {
                    const double arg = (gk[ig] * R_centre[wannier_index]) * ModuleBase::TWO_PI;
                    sf[ig] = std::complex<double>(cos(arg), -sin(arg));

                    trial_orbitals_k(wannier_index, ig) = sf[ig] * (bs2 * orb_s[ig]  + tmp_bs2 * orb_px[ig]);
                }
            }

            if (L[wannier_index] == -2)
            {
                if (m[wannier_index] == 0 || m[wannier_index] == 1)
                {
                    double tmp_bs2 = bs2;
                    if (m[wannier_index] == 1) tmp_bs2 = -bs2;

                    std::vector<std::complex<double>> orb_s(npw);
                    std::vector<std::complex<double>> orb_px(npw);
                    std::vector<std::complex<double>> orb_py(npw);

                    get_trial_orbitals_lm_k(0, 0, ylm, gk.data(), npw, radial_c, orb_s.data());
                    get_trial_orbitals_lm_k(1, 1, ylm, gk.data(), npw, radial_c + nqx, orb_px.data());
                    get_trial_orbitals_lm_k(1, 2, ylm, gk.data(), npw, radial_c + nqx, orb_py.data());

                    for (int ig = 0; ig < npw; ig++)
                    {
                        const double arg = (gk[ig] * R_centre[wannier_index]) * ModuleBase::TWO_PI;
                        sf[ig] = std::complex<double>(cos(arg), -sin(arg));

                        trial_orbitals_k(wannier_index, ig)
                            = sf[ig] * (bs3 * orb_s[ig] - bs6 * orb_px[ig] + tmp_bs2 * orb_py[ig]);
                    }
                }
                else if (m[wannier_index] == 2)
                {
                    std::vector<std::complex<double>> orb_s(npw);
                    std::vector<std::complex<double>> orb_px(npw);
                    get_trial_orbitals_lm_k(0, 0, ylm, gk.data(), npw, radial_c, orb_s.data());
                    get_trial_orbitals_lm_k(1, 1, ylm, gk.data(), npw, radial_c + nqx, orb_px.data());

                    for (int ig = 0; ig < npw; ig++)
                    {
                        const double arg = (gk[ig] * R_centre[wannier_index]) * ModuleBase::TWO_PI;
                        sf[ig] = std::complex<double>(cos(arg), -sin(arg));

                        trial_orbitals_k(wannier_index, ig) = sf[ig] * (bs3 * orb_s[ig] + 2.0 * bs6 * orb_px[ig]);
                    }
                }
            }

            if (L[wannier_index] == -3)
            {
                double m_px = 1.0;
                double m_py = 1.0;
                double m_pz = 1.0;

                if (m[wannier_index] == 1)
                {
                    m_py = -1.0;
                    m_pz = -1.0;
                }
                else if (m[wannier_index] == 2)
                {
                    m_px = -1.0;
                    m_pz = -1.0;
                }
                else if (m[wannier_index] == 3)
                {
                    m_px = -1.0;
                    m_py = -1.0;
                }

                std::vector<std::complex<double>> orb_s(npw);
                std::vector<std::complex<double>> orb_px(npw);
                std::vector<std::complex<double>> orb_py(npw);
                std::vector<std::complex<double>> orb_pz(npw);

                get_trial_orbitals_lm_k(0, 0, ylm, gk.data(), npw, radial_c, orb_s.data());
                get_trial_orbitals_lm_k(1, 1, ylm, gk.data(), npw, radial_c + nqx, orb_px.data());
                get_trial_orbitals_lm_k(1, 2, ylm, gk.data(), npw, radial_c + nqx, orb_py.data());
                get_trial_orbitals_lm_k(1, 0, ylm, gk.data(), npw, radial_c + nqx, orb_pz.data());

                for (int ig = 0; ig < npw; ig++)
                {
                    const double arg = (gk[ig] * R_centre[wannier_index]) * ModuleBase::TWO_PI;
                    sf[ig] = std::complex<double>(cos(arg), -sin(arg));

                    trial_orbitals_k(wannier_index, ig)
                            = sf[ig] * 0.5 * (orb_s[ig] + m_px * orb_px[ig] + m_py * orb_py[ig] + m_pz * orb_pz[ig]);
                }
            }
            
            if (L[wannier_index] == -4)
            {
                if (m[wannier_index] == 0 || m[wannier_index] == 1)
                {
                    double tmp_bs2 = bs2;
                    if (m[wannier_index] == 1) tmp_bs2 = -bs2;

                    std::vector<std::complex<double>> orb_s(npw);
                    std::vector<std::complex<double>> orb_px(npw);
                    std::vector<std::complex<double>> orb_py(npw);

                    get_trial_orbitals_lm_k(0, 0, ylm, gk.data(), npw, radial_c, orb_s.data());
                    get_trial_orbitals_lm_k(1, 1, ylm, gk.data(), npw, radial_c + nqx, orb_px.data());
                    get_trial_orbitals_lm_k(1, 2, ylm, gk.data(), npw, radial_c + nqx, orb_py.data());

                    for (int ig = 0; ig < npw; ig++)
                    {
                        const double arg = (gk[ig] * R_centre[wannier_index]) * ModuleBase::TWO_PI;
                        sf[ig] = std::complex<double>(cos(arg), -sin(arg));

                        trial_orbitals_k(wannier_index, ig)
                            = sf[ig] * (bs3 * orb_s[ig] - bs6 * orb_px[ig] + tmp_bs2 * orb_py[ig]);
                    }
                }
                else if (m[wannier_index] == 2)
                {
                    std::vector<std::complex<double>> orb_s(npw);
                    std::vector<std::complex<double>> orb_px(npw);

                    get_trial_orbitals_lm_k(0, 0, ylm, gk.data(), npw, radial_c, orb_s.data());
                    get_trial_orbitals_lm_k(1, 1, ylm, gk.data(), npw, radial_c + nqx, orb_px.data());

                    for (int ig = 0; ig < npw; ig++)
                    {
                        const double arg = (gk[ig] * R_centre[wannier_index]) * ModuleBase::TWO_PI;
                        sf[ig] = std::complex<double>(cos(arg), -sin(arg));

                        trial_orbitals_k(wannier_index, ig) = sf[ig] * (bs3 * orb_s[ig] + 2.0 * bs6 * orb_px[ig]);
                    }
                }
                else if (m[wannier_index] == 3 || m[wannier_index] == 4)
                {
                    double m_pz = 1.0;
                    if (m[wannier_index] == 4) m_pz = -1.0;

                    std::vector<std::complex<double>> orb_pz(npw);
                    std::vector<std::complex<double>> orb_dz2(npw);

                    get_trial_orbitals_lm_k(1, 0, ylm, gk.data(), npw, radial_c + nqx, orb_pz.data());
                    get_trial_orbitals_lm_k(2, 0, ylm, gk.data(), npw, radial_c + 2 * nqx, orb_dz2.data());

                    for (int ig = 0; ig < npw; ig++)
                    {
                        const double arg = (gk[ig] * R_centre[wannier_index]) * ModuleBase::TWO_PI;
                        sf[ig] = std::complex<double>(cos(arg), -sin(arg));

                        trial_orbitals_k(wannier_index, ig) = sf[ig] * bs2 * (m_pz * orb_pz[ig] + orb_dz2[ig]);
                    }
                }
            }

            if (L[wannier_index] == -5)
            {
                if (m[wannier_index] == 0 || m[wannier_index] == 1)
                {
                    double tmp_bs2 = -bs2;
                    double tmp_bs12 = -bs12;
                    double tmp_d = 0.5;

                    if (m[wannier_index] == 1)
                    {
                        tmp_bs2 = bs2;
                    }

                    std::vector<std::complex<double>> orb_s(npw);
                    std::vector<std::complex<double>> orb_px(npw);
                    std::vector<std::complex<double>> orb_dz2(npw);
                    std::vector<std::complex<double>> orb_dx2_y2(npw);

                    get_trial_orbitals_lm_k(0, 0, ylm, gk.data(), npw, radial_c, orb_s.data());
                    get_trial_orbitals_lm_k(1, 1, ylm, gk.data(), npw, radial_c + nqx, orb_px.data());
                    get_trial_orbitals_lm_k(2, 0, ylm, gk.data(), npw, radial_c + 2 * nqx, orb_dz2.data());
                    get_trial_orbitals_lm_k(2, 3, ylm, gk.data(), npw, radial_c + 2 * nqx, orb_dx2_y2.data());

                    for (int ig = 0; ig < npw; ig++)
                    {
                        const double arg = (gk[ig] * R_centre[wannier_index]) * ModuleBase::TWO_PI;
                        sf[ig] = std::complex<double>(cos(arg), -sin(arg));

                        trial_orbitals_k(wannier_index, ig)
                            = sf[ig] * (bs6 * orb_s[ig] + tmp_bs2 * orb_px[ig] + tmp_bs12 * orb_dz2[ig]
                                        + tmp_d * orb_dx2_y2[ig]);
                    }
                }
                else if (m[wannier_index] == 2 || m[wannier_index] == 3)
                {
                    double tmp_bs2 = -bs2;
                    double tmp_bs12 = -bs12;
                    double tmp_d = -0.5;

                    if (m[wannier_index] == 3)
                    {
                        tmp_bs2 = bs2;
                    }

                    std::vector<std::complex<double>> orb_s(npw);
                    std::vector<std::complex<double>> orb_py(npw);
                    std::vector<std::complex<double>> orb_dz2(npw);
                    std::vector<std::complex<double>> orb_dx2_y2(npw);

                    get_trial_orbitals_lm_k(0, 0, ylm, gk.data(), npw, radial_c, orb_s.data());
                    get_trial_orbitals_lm_k(1, 2, ylm, gk.data(), npw, radial_c + nqx, orb_py.data());
                    get_trial_orbitals_lm_k(2, 0, ylm, gk.data(), npw, radial_c + 2 * nqx, orb_dz2.data());
                    get_trial_orbitals_lm_k(2, 3, ylm, gk.data(), npw, radial_c + 2 * nqx, orb_dx2_y2.data());

                    for (int ig = 0; ig < npw; ig++)
                    {
                        const double arg = (gk[ig] * R_centre[wannier_index]) * ModuleBase::TWO_PI;
                        sf[ig] = std::complex<double>(cos(arg), -sin(arg));

                        trial_orbitals_k(wannier_index, ig)
                            = sf[ig] * (bs6 * orb_s[ig] + tmp_bs2 * orb_py[ig] + tmp_bs12 * orb_dz2[ig]
                                        + tmp_d * orb_dx2_y2[ig]);
                    }
                }
                else if (m[wannier_index] == 4 || m[wannier_index] == 5)
                {
                    double tmp_pz = -1.0;

                    if (m[wannier_index] == 5) tmp_pz = 1.0;

                    std::vector<std::complex<double>> orb_s(npw);
                    std::vector<std::complex<double>> orb_pz(npw);
                    std::vector<std::complex<double>> orb_dz2(npw);

                    get_trial_orbitals_lm_k(0, 0, ylm, gk.data(), npw, radial_c, orb_s.data());
                    get_trial_orbitals_lm_k(1, 0, ylm, gk.data(), npw, radial_c + nqx, orb_pz.data());
                    get_trial_orbitals_lm_k(2, 0, ylm, gk.data(), npw, radial_c + 2 * nqx, orb_dz2.data());

                    for (int ig = 0; ig < npw; ig++)
                    {
                        const double arg = (gk[ig] * R_centre[wannier_index]) * ModuleBase::TWO_PI;
                        sf[ig] = std::complex<double>(cos(arg), -sin(arg));

                        trial_orbitals_k(wannier_index, ig)
                            = sf[ig] * (bs6 * orb_s[ig] + tmp_pz * bs2 * orb_pz[ig] + bs3 * orb_dz2[ig]);
                    }
                }
            }

        }
    }

    for (int wannier_index = 0; wannier_index < num_wannier; wannier_index++)
    {
        std::complex<double> anorm(0.0, 0.0);
        for (int ig = 0; ig < npw; ig++)
        {
            anorm += conj(trial_orbitals_k(wannier_index, ig)) * trial_orbitals_k(wannier_index, ig);
        }

        Parallel_Reduce::reduce_all(anorm);
        for (int ig = 0; ig < npw; ig++)
        {
            trial_orbitals_k(wannier_index, ig) = trial_orbitals_k(wannier_index, ig) / sqrt(anorm);
        }
    }
}

void toW90_PW::get_trial_orbitals_lm_k(
    const int &orbital_L,
    const int &orbital_m,
    const ModuleBase::matrix &ylm,
    const ModuleBase::Vector3<double> *gk,
    const int &npw,
    double *radial_in_q_single,
    std::complex<double> *orbital_in_G_single
)
{
    const double tpiba = *this->tpiba;
    for (int ig = 0; ig < npw; ig++)
    {
        orbital_in_G_single[ig] = ModuleBase::PolyInt::Polynomial_Interpolation(
            radial_in_q_single, nqx_, dq_, gk[ig].norm() * tpiba);
    }

    std::complex<double> lphase = pow(ModuleBase::NEG_IMAG_UNIT, orbital_L);
    int index = orbital_L * orbital_L + orbital_m;
    for (int ig = 0; ig < npw; ig++)
    {
        orbital_in_G_single[ig] = lphase * ylm(index, ig) * orbital_in_G_single[ig];
    }

    return;
}

void toW90_PW::integral(
    const int meshr,
    const double *psir,
    const double *r,
    const double *rab,
    const int &l,
    double *table
)
{
    const double pref = ModuleBase::FOUR_PI / sqrt(*this->omega);

    std::vector<double> inner_part(meshr);
    for (int ir = 0; ir < meshr; ir++)
    {
        inner_part[ir] = psir[ir] * psir[ir];
    }

    double unit = 0.0;
    ModuleBase::Integral::Simpson_Integral(meshr, inner_part.data(), rab, unit);

    std::vector<double> aux(meshr);
    std::vector<double> vchi(meshr);
    for (int iq = 0; iq < nqx_; iq++)
    {
        const double q = dq_ * iq;
        ModuleBase::Sphbes::Spherical_Bessel(meshr, r, q, l, aux.data());
        for (int ir = 0; ir < meshr; ir++)
        {
            vchi[ir] = psir[ir] * aux[ir] * r[ir];
        }

        double vqint = 0.0;
        ModuleBase::Integral::Simpson_Integral(meshr, vchi.data(), rab, vqint);

        table[iq] = vqint * pref;
    }
    return;
}
