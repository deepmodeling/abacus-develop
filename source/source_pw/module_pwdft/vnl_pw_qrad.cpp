#include "vnl_pw.h"

#include "source_io/module_parameter/parameter.h"
#include "source_base/global_function.h"
#include "source_base/global_variable.h"
#include "source_base/math_polyint.h"
#include "source_base/math_sphbes.h"
#include "source_base/math_integral.h"

#include <cassert>
#include <cmath>

/**
 * @file vnl_pw_qrad.cpp
 * @brief Compute the radial Fourier transform of the USPP augmentation
 *        charge Q(r) and its interpolation table qrad.
 *
 * This file contains:
 * - compute_qrad(): build the qrad interpolation table
 * - radial_fft_q(): interpolate qrad on-the-fly for given G-vectors
 * - Explicit template instantiations for CPU float/double
 */

/**
 * @brief Build the radial Fourier transform interpolation table qrad.
 *
 * For each Vanderbilt USPP pseudopotential, compute
 * qrad(it, l, ijv, iq) = (4*pi/omega) * integral[ r^2 * j_l(q*r) * Q_{ij}(r) dr ]
 * using Simpson integration on the pseudopotential radial grid.
 *
 * @param cell unit cell (read-only, but passed by reference for compatibility)
 */
void pseudopot_cell_vnl::compute_qrad(UnitCell& cell)
{
    const double pref = ModuleBase::FOUR_PI / cell.omega;

    for (int it = 0; it < cell.ntype; it++)
    {
        Atom_pseudo* upf = &cell.atoms[it].ncpp;
        if (upf->tvanp)
        {
            const int nbeta = upf->nbeta;
            int kkbeta = upf->kkbeta;
            if ((kkbeta % 2 == 0) && kkbeta > 0)
            {
                kkbeta--;
            }
            std::vector<double> aux(kkbeta);
            std::vector<double> besr(kkbeta);

            for (int l = 0; l < upf->nqlc; l++)
            {
                for (int iq = 0; iq < PARAM.globalv.nqxq; iq++)
                {
                    const double q = iq * PARAM.globalv.dq;
                    // here we compute the spherical bessel function for each q_i
                    ModuleBase::Sphbes::Spherical_Bessel(kkbeta, upf->r.data(), q, l, besr.data());
                    for (int nb = 0; nb < nbeta; nb++)
                    {
                        // the Q are symmetric with respect to indices nb and mb
                        for (int mb = nb; mb < nbeta; mb++)
                        {
                            const int ijv = mb * (mb + 1) / 2 + nb;
                            if ((l >= std::abs(upf->lll[nb] - upf->lll[mb])) && (l <= (upf->lll[nb] + upf->lll[mb]))
                                && ((l + upf->lll[nb] + upf->lll[mb]) % 2 == 0))
                            {
                                for (int ir = 0; ir < kkbeta; ir++)
                                {
                                    aux[ir] = besr[ir] * upf->qfuncl(l, ijv, ir);
                                }
                                // then we integrate with all the Q functions
                                double vqint = 0.0;
                                ModuleBase::Integral::Simpson_Integral(kkbeta, aux.data(), upf->rab.data(), vqint);
                                qrad(it, l, ijv, iq) = vqint * pref;
                            }
                        }
                    }
                }
            }
        }
    }
}

/**
 * @brief Interpolate the radial Q-function on given G-vectors (CPU matrix version).
 *
 * This is the CPU-only version that works with ModuleBase::matrix for Ylm.
 * It is kept for backward compatibility with older callers.
 *
 * @param ng     number of G-vectors
 * @param ih     first beta function index
 * @param jh     second beta function index
 * @param itype  atom type index
 * @param qnorm  |G| values (size ng)
 * @param ylm    real spherical harmonics matrix
 * @param qg     output Q(G) array (size ng)
 */
void pseudopot_cell_vnl::radial_fft_q(const int ng,
                                      const int ih,
                                      const int jh,
                                      const int itype,
                                      const double* qnorm,
                                      const ModuleBase::matrix ylm,
                                      std::complex<double>* qg) const
{
    // computes the indices which correspond to ih,jh
    const int nb = indv(itype, ih);
    const int mb = indv(itype, jh);
    assert(nb < nbetam);
    assert(mb < nbetam);
    int ijv = 0;
    if (nb >= mb)
    {
        ijv = nb * (nb + 1) / 2 + mb;
    }
    else
    {
        ijv = mb * (mb + 1) / 2 + nb;
    }
    const int ivl = nhtolm(itype, ih);
    const int jvl = nhtolm(itype, jh);

    for (int ig = 0; ig < ng; ig++)
    {
        qg[ig] = {0, 0};
    }
    // makes the sum over the non zero LM
    int l = -1;
    std::complex<double> pref(0.0, 0.0);
    for (int lm = 0; lm < this->lpx(ivl, jvl); lm++)
    {
        int lp = this->lpl(ivl, jvl, lm);
        assert(lp >= 0);
        assert(lp < 49);
        if (lp == 0)
        {
            l = 0;
        }
        else if (lp < 4)
        {
            l = 1;
        }
        else if (lp < 9)
        {
            l = 2;
        }
        else if (lp < 16)
        {
            l = 3;
        }
        else if (lp < 25)
        {
            l = 4;
        }
        else if (lp < 36)
        {
            l = 5;
        }
        else
        {
            l = 6;
        }
        pref = pow(ModuleBase::NEG_IMAG_UNIT, l) * this->ap(lp, ivl, jvl);

        double qm1 = -1.0; // any number smaller than qnorm
        double work = 0.0;
        for (int ig = 0; ig < ng; ig++)
        {
            // calculate quantites depending on the module of G only when needed
            if (std::abs(qnorm[ig] - qm1) > 1e-6)
            {
                work = ModuleBase::PolyInt::Polynomial_Interpolation(this->qrad,
                                                                     itype,
                                                                     l,
                                                                     ijv,
                                                                     PARAM.globalv.nqxq,
                                                                     PARAM.globalv.dq,
                                                                     qnorm[ig]);
                qm1 = qnorm[ig];
            }
            qg[ig] += pref * work * ylm(lp, ig);
        }
    }
}

/**
 * @brief Interpolate the radial Q-function on given G-vectors (CPU raw-pointer template version).
 *
 * This version works with raw CPU pointers.
 * The angular momentum l is determined from the combined lm index.
 *
 * @param ctx    CPU device context
 * @param ng     number of G-vectors
 * @param ih     first beta function index
 * @param jh     second beta function index
 * @param itype  atom type index
 * @param qnorm  |G| values (size ng)
 * @param ylm    real spherical harmonics array (size ng * 49)
 * @param qg     output Q(G) array (size ng)
 */
template <typename FPTYPE, typename Device>
void pseudopot_cell_vnl::radial_fft_q(Device* ctx,
                                      const int ng,
                                      const int ih,
                                      const int jh,
                                      const int itype,
                                      const FPTYPE* qnorm,
                                      const FPTYPE* ylm,
                                      std::complex<FPTYPE>* qg) const
{
    using setmem_complex_op = base_device::memory::set_memory_op<std::complex<FPTYPE>, Device>;

    // computes the indices which correspond to ih,jh
    const int nb = indv(itype, ih);
    const int mb = indv(itype, jh);
    assert(nb < nbetam);
    assert(mb < nbetam);
    int ijv = 0;
    if (nb >= mb)
    {
        ijv = nb * (nb + 1) / 2 + mb;
    }
    else
    {
        ijv = mb * (mb + 1) / 2 + nb;
    }
    const int ivl = nhtolm(itype, ih);
    const int jvl = nhtolm(itype, jh);

    setmem_complex_op()(qg, 0, ng);

    // makes the sum over the non zero LM
    int l = -1;
    std::complex<FPTYPE> pref(0.0, 0.0);
    for (int lm = 0; lm < this->lpx(ivl, jvl); lm++)
    {
        int lp = this->lpl(ivl, jvl, lm);
        assert(lp >= 0);
        assert(lp < 49);
        if (lp == 0)
        {
            l = 0;
        }
        else if (lp < 4)
        {
            l = 1;
        }
        else if (lp < 9)
        {
            l = 2;
        }
        else if (lp < 16)
        {
            l = 3;
        }
        else if (lp < 25)
        {
            l = 4;
        }
        else if (lp < 36)
        {
            l = 5;
        }
        else
        {
            l = 6;
        }
        pref = static_cast<std::complex<FPTYPE>>(pow(ModuleBase::NEG_IMAG_UNIT, l) * this->ap(lp, ivl, jvl));

        double qm1 = -1.0; // any number smaller than qnorm
        double work = 0.0;
        for (int ig = 0; ig < ng; ig++)
        {
            const double qnorm_value = static_cast<double>(qnorm[ig]);
            if (std::abs(qnorm_value - qm1) > 1e-6)
            {
                work = ModuleBase::PolyInt::Polynomial_Interpolation(this->qrad,
                                                                     itype,
                                                                     l,
                                                                     ijv,
                                                                     PARAM.globalv.nqxq,
                                                                     PARAM.globalv.dq,
                                                                     qnorm_value);
                qm1 = qnorm_value;
            }
            qg[ig] += pref * static_cast<FPTYPE>(work) * ylm[lp * ng + ig];
        }
    }
}

// Explicit instantiations for CPU float/double precision.
// These must stay in the same translation unit as the template definition.
template void pseudopot_cell_vnl::radial_fft_q<float, base_device::DEVICE_CPU>(base_device::DEVICE_CPU*,
                                                                               const int,
                                                                               const int,
                                                                               const int,
                                                                               const int,
                                                                               const float*,
                                                                               const float*,
                                                                               std::complex<float>*) const;
template void pseudopot_cell_vnl::radial_fft_q<double, base_device::DEVICE_CPU>(base_device::DEVICE_CPU*,
                                                                                const int,
                                                                                const int,
                                                                                const int,
                                                                                const int,
                                                                                const double*,
                                                                                const double*,
                                                                                std::complex<double>*) const;
