#include "module_base/libm/libm.h"
#include "module_parameter/parameter.h"
#include "module_base/math_polyint.h"
#include "module_base/math_ylmreal.h"
#include "module_base/timer.h"
#include "module_elecstate/elecstate_pw.h"
#include "module_hamilt_pw/hamilt_pwdft/nonlocal_maths.hpp"
#include "stress_pw.h"

// computes the part of the crystal stress which is due
// to the dependence of the Q function on the atomic position in PW base
template <typename FPTYPE, typename Device>
void Stress_PW<FPTYPE, Device>::stress_us(ModuleBase::matrix& sigma,
                                          ModulePW::PW_Basis* rho_basis,
                                          const pseudopot_cell_vnl& nlpp,
                                          const UnitCell& ucell)
{
    ModuleBase::TITLE("Stress_Func", "stress_us");
    ModuleBase::timer::tick("Stress_Func", "stress_us");

    const int npw = rho_basis->npw;
    const int nh_tot = nlpp.nhm * (nlpp.nhm + 1) / 2;
    const std::complex<double> fac = ModuleBase::NEG_IMAG_UNIT * ucell.tpiba;
    const std::complex<double> ci_tpi = ModuleBase::IMAG_UNIT * ModuleBase::TWO_PI;
    double* becsum = static_cast<const elecstate::ElecStatePW<std::complex<FPTYPE>, Device>*>(this->pelec)->becsum;

    ModuleBase::matrix stressus(3, 3);

    ModuleBase::matrix veff = this->pelec->pot->get_effective_v();
    ModuleBase::ComplexMatrix vg(PARAM.inp.nspin, npw);
    // fourier transform of the total effective potential
    for (int is = 0; is < PARAM.inp.nspin; is++)
    {
        rho_basis->real2recip(&veff.c[is * veff.nc], &vg(is, 0));
    }

    ModuleBase::matrix ylmk0(nlpp.lmaxq * nlpp.lmaxq, npw);
    ModuleBase::YlmReal::Ylm_Real(nlpp.lmaxq * nlpp.lmaxq, npw, rho_basis->gcar, ylmk0);

    // double* qnorm = new double[npw];
    std::vector<double> qnorm_vec(npw);
    double* qnorm = qnorm_vec.data();
    for (int ig = 0; ig < npw; ig++)
    {
        qnorm[ig] = rho_basis->gcar[ig].norm() * ucell.tpiba;
    }

    // here we compute the integral Q*V for each atom,
    //      I = sum_G G_a exp(-iR.G) Q_nm v^*
    // (no contribution from G=0)
    ModuleBase::matrix dylmk0(nlpp.lmaxq * nlpp.lmaxq, npw);
    for (int ipol = 0; ipol < 3; ipol++)
    {
        double* gcar_ptr = reinterpret_cast<double*>(rho_basis->gcar);
        hamilt::Nonlocal_maths<FPTYPE, Device>::dylmr2(nlpp.lmaxq * nlpp.lmaxq,
                                                       npw,
                                                       gcar_ptr,
                                                       dylmk0.c,
                                                       ipol);
        for (int it = 0; it < ucell.ntype; it++)
        {
            Atom* atom = &ucell.atoms[it];
            if (atom->ncpp.tvanp)
            {
                // nij = max number of (ih,jh) pairs per atom type nt
                // qgm contains derivatives of the Fourier transform of the Q function
                const int nij = atom->ncpp.nh * (atom->ncpp.nh + 1) / 2;
                ModuleBase::ComplexMatrix qgm(nij, npw);
                ModuleBase::matrix tbecsum(PARAM.inp.nspin, nij);

                // Compute and store derivatives of Q(G) for this atomic species
                // (without structure factor)
                int ijh = 0;
                for (int ih = 0; ih < atom->ncpp.nh; ih++)
                {
                    for (int jh = ih; jh < atom->ncpp.nh; jh++)
                    {
                        this->dqvan2(nlpp,
                                     ih,
                                     jh,
                                     it,
                                     ipol,
                                     npw,
                                     rho_basis->gcar,
                                     qnorm,
                                     ucell.tpiba,
                                     ylmk0,
                                     dylmk0,
                                     &qgm(ijh, 0));
                        ijh++;
                    }
                }
                double* qgm_data = reinterpret_cast<double*>(qgm.c);

                for (int ia = 0; ia < atom->na; ia++)
                {
                    const int iat = ucell.itia2iat(it, ia);
                    for (int is = 0; is < PARAM.inp.nspin; is++)
                    {
                        for (int ij = 0; ij < nij; ij++)
                        {
                            tbecsum(is, ij) = becsum[is * ucell.nat * nh_tot + iat * nh_tot + ij];
                        }
                    }

                    ModuleBase::ComplexMatrix aux2(PARAM.inp.nspin, npw);
                    double* aux2_data = reinterpret_cast<double*>(aux2.c);

                    const char transa = 'N';
                    const char transb = 'N';
                    const int dim = 2 * npw;
                    const double one = 1;
                    const double zero = 0;
                    dgemm_(&transa,
                           &transb,
                           &dim,
                           &PARAM.inp.nspin,
                           &nij,
                           &one,
                           qgm_data,
                           &dim,
                           tbecsum.c,
                           &nij,
                           &zero,
                           aux2_data,
                           &dim);

                    for (int is = 0; is < PARAM.inp.nspin; is++)
                    {
                        for (int ig = 0; ig < npw; ig++)
                        {
                            aux2(is, ig) *= conj(vg(is, ig));
                        }
                    }

                    ModuleBase::ComplexMatrix aux1(3, npw);
                    double* aux1_data = reinterpret_cast<double*>(aux1.c);
                    for (int ig = 0; ig < npw; ig++)
                    {
                        double arg = rho_basis->gcar[ig] * atom->tau[ia];
                        std::complex<double> cfac = ucell.tpiba * ModuleBase::libm::exp(ci_tpi * arg);
                        for (int ipol = 0; ipol < 3; ipol++)
                        {
                            aux1(ipol, ig) = cfac * rho_basis->gcar[ig][ipol];
                        }
                    }

                    ModuleBase::matrix fac(PARAM.inp.nspin, 3);
                    const char transc = 'T';
                    const int three = 3;
                    dgemm_(&transc,
                           &transb,
                           &three,
                           &PARAM.inp.nspin,
                           &dim,
                           &one,
                           aux1_data,
                           &dim,
                           aux2_data,
                           &dim,
                           &zero,
                           fac.c,
                           &three);

                    for (int is = 0; is < PARAM.inp.nspin; is++)
                    {
                        for (int jpol = 0; jpol < 3; jpol++)
                        {
                            stressus(jpol, ipol) += fac(is, jpol);
                        }
                    }
                }
            }
        }
    }

    Parallel_Reduce::reduce_all(stressus.c, stressus.nr * stressus.nc);
    for (int l = 0; l < 3; l++)
    {
        for (int m = l; m < 3; m++)
        {
            stressus(m, l) = stressus(l, m);
        }
    }
    sigma += stressus;

    ModuleBase::timer::tick("Stress_Func", "stress_us");
    return;
}

template <typename FPTYPE, typename Device>
void Stress_Func<FPTYPE, Device>::dqvan2(const pseudopot_cell_vnl& nlpp,
                                         const int ih,
                                         const int jh,
                                         const int itype,
                                         const int ipol,
                                         const int ng,
                                         const ModuleBase::Vector3<FPTYPE>* g,
                                         const FPTYPE* qnorm,
                                         const FPTYPE& tpiba,
                                         const ModuleBase::matrix& ylmk0,
                                         const ModuleBase::matrix& dylmk0,
                                         std::complex<FPTYPE>* dqg)
{
    if (PARAM.inp.test_pp) {
        ModuleBase::TITLE("Stress", "dqvan2");
}

    // computes the indices which correspond to ih,jh
    const int nb = nlpp.indv(itype, ih);
    const int mb = nlpp.indv(itype, jh);
    assert(nb < nlpp.nbetam);
    assert(mb < nlpp.nbetam);
    int ijv = 0;
    if (nb >= mb)
    {
        ijv = nb * (nb + 1) / 2 + mb;
    }
    else
    {
        ijv = mb * (mb + 1) / 2 + nb;
    }
    const int ivl = nlpp.nhtolm(itype, ih);
    const int jvl = nlpp.nhtolm(itype, jh);

    for (int ig = 0; ig < ng; ig++)
    {
        dqg[ig] = {0, 0};
    }

    // make the sum over the non zero LM
    int l = -1;
    std::complex<double> pref(0.0, 0.0);
    for (int lm = 0; lm < nlpp.lpx(ivl, jvl); lm++)
    {
        int lp = nlpp.lpl(ivl, jvl, lm);
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
        pref = pow(ModuleBase::NEG_IMAG_UNIT, l) * nlpp.ap(lp, ivl, jvl);

        double qm1 = -1.0; // any number smaller than qnorm
        double work = 0.0, work1 = 0.0;
        for (int ig = 0; ig < ng; ig++)
        {
            if (std::abs(qnorm[ig] - qm1) > 1e-6)
            {
                work = ModuleBase::PolyInt::Polynomial_Interpolation(nlpp.qrad,
                                                                     itype,
                                                                     l,
                                                                     ijv,
                                                                     PARAM.globalv.nqxq,
                                                                     PARAM.globalv.dq,
                                                                     qnorm[ig]);
                work1 = this->Polynomial_Interpolation_nl(nlpp.qrad, itype, l, ijv, PARAM.globalv.dq, qnorm[ig]);
                qm1 = qnorm[ig];
            }
            dqg[ig] += pref * work * dylmk0(lp, ig) / tpiba;
            if (qnorm[ig] > 1e-9)
            {
                dqg[ig] += pref * work1 * ylmk0(lp, ig) * tpiba * g[ig][ipol] / qnorm[ig];
            }
        }
    }
}

template class Stress_PW<double, base_device::DEVICE_CPU>;
#if ((defined __CUDA) || (defined __ROCM))
template class Stress_PW<double, base_device::DEVICE_GPU>;
#endif