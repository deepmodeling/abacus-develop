#include "vnl_pw.h"

#include "source_io/module_parameter/parameter.h"
#include "source_base/clebsch_gordan_coeff.h"
#include "source_base/global_function.h"
#include "source_base/global_variable.h"
#include "source_base/math_integral.h"
#include "source_base/math_polyint.h"
#include "source_base/math_sphbes.h"
#include "source_base/math_ylmreal.h"
#include "source_base/parallel_reduce.h"
#include "source_base/timer.h"

#include <cmath>
#include <vector>

/**
 * @file vnl_pw_init_vnl.cpp
 * @brief Non-local pseudopotential init_vnl: build D terms, qq coefficients
 *        and the beta-function interpolation table.
 *
 * Responsibilities:
 * - set up nhtol/nhtolm/nhtoj/indv index tables
 * - initialize dvan / dvan_so (KB / USPP, with or without SOC)
 * - compute Clebsch-Gordan coefficients and qrad table for USPP
 * - compute qq_nt / qq_so (G=0 Q components) and broadcast via MPI
 * - fill tab interpolation table for beta projectors
 * - stage double/float data to device memory when GPU is enabled
 */

void pseudopot_cell_vnl::init_vnl(UnitCell& cell, const ModulePW::PW_Basis* rho_basis)
{
    ModuleBase::TITLE("pseudopot_cell_vnl", "init_vnl");
    ModuleBase::timer::start("ppcell_vnl", "init_vnl");

    this->omega_old = cell.omega;

    // from init_us_1
    //    a) For each non vanderbilt pseudopotential it computes the D and
    //       the betar in the same form of the Vanderbilt pseudopotential.
    //    b) It computes the indices indv which establish the correspondence
    //       nh <-> beta in the atom
    //    c) It computes the indices nhtol which establish the correspondence
    //       nh <-> angular momentum of the beta function
    //    d) It computes the indices nhtolm which establish the correspondence
    //       nh <-> combined (l,m) index for the beta function.

    // the following prevents an out-of-bound error: upf(nt)%nqlc=2*lmax+1
    // but in some versions of the PP files lmax is not set to the maximum
    // l of the beta functions but includes the l of the local potential
    for (int it = 0; it < cell.ntype; it++)
    {
        if (cell.atoms[it].ncpp.tvanp)
        {
            cell.atoms[it].ncpp.nqlc = std::min(cell.atoms[it].ncpp.nqlc, lmaxq);
            if (cell.atoms[it].ncpp.nqlc < 0) {
                cell.atoms[it].ncpp.nqlc = 0;
}
        }
    }

    // In the spin-orbit case we need the unitary matrix u which rotates the
    // real spherical harmonics and yields the complex ones.
    soc.fcoef.create(cell.ntype, this->nhm, this->nhm);
    if (PARAM.inp.lspinorb)
    {
        soc.rot_ylm(this->lmaxkb);
    }

    // For each pseudopotential we initialize the indices nhtol, nhtolm,
    // nhtoj, indv, and if the pseudopotential is of KB type we initialize
    // the atomic D terms
    this->dvan.zero_out();
    this->dvan_so.zero_out(); // added by zhengdy-soc
    this->indv_ijkb0.assign(cell.nat, 0);
    int ijkb0 = 0;
    for (int it = 0; it < cell.ntype; it++)
    {
        int BetaIndex = 0;
        const int Nprojectors = cell.atoms[it].ncpp.nh;
        for (int ib = 0; ib < cell.atoms[it].ncpp.nbeta; ib++)
        {
            const int l = cell.atoms[it].ncpp.lll[ib];
            const double j = cell.atoms[it].ncpp.jjj[ib];
            for (int m = 0; m < 2 * l + 1; m++)
            {
                this->nhtol(it, BetaIndex) = l;
                this->nhtolm(it, BetaIndex) = l * l + m;
                this->nhtoj(it, BetaIndex) = j;
                this->indv(it, BetaIndex) = ib;
                ++BetaIndex;
            }
        }

        // ijtoh map augmentation channel indexes ih and jh to composite
        // "triangular" index ijh
        for (int ih1 = 0; ih1 < nhm; ih1++)
        {
            for (int ih2 = ih1; ih2 < nhm; ih2++)
            {
                this->ijtoh(it, ih1, ih2) = -1;
                this->ijtoh(it, ih2, ih1) = -1;
            }
        }
        int ijv = 0;
        for (int ih1 = 0; ih1 < Nprojectors; ih1++)
        {
            for (int ih2 = ih1; ih2 < Nprojectors; ih2++)
            {
                ijv++; // liuyu I'm not sure from 0 or 1
                this->ijtoh(it, ih1, ih2) = ijv;
                this->ijtoh(it, ih2, ih1) = ijv;
            }
        }

        // ijkb0 points to the first beta "in the solid" for atom ia
        // i.e. ijkb0,.. ijkb0+nh(ityp(ia))-1 are the nh beta functions of
        // atom ia in the global list of beta functions (ijkb0=0 for ia=0)
        for (int ia = 0; ia < cell.nat; ia++)
        {
            if (it == cell.iat2it[ia])
            {
                this->indv_ijkb0[ia] = ijkb0;
                ijkb0 += Nprojectors;
            }
        }

        // From now on the only difference between KB and US pseudopotentials
        // is in the presence of the q and Q functions.
        // Here we initialize the D of the solid
        if (cell.atoms[it].ncpp.has_so)
        {
            for (int ip = 0; ip < Nprojectors; ip++)
            {
                const int l1 = this->nhtol(it, ip);
                const double j1 = this->nhtoj(it, ip);
                const int m1 = this->nhtolm(it, ip) - l1 * l1;
                // const int v1 = static_cast<int>( indv(it, ip ) );
                for (int ip2 = 0; ip2 < Nprojectors; ip2++)
                {
                    const int l2 = this->nhtol(it, ip2);
                    const double j2 = this->nhtoj(it, ip2);
                    const int m2 = this->nhtolm(it, ip2) - l2 * l2;
                    // const int v2 = static_cast<int>( indv(it, ip2 ) );
                    if (l1 == l2 && fabs(j1 - j2) < 1e-7)
                    {
                        for (int is1 = 0; is1 < 2; is1++)
                        {
                            for (int is2 = 0; is2 < 2; is2++)
                            {
                                soc.set_fcoef(l1, l2, is1, is2, m1, m2, j1, j2, it, ip, ip2);
                            }
                        }
                    }
                }
            }
            //
            //   and calculate the bare coefficients
            //
            for (int ip = 0; ip < Nprojectors; ++ip)
            {
                const int ir = static_cast<int>(indv(it, ip));
                for (int ip2 = 0; ip2 < Nprojectors; ++ip2)
                {
                    const int is = static_cast<int>(indv(it, ip2));
                    int ijs = 0;
                    for (int is1 = 0; is1 < 2; ++is1)
                    {
                        for (int is2 = 0; is2 < 2; ++is2)
                        {
                            this->dvan_so(ijs, it, ip, ip2)
                                = cell.atoms[it].ncpp.dion(ir, is) * soc.fcoef(it, is1, is2, ip, ip2);
                            ++ijs;
                            if (ir != is) {
                                soc.fcoef(it, is1, is2, ip, ip2) = std::complex<double>(0.0, 0.0);
}
                        }
                    }
                }
            }
        }
        else {
            for (int ip = 0; ip < Nprojectors; ip++)
            {
                for (int ip2 = 0; ip2 < Nprojectors; ip2++)
                {
                    if (this->nhtol(it, ip) == nhtol(it, ip2) && this->nhtolm(it, ip) == nhtolm(it, ip2))
                    {
                        const int ir = static_cast<int>(indv(it, ip));
                        const int is = static_cast<int>(indv(it, ip2));
                        if (PARAM.inp.lspinorb)
                        {
                            this->dvan_so(0, it, ip, ip2) = cell.atoms[it].ncpp.dion(ir, is);
                            this->dvan_so(3, it, ip, ip2) = cell.atoms[it].ncpp.dion(ir, is);
                        }
                        else
                        {
                            this->dvan(it, ip, ip2) = cell.atoms[it].ncpp.dion(ir, is);
                        }
                    }
                }
            }
}
    }

    // e) It computes the coefficients c_{LM}^{nm} which relates the
    //    spherical harmonics in the Q expansion
    // f) It computes the interpolation table "qrad" for Q(G)
    // g) It computes the qq terms which define the S matrix.

    // compute Clebsch-Gordan coefficients
    if (PARAM.globalv.use_uspp)
    {
        ModuleBase::Clebsch_Gordan::clebsch_gordan(lmaxkb + 1, ap, lpx, lpl);
    }

    // here for the US types we compute the Fourier transform of the Q functions.
    if (lmaxq > 0)
    {
        this->compute_qrad(cell);
    }

    // compute the qq coefficients by integrating the Q.
    // The qq are the g=0 components of Q
    if (rho_basis->ig_gge0 >= 0)
    {
        ModuleBase::matrix ylmk0(lmaxq * lmaxq, 1);
        const double qnorm = 0.0; // only G=0 term
        std::complex<double> qgm(0.0, 0.0);
        ModuleBase::YlmReal::Ylm_Real(lmaxq * lmaxq, 1, &(rho_basis->gcar[rho_basis->ig_gge0]), ylmk0);
        for (int it = 0; it < cell.ntype; it++)
        {
            Atom_pseudo* upf = &cell.atoms[it].ncpp;
            if (upf->tvanp)
            {
                if (upf->has_so)
                {
                    for (int ih = 0; ih < upf->nh; ih++)
                    {
                        for (int jh = 0; jh < upf->nh; jh++)
                        {
                            this->radial_fft_q(1, ih, jh, it, &qnorm, ylmk0, &qgm);
                            this->qq_nt(it, ih, jh) = cell.omega * qgm.real();
                            for (int kh = 0; kh < upf->nh; kh++)
                            {
                                for (int lh = 0; lh < upf->nh; lh++)
                                {
                                    int ijs = 0;
                                    for (int is1 = 0; is1 < 2; ++is1)
                                    {
                                        for (int is2 = 0; is2 < 2; ++is2)
                                        {
                                            for (int is = 0; is < 2; is++)
                                            {
                                                this->qq_so(it, ijs, kh, lh) += cell.omega * qgm.real()
                                                                                * soc.fcoef(it, is1, is, kh, ih)
                                                                                * soc.fcoef(it, is, is2, jh, lh);
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                else
                {
                    for (int ih = 0; ih < upf->nh; ih++)
                    {
                        for (int jh = ih; jh < upf->nh; jh++)
                        {
                            this->radial_fft_q(1, ih, jh, it, &qnorm, ylmk0, &qgm);
                            if (PARAM.inp.lspinorb)
                            {
                                this->qq_so(it, 0, ih, jh) = cell.omega * qgm.real();
                                this->qq_so(it, 0, jh, ih) = this->qq_so(it, 0, ih, jh);
                                this->qq_so(it, 3, ih, jh) = this->qq_so(it, 0, ih, jh);
                                this->qq_so(it, 3, jh, ih) = this->qq_so(it, 0, ih, jh);
                            }
                            this->qq_nt(it, ih, jh) = cell.omega * qgm.real();
                            this->qq_nt(it, jh, ih) = cell.omega * qgm.real();
                        }
                    }
                }
            }
        }
    }

#ifdef __MPI
    Parallel_Reduce::reduce_pool(this->qq_nt.ptr, this->qq_nt.getSize());
    Parallel_Reduce::reduce_pool(this->qq_so.ptr, this->qq_so.getSize());
#endif

    // set the atomic specific qq_at matrices
    for (int ia = 0; ia < cell.nat; ia++)
    {
        int it = cell.iat2it[ia];
        for (int ih = 0; ih < nhm; ih++)
        {
            for (int jh = ih; jh < nhm; jh++)
            {
                this->qq_at(ia, ih, jh) = qq_nt(it, ih, jh);
                this->qq_at(ia, jh, ih) = qq_nt(it, jh, ih);
            }
        }
    }

    // h) It fills the interpolation table for the beta functions
    /**********************************************************
    // He Lixin: this block is used for non-local potential
    // fill the interpolation table tab
    ************************************************************/

    const double pref = ModuleBase::FOUR_PI / sqrt(cell.omega);
    this->tab.zero_out();
    GlobalV::ofs_running << "\n Init Non-Local PseudoPotential table : ";
    for (int it = 0; it < cell.ntype; it++)
    {
        const int nbeta = cell.atoms[it].ncpp.nbeta;
        int kkbeta = cell.atoms[it].ncpp.kkbeta;

        // mohan modify 2008-3-31
        // mohan add kkbeta>0 2009-2-27
        if ((kkbeta % 2 == 0) && kkbeta > 0)
        {
            kkbeta--;
        }

        std::vector<double> jl(kkbeta);
        std::vector<double> aux(kkbeta);
        for (int ib = 0; ib < nbeta; ib++)
        {
            const int l = cell.atoms[it].ncpp.lll[ib];
            for (int iq = 0; iq < PARAM.globalv.nqx; iq++)
            {
                const double q = iq * PARAM.globalv.dq;
                ModuleBase::Sphbes::Spherical_Bessel(kkbeta, cell.atoms[it].ncpp.r.data(), q, l, jl.data());
                for (int ir = 0; ir < kkbeta; ir++)
                {
                    aux[ir] = cell.atoms[it].ncpp.betar(ib, ir) * jl[ir] * cell.atoms[it].ncpp.r[ir];
                }
                double vqint=0.0;
                ModuleBase::Integral::Simpson_Integral(kkbeta, aux.data(), cell.atoms[it].ncpp.rab.data(), vqint);
                this->tab(it, ib, iq) = vqint * pref;
            }
        }
    }
    if (this->use_gpu_)
    {
        if (PARAM.globalv.has_float_data)
        {
            castmem_d2s_h2d_op()(this->s_indv, this->indv.c, this->indv.nr * this->indv.nc);
            castmem_d2s_h2d_op()(this->s_nhtol, this->nhtol.c, this->nhtol.nr * this->nhtol.nc);
            castmem_d2s_h2d_op()(this->s_nhtolm, this->nhtolm.c, this->nhtolm.nr * this->nhtolm.nc);
            castmem_d2s_h2d_op()(this->s_tab, this->tab.ptr, this->tab.getSize());
            if (this->nhm > 0)
            {
                castmem_d2s_h2d_op()(this->s_qq_nt, this->qq_nt.ptr, this->qq_nt.getSize());
                castmem_z2c_h2d_op()(this->c_qq_so, this->qq_so.ptr, this->qq_so.getSize());
            }
        }
        if (PARAM.globalv.has_double_data && this->nhm > 0)
        {
            syncmem_z2z_h2d_op()(this->z_qq_so, this->qq_so.ptr, this->qq_so.getSize());
        }
        // Even when the single precision flag is enabled,
        // these variables are utilized in the Force/Stress calculation as well.
        // modified by denghuilu at 2023-05-15
        syncmem_d2d_h2d_op()(this->d_indv, this->indv.c, this->indv.nr * this->indv.nc);
        syncmem_d2d_h2d_op()(this->d_nhtol, this->nhtol.c, this->nhtol.nr * this->nhtol.nc);
        syncmem_d2d_h2d_op()(this->d_nhtolm, this->nhtolm.c, this->nhtolm.nr * this->nhtolm.nc);
        syncmem_d2d_h2d_op()(this->d_tab, this->tab.ptr, this->tab.getSize());
        if (this->nhm > 0)
        {
            syncmem_d2d_h2d_op()(this->d_qq_nt, this->qq_nt.ptr, this->qq_nt.getSize());
        }
    }
    else
    {
        if (PARAM.globalv.has_float_data)
        {
            castmem_d2s_h2h_op()(this->s_indv, this->indv.c, this->indv.nr * this->indv.nc);
            castmem_d2s_h2h_op()(this->s_nhtol, this->nhtol.c, this->nhtol.nr * this->nhtol.nc);
            castmem_d2s_h2h_op()(this->s_nhtolm, this->nhtolm.c, this->nhtolm.nr * this->nhtolm.nc);
            castmem_d2s_h2h_op()(this->s_tab, this->tab.ptr, this->tab.getSize());
            if (this->nhm > 0)
            {
                castmem_d2s_h2h_op()(this->s_qq_nt, this->qq_nt.ptr, this->qq_nt.getSize());
                castmem_z2c_h2h_op()(this->c_qq_so, this->qq_so.ptr, this->qq_so.getSize());
            }
        }
        // There's no need to synchronize double precision pointers while in a CPU environment.
    }
    ModuleBase::timer::end("ppcell_vnl", "init_vnl");
    GlobalV::ofs_running << "\n Init Non-Local-Pseudopotential done." << std::endl;
    return;
}
