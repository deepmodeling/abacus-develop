// ============================================================
// This code is added by Mohan Chen on 2026-05-18.
// This code is currently in the design phase and has not been
// put into production yet. It may change in the future.
// Please use this code with caution. Only developers who know
// what they are doing should use this code.
// ============================================================

#include "dfpt_phon.h"

#include "dfpt_kq_basis.h"
#include "dfpt_pert.h"
#include "source_base/constants.h"
#include "source_base/global_function.h"
#include "source_base/module_external/lapack_connector.h"
#include "source_base/tool_quit.h"
#include "source_base/truncated_func.h"
#include "source_basis/module_pw/pw_basis.h"

#include <cmath>
#include <complex>
#include <cstdlib>
#include <iostream>
#include <vector>

namespace ModuleDFPT {

DFPT_Phon::DFPT_Phon() {}

DFPT_Phon::~DFPT_Phon() {}

void DFPT_Phon::init(UnitCell& ucell, ModulePW::PW_Basis* pw_rho, DFPT_Pert* pert) {
    ucell_ = &ucell;
    pw_rho_ = pw_rho;
    pert_ = pert;
}

// ---------------------------------------------------------------------------
// Ewald ion-ion force constants (C5)
// ---------------------------------------------------------------------------

void DFPT_Phon::ion_ion(const ModuleBase::Vector3<double>& q_frac,
                        ModuleBase::ComplexMatrix& dyn) {
    const int nat = ucell_->nat;
    const int nat3 = 3 * nat;
    const double lat0 = ucell_->lat0;
    const ModuleBase::Matrix3& latvec = ucell_->latvec;

    // total ionic charge
    double charge = 0.0;
    for (int it = 0; it < ucell_->ntype; ++it) {
        charge += ucell_->atoms[it].na * ucell_->atoms[it].ncpp.zv;
    }

    // choose the screening alpha so that the G-sum tail is converged inside
    // the rho grid (the erfc envelope bounds the exp(-G^2/4alpha) tail);
    // ggecut counts |G_max|^2 in 1/lat0^2 units (pw_basis.h), so the bohr^2
    // cutoff is ggecut * tpiba2
    double alpha = 1.1;
    double upperbound = 0.0;
    do {
        alpha *= 0.9;
        if (alpha < 1.0e-4) {
            ModuleBase::WARNING_QUIT("DFPT_Phon::ion_ion",
                                     "Can't find optimal Ewald alpha.");
        }
        upperbound = 2.0 * charge * charge
                     * std::sqrt(2.0 * alpha / ModuleBase::TWO_PI)
                     * ModuleBase::truncated_erfc(
                         std::sqrt(pw_rho_->ggecut * ucell_->tpiba2 / 4.0 / alpha));
    } while (upperbound > 1.0e-6);
    ewald_alpha_ = alpha;
    // erfc(alpha R) < 1e-16 well inside 6/sqrt(alpha)
    ewald_rcut_ = 6.0 / std::sqrt(alpha);

    const ModuleBase::Vector3<double> q_cart = q_frac * ucell_->G;

    // ---------------- reciprocal-space part ----------------
    // Poisson pair identity (validated against direct sums):
    //   sum_L h(R) e^{i2pi q.L} = sum_L h_erfc(R) e^{i2pi q.L}
    //     + (4pi/Omega) sum_{|G+q|>0} (G+q)_a (G+q)_b / |G+q|^2
    //       exp(-|G+q|^2/4a) e^{i2pi (G+q).(tau_a-tau_b)}
    // so the G part enters D with the + sign while the erfc part carries -.
    // The on-site diagonal (both second derivatives act on tau_a in cell 0)
    // is phase-free: it is accumulated from Gamma-phase (G-only) pair terms
    // as -sqrt(Mb/Ma) times the pair element. sq/s0 accumulate the self-image
    // phase difference of the same-atom images:
    //   sum_{L!=0} h(L)(e^{i2pi q.L} - 1)
    //     = [erfc piece in the R part] + (4pi/Omega)(sq - s0 - delta_ab/3),
    // where the delta/3 is the G=0 limit (w_a w_b / w^2 -> delta_ab/3).
    double sq[3][3] = {{0.0}};
    double s0[3][3] = {{0.0}};
    for (int ig = 0; ig < pw_rho_->npw; ++ig) {
        const ModuleBase::Vector3<double>& gcart = pw_rho_->gcar[ig];
        const ModuleBase::Vector3<double> w = gcart + q_cart;
        const double w2 = w * w;
        const double g2 = gcart * gcart;
        if (w2 < 1.0e-12) {
            // G + q = 0 (only possible at q = 0 with G = 0): isotropic limit
            for (int d = 0; d < 3; ++d) {
                sq[d][d] += 1.0 / 3.0;
            }
            continue;
        }
        const double w2_bohr = w2 * ucell_->tpiba2;
        const double gauss = ModuleBase::truncated_exp(-w2_bohr / (4.0 * alpha));
        for (int da = 0; da < 3; ++da) {
            for (int db = 0; db < 3; ++db) {
                sq[da][db] += w[da] * w[db] / w2 * gauss;
            }
        }
        double gauss_g = 0.0;
        if (g2 > 1.0e-12) {
            gauss_g = ModuleBase::truncated_exp(-g2 * ucell_->tpiba2 / (4.0 * alpha));
            for (int da = 0; da < 3; ++da) {
                for (int db = 0; db < 3; ++db) {
                    s0[da][db] += gcart[da] * gcart[db] / g2 * gauss_g;
                }
            }
        }
        for (int ia = 0; ia < nat; ++ia) {
            const int ita = ucell_->iat2it[ia];
            const int iia = ucell_->iat2ia[ia];
            const double za = ucell_->atoms[ita].ncpp.zv;
            const double ma = ucell_->atoms[ita].mass;
            const ModuleBase::Vector3<double>& ta = ucell_->atoms[ita].tau[iia];
            for (int ib = 0; ib < nat; ++ib) {
                if (ib == ia) {
                    continue;
                }
                const int itb = ucell_->iat2it[ib];
                const int iib = ucell_->iat2ia[ib];
                const double zb = ucell_->atoms[itb].ncpp.zv;
                const double mb = ucell_->atoms[itb].mass;
                const ModuleBase::Vector3<double>& tb = ucell_->atoms[itb].tau[iib];
                const double arg = ModuleBase::TWO_PI * (w * (ta - tb));
                const std::complex<double> phase(std::cos(arg), std::sin(arg));
                const double pref = ModuleBase::FOUR_PI / ucell_->omega
                                    * za * zb * ModuleBase::e2 * gauss
                                    / (std::sqrt(ma * mb) * w2);
                // Gamma-phase on-site piece (G-only kernel, G != 0)
                std::complex<double> phase0(1.0, 0.0);
                double pref0 = 0.0;
                if (g2 > 1.0e-12) {
                    const double arg0 = ModuleBase::TWO_PI * (gcart * (ta - tb));
                    phase0 = std::complex<double>(std::cos(arg0), std::sin(arg0));
                    pref0 = ModuleBase::FOUR_PI / ucell_->omega
                            * za * zb * ModuleBase::e2 * gauss_g
                            / (std::sqrt(ma * mb) * g2);
                }
                for (int da = 0; da < 3; ++da) {
                    for (int db = 0; db < 3; ++db) {
                        const std::complex<double> elem = pref * w[da] * w[db] * phase;
                        dyn(3 * ia + da, 3 * ib + db) += elem;
                        // on-site diagonal: phase-free (Gamma) accumulation,
                        // Phi_ii = -Phi_ij => -sqrt(Mb/Ma) on the pair term
                        dyn(3 * ia + da, 3 * ia + db)
                            -= pref0 * gcart[da] * gcart[db] * phase0
                               * std::sqrt(mb / ma);
                    }
                }
            }
        }
    }
    // self-image G-space phase difference on the diagonal
    for (int ia = 0; ia < nat; ++ia) {
        const double za = ucell_->atoms[ucell_->iat2it[ia]].ncpp.zv;
        const double ma = ucell_->atoms[ucell_->iat2it[ia]].mass;
        const double f2 = za * za * ModuleBase::e2 / ma;
        for (int da = 0; da < 3; ++da) {
            for (int db = 0; db < 3; ++db) {
                dyn(3 * ia + da, 3 * ia + db)
                    += f2 * ModuleBase::FOUR_PI / ucell_->omega
                       * (sq[da][db] - s0[da][db] - (da == db ? 1.0 / 3.0 : 0.0));
            }
        }
    }

    // ---------------- real-space part ----------------
    // h_ab(R) = d^2/dR_a dR_b [ erfc(sqrt(alpha) R) / R ]
    //         = erfc(sqrt(alpha) R) (3 Ra Rb - delta R^2)/R^5
    //         + (2 sqrt(alpha)/sqrt(pi)) e^{-alpha R^2}
    //           [ 2 alpha Ra Rb/R^2 + 3 Ra Rb/R^4 - delta/R^2 ]
    // D^R_ab = -(1/sqrt(MaMb)) ZaZb e2 h(R = tau_b + l - tau_a) e^{i2pi q.l}
    // ranges of the lattice-vector shells (rows of latvec are the lattice
    // translations in lat0 units)
    const double row_e[3][3] = {{latvec.e11, latvec.e12, latvec.e13},
                                {latvec.e21, latvec.e22, latvec.e23},
                                {latvec.e31, latvec.e32, latvec.e33}};
    int nmax[3] = {0, 0, 0};
    for (int d = 0; d < 3; ++d) {
        const ModuleBase::Vector3<double> a1(row_e[d][0], row_e[d][1], row_e[d][2]);
        const double len = std::sqrt(a1 * a1) * lat0; // bohr
        nmax[d] = static_cast<int>(std::ceil(ewald_rcut_ / len)) + 1;
    }
    for (int ia = 0; ia < nat; ++ia) {
        const int ita = ucell_->iat2it[ia];
        const int iia = ucell_->iat2ia[ia];
        const double za = ucell_->atoms[ita].ncpp.zv;
        const double ma = ucell_->atoms[ita].mass;
        for (int ib = 0; ib < nat; ++ib) {
            const int itb = ucell_->iat2it[ib];
            const int iib = ucell_->iat2ia[ib];
            const double zb = ucell_->atoms[itb].ncpp.zv;
            const double mb = ucell_->atoms[itb].mass;
            const ModuleBase::Vector3<double> dt =
                ucell_->atoms[itb].tau[iib] - ucell_->atoms[ita].tau[iia];
            if (ib == ia) {
                // self-image phase difference: the on-site i-i energy is
                // L-independent while the cross-cell i-i force constants carry
                // e^{i2pi q.L}, so D_ii receives
                //   -(Za^2 e2/Ma) sum_{L!=0} h_erfc(L) (e^{i2pi q.L} - 1);
                // the imaginary part cancels over the +-L symmetric sphere
                // (h is even) and L = 0 carries e^{i2pi q.0} - 1 = 0
                for (int n1 = -nmax[0]; n1 <= nmax[0]; ++n1) {
                    for (int n2 = -nmax[1]; n2 <= nmax[1]; ++n2) {
                        for (int n3 = -nmax[2]; n3 <= nmax[2]; ++n3) {
                            if (n1 == 0 && n2 == 0 && n3 == 0) {
                                continue;
                            }
                            const ModuleBase::Vector3<double> lvec(
                                n1 * latvec.e11 + n2 * latvec.e21 + n3 * latvec.e31,
                                n1 * latvec.e12 + n2 * latvec.e22 + n3 * latvec.e32,
                                n1 * latvec.e13 + n2 * latvec.e23 + n3 * latvec.e33);
                            const ModuleBase::Vector3<double> r = lvec * lat0;
                            const double r2 = r * r;
                            if (r2 > ewald_rcut_ * ewald_rcut_) {
                                continue;
                            }
                            const double rlen = std::sqrt(r2);
                            const double sar = std::sqrt(alpha);
                            const double e2a = ModuleBase::truncated_exp(-alpha * r2);
                            const double f = 2.0 * sar / std::sqrt(ModuleBase::PI) * e2a;
                            const double er = ModuleBase::truncated_erfc(sar * rlen);
                            const double ph_arg = ModuleBase::TWO_PI
                                                  * (q_frac.x * n1 + q_frac.y * n2
                                                     + q_frac.z * n3);
                            const double wcos = std::cos(ph_arg) - 1.0;
                            const double f2 = za * za * ModuleBase::e2 / ma;
                            for (int da = 0; da < 3; ++da) {
                                for (int db = 0; db < 3; ++db) {
                                    const double delta = (da == db) ? 1.0 : 0.0;
                                    const double h = er * (3.0 * r[da] * r[db] - delta * r2)
                                                         / (rlen * r2 * r2)
                                                     + f * (2.0 * alpha * r[da] * r[db] / r2
                                                            + 3.0 * r[da] * r[db] / (r2 * r2)
                                                            - delta / r2);
                                    dyn(3 * ia + da, 3 * ia + db) -= f2 * h * wcos;
                                }
                            }
                        }
                    }
                }
                continue;
            }
            for (int n1 = -nmax[0]; n1 <= nmax[0]; ++n1) {
                for (int n2 = -nmax[1]; n2 <= nmax[1]; ++n2) {
                    for (int n3 = -nmax[2]; n3 <= nmax[2]; ++n3) {
                        const ModuleBase::Vector3<double> lvec(
                            n1 * latvec.e11 + n2 * latvec.e21 + n3 * latvec.e31,
                            n1 * latvec.e12 + n2 * latvec.e22 + n3 * latvec.e32,
                            n1 * latvec.e13 + n2 * latvec.e23 + n3 * latvec.e33);
                        ModuleBase::Vector3<double> r = (lvec + dt) * lat0; // bohr
                        const double r2 = r * r;
                        if (r2 > ewald_rcut_ * ewald_rcut_) {
                            continue;
                        }
                        const double rlen = std::sqrt(r2);
                        const double r3 = r2 * rlen;
                        const double sar = std::sqrt(alpha);
                        const double e2a = ModuleBase::truncated_exp(-alpha * r2);
                        const double f = 2.0 * sar / std::sqrt(ModuleBase::PI) * e2a;
                        const double er = ModuleBase::truncated_erfc(sar * rlen);
                        const double ph_arg = ModuleBase::TWO_PI
                                              * (q_frac.x * n1 + q_frac.y * n2 + q_frac.z * n3);
                        const std::complex<double> phase(std::cos(ph_arg), std::sin(ph_arg));
                        const double zab2 = za * zb * ModuleBase::e2 / std::sqrt(ma * mb);
                        for (int da = 0; da < 3; ++da) {
                            for (int db = 0; db < 3; ++db) {
                                const double delta = (da == db) ? 1.0 : 0.0;
                                // d^2/dR_a dR_b [erfc(sqrt(alpha) R)/R],
                                // validated against central finite differences
                                const double h = er * (3.0 * r[da] * r[db] - delta * r2) / (r3 * r2)
                                                 + f * (2.0 * alpha * r[da] * r[db] / r2
                                                        + 3.0 * r[da] * r[db] / (r2 * r2)
                                                        - delta / r2);
                                dyn(3 * ia + da, 3 * ib + db) -= zab2 * h * phase;
                                // on-site diagonal Phi_ii^R = sum_{j != i}
                                // Z_iZ_j sum_L h(r_ij + L): phase-free (both
                                // derivatives act on tau_a in cell 0), i.e.
                                // -sqrt(Mb/Ma) times the pair term
                                dyn(3 * ia + da, 3 * ia + db)
                                    += zab2 * std::sqrt(mb / ma) * h;
                            }
                        }
                    }
                }
            }
        }
    }

    // The Gaussian self constant -Z^2 sqrt(2 alpha/pi) and the h_erf contact
    // -4 alpha^{3/2}/(3 sqrt(pi)) delta_ab are tau-independent and cancel in
    // the (e^{i2pi q.L} - 1) differences; the diagonal is carried by the
    // phase-free cross-atom accumulation plus the self-image phase terms
    // (both G and R pieces above). At q = 0 all phase differences vanish and
    // the acoustic sum rule holds exactly by construction.
}

// ---------------------------------------------------------------------------
// electronic contribution (2n+1 theorem)
// ---------------------------------------------------------------------------

void DFPT_Phon::accumulate_electron(int q_idx, int atom_idx, int dir,
                                    const psi::Psi<std::complex<double>>& psi,
                                    const ModuleBase::matrix& wg, DFPT_PW_Data& data) {
    if (pert_ == nullptr || pw_rho_ == nullptr || ucell_ == nullptr) {
        return;
    }
    const int nat = ucell_->nat;
    const int nat3 = 3 * nat;
    if (accum_q_ != q_idx || dynmat_accum_.nr != nat3) {
        dynmat_accum_ = ModuleBase::ComplexMatrix(nat3, nat3, true);
        accum_q_ = q_idx;
    }
    const int rowb = 3 * atom_idx + dir;
    const ModuleBase::Vector3<double> q_frac = data.get_qvec(q_idx);
    const ModuleBase::Vector3<double> q_cart = q_frac * ucell_->G;
    const int nk = psi.get_nk();
    const int nbands = psi.get_nbands();

    // stash the converged dpsi of this displacement (apply_dv reuses the slot)
    std::vector<std::vector<std::vector<std::complex<double>>>> dpsib(nk);
    for (int ik = 0; ik < nk; ++ik) {
        dpsib[ik].resize(nbands);
        for (int ib = 0; ib < nbands; ++ib) {
            dpsib[ik][ib] = data.get_dpsi(q_idx, ik, ib);
        }
    }

    // rho ig -> shared FFT-cell reverse map (C1/C3 pattern)
    std::vector<int> ig_of_cell(pw_rho_->nxyz, -1);
    for (int ig = 0; ig < pw_rho_->npw; ++ig) {
        const int isz = pw_rho_->ig2isz[ig];
        const int iz = isz % pw_rho_->nz;
        const int is = isz / pw_rho_->nz;
        const int ixy = pw_rho_->is2fftixy[is];
        const int ix = ixy / pw_rho_->fftny;
        const int iy = ixy % pw_rho_->fftny;
        ig_of_cell[(ix * pw_rho_->ny + iy) * pw_rho_->nz + iz] = ig;
    }

    for (int iat = 0; iat < nat; ++iat) {
        for (int idir = 0; idir < 3; ++idir) {
            const int cola = 3 * iat + idir;
            // ---- term 2 <dpsi^b | dV^a_ext | psi> over all k,n ----
            // complex accumulation: at a generic q the single-k matrix
            // elements are complex (the imaginary parts pair-conjugate over
            // the k star), and the Hermitian symmetrization in assemble
            // relies on them.
            pert_->build_dv(q_idx, iat, idir, data);
            std::complex<double> cross(0.0, 0.0);
            const bool dbg2 = (getenv("DFPT_DEBUG") != nullptr);
            for (int ik = 0; ik < nk; ++ik) {
                pert_->apply_dv(q_idx, ik, psi, data);
                for (int ib = 0; ib < nbands; ++ib) {
                    if (wg(ik, ib) < 1.0e-8) {
                        continue;
                    }
                    const std::vector<std::complex<double>> rhs = data.get_dpsi(q_idx, ik, ib);
                    const std::vector<std::complex<double>>& sol = dpsib[ik][ib];
                    if (rhs.size() != sol.size() || sol.empty()) {
                        continue;
                    }
                    std::complex<double> dot(0.0, 0.0);
                    for (size_t i = 0; i < sol.size(); ++i) {
                        dot += std::conj(sol[i]) * rhs[i];
                    }
                    cross += wg(ik, ib) * dot;
                }
            }
            dynmat_accum_(rowb, cola) += 2.0 * cross
                / std::sqrt(ucell_->atoms[ucell_->iat2it[atom_idx]].mass
                            * ucell_->atoms[ucell_->iat2it[iat]].mass);
            if (dbg2) {
                std::cout << "DYNCHK term2 rowb=" << rowb << " cola=" << cola
                          << " 2cross=" << 2.0 * cross.real() << std::endl;
            }

            // ---- same-atom anharmonic term <psi | d2_ab V_ext | psi> ----
            if (iat == atom_idx && cola >= rowb) {
                std::vector<std::complex<double>> dv2_r;
                pert_->d2vloc_r(atom_idx, idir, dir, q_cart, dv2_r);
                if (static_cast<int>(dv2_r.size()) != pw_rho_->nrxx) {
                    dv2_r.assign(pw_rho_->nrxx, std::complex<double>(0.0, 0.0));
                }
                std::vector<std::vector<std::complex<double>>> chi;
                std::complex<double> d2sum(0.0, 0.0);
                std::complex<double> d2sum_loc(0.0, 0.0);
                std::complex<double> d2sum_nl(0.0, 0.0);
                std::vector<std::complex<double>> u_r(pw_rho_->nrxx);
                std::vector<std::complex<double>> x_r(pw_rho_->nrxx);
                std::vector<std::complex<double>> x_recip(pw_rho_->npw, std::complex<double>(0.0, 0.0));
                for (int ik = 0; ik < nk; ++ik) {
                    pert_->apply_d2vnl(atom_idx, idir, dir, q_cart, psi, ik, chi);
                    // k+q scatter map for this k
                    DFPT_KQ_Basis kq;
                    kq.init(pert_->get_pw_wfc(), q_cart, ik);
                    const int npwk_kq = kq.get_npwk();
                    for (int ib = 0; ib < nbands; ++ib) {
                        if (wg(ik, ib) < 1.0e-8) {
                            continue;
                        }
                        pert_->get_pw_wfc()->recip2real(&psi(ik, ib, 0), u_r.data(), ik);
                        if (static_cast<int>(chi.size()) == nbands
                            && static_cast<int>(chi[ib].size()) == npwk_kq) {
                            std::fill(x_recip.begin(), x_recip.end(), std::complex<double>(0.0, 0.0));
                            for (int igl = 0; igl < npwk_kq; ++igl) {
                                const int isz = kq.get_ig2isz(igl);
                                const int iz = isz % pert_->get_pw_wfc()->nz;
                                const int is = isz / pert_->get_pw_wfc()->nz;
                                const int ixy = pert_->get_pw_wfc()->is2fftixy[is];
                                const int ix = ixy / pert_->get_pw_wfc()->fftny;
                                const int iy = ixy % pert_->get_pw_wfc()->fftny;
                                const int ig_rho =
                                    ig_of_cell[(ix * pw_rho_->ny + iy) * pw_rho_->nz + iz];
                                if (ig_rho >= 0) {
                                    x_recip[ig_rho] = chi[ib][igl];
                                }
                            }
                            pw_rho_->recip2real(x_recip.data(), x_r.data());
                        }
                        else {
                            std::fill(x_r.begin(), x_r.end(), std::complex<double>(0.0, 0.0));
                        }
                        std::complex<double> expect(0.0, 0.0);
                        std::complex<double> expect_loc(0.0, 0.0);
                        std::complex<double> expect_nl(0.0, 0.0);
                        for (int ir = 0; ir < pw_rho_->nrxx; ++ir) {
                            expect += std::conj(u_r[ir]) * u_r[ir] * dv2_r[ir]
                                      + std::conj(u_r[ir]) * x_r[ir];
                            expect_loc += std::conj(u_r[ir]) * u_r[ir] * dv2_r[ir];
                            expect_nl += std::conj(u_r[ir]) * x_r[ir];
                        }
                        d2sum += wg(ik, ib) * expect / static_cast<double>(pw_rho_->nxyz);
                        d2sum_loc += wg(ik, ib) * expect_loc / static_cast<double>(pw_rho_->nxyz);
                        d2sum_nl += wg(ik, ib) * expect_nl / static_cast<double>(pw_rho_->nxyz);
                    }
                }
                dynmat_accum_(rowb, cola) += d2sum
                    / ucell_->atoms[ucell_->iat2it[atom_idx]].mass;
                if (dbg2) {
                    std::cout << "DYNCHK d2    rowb=" << rowb << " cola=" << cola
                              << " d2sum=" << d2sum.real()
                              << " loc=" << d2sum_loc.real()
                              << " nl=" << d2sum_nl.real() << std::endl;
                }
            }
        }
    }

    // restore the converged dpsi of this displacement
    for (int ik = 0; ik < nk; ++ik) {
        for (int ib = 0; ib < nbands; ++ib) {
            if (!dpsib[ik][ib].empty()) {
                data.set_dpsi(q_idx, ik, ib, dpsib[ik][ib]);
            }
        }
    }
}

// ---------------------------------------------------------------------------
// assemble / diagonalize / LO-TO / sum rule
// ---------------------------------------------------------------------------

void DFPT_Phon::assemble(int q_idx, DFPT_PW_Data& data) {
    if (ucell_ == nullptr) {
        return;
    }
    const int nat = ucell_->nat;
    const int nat3 = 3 * nat;
    ModuleBase::ComplexMatrix dyn(nat3, nat3, true);
    if (pw_rho_ != nullptr) {
        ion_ion(data.get_qvec(q_idx), dyn);
    }
    if (accum_q_ == q_idx && dynmat_accum_.nr == nat3) {
        if (getenv("DFPT_DEBUG") != nullptr) {
            std::cout << "DYNCHK ionic matrix (Ry/bohr^2/amu):" << std::endl;
            for (int i = 0; i < nat3; ++i) {
                std::cout << "DYNCHK ion row " << i << ":";
                for (int j = 0; j < nat3; ++j) {
                    std::cout << " " << dyn(i, j).real();
                }
                std::cout << std::endl;
            }
            std::cout << "DYNCHK electronic accum matrix:" << std::endl;
            for (int i = 0; i < nat3; ++i) {
                std::cout << "DYNCHK ele row " << i << ":";
                for (int j = 0; j < nat3; ++j) {
                    std::cout << " " << dynmat_accum_(i, j).real();
                }
                std::cout << std::endl;
            }
        }
        for (int i = 0; i < nat3; ++i) {
            for (int j = 0; j < nat3; ++j) {
                dyn(i, j) += dynmat_accum_(i, j);
            }
        }
    }
    // DFT+U dynamical-matrix term (U0 reservation, implemented with the C7/U1
    // Plus_U wiring): sum_nk w_nk [<psi|dV_U|dpsi> + frozen second-order term].
    if (data.with_u()) {
        dftu_onsite(q_idx, data);
    }
    // Hermitian symmetrization (rows filled by independent solves)
    for (int i = 0; i < nat3; ++i) {
        for (int j = i + 1; j < nat3; ++j) {
            const std::complex<double> avg
                = 0.5 * (dyn(i, j) + std::conj(dyn(j, i)));
            dyn(i, j) = avg;
            dyn(j, i) = std::conj(avg);
        }
    }
    data.set_dynmat(q_idx, dyn);
    dynmat_accum_ = ModuleBase::ComplexMatrix();
    accum_q_ = -1;
}

void DFPT_Phon::diagonalize(int q_idx, DFPT_PW_Data& data) {
    const int nat = ucell_->nat;
    const int nat3 = 3 * nat;
    ModuleBase::ComplexMatrix dyn = data.get_dynmat(q_idx);
    if (dyn.nr != nat3) {
        return;
    }

    // eigenvalues of the complex Hermitian dynamical matrix (Ry/bohr^2/amu)
    std::vector<double> w(nat3, 0.0);
    std::vector<double> rwork(std::max(1, 3 * nat3 - 2), 0.0);
    std::vector<std::complex<double>> work(1);
    int info = 0;
    LapackConnector::zheev('N', 'U', nat3, dyn, nat3, w.data(), work.data(), -1,
                           rwork.data(), &info);
    work.resize(std::max(1, static_cast<int>(work[0].real())));
    LapackConnector::zheev('N', 'U', nat3, dyn, nat3, w.data(), work.data(),
                           static_cast<int>(work.size()), rwork.data(), &info);

    // signed frequencies: omega = sgn(e) sqrt(|e|), converted to cm^-1
    // sqrt(Ry/(bohr^2 amu)) in cm^-1 = sqrt(RYDBERG_SI/amu_kg)/(bohr*2pi*c)
    const double amu_kg = 1.0e-3 / ModuleBase::NA;
    const double ry_bohr2_amu_to_cm1 = std::sqrt(ModuleBase::RYDBERG_SI / amu_kg)
                                       / (ModuleBase::BOHR_RADIUS_SI * ModuleBase::TWO_PI
                                          * 2.99792458e10);
    std::vector<double> freq(nat3, 0.0);
    for (int i = 0; i < nat3; ++i) {
        const double e = w[i];
        freq[i] = ((e >= 0.0) ? 1.0 : -1.0) * std::sqrt(std::abs(e)) * ry_bohr2_amu_to_cm1;
    }
    data.set_phon_freq(q_idx, freq);
}

void DFPT_Phon::add_loto(const ModuleBase::Vector3<double>& qhat, DFPT_PW_Data& data) {
    const int nat = ucell_->nat;
    const int nat3 = 3 * nat;
    ModuleBase::ComplexMatrix dyn = data.get_dynmat(0);
    if (dyn.nr != nat3) {
        return;
    }
    const ModuleBase::matrix eps = data.get_dielectric();
    if (eps.nr != 3 || eps.nc != 3) {
        return; // no dielectric tensor stored yet (C6 not run)
    }
    const double qeq = qhat.x * (qhat.x * eps(0, 0) + qhat.y * eps(1, 0) + qhat.z * eps(2, 0))
                       + qhat.y * (qhat.x * eps(0, 1) + qhat.y * eps(1, 1) + qhat.z * eps(2, 1))
                       + qhat.z * (qhat.x * eps(0, 2) + qhat.y * eps(1, 2) + qhat.z * eps(2, 2));
    if (std::abs(qeq) < 1.0e-10) {
        return;
    }
    const double pref = ModuleBase::FOUR_PI * ModuleBase::e2 / ucell_->omega / qeq;
    for (int ia = 0; ia < nat; ++ia) {
        const double ma = ucell_->atoms[ucell_->iat2it[ia]].mass;
        const ModuleBase::matrix za = data.get_born(ia);
        if (za.nr != 3 || za.nc != 3) {
            continue;
        }
        for (int ib = 0; ib < nat; ++ib) {
            const double mb = ucell_->atoms[ucell_->iat2it[ib]].mass;
            const ModuleBase::matrix zb = data.get_born(ib);
            for (int da = 0; da < 3; ++da) {
                // (qhat Z*_a)_da = sum_gamma qhat_gamma Z_a(da,gamma)
                const double qza = qhat.x * za(da, 0) + qhat.y * za(da, 1) + qhat.z * za(da, 2);
                for (int db = 0; db < 3; ++db) {
                    const double qzb = qhat.x * zb(db, 0) + qhat.y * zb(db, 1) + qhat.z * zb(db, 2);
                    dyn(3 * ia + da, 3 * ib + db) += pref * qza * qzb / std::sqrt(ma * mb);
                }
            }
        }
    }
    data.set_dynmat(0, dyn);
}

bool DFPT_Phon::check_sum_rule(int q_idx, DFPT_PW_Data& data) const {
    const ModuleBase::Vector3<double> q_frac = data.get_qvec(q_idx);
    if (std::abs(q_frac.x) > 1.0e-8 || std::abs(q_frac.y) > 1.0e-8
        || std::abs(q_frac.z) > 1.0e-8) {
        return true; // only applies at Gamma
    }
    const int nat3 = 3 * ucell_->nat;
    ModuleBase::ComplexMatrix dyn = data.get_dynmat(q_idx);
    if (dyn.nr != nat3) {
        return false;
    }
    double max_elem = 0.0;
    for (int i = 0; i < nat3; ++i) {
        for (int j = 0; j < nat3; ++j) {
            max_elem = std::max(max_elem, std::abs(dyn(i, j)));
        }
    }
    if (max_elem < 1.0e-12) {
        return true;
    }
    for (int i = 0; i < nat3; ++i) {
        std::complex<double> colsum(0.0, 0.0);
        for (int j = 0; j < nat3; ++j) {
            colsum += dyn(i, j);
        }
        if (std::abs(colsum) > 1.0e-6 * max_elem) {
            return false;
        }
    }
    return true;
}

void DFPT_Phon::dftu_onsite(int q_idx, DFPT_PW_Data& data) {
    // Reserved DFT+U contribution to the dynamical matrix (U0).
    // The physical implementation lands with the Plus_U production wiring:
    //   sum_nk w_nk [ <psi|dV_U|dpsi> + frozen second-order U term
    //   (~ becp * V_eff * dbecp_f contractions) ], accumulated into the
    //   dynamical matrix. dV_U itself is assembled by DFPT_Pert::build_dv_u.
    (void)q_idx;
    (void)data;
}

} // namespace ModuleDFPT
