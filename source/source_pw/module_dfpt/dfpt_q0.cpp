// ============================================================
// This code is added by Mohan Chen on 2026-05-18.
// This code is currently in design phase and has not been
// put into production yet. It may change in the future.
// Please use this code with caution. Only developers who know
// what they are doing should use this code.
// ============================================================

#include "dfpt_q0.h"

#include "dfpt_pert.h"
#include "source_base/constants.h"
#include "source_base/global_function.h"

#include <cmath>
#include <complex>
#include <vector>

namespace ModuleDFPT {

DFPT_Q0::DFPT_Q0() {}

DFPT_Q0::~DFPT_Q0() {}

void DFPT_Q0::init(UnitCell& ucell, ModulePW::PW_Basis* pw_rho,
                   ModulePW::PW_Basis_K* pw_wfc, DFPT_Pert* pert) {
    ucell_ = &ucell;
    pw_rho_ = pw_rho;
    pw_wfc_ = pw_wfc;
    pert_ = pert;
}

void DFPT_Q0::pos_matrix(const psi::Psi<std::complex<double>>& psi,
                         const ModuleBase::matrix& eig,
                         std::vector<std::vector<std::vector<ModuleBase::Vector3<std::complex<double>>>>>& r_mat) {
    const int nk = psi.get_nk();
    const int nbands = psi.get_nbands();
    r_mat.assign(nk,
                 std::vector<std::vector<ModuleBase::Vector3<std::complex<double>>>>(
                     nbands,
                     std::vector<ModuleBase::Vector3<std::complex<double>>>(
                         nbands, ModuleBase::Vector3<std::complex<double>>(0.0, 0.0, 0.0))));
    if (pw_wfc_ == nullptr || ucell_ == nullptr || pert_ == nullptr) {
        return;
    }
    const double tpiba = ucell_->tpiba;
    const double tpiba2 = tpiba * tpiba;
    for (int ik = 0; ik < nk; ++ik) {
        const int npwk = pw_wfc_->npwk[ik];
        std::vector<ModuleBase::Vector3<double>> gk(npwk);
        for (int ig = 0; ig < npwk; ++ig) {
            gk[ig] = pw_wfc_->getgpluskcar(ik, ig);
        }
        // velocity operator dH/dk matrix elements, with the k derivative in
        // the same dimensionless 2*pi/lat0 units build_vkb_dk uses:
        //   p^d_{mn} = <u_m| 2 tpiba^2 (k+G)_d + dV_nl/dk_d |u_n>
        // V_loc is k-independent; the DFT+U commutator is the U0 reservation.
        std::vector<std::vector<ModuleBase::Vector3<std::complex<double>>>> p_mat(
            nbands,
            std::vector<ModuleBase::Vector3<std::complex<double>>>(
                nbands, ModuleBase::Vector3<std::complex<double>>(0.0, 0.0, 0.0)));
        // diagonal kinetic part: T = tpiba^2 |k+G|^2 (Ry a.u.)
        for (int m = 0; m < nbands; ++m) {
            for (int n = 0; n < nbands; ++n) {
                std::complex<double> dot[3] = {std::complex<double>(0.0, 0.0),
                                               std::complex<double>(0.0, 0.0),
                                               std::complex<double>(0.0, 0.0)};
                for (int ig = 0; ig < npwk; ++ig) {
                    const std::complex<double> cc =
                        std::conj(psi(ik, m, ig)) * psi(ik, n, ig);
                    for (int d = 0; d < 3; ++d) {
                        dot[d] += 2.0 * tpiba2 * gk[ig][d] * cc;
                    }
                }
                for (int d = 0; d < 3; ++d) {
                    p_mat[m][n][d] = dot[d];
                }
            }
        }
        // nonlocal derivative part: dV_nl/dk_d = sum_{mu,nu} (|dvkb_mu> D_{mu,nu} <vkb_nu|
        //                                              + |vkb_mu> D_{mu,nu} <dvkb_nu|)
        // with D_{mu,nu} = dion(ib_mu, ib_nu) delta_{m_mu, m_nu} (dVnl_dtau layout).
        for (int it = 0; it < ucell_->ntype; ++it) {
            const pseudo& ncpp = ucell_->atoms[it].ncpp;
            const int nh = ncpp.nh;
            if (nh == 0) {
                continue;
            }
            if (ncpp.tvanp || ncpp.has_so) {
                ModuleBase::WARNING_QUIT("DFPT_Q0::pos_matrix",
                                         "DFPT velocity operator is implemented for "
                                         "normal-conserving separable pseudopotentials only.");
            }
            // projector -> (radial beta index, m channel) table, matching build_vkb
            std::vector<int> mu_ib(nh, 0);
            std::vector<int> mu_m(nh, 0);
            int mu_idx = 0;
            for (int ib = 0; ib < ncpp.nbeta; ++ib) {
                const int l = ncpp.lll[ib];
                for (int m = 0; m < 2 * l + 1; ++m) {
                    if (mu_idx < nh) {
                        mu_ib[mu_idx] = ib;
                        mu_m[mu_idx] = m;
                    }
                    ++mu_idx;
                }
            }
            for (int ia = 0; ia < ucell_->atoms[it].na; ++ia) {
                std::vector<std::vector<std::complex<double>>> vkb;
                pert_->build_vkb(it, ia, gk, vkb);
                // becp_b[mu] = <vkb_mu|u_b> for all bands
                std::vector<std::vector<std::complex<double>>> becp(nbands);
                for (int b = 0; b < nbands; ++b) {
                    becp[b].assign(nh, std::complex<double>(0.0, 0.0));
                    for (int mu = 0; mu < nh; ++mu) {
                        for (int ig = 0; ig < npwk; ++ig) {
                            becp[b][mu] += std::conj(vkb[mu][ig]) * psi(ik, b, ig);
                        }
                    }
                }
                for (int d = 0; d < 3; ++d) {
                    std::vector<std::vector<std::complex<double>>> dvkb;
                    pert_->build_vkb_dk(it, ia, d, gk, vkb, dvkb);
                    // dbecp_b[mu] = <dvkb_mu|u_b>
                    std::vector<std::vector<std::complex<double>>> dbecp(nbands);
                    for (int b = 0; b < nbands; ++b) {
                        dbecp[b].assign(nh, std::complex<double>(0.0, 0.0));
                        for (int mu = 0; mu < nh; ++mu) {
                            for (int ig = 0; ig < npwk; ++ig) {
                                dbecp[b][mu] += std::conj(dvkb[mu][ig]) * psi(ik, b, ig);
                            }
                        }
                    }
                    // accumulate the two Hermitian-conjugate projector terms
                    for (int m = 0; m < nbands; ++m) {
                        for (int n = 0; n < nbands; ++n) {
                            std::complex<double> term(0.0, 0.0);
                            for (int mu = 0; mu < nh; ++mu) {
                                std::complex<double> out_m(0.0, 0.0);
                                std::complex<double> in_n(0.0, 0.0);
                                for (int nu = 0; nu < nh; ++nu) {
                                    if (mu_m[mu] != mu_m[nu]) {
                                        continue;
                                    }
                                    const double dij = ncpp.dion(mu_ib[mu], mu_ib[nu]);
                                    out_m += dij * becp[n][nu];
                                    in_n += dij * dbecp[n][nu];
                                }
                                // <u_m|dvkb_mu> D <vkb|u_n> + <u_m|vkb_mu> D <dvkb|u_n>
                                term += std::conj(dbecp[m][mu]) * out_m
                                        + std::conj(becp[m][mu]) * in_n;
                            }
                            p_mat[m][n][d] += term;
                        }
                    }
                }
            }
        }
        // velocity -> position: r = -i v / (tpiba (eps_m - eps_n)), r in bohr
        // (from [H, r] = -i dH/dk in Ry a.u.); degenerate pairs are skipped,
        // their gauge-dependent matrix elements carry no unique value.
        for (int m = 0; m < nbands; ++m) {
            for (int n = 0; n < nbands; ++n) {
                if (m == n) {
                    continue;
                }
                const double de = eig(ik, m) - eig(ik, n);
                if (std::abs(de) < 1.0e-8) {
                    continue;
                }
                for (int d = 0; d < 3; ++d) {
                    r_mat[ik][m][n][d] = std::complex<double>(0.0, -1.0) * p_mat[m][n][d]
                                          / (tpiba * de);
                }
            }
        }
    }
}

void DFPT_Q0::compute_eps(const psi::Psi<std::complex<double>>& psi,
                          const ModuleBase::matrix& wg,
                          const ModuleBase::matrix& eig, DFPT_PW_Data& data) {
    if (ucell_ == nullptr) {
        return;
    }
    std::vector<std::vector<std::vector<ModuleBase::Vector3<std::complex<double>>>>> r_mat;
    pos_matrix(psi, eig, r_mat);
    const int nk = psi.get_nk();
    const int nbands = psi.get_nbands();
    ModuleBase::matrix eps(3, 3, true);
    for (int a = 0; a < 3; ++a) {
        for (int b = 0; b < 3; ++b) {
            double chi = 0.0;
            for (int ik = 0; ik < nk; ++ik) {
                for (int v = 0; v < nbands; ++v) {
                    if (!dfpt_band_occupied(wg, ik, v)) {
                        continue; // empty
                    }
                    for (int c = 0; c < nbands; ++c) {
                        if (dfpt_band_occupied(wg, ik, c)) {
                            continue; // occupied
                        }
                        const double de = eig(ik, c) - eig(ik, v);
                        if (std::abs(de) < 1.0e-8) {
                            continue;
                        }
                        chi += wg(ik, v)
                               * (r_mat[ik][v][c][a] * r_mat[ik][c][v][b]).real() / de;
                    }
                }
            }
            eps(a, b) = ((a == b) ? 1.0 : 0.0)
                        + 8.0 * ModuleBase::PI / ucell_->omega * chi / nk;
        }
    }
    data.set_dielectric(eps);
}

void DFPT_Q0::compute_born(const psi::Psi<std::complex<double>>& psi,
                           const ModuleBase::matrix& wg,
                           const ModuleBase::matrix& eig, DFPT_PW_Data& data) {
    if (ucell_ == nullptr || pert_ == nullptr) {
        return;
    }
    std::vector<std::vector<std::vector<ModuleBase::Vector3<std::complex<double>>>>> r_mat;
    pos_matrix(psi, eig, r_mat);
    const int nk = psi.get_nk();
    const int nbands = psi.get_nbands();
    const int nat = ucell_->nat;

    // stash the q=0 dpsi slots (apply_dv reuses them, phon backup pattern)
    std::vector<std::vector<std::vector<std::complex<double>>>> dpsib(nk);
    for (int ik = 0; ik < nk; ++ik) {
        dpsib[ik].resize(nbands);
        for (int ib = 0; ib < nbands; ++ib) {
            dpsib[ik][ib] = data.get_dpsi(0, ik, ib);
        }
    }

    for (int iat = 0; iat < nat; ++iat) {
        ModuleBase::matrix zstar(3, 3, true);
        for (int idir = 0; idir < 3; ++idir) {
            // dV matrix elements at q = 0 through the C1 path; apply_dv
            // delivers dV|u_v> on the k+q = k basis for every k.
            pert_->build_dv(0, iat, idir, data);
            std::vector<double> acc(3, 0.0);
            for (int ik = 0; ik < nk; ++ik) {
                pert_->apply_dv(0, ik, psi, data);
                for (int v = 0; v < nbands; ++v) {
                    if (!dfpt_band_occupied(wg, ik, v)) {
                        continue; // empty
                    }
                    const std::vector<std::complex<double>> rhs =
                        data.get_dpsi(0, ik, v);
                    if (rhs.empty()) {
                        continue;
                    }
                    for (int m = 0; m < nbands; ++m) {
                        const double de = eig(ik, m) - eig(ik, v);
                        if (std::abs(de) < 1.0e-8) {
                            continue; // m == v or degenerate partner
                        }
                        std::complex<double> dv_mv(0.0, 0.0);
                        for (size_t ig = 0; ig < rhs.size(); ++ig) {
                            dv_mv += std::conj(psi(ik, m, ig)) * rhs[ig];
                        }
                        // <u_v|dV|u_m> = conj(dv_mv), multiplied from the
                        // right by <u_m|r_a|u_v> (Gonze-Lee ordering)
                        for (int a = 0; a < 3; ++a) {
                            acc[a] += wg(ik, v)
                                      * (std::conj(dv_mv) * r_mat[ik][m][v][a]).real() / de;
                        }
                    }
                }
            }
            for (int a = 0; a < 3; ++a) {
                zstar(a, idir) = -4.0 / nk * acc[a];
            }
        }
        // ionic rigid-ion charge on the diagonal (a == b directions)
        const int it = ucell_->iat2it[iat];
        const double zion = ucell_->atoms[it].ncpp.zv;
        for (int d = 0; d < 3; ++d) {
            zstar(d, d) += zion;
        }
        data.set_born(iat, zstar);
    }

    // restore the stashed q=0 dpsi
    for (int ik = 0; ik < nk; ++ik) {
        for (int ib = 0; ib < nbands; ++ib) {
            if (!dpsib[ik][ib].empty()) {
                data.set_dpsi(0, ik, ib, dpsib[ik][ib]);
            }
        }
    }
}

void DFPT_Q0::compute_q0_response(DFPT_PW_Data& data) {
    // DFT+U reservation (U0): V_U is nonlocal (onsite projector), so the
    // position operator does NOT commute with the DFT+U potential. The
    // [r, V_U] commutator term must be handled separately in addition to
    // the occupation-matrix response (docc) when u_active() runs; this is
    // the hardest DFT+U piece and is deferred with the Plus_U wiring.
    (void)data;
}

} // namespace ModuleDFPT
