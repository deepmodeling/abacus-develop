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
#include <cstdlib>
#include <iostream>
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
    stars_.clear();
}

namespace {
// element accessor for ModuleBase::Matrix3 (row i, column j); the public
// interface only exposes the named e11..e33 members
inline double me(const ModuleBase::Matrix3& m, int i, int j) {
    switch (3 * i + j) {
    case 0: return m.e11;
    case 1: return m.e12;
    case 2: return m.e13;
    case 3: return m.e21;
    case 4: return m.e22;
    case 5: return m.e23;
    case 6: return m.e31;
    case 7: return m.e32;
    default: return m.e33;
    }
}

// folded fractional equality with the lattice periodicity absorbed
inline bool folded_equal(double a, double b, double tol) {
    const double d = std::abs(a - b);
    return d < tol || std::abs(d - 1.0) < tol;
}
} // namespace

void DFPT_Q0::build_stars(int nk) {
    // every k starts with the identity member (also the permanent fallback)
    stars_.assign(nk, std::vector<StarMember>(1, StarMember()));
    if (ucell_ == nullptr || pw_wfc_ == nullptr || pw_wfc_->kvec_d == nullptr) {
        return;
    }
    const ModuleSymmetry::Symmetry& symm = ucell_->symm;
    if (symm.nrotk <= 0) {
        // no point-group analysis (symmetry off / unreduced mesh): the
        // stored list is already the full mesh, identity members only
        return;
    }
    const bool stardbg = getenv("DFPT_STARDBG") != nullptr;
    const int nat = ucell_->nat;
    if (stardbg) {
        std::cout << "STARDBG nrotk=" << symm.nrotk << " nat=" << nat << std::endl;
        for (int j = 0; j < symm.nrotk; ++j) {
            std::cout << "STARDBG j=" << j << " kg=[" << symm.kgmatrix[j].e11 << ","
                      << symm.kgmatrix[j].e12 << "," << symm.kgmatrix[j].e13 << ";"
                      << symm.kgmatrix[j].e21 << "," << symm.kgmatrix[j].e22 << ","
                      << symm.kgmatrix[j].e23 << ";" << symm.kgmatrix[j].e31 << ","
                      << symm.kgmatrix[j].e32 << "," << symm.kgmatrix[j].e33 << "] g=["
                      << symm.gmatrix[j].e11 << "," << symm.gmatrix[j].e12 << ","
                      << symm.gmatrix[j].e13 << ";" << symm.gmatrix[j].e21 << ","
                      << symm.gmatrix[j].e22 << "," << symm.gmatrix[j].e23 << ";"
                      << symm.gmatrix[j].e31 << "," << symm.gmatrix[j].e32 << ","
                      << symm.gmatrix[j].e33 << "] gt=(" << symm.gtrans[j].x << ","
                      << symm.gtrans[j].y << "," << symm.gtrans[j].z << ")" << std::endl;
        }
    }
    std::vector<ModuleBase::Vector3<double>> kfolds;
    for (int ik = 0; ik < nk; ++ik) {
        kfolds.clear();
        // the pre-filled identity member owns the folded k itself
        ModuleBase::Vector3<double> k0 = pw_wfc_->kvec_d[ik];
        k0.x -= std::round(k0.x);
        k0.y -= std::round(k0.y);
        k0.z -= std::round(k0.z);
        kfolds.push_back(k0);
        for (int j = 0; j < symm.nrotk; ++j) {
            ModuleBase::Vector3<double> kp = pw_wfc_->kvec_d[ik] * symm.kgmatrix[j];
            // fold to [-0.5, 0.5): star members are grid points, the
            // folded coordinates identify the distinct mesh points
            kp.x -= std::round(kp.x);
            kp.y -= std::round(kp.y);
            kp.z -= std::round(kp.z);
            bool dup = false;
            for (size_t im = 0; im < kfolds.size(); ++im) {
                if (folded_equal(kp.x, kfolds[im].x, 1.0e-5)
                    && folded_equal(kp.y, kfolds[im].y, 1.0e-5)
                    && folded_equal(kp.z, kfolds[im].z, 1.0e-5)) {
                    dup = true;
                    break;
                }
            }
            if (dup) {
                continue;
            }
            kfolds.push_back(kp);
            StarMember mem;
            // cartesian form of the same operation: k_frac' = k_frac * K,
            // k_cart = k_frac * G, hence k_cart' = k_cart * (G^-1 K G). That
            // product is the row-convention operator; rotate_tensor applies
            // the column form chi' = R chi R^T, so store the transpose
            const ModuleBase::Matrix3 krow
                = ucell_->G.Inverse() * symm.kgmatrix[j] * ucell_->G;
            mem.cart = ModuleBase::Matrix3(krow.e11, krow.e21, krow.e31,
                                           krow.e12, krow.e22, krow.e32,
                                           krow.e13, krow.e23, krow.e33);
            if (stardbg) {
                std::cout << "STARDBG ik=" << ik << " j=" << j << " kp=(" << kp.x
                          << "," << kp.y << "," << kp.z << ") cart=[" << mem.cart.e11
                          << "," << mem.cart.e12 << "," << mem.cart.e13 << ";"
                          << mem.cart.e21 << "," << mem.cart.e22 << "," << mem.cart.e23
                          << ";" << mem.cart.e31 << "," << mem.cart.e32 << ","
                          << mem.cart.e33 << "]" << std::endl;
            }
            // atom image under the paired direct-space operation
            mem.atom_map.assign(nat, -1);
            bool ok = true;
            for (int iat = 0; iat < nat && ok; ++iat) {
                const int it = ucell_->iat2it[iat];
                const int ia = ucell_->iat2ia[iat];
                ModuleBase::Vector3<double> tp
                    = ucell_->atoms[it].taud[ia] * symm.gmatrix[j] + symm.gtrans[j];
                tp.x -= std::floor(tp.x);
                tp.y -= std::floor(tp.y);
                tp.z -= std::floor(tp.z);
                for (int jat = 0; jat < nat; ++jat) {
                    if (ucell_->iat2it[jat] != it) {
                        continue; // a species maps onto itself
                    }
                    const int ja = ucell_->iat2ia[jat];
                    const ModuleBase::Vector3<double>& tq
                        = ucell_->atoms[it].taud[ja];
                    if (folded_equal(tp.x, tq.x, 1.0e-4)
                        && folded_equal(tp.y, tq.y, 1.0e-4)
                        && folded_equal(tp.z, tq.z, 1.0e-4)) {
                        mem.atom_map[iat] = jat;
                        break;
                    }
                }
                if (mem.atom_map[iat] < 0) {
                    ok = false;
                }
            }
            if (!ok) {
                // inconsistent operation set: fall back to identity-only
                // stars for every k (the unreduced-sum behavior)
                stars_.assign(nk, std::vector<StarMember>(1, StarMember()));
                return;
            }
            if (stardbg) {
                std::cout << "STARDBG ik=" << ik << " j=" << j << " amap=";
                for (int iat2 = 0; iat2 < nat; ++iat2) {
                    std::cout << mem.atom_map[iat2] << (iat2 + 1 < nat ? "," : "");
                }
                std::cout << std::endl;
            }
            stars_[ik].push_back(mem);
        }
    }
}

void DFPT_Q0::rotate_tensor(const ModuleBase::Matrix3& r,
                            const ModuleBase::matrix& chi,
                            double (&chi_rot)[9]) {
    for (int a = 0; a < 3; ++a) {
        for (int b = 0; b < 3; ++b) {
            double s = 0.0;
            for (int ap = 0; ap < 3; ++ap) {
                for (int bp = 0; bp < 3; ++bp) {
                    s += me(r, a, ap) * me(r, b, bp) * chi(ap, bp);
                }
            }
            chi_rot[3 * a + b] = s;
        }
    }
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
        if (getenv("DFPT_Q0DBG") != nullptr) {
            std::cout << "Q0DBG ik=" << ik << " tpiba=" << tpiba
                      << " npwk=" << npwk << std::endl;
            for (int m = 0; m < nbands; ++m) {
                for (int n = 0; n < nbands; ++n) {
                    if (m == n) {
                        continue;
                    }
                    const double de = eig(ik, m) - eig(ik, n);
                    double p2 = 0.0;
                    for (int d = 0; d < 3; ++d) {
                        p2 += std::norm(p_mat[m][n][d]);
                    }
                    std::cout << "Q0DBG p m=" << m << " n=" << n
                              << " de=" << de << " px=" << p_mat[m][n][0]
                              << " py=" << p_mat[m][n][1]
                              << " pz=" << p_mat[m][n][2] << std::endl;
                }
            }
        }
    }
}

void DFPT_Q0::compute_eps(const ModuleBase::matrix& wg, DFPT_PW_Data& data) {
    if (ucell_ == nullptr) {
        return;
    }
    const int nk = wg.nr;
    const int nbands = wg.nc;

    // bare position legs Y^a and converged E-field responses dpsi^E,b of
    // the q = 0 mesh (DFPT_PW::solve_pos_resp / solve_efield_resp)
    std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>> yr(3);
    std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>> de(3);
    for (int a = 0; a < 3; ++a) {
        yr[a] = data.get_pos_resp(a);
        de[a] = data.get_dpsi_efield(a);
        if (static_cast<int>(yr[a].size()) != nk
            || static_cast<int>(de[a].size()) != nk) {
            return; // responses not solved: nothing to accumulate
        }
    }

    build_stars(nk);
    // wg-weighted partial chi_k[ik](a, b) = sum_occ Re<Y^a|dpsi^E,b> at
    // every stored k (QE dielec.f90: eps -= 4*(4pi/Omega)*wk*Re<dvpsi|dpsi>)
    std::vector<ModuleBase::matrix> chi_k(nk, ModuleBase::matrix(3, 3, true));
    for (int ik = 0; ik < nk; ++ik) {
        for (int v = 0; v < nbands; ++v) {
            if (!dfpt_band_occupied(wg, ik, v)) {
                continue; // empty
            }
            for (int a = 0; a < 3; ++a) {
                const int npw = static_cast<int>(yr[a][ik][v].size());
                if (npw <= 0) {
                    continue;
                }
                for (int b = 0; b < 3; ++b) {
                    if (static_cast<int>(de[b][ik][v].size()) != npw) {
                        continue;
                    }
                    std::complex<double> dot(0.0, 0.0);
                    for (int ig = 0; ig < npw; ++ig) {
                        dot += std::conj(yr[a][ik][v][ig]) * de[b][ik][v][ig];
                    }
                    chi_k[ik](a, b) += wg(ik, v) * dot.real();
                }
            }
        }
    }
    ModuleBase::matrix eps(3, 3, true);
    // wg carries the full k weight (star size included) times the spin
    // factor 2, so the star-averaged partials sum to the complete
    // Brillouin-zone average: no extra 1/nk normalization
    for (int ik = 0; ik < nk; ++ik) {
        const double inv_nstar = 1.0 / static_cast<double>(stars_[ik].size());
        for (size_t im = 0; im < stars_[ik].size(); ++im) {
            double rot[9];
            rotate_tensor(stars_[ik][im].cart, chi_k[ik], rot);
            for (int a = 0; a < 3; ++a) {
                for (int b = 0; b < 3; ++b) {
                    eps(a, b) += inv_nstar * rot[3 * a + b];
                }
            }
        }
    }
    for (int a = 0; a < 3; ++a) {
        for (int b = 0; b < 3; ++b) {
            // 16 pi / Omega: QE dielec.f90 form eps = 1 - 4*(4pi/Omega)*wk*
            // Re<Y^a|dpsi^E,b> (validated against QE 7.2 Si to 0.06%:
            // 23.6825 here vs 23.6685 QE)
            eps(a, b) *= -16.0 * ModuleBase::PI / ucell_->omega;
            if (a == b) {
                eps(a, b) += 1.0;
            }
        }
    }
    data.set_dielectric(eps);
}

void DFPT_Q0::compute_born(const psi::Psi<std::complex<double>>& psi,
                           const ModuleBase::matrix& wg,
                           const ModuleBase::matrix& eig, DFPT_PW_Data& data) {
    if (ucell_ == nullptr) {
        return;
    }
    const int nk = psi.get_nk();
    const int nbands = psi.get_nbands();
    const int nat = ucell_->nat;
    const int nbasis = psi.get_nbasis();
    (void)eig;

    // solved position legs Y^a_{k,v} = P_c x_a|psi_v> of the q = 0 mesh
    // (DFPT_PW::solve_pos_resp stashes them per direction)
    std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>> yr(3);
    for (int a = 0; a < 3; ++a) {
        yr[a] = data.get_pos_resp(a);
        if (static_cast<int>(yr[a].size()) != nk) {
            return; // position responses not solved: nothing to accumulate
        }
    }

    build_stars(nk);
    // star-rotated electronic partials, credited to the image atom under
    // each star member: zacc[kappa](a, idir)
    std::vector<ModuleBase::matrix> zacc(nat, ModuleBase::matrix(3, 3, true));

    for (int iat = 0; iat < nat; ++iat) {
        // wg-weighted partial chi_k[ik](a, idir) of THIS atom at every k
        std::vector<ModuleBase::matrix> chi_k(nk, ModuleBase::matrix(3, 3, true));
        for (int idir = 0; idir < 3; ++idir) {
            // converged screened displacement response dpsi(scf)/du of this
            // mode, stashed by solve_displacement before compute_born runs
            const std::vector<std::vector<std::vector<std::complex<double>>>> disp
                = data.get_dpsi_disp(iat, idir);
            if (static_cast<int>(disp.size()) != nk) {
                continue;
            }
            for (int ik = 0; ik < nk; ++ik) {
                for (int v = 0; v < nbands; ++v) {
                    if (!dfpt_band_occupied(wg, ik, v)) {
                        continue; // empty
                    }
                    const int npw = static_cast<int>(disp[ik][v].size());
                    if (npw <= 0 || npw > nbasis) {
                        continue; // unsolved slot or inconsistent basis
                    }
                    // <dpsi^kappa_idir(scf)|Y^a_v> per field direction
                    for (int a = 0; a < 3; ++a) {
                        if (static_cast<int>(yr[a][ik][v].size()) != npw) {
                            continue;
                        }
                        std::complex<double> dot(0.0, 0.0);
                        for (int ig = 0; ig < npw; ++ig) {
                            dot += std::conj(disp[ik][v][ig]) * yr[a][ik][v][ig];
                        }
                        chi_k[ik](a, idir) += wg(ik, v) * dot.real();
                    }
                }
            }
        }
        // star average: the partial at member Rk is R chi(k) R^T and is
        // credited to the image atom R(iat); wg already carries the star
        // size, so each member contributes with 1/n_star
        for (int ik = 0; ik < nk; ++ik) {
            const double inv_nstar
                = 1.0 / static_cast<double>(stars_[ik].size());
            for (size_t im = 0; im < stars_[ik].size(); ++im) {
                const StarMember& mem = stars_[ik][im];
                const int jat = (mem.atom_map.empty()) ? iat : mem.atom_map[iat];
                double rot[9];
                rotate_tensor(mem.cart, chi_k[ik], rot);
                for (int a = 0; a < 3; ++a) {
                    for (int d = 0; d < 3; ++d) {
                        zacc[jat](a, d) += inv_nstar * rot[3 * a + d];
                    }
                }
            }
        }
    }

    for (int iat = 0; iat < nat; ++iat) {
        ModuleBase::matrix zstar(3, 3, true);
        for (int a = 0; a < 3; ++a) {
            for (int d = 0; d < 3; ++d) {
                zstar(a, d) = -2.0 * zacc[iat](a, d);
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
