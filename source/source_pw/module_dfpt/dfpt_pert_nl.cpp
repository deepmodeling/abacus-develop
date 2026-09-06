// Nonlocal first- and second-order potentials of DFPT_Pert (normal-
// conserving separable case), split out of dfpt_pert.cpp: dVnl_dtau
// and apply_d2vnl with the shared projector table, D contraction and
// projector-sum helpers. All formulas are moved verbatim from the
// original bodies; the always-true include_middle knob of
// apply_d2vnl is dropped (its q-independence is established QE ground
// truth, see the comment inside apply_d2vnl).

#include "dfpt_pert.h"

#include "source_base/constants.h"
#include "source_base/tool_quit.h"
#include "source_cell/atom_pseudo.h"

#include <cmath>
#include <complex>
#include <vector>

#include "source_base/timer.h"
#include "source_base/tool_title.h"

namespace ModuleDFPT
{

namespace
{

/// projector -> (radial beta index, m channel) table matching build_vkb
void nl_projector_table(const pseudo& ncpp, int nh, std::vector<int>& mu_ib, std::vector<int>& mu_m)
{
    mu_ib.assign(nh, 0);
    mu_m.assign(nh, 0);
    int mu_idx = 0;
    for (int ib = 0; ib < ncpp.nbeta; ++ib)
    {
        const int l = ncpp.lll[ib];
        for (int m = 0; m < 2 * l + 1; ++m)
        {
            if (mu_idx < nh)
            {
                mu_ib[mu_idx] = ib;
                mu_m[mu_idx] = m;
            }
            ++mu_idx;
        }
    }
}

/// cartesian list of the k (+q) plane waves of one k point (2*pi/lat0 units)
std::vector<ModuleBase::Vector3<double>> nl_gk_list(const ModulePW::PW_Basis_K& pw_wfc, int k_idx, int npwk)
{
    std::vector<ModuleBase::Vector3<double>> gk(npwk);
    for (int ig = 0; ig < npwk; ++ig)
    {
        gk[ig] = pw_wfc.getgpluskcar(k_idx, ig);
    }
    return gk;
}

/// becp_nu(k) = sum_G' conj(vkb_in[nu][G']) psi(G') for one band
void nl_becp(int npwk,
             const std::vector<std::vector<std::complex<double>>>& vkb_in,
             const std::complex<double>* psi_in,
             std::vector<std::complex<double>>& becp)
{
    const int nh = static_cast<int>(vkb_in.size());
    becp.assign(nh, std::complex<double>(0.0, 0.0));
    for (int nu = 0; nu < nh; ++nu)
    {
        for (int ig = 0; ig < npwk; ++ig)
        {
            becp[nu] += std::conj(vkb_in[nu][ig]) * psi_in[ig];
        }
    }
}

/// becp of the momentum-weighted band i (k+G')_dir |psi> (term B carrier)
void nl_becp_dpsi(int npwk,
                  int dir,
                  double tpiba,
                  const std::vector<ModuleBase::Vector3<double>>& gk_in,
                  const std::vector<std::vector<std::complex<double>>>& vkb_in,
                  const std::complex<double>* psi_in,
                  std::vector<std::complex<double>>& becp_dpsi)
{
    const int nh = static_cast<int>(vkb_in.size());
    becp_dpsi.assign(nh, std::complex<double>(0.0, 0.0));
    for (int nu = 0; nu < nh; ++nu)
    {
        for (int ig = 0; ig < npwk; ++ig)
        {
            const std::complex<double> dpsi_ig
                = std::complex<double>(0.0, 1.0) * (tpiba * gk_in[ig][dir]) * psi_in[ig];
            becp_dpsi[nu] += std::conj(vkb_in[nu][ig]) * dpsi_ig;
        }
    }
}

/// dcbecp = D * becp with D_{mu,nu} = dion(ib_mu, ib_nu) delta_{m_mu, m_nu}
void nl_d_contract(const pseudo& ncpp,
                   const std::vector<int>& mu_ib,
                   const std::vector<int>& mu_m,
                   const std::vector<std::complex<double>>& becp,
                   std::vector<std::complex<double>>& dcbecp)
{
    const int nh = static_cast<int>(becp.size());
    dcbecp.assign(nh, std::complex<double>(0.0, 0.0));
    for (int mu = 0; mu < nh; ++mu)
    {
        for (int nu = 0; nu < nh; ++nu)
        {
            if (mu_m[mu] != mu_m[nu])
            {
                continue;
            }
            dcbecp[mu] += ncpp.dion(mu_ib[mu], mu_ib[nu]) * becp[nu];
        }
    }
}

/// D contraction of the four becp variants with the same-m selection rule
/// shared with dVnl_dtau; dx rows are 0 = d0 (plain), 1 = da_ (k+G')_da,
/// 2 = db_ (k+G')_db, 3 = dab (k+G')_da (k+G')_db
void nl_d_contract_x(const pseudo& ncpp,
                     const std::vector<int>& mu_ib,
                     const std::vector<int>& mu_m,
                     const std::vector<std::vector<std::complex<double>>>& becp_x,
                     std::vector<std::vector<std::complex<double>>>& dx)
{
    const int nh = static_cast<int>(becp_x[0].size());
    dx.assign(4, std::vector<std::complex<double>>(nh, std::complex<double>(0.0, 0.0)));
    for (int mu = 0; mu < nh; ++mu)
    {
        for (int nu = 0; nu < nh; ++nu)
        {
            if (mu_m[mu] != mu_m[nu])
            {
                continue;
            }
            const double dij = ncpp.dion(mu_ib[mu], mu_ib[nu]);
            for (int iw = 0; iw < 4; ++iw)
            {
                dx[iw][mu] += dij * becp_x[iw][nu];
            }
        }
    }
}

/// (Vnl |carrier>)_igl = sum_mu vkb_out[mu][igl] coeff[mu] on the k+q basis
std::complex<double> nl_sum_projectors(const std::vector<std::vector<std::complex<double>>>& vkb_out,
                                       const std::vector<std::complex<double>>& coeff,
                                       int nh,
                                       int igl)
{
    std::complex<double> vnlpsi(0.0, 0.0);
    for (int mu = 0; mu < nh; ++mu)
    {
        vnlpsi += vkb_out[mu][igl] * coeff[mu];
    }
    return vnlpsi;
}

} // namespace

void DFPT_Pert::dVnl_dtau(int atom_idx,
                          int dir,
                          const ModuleBase::Vector3<double>& q_cart,
                          const psi::Psi<std::complex<double>>& psi,
                          int k_idx,
                          std::vector<std::vector<std::complex<double>>>& dv_psi)
{
    ModuleBase::TITLE("DFPT_Pert", "dVnl_dtau");
    ModuleBase::timer::start("DFPT_Pert", "dVnl_dtau");
    int it = 0;
    int ia = 0;
    atom_index(atom_idx, it, ia);
    if (ia < 0)
    {
        ModuleBase::timer::end("DFPT_Pert", "dVnl_dtau");
        return;
    }
    const pseudo& ncpp = ucell_->atoms[it].ncpp;
    if (ncpp.tvanp || ncpp.has_so)
    {
        // the separable NC path documented in C1; ultrasoft and spin-orbit
        // projectors are deferred (their D and augmentation have |k+q| shifts
        // that need the USPP machinery).
        ModuleBase::WARNING_QUIT("DFPT_Pert::dVnl_dtau",
                                 "DFPT nonlocal first-order potential is implemented "
                                 "for normal-conserving separable pseudopotentials only.");
    }
    const int nh = ncpp.nh;

    // projector -> (radial beta index, m channel) table, matching build_vkb.
    std::vector<int> mu_ib;
    std::vector<int> mu_m;
    nl_projector_table(ncpp, nh, mu_ib, mu_m);

    // incoming k basis: G = k + G' (pw_wfc k-basis index)
    const int npwk = pw_wfc_->npwk[k_idx];
    const std::vector<ModuleBase::Vector3<double>> gk_in = nl_gk_list(*pw_wfc_, k_idx, npwk);
    std::vector<std::vector<std::complex<double>>> vkb_in;
    build_vkb(it, ia, gk_in, vkb_in);

    // outgoing k+q basis
    DFPT_KQ_Basis kq;
    kq.init(pw_wfc_, pw_rho_, q_cart, k_idx);
    const int npwk_kq = kq.get_npwk();
    std::vector<ModuleBase::Vector3<double>> gk_out(npwk_kq);
    for (int igl = 0; igl < npwk_kq; ++igl)
    {
        gk_out[igl] = kq.get_gpluskq(igl);
    }
    std::vector<std::vector<std::complex<double>>> vkb_out;
    build_vkb(it, ia, gk_out, vkb_out);

    const int nbands = psi.get_nbands();
    dv_psi.assign(nbands, std::vector<std::complex<double>>(npwk_kq, std::complex<double>(0.0, 0.0)));

    for (int iband = 0; iband < nbands; ++iband)
    {
        // becp_nu(k) = sum_G' conj(vkb_in[nu][G']) psi(G')
        std::vector<std::complex<double>> becp;
        nl_becp(npwk, vkb_in, &psi(k_idx, iband, 0), becp);
        // dcbecp = D * becp with D_{mu,nu} = dion(ib_mu, ib_nu) delta_{m_mu, m_nu}
        std::vector<std::complex<double>> dcbecp;
        nl_d_contract(ncpp, mu_ib, mu_m, becp, dcbecp);
        // term A: i (k+q+G'')_dir * (Vnl |psi>) on the k+q basis
        std::vector<std::complex<double>> term_a(npwk_kq, std::complex<double>(0.0, 0.0));
        for (int igl = 0; igl < npwk_kq; ++igl)
        {
            term_a[igl] = std::complex<double>(0.0, 1.0) * (ucell_->tpiba * gk_out[igl][dir])
                          * nl_sum_projectors(vkb_out, dcbecp, nh, igl);
        }
        // term B: Vnl [i (k+G')_dir |psi>]
        std::vector<std::complex<double>> becp_dpsi;
        nl_becp_dpsi(npwk, dir, ucell_->tpiba, gk_in, vkb_in, &psi(k_idx, iband, 0), becp_dpsi);
        std::vector<std::complex<double>> dcbecp_dpsi;
        nl_d_contract(ncpp, mu_ib, mu_m, becp_dpsi, dcbecp_dpsi);
        for (int igl = 0; igl < npwk_kq; ++igl)
        {
            // GS exp(-2pi gk.tau) projector convention: dVnl/dtau_dir
            // |psi> = -i (k+q+G'')_dir (Vnl|psi>) + Vnl[i (k+G')_dir |psi>]
            dv_psi[iband][igl] = nl_sum_projectors(vkb_out, dcbecp_dpsi, nh, igl) - term_a[igl];
        }
    }
    ModuleBase::timer::end("DFPT_Pert", "dVnl_dtau");
}

void DFPT_Pert::apply_d2vnl(int atom_idx,
                            int da,
                            int db,
                            const ModuleBase::Vector3<double>& q_eff,
                            const psi::Psi<std::complex<double>>& psi,
                            int k_idx,
                            std::vector<std::vector<std::complex<double>>>& d2v_psi) const
{
    ModuleBase::TITLE("DFPT_Pert", "apply_d2vnl");
    ModuleBase::timer::start("DFPT_Pert", "apply_d2vnl");
    int it = 0;
    int ia = 0;
    atom_index(atom_idx, it, ia);
    if (ia < 0)
    {
        ModuleBase::timer::end("DFPT_Pert", "apply_d2vnl");
        return;
    }
    const pseudo& ncpp = ucell_->atoms[it].ncpp;
    if (ncpp.tvanp || ncpp.has_so)
    {
        ModuleBase::WARNING_QUIT("DFPT_Pert::apply_d2vnl",
                                 "DFPT second-order nonlocal potential is implemented "
                                 "for normal-conserving separable pseudopotentials only.");
    }
    const int nh = ncpp.nh;
    const int nbands = psi.get_nbands();

    // projector -> (radial index, m channel) table, matching build_vkb
    std::vector<int> mu_ib;
    std::vector<int> mu_m;
    nl_projector_table(ncpp, nh, mu_ib, mu_m);

    // incoming k basis and outgoing k+q basis projectors (same atom)
    const int npwk = pw_wfc_->npwk[k_idx];
    const std::vector<ModuleBase::Vector3<double>> gk_in = nl_gk_list(*pw_wfc_, k_idx, npwk);
    std::vector<std::vector<std::complex<double>>> vkb_in;
    build_vkb(it, ia, gk_in, vkb_in);
    DFPT_KQ_Basis kq;
    kq.init(pw_wfc_, pw_rho_, q_eff, k_idx);
    const int npwk_kq = kq.get_npwk();
    std::vector<ModuleBase::Vector3<double>> gk_out(npwk_kq);
    for (int igl = 0; igl < npwk_kq; ++igl)
    {
        gk_out[igl] = kq.get_gpluskq(igl);
    }
    std::vector<std::vector<std::complex<double>>> vkb_out;
    build_vkb(it, ia, gk_out, vkb_out);

    d2v_psi.assign(nbands, std::vector<std::complex<double>>(npwk_kq, std::complex<double>(0.0, 0.0)));
    for (int iband = 0; iband < nbands; ++iband)
    {
        // becp and its (k+G')-weighted variants: becp_x = sum x(G') |beta><psi|
        // rows 0 = plain, 1 = (k+G')_da, 2 = (k+G')_db, 3 = (k+G')_da (k+G')_db
        std::vector<std::vector<std::complex<double>>> becp_x(
            4, std::vector<std::complex<double>>(nh, std::complex<double>(0.0, 0.0)));
        for (int nu = 0; nu < nh; ++nu)
        {
            for (int ig = 0; ig < npwk; ++ig)
            {
                const std::complex<double> vc = std::conj(vkb_in[nu][ig]) * psi(k_idx, iband, ig);
                const double kp_da = ucell_->tpiba * gk_in[ig][da];
                const double kp_db = ucell_->tpiba * gk_in[ig][db];
                becp_x[0][nu] += vc;
                becp_x[1][nu] += kp_da * vc;
                becp_x[2][nu] += kp_db * vc;
                becp_x[3][nu] += kp_da * kp_db * vc;
            }
        }
        // D contraction with the same-m selection rule as dVnl_dtau
        std::vector<std::vector<std::complex<double>>> dx;
        nl_d_contract_x(ncpp, mu_ib, mu_m, becp_x, dx);
        // chi(G'') = sum_mu vkb_out,mu [ -kq_da kq_db d0 - dab
        //                              + kq_da db_ + kq_db da_ ]_mu
        // QE ground truth (dynmat_us.f90 + phq_init.f90): the KB second-order
        // term pairs gammap (integer-G (k+G)_da(k+G)_db derivative of beta)
        // with becp1 = <beta_k|psi_k> and the same-atom alphap_a* alphap_b
        // middle product; everything is built at k with integer-G momentum
        // factors, so the caller passes q_eff = 0 and the kernel is
        // q-independent for every q.
        for (int igl = 0; igl < npwk_kq; ++igl)
        {
            const double kq_da = ucell_->tpiba * gk_out[igl][da];
            const double kq_db = ucell_->tpiba * gk_out[igl][db];
            std::complex<double> chi(0.0, 0.0);
            for (int mu = 0; mu < nh; ++mu)
            {
                chi += vkb_out[mu][igl] * (-kq_da * kq_db * dx[0][mu] - dx[3][mu]);
                chi += vkb_out[mu][igl] * (kq_da * dx[2][mu] + kq_db * dx[1][mu]);
            }
            d2v_psi[iband][igl] = chi;
        }
    }
    ModuleBase::timer::end("DFPT_Pert", "apply_d2vnl");
}

} // namespace ModuleDFPT
