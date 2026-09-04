#include "dfpt_q0.h"

#include "dfpt_pert.h"
#include "source_base/global_function.h"
#include "source_base/timer.h"
#include "source_base/tool_title.h"

#include <cmath>
#include <complex>
#include <vector>

namespace ModuleDFPT
{

namespace
{

/// [m][n] velocity matrix elements of one k point
typedef std::vector<std::vector<ModuleBase::Vector3<std::complex<double>>>> VelocityMat;

/// diagonal kinetic velocity, p^d_{mn} = <u_m| 2 tpiba^2 (k+G)_d |u_n>,
/// with the k derivative in the same dimensionless 2*pi/lat0 units
/// build_vkb_dk uses (T = tpiba^2 |k+G|^2 in Ry a.u.)
void kinetic_velocity(int nbands,
                      int npwk,
                      double tpiba2,
                      const std::vector<ModuleBase::Vector3<double>>& gk,
                      const psi::Psi<std::complex<double>>& psi,
                      int ik,
                      VelocityMat& p_mat)
{
    for (int m = 0; m < nbands; ++m)
    {
        for (int n = 0; n < nbands; ++n)
        {
            std::complex<double> dot[3]
                = {std::complex<double>(0.0, 0.0), std::complex<double>(0.0, 0.0), std::complex<double>(0.0, 0.0)};
            for (int ig = 0; ig < npwk; ++ig)
            {
                const std::complex<double> cc = std::conj(psi(ik, m, ig)) * psi(ik, n, ig);
                for (int d = 0; d < 3; ++d)
                {
                    dot[d] += 2.0 * tpiba2 * gk[ig][d] * cc;
                }
            }
            for (int d = 0; d < 3; ++d)
            {
                p_mat[m][n][d] = dot[d];
            }
        }
    }
}

/// projector -> (radial beta index, m channel) table, matching build_vkb
void build_projector_table(const pseudo& ncpp, std::vector<int>& mu_ib, std::vector<int>& mu_m)
{
    const int nh = ncpp.nh;
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

/// becp_b[mu] = <wfc_mu|u_b> for all bands (wfc = vkb or dvkb)
void project_bands(int nbands,
                   int nh,
                   int npwk,
                   const std::vector<std::vector<std::complex<double>>>& wfc,
                   const psi::Psi<std::complex<double>>& psi,
                   int ik,
                   std::vector<std::vector<std::complex<double>>>& becp)
{
    const std::complex<double> zero(0.0, 0.0);
    becp.assign(nbands, std::vector<std::complex<double>>(nh, zero));
    for (int b = 0; b < nbands; ++b)
    {
        for (int mu = 0; mu < nh; ++mu)
        {
            for (int ig = 0; ig < npwk; ++ig)
            {
                becp[b][mu] += std::conj(wfc[mu][ig]) * psi(ik, b, ig);
            }
        }
    }
}

/// accumulate the two Hermitian-conjugate projector terms of one k
/// derivative direction into p_dir[m][n]:
///   <u_m|dvkb_mu> D <vkb|u_n> + <u_m|vkb_mu> D <dvkb|u_n>
/// with D_{mu,nu} = dion(ib_mu, ib_nu) delta_{m_mu, m_nu} (dVnl_dtau layout)
void add_nonlocal_dk(int nh,
                     const pseudo& ncpp,
                     const std::vector<int>& mu_ib,
                     const std::vector<int>& mu_m,
                     const std::vector<std::vector<std::complex<double>>>& becp,
                     const std::vector<std::vector<std::complex<double>>>& dbecp,
                     std::vector<std::vector<std::complex<double>>>& p_dir)
{
    const int nbands = static_cast<int>(p_dir.size());
    for (int m = 0; m < nbands; ++m)
    {
        for (int n = 0; n < static_cast<int>(p_dir[m].size()); ++n)
        {
            std::complex<double> term(0.0, 0.0);
            for (int mu = 0; mu < nh; ++mu)
            {
                std::complex<double> out_m(0.0, 0.0);
                std::complex<double> in_n(0.0, 0.0);
                for (int nu = 0; nu < nh; ++nu)
                {
                    if (mu_m[mu] != mu_m[nu])
                    {
                        continue; // a radial m channel maps onto itself
                    }
                    const double dij = ncpp.dion(mu_ib[mu], mu_ib[nu]);
                    out_m += dij * becp[n][nu];
                    in_n += dij * dbecp[n][nu];
                }
                term += std::conj(dbecp[m][mu]) * out_m + std::conj(becp[m][mu]) * in_n;
            }
            p_dir[m][n] += term;
        }
    }
}

/// nonlocal derivative part of the velocity operator:
///   dV_nl/dk_d = sum_{mu,nu} (|dvkb_mu> D_{mu,nu} <vkb_nu|
///                             + |vkb_mu> D_{mu,nu} <dvkb_nu|)
/// V_loc is k-independent; the DFT+U commutator is the U0 reservation.
void nonlocal_velocity(const UnitCell& ucell,
                       const DFPT_Pert& pert,
                       const std::vector<ModuleBase::Vector3<double>>& gk,
                       const psi::Psi<std::complex<double>>& psi,
                       int ik,
                       int npwk,
                       VelocityMat& p_mat)
{
    const int nbands = psi.get_nbands();
    for (int it = 0; it < ucell.ntype; ++it)
    {
        const pseudo& ncpp = ucell.atoms[it].ncpp;
        const int nh = ncpp.nh;
        if (nh == 0)
        {
            continue;
        }
        if (ncpp.tvanp || ncpp.has_so)
        {
            ModuleBase::WARNING_QUIT("DFPT_Q0::pos_matrix",
                                     "DFPT velocity operator is implemented for "
                                     "normal-conserving separable pseudopotentials only.");
        }
        std::vector<int> mu_ib;
        std::vector<int> mu_m;
        build_projector_table(ncpp, mu_ib, mu_m);
        for (int ia = 0; ia < ucell.atoms[it].na; ++ia)
        {
            std::vector<std::vector<std::complex<double>>> vkb;
            pert.build_vkb(it, ia, gk, vkb);
            // becp_b[mu] = <vkb_mu|u_b> for all bands
            std::vector<std::vector<std::complex<double>>> becp;
            project_bands(nbands, nh, npwk, vkb, psi, ik, becp);
            for (int d = 0; d < 3; ++d)
            {
                std::vector<std::vector<std::complex<double>>> dvkb;
                pert.build_vkb_dk(it, ia, d, gk, vkb, dvkb);
                // dbecp_b[mu] = <dvkb_mu|u_b>
                std::vector<std::vector<std::complex<double>>> dbecp;
                project_bands(nbands, nh, npwk, dvkb, psi, ik, dbecp);
                const std::complex<double> zero(0.0, 0.0);
                std::vector<std::vector<std::complex<double>>> p_dir(
                    nbands, std::vector<std::complex<double>>(nbands, zero));
                add_nonlocal_dk(nh, ncpp, mu_ib, mu_m, becp, dbecp, p_dir);
                for (int m = 0; m < nbands; ++m)
                {
                    for (int n = 0; n < nbands; ++n)
                    {
                        p_mat[m][n][d] += p_dir[m][n];
                    }
                }
            }
        }
    }
}

/// velocity -> position: r = -i v / (tpiba (eps_m - eps_n)), r in bohr
/// (from [H, r] = -i dH/dk in Ry a.u.); degenerate pairs are skipped,
/// their gauge-dependent matrix elements carry no unique value.
void velocity_to_position(int nbands,
                          double tpiba,
                          const ModuleBase::matrix& eig,
                          int ik,
                          const VelocityMat& p_mat,
                          VelocityMat& r_mat_ik)
{
    const double degen_tol = 1.0e-8; ///< empirical parameter: eigenvalue gap (Ry) for the degenerate-pair skip
    for (int m = 0; m < nbands; ++m)
    {
        for (int n = 0; n < nbands; ++n)
        {
            if (m == n)
            {
                continue;
            }
            const double de = eig(ik, m) - eig(ik, n);
            if (std::abs(de) < degen_tol)
            {
                continue;
            }
            for (int d = 0; d < 3; ++d)
            {
                r_mat_ik[m][n][d] = std::complex<double>(0.0, -1.0) * p_mat[m][n][d] / (tpiba * de);
            }
        }
    }
}

} // namespace

void DFPT_Q0::pos_matrix(const psi::Psi<std::complex<double>>& psi,
                         const ModuleBase::matrix& eig,
                         std::vector<std::vector<std::vector<ModuleBase::Vector3<std::complex<double>>>>>& r_mat)
{
    ModuleBase::TITLE("DFPT_Q0", "pos_matrix");
    ModuleBase::timer::start("DFPT_Q0", "pos_matrix");
    const int nk = psi.get_nk();
    const int nbands = psi.get_nbands();
    r_mat.assign(nk,
                 std::vector<std::vector<ModuleBase::Vector3<std::complex<double>>>>(
                     nbands,
                     std::vector<ModuleBase::Vector3<std::complex<double>>>(
                         nbands,
                         ModuleBase::Vector3<std::complex<double>>(0.0, 0.0, 0.0))));
    if (pw_wfc_ == nullptr || ucell_ == nullptr || pert_ == nullptr)
    {
        ModuleBase::timer::end("DFPT_Q0", "pos_matrix");
        return;
    }
    const double tpiba = ucell_->tpiba;
    const double tpiba2 = tpiba * tpiba;
    for (int ik = 0; ik < nk; ++ik)
    {
        const int npwk = pw_wfc_->npwk[ik];
        std::vector<ModuleBase::Vector3<double>> gk(npwk);
        for (int ig = 0; ig < npwk; ++ig)
        {
            gk[ig] = pw_wfc_->getgpluskcar(ik, ig);
        }
        // velocity operator dH/dk matrix elements (kinetic + nonlocal parts)
        VelocityMat p_mat(nbands,
                          std::vector<ModuleBase::Vector3<std::complex<double>>>(
                              nbands,
                              ModuleBase::Vector3<std::complex<double>>(0.0, 0.0, 0.0)));
        kinetic_velocity(nbands, npwk, tpiba2, gk, psi, ik, p_mat);
        nonlocal_velocity(*ucell_, *pert_, gk, psi, ik, npwk, p_mat);
        velocity_to_position(nbands, tpiba, eig, ik, p_mat, r_mat[ik]);
    }
    ModuleBase::timer::end("DFPT_Q0", "pos_matrix");
}

} // namespace ModuleDFPT
