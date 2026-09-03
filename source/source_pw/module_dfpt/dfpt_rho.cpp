#include "dfpt_rho.h"

#include "dfpt_kq_basis.h"
#include "source_base/constants.h"
#include "source_base/global_function.h"
#include "source_base/timer.h"
#include "source_base/tool_title.h"
#include "source_base/module_mixing/plain_mixing.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <functional>
#include <vector>

namespace ModuleDFPT
{

DFPT_Rho::DFPT_Rho()
{
}

// Defined here (after plain_mixing.h is included) so the unique_ptr deleter
// sees the complete Plain_Mixing type. The body is implicitly generated.
DFPT_Rho::~DFPT_Rho() = default;

void DFPT_Rho::init(const Config& cfg)
{
    ModuleBase::TITLE("DFPT_Rho", "init");
    ModuleBase::timer::start("DFPT_Rho", "init");
    nspin_ = cfg.nspin;
    nrxx_ = cfg.nrxx;
    pw_rho_ = cfg.pw_rho;
    pw_wfc_ = cfg.pw_wfc;
    recip_matrix_ = cfg.recip_matrix;
    mix_beta_ = cfg.mix_beta;
    mix_type_ = cfg.mix_type;
    kerker_a2_ = cfg.kerker_a2;
    if (cfg.mix_type != "plain" && cfg.mix_type != "kerker")
    {
        ModuleBase::WARNING_QUIT("DFPT_Rho", "unsupported mix_type, expected plain or kerker");
    }
    // make_unique is C++14; C++11 uses new directly through reset()
    mixer_.reset(new Base_Mixing::Plain_Mixing(mix_beta_));
    ModuleBase::timer::end("DFPT_Rho", "init");
}

bool DFPT_Rho::is_gamma_q_(const ModuleBase::Vector3<double>& q)
{
    const double gamma_tol = 1.0e-10; ///< empirical parameter: fractional-q tolerance for the Gamma test
    return std::abs(q.x) < gamma_tol && std::abs(q.y) < gamma_tol && std::abs(q.z) < gamma_tol;
}

void DFPT_Rho::add_band_(int ik,
                         int ib,
                         double w,
                         const std::complex<double>* c_ptr,
                         const std::vector<std::complex<double>>& dpsi,
                         const DFPT_KQ_Basis& kq,
                         std::vector<std::complex<double>>& a_r)
{
    const int npw_kq = kq.get_npwk();
    // k+q G index -> rho-grid ig (both bases share the FFT cell)
    std::vector<int> kq2rho(npw_kq, -1);
    for (int igl = 0; igl < npw_kq; ++igl)
    {
        kq2rho[igl] = kq.get_ig_rho(igl);
    }
    // periodic part u_nk(r) on the shared grid (phase-free FFT)
    std::vector<std::complex<double>> u_r(pw_rho_->nrxx);
    pw_wfc_->recip2real(c_ptr, u_r.data(), ik);
    // periodic part du_nk(r): scatter the k+q coefficients onto the
    // rho grid and transform (same convention, so the product is
    // consistent with the u transform)
    std::vector<std::complex<double>> d_recip(pw_rho_->npw, std::complex<double>(0.0, 0.0));
    const int nd = std::min(npw_kq, static_cast<int>(dpsi.size()));
    for (int igl = 0; igl < nd; ++igl)
    {
        if (kq2rho[igl] >= 0)
        {
            d_recip[kq2rho[igl]] = dpsi[igl];
        }
    }
    std::vector<std::complex<double>> d_r(pw_rho_->nrxx);
    pw_rho_->recip2real(d_recip.data(), d_r.data());
    // same normalization as the GS density accumulation
    // (elecstate_pw.cpp rhoBandK: w1 = wg / omega), including the
    // spin factor 2: QE incdrhoscf uses wgt = 2 * weight / omega
    // at every q (the factor 2 is the spin degeneracy, not a
    // Hermitian completion)
    const double w1 = 2.0 * w / pw_rho_->omega;
    for (int ir = 0; ir < pw_rho_->nrxx; ++ir)
    {
        a_r[ir] += w1 * std::conj(u_r[ir]) * d_r[ir];
    }
}

void DFPT_Rho::make_drho_r_(const std::vector<std::complex<double>>& drho_g,
                            const ModuleBase::Vector3<double>& q_frac,
                            std::vector<double>& drho_r) const
{
    // real-space manifest density: at q = 0 the completed coefficients are
    // already the full (real) response; away from q = 0 the manifest is the
    // real combination 2 Re[e^{i q r} A(r)] of the one-sided amplitude
    std::vector<std::complex<double>> a_clean(pw_rho_->nrxx);
    pw_rho_->recip2real(drho_g.data(), a_clean.data());
    drho_r.resize(pw_rho_->nrxx);
    if (is_gamma_q_(q_frac))
    {
        for (int ir = 0; ir < pw_rho_->nrxx; ++ir)
        {
            drho_r[ir] = a_clean[ir].real();
        }
    }
    else
    {
        for (int ix = 0; ix < pw_rho_->nx; ++ix)
        {
            for (int iy = 0; iy < pw_rho_->ny; ++iy)
            {
                for (int iz = 0; iz < pw_rho_->nz; ++iz)
                {
                    const int ir = (ix * pw_rho_->ny + iy) * pw_rho_->nz + iz;
                    const double theta
                        = ModuleBase::TWO_PI
                          * (q_frac.x * ix / pw_rho_->nx + q_frac.y * iy / pw_rho_->ny + q_frac.z * iz / pw_rho_->nz);
                    drho_r[ir] = 2.0 * (a_clean[ir].real() * std::cos(theta) - a_clean[ir].imag() * std::sin(theta));
                }
            }
        }
    }
}

void DFPT_Rho::zero_neg_q_(const ModuleBase::Vector3<double>& q_cart,
                           std::vector<std::complex<double>>& drho_g) const
{
    // charge conservation: the Delta = -q harmonic (G+q = 0 component of the
    // response density) must vanish whenever -q falls on a reciprocal
    // lattice vector; for a generic q inside the cell this never triggers
    const ModuleBase::Vector3<double> mq_cart(-q_cart.x, -q_cart.y, -q_cart.z);
    const ModuleBase::Vector3<double> mfrac = mq_cart * recip_matrix_.Inverse();
    const double mr[3] = {std::round(mfrac.x), std::round(mfrac.y), std::round(mfrac.z)};
    const double gvec_match_tol = 1.0e-6; ///< empirical parameter: folded -q reciprocal-vector match tolerance
    if (std::abs(mfrac.x - mr[0]) < gvec_match_tol && std::abs(mfrac.y - mr[1]) < gvec_match_tol
        && std::abs(mfrac.z - mr[2]) < gvec_match_tol)
    {
        // locate the rho-grid G equal to -q through its FFT cell
        const int cix = (static_cast<int>(mr[0]) % pw_rho_->nx + pw_rho_->nx) % pw_rho_->nx;
        const int ciy = (static_cast<int>(mr[1]) % pw_rho_->ny + pw_rho_->ny) % pw_rho_->ny;
        const int ciz = (static_cast<int>(mr[2]) % pw_rho_->nz + pw_rho_->nz) % pw_rho_->nz;
        for (int ig = 0; ig < pw_rho_->npw; ++ig)
        {
            const int isz = pw_rho_->ig2isz[ig];
            const int iz = isz % pw_rho_->nz;
            const int is = isz / pw_rho_->nz;
            const int ixy = pw_rho_->is2fftixy[is];
            const int ix = ixy / pw_rho_->fftny;
            const int iy = ixy % pw_rho_->fftny;
            if (ix == cix && iy == ciy && iz == ciz)
            {
                drho_g[ig] = std::complex<double>(0.0, 0.0);
                break;
            }
        }
    }
}

void DFPT_Rho::compute_drho(const psi::Psi<std::complex<double>>& psi,
                            const ModuleBase::matrix& wg,
                            int q_idx,
                            DFPT_PW_Data& data)
{
    ModuleBase::TITLE("DFPT_Rho", "compute_drho");
    ModuleBase::timer::start("DFPT_Rho", "compute_drho");
    if (pw_rho_ == nullptr || pw_wfc_ == nullptr)
    {
        ModuleBase::timer::end("DFPT_Rho", "compute_drho");
        return;
    }
    if (nspin_ != 1)
    {
        ModuleBase::WARNING_QUIT("DFPT_Rho", "only nspin = 1 is supported in the design phase");
    }
    const int nk = psi.get_nk();
    const int nbands = psi.get_nbands();
    const ModuleBase::Vector3<double> q_frac = data.get_qvec(q_idx);
    const ModuleBase::Vector3<double> q_cart = q_frac * recip_matrix_;

    std::vector<std::complex<double>> a_r(pw_rho_->nrxx, std::complex<double>(0.0, 0.0));
    DFPT_KQ_Basis kq;
    for (int ik = 0; ik < nk; ++ik)
    {
        kq.init(pw_wfc_, pw_rho_, q_cart, ik);
        for (int ib = 0; ib < nbands; ++ib)
        {
            if (!dfpt_band_occupied(wg, ik, ib))
            {
                continue; // unoccupied band: no contribution to the density
            }
            const double w = wg(ik, ib);
            const std::vector<std::complex<double>> dpsi = data.get_dpsi(q_idx, ik, ib);
            add_band_(ik, ib, w, &psi(ik, ib, 0), dpsi, kq, a_r);
        }
    }

    // Hermitian completion at q = 0: the band loop above stores only the
    // u_n^* du_n piece (spin factor 2 already included in w1); the physical
    // (real) response density at q = 0 also contains the du_n u_n^* piece,
    // which coincides with the conjugate of the stored one, so keep only
    // the real part of the amplitude before the FFT. The resulting
    // coefficients are exactly Hermitian on the sphere, including
    // one-sided sticks whose -G falls outside it. Away from q = 0 the +q
    // harmonic of the response is exactly the one-sided object and no
    // completion applies.
    if (is_gamma_q_(q_frac))
    {
        for (int ir = 0; ir < pw_rho_->nrxx; ++ir)
        {
            a_r[ir] = std::complex<double>(a_r[ir].real(), 0.0);
        }
    }

    // q-shifted coefficients A_Delta on the rho grid
    std::vector<std::complex<double>> drho_g(pw_rho_->npw);
    pw_rho_->real2recip(a_r.data(), drho_g.data());

    // charge conservation: zero the Delta = -q harmonic
    zero_neg_q_(q_cart, drho_g);
    data.set_drho_g(q_idx, 0, drho_g);

    // real-space manifest density
    std::vector<double> drho_r;
    make_drho_r_(drho_g, q_frac, drho_r);
    data.set_drho_r(q_idx, 0, drho_r);

    // remember the freshly computed output for the mixing step
    if (q_idx >= static_cast<int>(drho_out_.size()))
    {
        drho_out_.resize(q_idx + 1);
    }
    drho_out_[q_idx].assign(1, drho_g);
    ModuleBase::timer::end("DFPT_Rho", "compute_drho");
}

void DFPT_Rho::cal_docc(const psi::Psi<std::complex<double>>& psi,
                        const ModuleBase::matrix& wg,
                        int q_idx,
                        DFPT_PW_Data& data)
{
    ModuleBase::TITLE("DFPT_Rho", "cal_docc");
    ModuleBase::timer::start("DFPT_Rho", "cal_docc");
    // Reserved first-order occupation matrix (docc) for DFT+U (U0).
    // The physical cross terms need the beta projectors at both k and k+q
    // (a PW-side adapter of the build_vkb machinery); they land together
    // with the Plus_U production wiring in the C7/U1 window, when dpsi and
    // a usable Plus_U provider coexist. Pure-PW runs keep u_active() false
    // and never reach this accumulation:
    //   cross term: Re(becp(k+q, dpsi) * becp(k, psi))      (response)
    //   frozen term: becp(k, psi) * dbecp_f(k, psi)         (GS k)
    if (!data.with_u())
    {
        ModuleBase::timer::end("DFPT_Rho", "cal_docc");
        return;
    }
    (void)psi;
    (void)wg;
    (void)q_idx;
    (void)data;
    ModuleBase::timer::end("DFPT_Rho", "cal_docc");
}

void DFPT_Rho::reset_mixing(int q_idx)
{
    ModuleBase::TITLE("DFPT_Rho", "reset_mixing");
    ModuleBase::timer::start("DFPT_Rho", "reset_mixing");
    if (q_idx < 0)
    {
        ModuleBase::timer::end("DFPT_Rho", "reset_mixing");
        return;
    }
    if (q_idx < static_cast<int>(drho_in_.size()))
    {
        drho_in_[q_idx].clear();
    }
    if (q_idx < static_cast<int>(residual_.size()))
    {
        residual_[q_idx] = 0.0;
    }
    ModuleBase::timer::end("DFPT_Rho", "reset_mixing");
}

void DFPT_Rho::mix_drho(int q_idx, DFPT_PW_Data& data)
{
    ModuleBase::TITLE("DFPT_Rho", "mix_drho");
    ModuleBase::timer::start("DFPT_Rho", "mix_drho");
    if (mixer_ == nullptr || pw_rho_ == nullptr)
    {
        ModuleBase::timer::end("DFPT_Rho", "mix_drho");
        return;
    }
    const std::vector<std::complex<double>> out = data.get_drho_g(q_idx, 0);
    if (out.empty() || static_cast<int>(out.size()) != pw_rho_->npw)
    {
        ModuleBase::timer::end("DFPT_Rho", "mix_drho");
        return;
    }
    const int npw = pw_rho_->npw;
    if (q_idx >= static_cast<int>(drho_in_.size()))
    {
        drho_in_.resize(q_idx + 1);
        residual_.resize(q_idx + 1, 0.0);
    }
    // first iteration starts from a zero input density
    if (drho_in_[q_idx].empty())
    {
        drho_in_[q_idx].assign(1, std::vector<std::complex<double>>(npw, std::complex<double>(0.0, 0.0)));
    }
    const std::vector<std::complex<double>>& rin = drho_in_[q_idx][0];
    std::vector<std::complex<double>> mixed(npw);
    const double w2_floor = 1.0e-12; ///< empirical parameter: |G+q|^2 zero-shell guard (1/lat0^2, Kerker freeze)
    // the fractional q is needed both by the Kerker screen and by the
    // real-space manifest below; the q-shifted |G+q| convention matches
    // v_hartree_q (gcar + q_frac * recip, 1/lat0^2 units)
    const ModuleBase::Vector3<double> q_frac = data.get_qvec(q_idx);
    if (mix_type_ == "kerker")
    {
        const ModuleBase::Vector3<double> q_cart = q_frac * recip_matrix_;
        std::vector<std::complex<double>> rin_s(npw);
        std::vector<std::complex<double>> out_s(npw);
        std::vector<std::complex<double>> mixed_s(npw);
        for (int ig = 0; ig < npw; ++ig)
        {
            const ModuleBase::Vector3<double> w = pw_rho_->gcar[ig] + q_cart;
            const double w2 = w * w;
            // |G+q| = 0 harmonic: f = 0, frozen at rin (that harmonic is
            // dropped by compute_drho, so both inputs are zero there)
            const double f = (w2 < w2_floor) ? 0.0 : w2 / (w2 + kerker_a2_);
            rin_s[ig] = f * rin[ig];
            out_s[ig] = f * out[ig];
        }
        mixer_->plain_mix(mixed_s.data(),
                          rin_s.data(),
                          out_s.data(),
                          npw,
                          std::function<void(std::complex<double>*)>());
        // add back the screened-out part: mixed = rin + beta f (out - rin),
        // i.e. a plain mix with the per-shell coefficient beta f_g while
        // the stored density stays physical (not screen-scaled)
        for (int ig = 0; ig < npw; ++ig)
        {
            mixed[ig] = rin[ig] + (mixed_s[ig] - rin_s[ig]);
        }
    }
    else
    {
        mixer_->plain_mix(mixed.data(), rin.data(), out.data(), npw, std::function<void(std::complex<double>*)>());
    }
    // relative residual ||out - in|| / ||out||
    double dn2 = 0.0;
    double o2 = 0.0;
    for (int ig = 0; ig < npw; ++ig)
    {
        dn2 += std::norm(out[ig] - rin[ig]);
        o2 += std::norm(out[ig]);
    }
    residual_[q_idx] = (o2 > 0.0) ? std::sqrt(dn2 / o2) : 0.0;
    drho_in_[q_idx][0] = mixed;
    data.set_drho_g(q_idx, 0, mixed);

    // rebuild the real-space manifest from the mixed coefficients (q = 0:
    // completed coefficients are the full real response; otherwise the
    // one-sided 2 Re[e^{i q r} A(r)] manifest)
    std::vector<double> drho_r;
    make_drho_r_(mixed, q_frac, drho_r);
    data.set_drho_r(q_idx, 0, drho_r);
    ModuleBase::timer::end("DFPT_Rho", "mix_drho");
}

double DFPT_Rho::get_residual(int q_idx, DFPT_PW_Data& data) const
{
    ModuleBase::TITLE("DFPT_Rho", "get_residual");
    ModuleBase::timer::start("DFPT_Rho", "get_residual");
    (void)data;
    if (q_idx < 0 || q_idx >= static_cast<int>(residual_.size()))
    {
        ModuleBase::timer::end("DFPT_Rho", "get_residual");
        return 0.0;
    }
    ModuleBase::timer::end("DFPT_Rho", "get_residual");
    return residual_[q_idx];
}

void DFPT_Rho::v_hartree_q(const ModuleBase::Vector3<double>& q_cart,
                           const std::vector<std::complex<double>>& drho_g,
                           std::vector<std::complex<double>>& dv_ha_g) const
{
    ModuleBase::TITLE("DFPT_Rho", "v_hartree_q");
    ModuleBase::timer::start("DFPT_Rho", "v_hartree_q");
    if (pw_rho_ == nullptr)
    {
        dv_ha_g.clear();
        ModuleBase::timer::end("DFPT_Rho", "v_hartree_q");
        return;
    }
    const int npw = pw_rho_->npw;
    if (static_cast<int>(drho_g.size()) != npw)
    {
        dv_ha_g.clear();
        ModuleBase::timer::end("DFPT_Rho", "v_hartree_q");
        return;
    }
    dv_ha_g.assign(npw, std::complex<double>(0.0, 0.0));
    const double w2_floor = 1.0e-12; ///< empirical parameter: |G+q|^2 zero-shell guard (1/lat0^2, Coulomb skip)
    for (int ig = 0; ig < npw; ++ig)
    {
        const ModuleBase::Vector3<double> w = pw_rho_->gcar[ig] + q_cart;
        const double w2_lat0 = w * w; // 1/lat0^2 units, like pw_rho_->gg
        // skip |G+q| = 0 (ig = -q): the q-shifted G=0 harmonic of the
        // Hartree kernel (v_hartree skips ig_gge0 the same way)
        if (w2_lat0 < w2_floor)
        {
            continue;
        }
        const double fac = ModuleBase::e2 * ModuleBase::FOUR_PI / (pw_rho_->tpiba2 * w2_lat0);
        dv_ha_g[ig] = fac * drho_g[ig];
    }
    ModuleBase::timer::end("DFPT_Rho", "v_hartree_q");
}

} // namespace ModuleDFPT
