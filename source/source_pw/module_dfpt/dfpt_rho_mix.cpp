// ============================================================
// Response-density mixing family of DFPT_Rho: the per-q input/residual reset,
// the plain / Kerker / Modified-Broyden mix_drho stack, and the
// relative-residual getter. The companion tests live in
// test_serial/dfpt_rho_serial_test.cpp; the response kernels
// (compute_drho / drhoc / v_hartree_q) stay in dfpt_rho.cpp.
// ============================================================

#include "dfpt_rho.h"
#include "dfpt_pw_data.h"

#include "source_base/timer.h"
#include "source_base/tool_title.h"
#include "source_base/module_mixing/plain_mixing.h"
#include "source_base/module_mixing/broyden_mixing.h"
#include "source_base/parallel_reduce.h"

#include <cmath>
#include <complex>
#include <functional>
#include <limits>
#include <vector>

namespace ModuleDFPT
{

namespace
{
/// relative defect ||out - in|| / ||out|| with convergence-safe guards:
/// a zero output norm only reports convergence when the defect itself is
/// zero (v_out = 1 - c v_in maps a nonzero input to zero far from
/// its fixed point), and non-finite norms never report convergence
double safe_relative_defect(double dn2, double o2)
{
    if (!std::isfinite(dn2) || !std::isfinite(o2))
    {
        return std::numeric_limits<double>::max();
    }
    if (o2 <= 0.0)
    {
        return (dn2 <= 0.0) ? 0.0 : std::numeric_limits<double>::max();
    }
    return std::sqrt(dn2 / o2);
}
} // namespace

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
    // Broyden history is per (q, displacement-cycle): drop it so the next
    // displacement starts fresh
    if (q_idx < static_cast<int>(broyden_mdata_.size()))
    {
        broyden_mdata_[q_idx].reset();
    }
    // the shared Broyden mixer must be reset on every valid cycle,
    // independently of whether this q already has a history entry: when a
    // new q is entered for the first time the mixer is still address-bound
    // to the previous q's Mixing_Data, and its first push would exit
    // through the one-object-one-data address guard (broyden_mixing.cpp)
    if (broyden_mixer_ != nullptr)
    {
        broyden_mixer_->reset();
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
    if (mix_type_ == "rho_broyden" && broyden_mixer_ != nullptr)
    {
        // Modified Broyden (Johnson 1988): maintain a per-q history of
        // residual differences, solve a small least-squares each step
        if (q_idx >= static_cast<int>(broyden_mdata_.size()))
        {
            broyden_mdata_.resize(q_idx + 1);
        }
        if (broyden_mdata_[q_idx] == nullptr)
        {
            broyden_mdata_[q_idx] = std::unique_ptr<Base_Mixing::Mixing_Data>(
                new Base_Mixing::Mixing_Data());
            broyden_mixer_->init_mixing_data(*broyden_mdata_[q_idx], npw, sizeof(std::complex<double>));
        }
        Base_Mixing::Mixing_Data& mdata = *broyden_mdata_[q_idx];
        // push the current in/out pair into the history, then solve for
        // the Broyden coefficients and mix
        auto screen = std::function<void(std::complex<double>*)>();
        // the Broyden push_data calls mix(data, data_in, F) unconditionally
        // (broyden_mixing.cpp:47, no nullptr guard): provide the plain
        // preliminary mix data = data_in + beta * F
        const double beta = mix_beta_;
        auto combine =
            std::function<void(std::complex<double>*, const std::complex<double>*, const std::complex<double>*)>(
                [npw, beta](std::complex<double>* data,
                            const std::complex<double>* data_in,
                            const std::complex<double>* F) {
                    for (int i = 0; i < npw; ++i)
                    {
                        data[i] = data_in[i] + beta * F[i];
                    }
                });
        broyden_mixer_->push_data(mdata,
                                  rin.data(),
                                  out.data(),
                                  screen,
                                  combine,
                                  true);
        // the inner product is a pool-global quantity: each rank only
        // holds its G-shell share, so the partial sum must be reduced
        // before the Gram system is assembled (otherwise every rank
        // solves a different least-squares problem)
        broyden_mixer_->cal_coef(mdata,
            std::function<double(std::complex<double>*, std::complex<double>*)>(
            [npw](std::complex<double>* a, std::complex<double>* b) {
                double s = 0.0;
                for (int i = 0; i < npw; ++i)
                {
                    s += (a[i] * std::conj(b[i])).real();
                }
#ifdef __MPI
                Parallel_Reduce::reduce_pool(s);
#endif
                return s;
            }));
        broyden_mixer_->mix_data(mdata, mixed.data());
    }
    else if (mix_type_ == "rho_kerker")
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
    // relative residual ||out - in|| / ||out||; the norms are pool-global
    // quantities (each rank sums its G shells only)
    double dn2 = 0.0;
    double o2 = 0.0;
    for (int ig = 0; ig < npw; ++ig)
    {
        dn2 += std::norm(out[ig] - rin[ig]);
        o2 += std::norm(out[ig]);
    }
#ifdef __MPI
    Parallel_Reduce::reduce_pool(dn2);
    Parallel_Reduce::reduce_pool(o2);
#endif
    residual_[q_idx] = safe_relative_defect(dn2, o2);
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

} // namespace ModuleDFPT
