#include "xc_rvv10_pw.h"

#include "source_base/parallel_reduce.h"
#include "source_basis/module_pw/pw_basis.h"
#include "xc_functional.h"

#include <cmath>
#include <complex>
#include <stdexcept>

namespace Rvv10
{
Evaluator::Evaluator(double b, double c) : model_(b, c)
{
}

Evaluation Evaluator::evaluate(const ModulePW::PW_Basis& pw,
                                     const std::vector<double>& total,
                                     const std::vector<double>& valence) const
{
    // Reject unsupported layouts before any collective FFT can start. The
    // existing PW transforms already handle pool-distributed real and G grids.
    if (pw.gamma_only || pw.get_device() != "cpu" || pw.get_precision() != "double" || pw.nrxx <= 0
        || pw.nxyz <= 0 || (pw.npw > 0 && (!pw.gcar || !pw.gg)) || !std::isfinite(pw.omega)
        || pw.omega <= 0 || !std::isfinite(pw.tpiba) || pw.tpiba <= 0
        || total.size() != static_cast<std::size_t>(pw.nrxx) || valence.size() != total.size())
    {
        throw std::invalid_argument("rVV10 evaluator requires an initialized CPU/double full-complex density "
                                    "basis and matching local density arrays");
    }
    for (int r = 0; r < pw.nrxx; ++r)
    {
        if (!std::isfinite(total[r]) || !std::isfinite(valence[r]))
            throw std::invalid_argument("rVV10 density contains a nonfinite value");
    }
    for (int g = 0; g < pw.npw; ++g)
    {
        if (!std::isfinite(pw.gg[g]) || pw.gg[g] < 0 || std::sqrt(pw.gg[g]) * pw.tpiba >= kernel_.maximum_g())
            throw std::out_of_range("rVV10 density basis exceeds the reciprocal kernel table");
    }

    typedef std::complex<double> Complex;
    std::vector<Complex> density_g(pw.npw);
    std::vector<ModuleBase::Vector3<double>> gradient(pw.nrxx);
    std::vector<ModuleBase::Vector3<double>> field(pw.nrxx);
    pw.real2recip(total.data(), density_g.data());
    XC_Functional::grad_rho(density_g.data(), gradient.data(), &pw, pw.tpiba);
    std::vector<LocalValues> local(pw.nrxx);
    for (int r = 0; r < pw.nrxx; ++r)
    {
        const auto& v = gradient[r];
        local[r] = model_.evaluate(total[r], v.x * v.x + v.y * v.y + v.z * v.z);
        field[r] = ModuleBase::Vector3<double>(0, 0, 0);
    }

    // Channel-major storage is transformed in place from theta(G) to u(G).
    // Regenerate local spline fields rather than storing 60 doubles per grid point.
    const std::size_t ng = pw.npw;
    std::vector<Complex> channels(channel_count * ng);
    std::vector<double> work(pw.nrxx);
    for (int a = 0; a < channel_count; ++a)
    {
        for (int r = 0; r < pw.nrxx; ++r)
            work[r] = basis_.channel(local[r], a).theta;
        pw.real2recip(work.data(), channels.data() + a * ng);
    }

    Evaluation result = {0.0, 0.0, std::vector<double>(pw.nrxx, model_.beta())};
    for (int g = 0; g < pw.npw; ++g)
    {
        const auto k = kernel_.values(std::sqrt(pw.gg[g]) * pw.tpiba);
        std::array<Complex, channel_count> theta;
        for (int a = 0; a < channel_count; ++a)
            theta[a] = channels[a * ng + g];
        for (int a = 0; a < channel_count; ++a)
        {
            Complex u = 0.0;
            for (int b = 0; b < channel_count; ++b)
                u += k[a * channel_count + b] * theta[b];
            result.energy += 0.5 * pw.omega * std::real(std::conj(theta[a]) * u);
            channels[a * ng + g] = u;
        }
    }
    for (int a = 0; a < channel_count; ++a)
    {
        pw.recip2real(channels.data() + a * ng, work.data());
        for (int r = 0; r < pw.nrxx; ++r)
        {
            const auto t = basis_.channel(local[r], a);
            result.potential[r] += work[r] * t.dtheta_dn;
            const double coefficient = 2 * work[r] * t.dtheta_dsigma;
            field[r].x += coefficient * gradient[r].x;
            field[r].y += coefficient * gradient[r].y;
            field[r].z += coefficient * gradient[r].z;
        }
    }
    XC_Functional::grad_dot(field.data(), work.data(), &pw, pw.tpiba);
    const double dv = pw.omega / pw.nxyz;
    for (int r = 0; r < pw.nrxx; ++r)
    {
        result.energy += dv * model_.beta() * total[r];
        result.potential[r] -= work[r];
        result.vtxc += dv * valence[r] * result.potential[r];
    }
    Parallel_Reduce::reduce_pool(result.energy);
    Parallel_Reduce::reduce_pool(result.vtxc);
    return result;
}
} // namespace Rvv10
