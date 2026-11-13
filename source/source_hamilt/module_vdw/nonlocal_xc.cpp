#include "nonlocal_xc.h"

#include "rvv10_workspace.h"
#include "source_base/constants.h"
#include "source_base/matrix.h"
#include "source_base/parallel_device.h"
#include "source_base/timer.h"
#include "source_base/tool_title.h"
#include "source_basis/module_pw/pw_basis.h"
#include "source_estate/module_charge/charge.h"
#include "source_hamilt/module_xc/xc_functional.h"

#include <cassert>
#include <complex>
#include <vector>

namespace vdw
{
namespace
{

/// Spin channels of chr.rho that add up to the valence density (for nspin = 4, rho[0] is the total density).
int density_channels(const Charge& chr)
{
    return (chr.nspin == 2) ? 2 : 1;
}

/// Total density of chr (valence plus core), its plane-wave gradient (a copy of total_gradient when that
/// is given) and sigma = |grad n|^2.
void total_density(const Charge& chr,
                   const std::vector<ModuleBase::Vector3<double>>* total_gradient,
                   std::vector<double>& rho,
                   std::vector<ModuleBase::Vector3<double>>& grad,
                   std::vector<double>& sigma)
{
    const ModulePW::PW_Basis& pw = *chr.rhopw;
    const int nrxx = chr.nrxx;
    const int nspin0 = density_channels(chr);
    const double* core_begin = chr.rho_core;
    const double* core_end = core_begin + nrxx;
    rho.assign(core_begin, core_end);
    for (int is = 0; is < nspin0; ++is)
    {
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
        for (int ir = 0; ir < nrxx; ++ir)
        {
            rho[ir] += chr.rho[is][ir];
        }
    }
    // the branch is collective (the transforms below are): the caller decides it for the whole pool
    if (total_gradient != nullptr)
    {
        assert(total_gradient->size() == static_cast<std::size_t>(nrxx));
        grad = *total_gradient;
    }
    else
    {
        std::vector<std::complex<double>> rhog(pw.npw);
        pw.real2recip(rho.data(), rhog.data());
        grad.resize(nrxx);
        XC_Functional::grad_rho(rhog.data(), grad.data(), &pw, pw.tpiba);
    }
    sigma.resize(nrxx);
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
    for (int ir = 0; ir < nrxx; ++ir)
    {
        sigma[ir] = grad[ir].norm2();
    }
}

/// rVV10 energy (Ha) of the total density of chr, with the density gradient (see total_density), dE/dn
/// (unless dedn is nullptr), dE/dsigma and, with stress != nullptr, the explicit strain derivative.
double evaluate_total_density(const XCFunctionalParameters& parameters,
                              const Charge& chr,
                              const std::vector<ModuleBase::Vector3<double>>* total_gradient,
                              std::vector<ModuleBase::Vector3<double>>& grad,
                              std::vector<double>* dedn,
                              std::vector<double>& dedsigma,
                              double* stress)
{
    std::vector<double> rho;
    std::vector<double> sigma;
    total_density(chr, total_gradient, rho, grad, sigma);
    if (dedn != nullptr)
    {
        dedn->resize(chr.nrxx);
    }
    dedsigma.resize(chr.nrxx);
    return Rvv10Workspace::instance().evaluate(*chr.rhopw, parameters, rho, sigma, dedn, dedsigma, stress);
}
} // namespace

void add_rvv10(const XCFunctionalParameters& parameters,
               const Charge& chr,
               double& etxc,
               double& vtxc,
               ModuleBase::matrix& v)
{
    add_rvv10(parameters, chr, nullptr, etxc, vtxc, v);
}

void add_rvv10(const XCFunctionalParameters& parameters,
               const Charge& chr,
               const std::vector<ModuleBase::Vector3<double>>* total_gradient,
               double& etxc,
               double& vtxc,
               ModuleBase::matrix& v)
{
    if (!parameters.rvv10)
    {
        return;
    }
    ModuleBase::TITLE("vdw", "add_rvv10");
    ModuleBase::timer::start("vdw", "add_rvv10");
    const ModulePW::PW_Basis& pw = *chr.rhopw;
    const int nrxx = chr.nrxx;
    std::vector<ModuleBase::Vector3<double>> grad;
    std::vector<double> dedn;
    std::vector<double> dedsigma;
    const double energy = evaluate_total_density(parameters, chr, total_gradient, grad, &dedn, dedsigma, nullptr);
    // v = dE/dn - div(2 dE/dsigma grad n), Hartree to Ry by e2
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
    for (int ir = 0; ir < nrxx; ++ir)
    {
        grad[ir] *= 2.0 * dedsigma[ir];
    }
    std::vector<double> divergence(nrxx);
    XC_Functional::grad_dot(grad.data(), divergence.data(), &pw, pw.tpiba);
    // for nspin = 4 the potential of the total density goes to v(0) only
    const int nspin0 = density_channels(chr);
    // this reduction and the one of the stress stay serial, so the result does not depend on the thread count
    double vtxc_nl = 0.0;
    for (int ir = 0; ir < nrxx; ++ir)
    {
        const double v_nl = ModuleBase::e2 * (dedn[ir] - divergence[ir]);
        for (int is = 0; is < nspin0; ++is)
        {
            v(is, ir) += v_nl;
            vtxc_nl += v_nl * chr.rho[is][ir];
        }
    }
#ifdef __MPI
    Parallel_Common::reduce_data(&vtxc_nl, 1, pw.pool_world);
#endif
    etxc += ModuleBase::e2 * energy;
    vtxc += vtxc_nl * pw.omega / pw.nxyz;
    ModuleBase::timer::end("vdw", "add_rvv10");
}

void add_rvv10_stress(const XCFunctionalParameters& parameters, const Charge& chr, ModuleBase::matrix& sigma)
{
    if (!parameters.rvv10)
    {
        return;
    }
    ModuleBase::TITLE("vdw", "add_rvv10_stress");
    ModuleBase::timer::start("vdw", "add_rvv10_stress");
    const ModulePW::PW_Basis& pw = *chr.rhopw;
    const int nrxx = chr.nrxx;
    std::vector<ModuleBase::Vector3<double>> grad;
    std::vector<double> dedsigma;
    double explicit_stress[rvv10_stress_components] = {0.0};
    const double energy = evaluate_total_density(parameters, chr, nullptr, grad, nullptr, dedsigma, explicit_stress);
    // integral of dE/dsigma d_a n d_b n, in units of the volume element
    double gradient_term[rvv10_stress_components] = {0.0};
    for (int ir = 0; ir < nrxx; ++ir)
    {
        const double g[3] = {grad[ir].x, grad[ir].y, grad[ir].z};
        for (int a = 0; a < 3; ++a)
        {
            for (int b = 0; b < 3; ++b)
            {
                gradient_term[3 * a + b] += dedsigma[ir] * g[a] * g[b];
            }
        }
    }
#ifdef __MPI
    Parallel_Common::reduce_data(gradient_term, rvv10_stress_components, pw.pool_world);
#endif
    // The full strain derivative is T - delta_ab Int(n dE/dn + 2 sigma dE/dsigma) - 2 Int dE/dsigma d_a n d_b n,
    // T the explicit part from libvdwxc. The XC diagonal -(etxc - vtxc)/omega in PW (in LCAO -etxc/omega plus the
    // Pulay term of the potential) already holds -e2 (E - Int(n dE/dn + 2 sigma dE/dsigma)) / omega; the rest is below.
    for (int a = 0; a < 3; ++a)
    {
        for (int b = 0; b < 3; ++b)
        {
            const double diagonal = (a == b) ? energy : 0.0;
            sigma(a, b) += -ModuleBase::e2 / pw.omega * (explicit_stress[3 * a + b] - diagonal)
                           + 2.0 * ModuleBase::e2 * gradient_term[3 * a + b] / pw.nxyz;
        }
    }
    ModuleBase::timer::end("vdw", "add_rvv10_stress");
}

void finalize_rvv10()
{
    Rvv10Workspace::release();
}

} // namespace vdw
