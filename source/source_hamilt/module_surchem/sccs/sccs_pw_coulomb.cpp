#include "sccs_pw_coulomb.h"

#include "source_base/constants.h"
#include "source_basis/module_pw/pw_basis.h"

#include <cmath>
#include <complex>
#include <stdexcept>

namespace ModuleSccs
{

std::vector<ModuleBase::Vector3<double>> periodic_gradient(
    const std::vector<double>& values,
    const ModulePW::PW_Basis& basis,
    const double tpiba)
{
    if (values.size() != static_cast<std::size_t>(basis.nrxx))
    {
        throw std::invalid_argument("SCCS scalar array does not match the local PW real-space grid");
    }
    if (!std::isfinite(tpiba) || tpiba <= 0.0 || basis.npw <= 0 || basis.nrxx <= 0
        || basis.gcar == nullptr)
    {
        throw std::invalid_argument("SCCS periodic gradient requires an initialized PW basis and positive tpiba");
    }

    std::vector<std::complex<double>> values_g(basis.npw);
    basis.real2recip(values.data(), values_g.data());
    std::vector<std::complex<double>> gradient_g(basis.npw);
    std::vector<double> gradient_r(basis.nrxx);
    std::vector<ModuleBase::Vector3<double>> gradient(basis.nrxx);
    for (int direction = 0; direction < 3; ++direction)
    {
        for (int ig = 0; ig < basis.npw; ++ig)
        {
            gradient_g[ig]
                = ModuleBase::IMAG_UNIT * tpiba * basis.gcar[ig][direction] * values_g[ig];
        }
        basis.recip2real(gradient_g.data(), gradient_r.data());
        for (int ir = 0; ir < basis.nrxx; ++ir)
        {
            gradient[ir][direction] = gradient_r[ir];
        }
    }
    return gradient;
}

PeriodicCoulombOperator::PeriodicCoulombOperator(const ModulePW::PW_Basis& basis,
                                                 const double tpiba)
    : basis_(basis), tpiba_(tpiba)
{
    if (!std::isfinite(tpiba_) || tpiba_ <= 0.0)
    {
        throw std::invalid_argument("SCCS periodic Coulomb operator requires positive finite tpiba");
    }
    if (basis_.npw <= 0 || basis_.nrxx <= 0 || basis_.gg == nullptr || basis_.gcar == nullptr)
    {
        throw std::invalid_argument("SCCS periodic Coulomb operator requires an initialized PW basis");
    }
}

void PeriodicCoulombOperator::apply(const std::vector<double>& charge,
                                    ElectrostaticField& field) const
{
    if (charge.size() != static_cast<std::size_t>(basis_.nrxx))
    {
        throw std::invalid_argument("SCCS charge array does not match the local PW real-space grid");
    }

    std::vector<std::complex<double>> charge_g(basis_.npw);
    std::vector<std::complex<double>> potential_g(basis_.npw);
    basis_.real2recip(charge.data(), charge_g.data());
    const double tpiba2 = tpiba_ * tpiba_;
    for (int ig = 0; ig < basis_.npw; ++ig)
    {
        if (ig == basis_.ig_gge0 || basis_.gg[ig] == 0.0)
        {
            potential_g[ig] = std::complex<double>();
        }
        else
        {
            potential_g[ig] = ModuleBase::FOUR_PI * charge_g[ig] / (tpiba2 * basis_.gg[ig]);
        }
    }

    field.potential.resize(basis_.nrxx);
    field.gradient.assign(basis_.nrxx, ModuleBase::Vector3<double>());
    basis_.recip2real(potential_g.data(), field.potential.data());

    std::vector<std::complex<double>> gradient_g(basis_.npw);
    std::vector<double> gradient_r(basis_.nrxx);
    for (int direction = 0; direction < 3; ++direction)
    {
        for (int ig = 0; ig < basis_.npw; ++ig)
        {
            gradient_g[ig] = ModuleBase::IMAG_UNIT * tpiba_ * basis_.gcar[ig][direction]
                             * potential_g[ig];
        }
        basis_.recip2real(gradient_g.data(), gradient_r.data());
        for (int ir = 0; ir < basis_.nrxx; ++ir)
        {
            field.gradient[ir][direction] = gradient_r[ir];
        }
    }
}

} // namespace ModuleSccs
