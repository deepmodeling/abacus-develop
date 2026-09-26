#include "sccs_pw_coulomb.h"

#include "source_base/constants.h"
#include "source_basis/module_pw/pw_basis.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <stdexcept>

namespace ModuleSccs
{

namespace
{

void adjoint_gradient_transform(
    const std::vector<ModuleBase::Vector3<double>>& field,
    const ModulePW::PW_Basis& basis,
    const double tpiba,
    const bool apply_coulomb,
    std::vector<double>& component,
    std::vector<std::complex<double>>& component_g,
    std::vector<std::complex<double>>& sum)
{
    if (field.size() != static_cast<std::size_t>(basis.nrxx)
        || !std::isfinite(tpiba) || tpiba <= 0.0)
    {
        throw std::invalid_argument("SCCS gradient adjoint requires a matching PW grid");
    }
    component.resize(basis.nrxx);
    component_g.resize(basis.npw);
    sum.resize(basis.npw);
    // Unlike the transformed component, the divergence is accumulated.
    const std::complex<double> zero;
    std::fill(sum.begin(), sum.end(), zero);
    for (int direction = 0; direction < 3; ++direction)
    {
        for (int ir = 0; ir < basis.nrxx; ++ir)
        {
            component[ir] = field[ir][direction];
        }
        basis.real2recip(component.data(), component_g.data());
        for (int ig = 0; ig < basis.npw; ++ig)
        {
            sum[ig] -= ModuleBase::IMAG_UNIT * tpiba * basis.gcar[ig][direction] * component_g[ig];
        }
    }
    if (apply_coulomb)
    {
        for (int ig = 0; ig < basis.npw; ++ig)
        {
            const double kernel = basis.gg[ig] == 0.0
                                      ? 0.0 : ModuleBase::FOUR_PI / (tpiba * tpiba * basis.gg[ig]);
            sum[ig] *= kernel;
        }
    }
    basis.recip2real(sum.data(), component.data());
}

} // namespace

std::vector<double> periodic_negative_divergence(
    const std::vector<ModuleBase::Vector3<double>>& field,
    const ModulePW::PW_Basis& basis,
    const double tpiba)
{
    std::vector<double> result;
    std::vector<std::complex<double>> component_g;
    std::vector<std::complex<double>> sum;
    adjoint_gradient_transform(field, basis, tpiba, false, result, component_g, sum);
    return result;
}

void PeriodicCoulombOperator::apply_gradient_adjoint(
    const std::vector<ModuleBase::Vector3<double>>& field,
    std::vector<double>& result) const
{
    adjoint_gradient_transform(field, basis_, tpiba_, true,
                               result, reciprocal_aux_, reciprocal_work_);
}

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

    reciprocal_work_.resize(basis_.npw);
    reciprocal_aux_.resize(basis_.npw);
    real_work_.resize(basis_.nrxx);
    // Convert charge to potential in place; the original Fourier charge is
    // not needed after multiplying by the Coulomb kernel.
    basis_.real2recip(charge.data(), reciprocal_work_.data());
    const double tpiba2 = tpiba_ * tpiba_;
    for (int ig = 0; ig < basis_.npw; ++ig)
    {
        if (ig == basis_.ig_gge0 || basis_.gg[ig] == 0.0)
        {
            reciprocal_work_[ig] = std::complex<double>();
        }
        else
        {
            reciprocal_work_[ig] = ModuleBase::FOUR_PI * reciprocal_work_[ig] / (tpiba2 * basis_.gg[ig]);
        }
    }

    field.potential.resize(basis_.nrxx);
    field.gradient.resize(basis_.nrxx);
    basis_.recip2real(reciprocal_work_.data(), field.potential.data());
    for (int direction = 0; direction < 3; ++direction)
    {
        for (int ig = 0; ig < basis_.npw; ++ig)
        {
            reciprocal_aux_[ig] = ModuleBase::IMAG_UNIT * tpiba_ * basis_.gcar[ig][direction]
                                  * reciprocal_work_[ig];
        }
        basis_.recip2real(reciprocal_aux_.data(), real_work_.data());
        for (int ir = 0; ir < basis_.nrxx; ++ir)
        {
            field.gradient[ir][direction] = real_work_[ir];
        }
    }
}

} // namespace ModuleSccs
