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

void counted_forward(const ModulePW::PW_Basis& basis,
                     const double* real,
                     std::complex<double>* reciprocal,
                     CoulombTransformCounts& counts)
{
    basis.real2recip(real, reciprocal);
    ++counts.forward_calls;
}

void counted_inverse(const ModulePW::PW_Basis& basis,
                     const std::complex<double>* reciprocal,
                     double* real,
                     CoulombTransformCounts& counts)
{
    basis.recip2real(reciprocal, real);
    ++counts.inverse_calls;
}

void adjoint_gradient_transform(
    const std::vector<ModuleBase::Vector3<double>>& field,
    const ModulePW::PW_Basis& basis,
    const double tpiba,
    std::vector<double>& component,
    std::vector<std::complex<double>>& component_g,
    std::vector<std::complex<double>>& sum,
    CoulombTransformCounts& counts)
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
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
        for (int ir = 0; ir < basis.nrxx; ++ir)
        {
            component[ir] = field[ir][direction];
        }
        counted_forward(basis, component.data(), component_g.data(), counts);
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
        for (int ig = 0; ig < basis.npw; ++ig)
        {
            sum[ig] -= ModuleBase::IMAG_UNIT * tpiba * basis.gcar[ig][direction] * component_g[ig];
        }
    }
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
    for (int ig = 0; ig < basis.npw; ++ig)
    {
        const double kernel = basis.gg[ig] == 0.0
                                  ? 0.0 : ModuleBase::FOUR_PI / (tpiba * tpiba * basis.gg[ig]);
        sum[ig] *= kernel;
    }
    counted_inverse(basis, sum.data(), component.data(), counts);
}

} // namespace

void PeriodicCoulombOperator::apply_gradient_adjoint(
    const std::vector<ModuleBase::Vector3<double>>& field,
    std::vector<double>& result) const
{
    adjoint_gradient_transform(field, basis_, tpiba_,
                               result, reciprocal_aux_, reciprocal_work_, counts_);
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
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
        for (int ig = 0; ig < basis.npw; ++ig)
        {
            gradient_g[ig]
                = ModuleBase::IMAG_UNIT * tpiba * basis.gcar[ig][direction] * values_g[ig];
        }
        basis.recip2real(gradient_g.data(), gradient_r.data());
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
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
    apply_impl(charge, &field.gradient, &field.potential);
}

void PeriodicCoulombOperator::apply_potential(
    const std::vector<double>& charge,
    std::vector<double>& potential) const
{
    apply_impl(charge, nullptr, &potential);
}

void PeriodicCoulombOperator::apply_gradient(
    const std::vector<double>& charge,
    std::vector<ModuleBase::Vector3<double>>& gradient) const
{
    apply_impl(charge, &gradient, nullptr);
}

void PeriodicCoulombOperator::apply_impl(
    const std::vector<double>& charge,
    std::vector<ModuleBase::Vector3<double>>* gradient,
    std::vector<double>* potential) const
{
    if (charge.size() != static_cast<std::size_t>(basis_.nrxx))
    {
        throw std::invalid_argument("SCCS charge array does not match the local PW real-space grid");
    }

    reciprocal_work_.resize(basis_.npw);
    // Convert charge to potential in place; the original Fourier charge is
    // not needed after multiplying by the Coulomb kernel.
    counted_forward(basis_, charge.data(), reciprocal_work_.data(), counts_);
    const double tpiba2 = tpiba_ * tpiba_;
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
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

    if (potential != nullptr)
    {
        potential->resize(basis_.nrxx);
        counted_inverse(basis_, reciprocal_work_.data(), potential->data(), counts_);
    }
    if (gradient != nullptr)
    {
        reciprocal_aux_.resize(basis_.npw);
        real_work_.resize(basis_.nrxx);
        gradient->resize(basis_.nrxx);
        for (int direction = 0; direction < 3; ++direction)
        {
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
            for (int ig = 0; ig < basis_.npw; ++ig)
            {
                reciprocal_aux_[ig] = ModuleBase::IMAG_UNIT * tpiba_ * basis_.gcar[ig][direction]
                                      * reciprocal_work_[ig];
            }
            counted_inverse(basis_, reciprocal_aux_.data(), real_work_.data(), counts_);
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
            for (int ir = 0; ir < basis_.nrxx; ++ir)
            {
                (*gradient)[ir][direction] = real_work_[ir];
            }
        }
    }
}

} // namespace ModuleSccs
