#include "ekinetic.h"
#include "operator_fs_utils.hpp"
#include "source_base/timer.h"

namespace hamilt
{

template <typename TK, typename TR>
void EKinetic<OperatorLCAO<TK, TR>>::cal_force_stress(const bool cal_force,
                                                          const bool cal_stress,
                                                          const HContainer<double>* dmR,
                                                          ModuleBase::matrix& force,
                                                          ModuleBase::matrix& stress)
{
    ModuleBase::TITLE("EKinetic", "cal_force_stress");
    ModuleBase::timer::start("EKinetic", "cal_force_stress");

    // Lambda function to calculate kinetic integral and its gradient
    auto integral_calc = [this](int T1, int L1, int N1, int M1,
                                 int T2, int L2, int N2, int M2,
                                 const ModuleBase::Vector3<double>& dtau,
                                 double* olm) {
        this->intor_->calculate(T1, L1, N1, M1, T2, L2, N2, M2,
                               dtau * this->ucell->lat0, &olm[0], &olm[1]);
    };

    // Use unified template with ForceSign=+1, StressSign=-1 for kinetic operator
    OperatorForceStress::cal_force_stress_2center<TK, TR, decltype(integral_calc), +1, -1>(
        cal_force, cal_stress, dmR, this->ucell, this->gridD,
        this->orb_cutoff_, this->basis_, this->layout_, integral_calc, force, stress);

    ModuleBase::timer::end("EKinetic", "cal_force_stress");
}

// explicit member function instantiations
template void EKinetic<OperatorLCAO<double, double>>::cal_force_stress(
    const bool cal_force, const bool cal_stress,
    const HContainer<double>* dmR,
    ModuleBase::matrix& force, ModuleBase::matrix& stress);
template void EKinetic<OperatorLCAO<std::complex<double>, double>>::cal_force_stress(
    const bool cal_force, const bool cal_stress,
    const HContainer<double>* dmR,
    ModuleBase::matrix& force, ModuleBase::matrix& stress);
template void EKinetic<OperatorLCAO<std::complex<double>, std::complex<double>>>::cal_force_stress(
    const bool cal_force, const bool cal_stress,
    const HContainer<double>* dmR,
    ModuleBase::matrix& force, ModuleBase::matrix& stress);

} // namespace hamilt
