#pragma once
#include "ekinetic.h"
#include "operator_force_stress_utils.hpp"
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
        this->orb_cutoff_, dmR->get_paraV(), integral_calc, force, stress);

    ModuleBase::timer::end("EKinetic", "cal_force_stress");
}

// Dummy implementations for cal_force_IJR and cal_stress_IJR
// These are not used in the simplified approach above
template <typename TK, typename TR>
void EKinetic<OperatorLCAO<TK, TR>>::cal_force_IJR(
    const int& iat1,
    const int& iat2,
    const Parallel_Orbitals* paraV,
    const std::unordered_map<int, std::vector<double>>& nlm1_all,
    const std::unordered_map<int, std::vector<double>>& nlm2_all,
    const hamilt::BaseMatrix<TR>* dmR_pointer,
    double* force1,
    double* force2)
{
    // Not used in current implementation
}

template <typename TK, typename TR>
void EKinetic<OperatorLCAO<TK, TR>>::cal_stress_IJR(
    const int& iat1,
    const int& iat2,
    const Parallel_Orbitals* paraV,
    const std::unordered_map<int, std::vector<double>>& nlm1_all,
    const std::unordered_map<int, std::vector<double>>& nlm2_all,
    const hamilt::BaseMatrix<TR>* dmR_pointer,
    const ModuleBase::Vector3<double>& dis1,
    const ModuleBase::Vector3<double>& dis2,
    double* stress)
{
    // Not used in current implementation
}

} // namespace hamilt
