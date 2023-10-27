#include "spin_constrain.h"

template <>
void SpinConstrain<double, psi::DEVICE_CPU>::cal_h_lambda(std::complex<double>* h_lambda,
                                                          const std::vector<std::complex<double>>& Sloc2,
                                                          bool column_major)
{
}

template <>
void SpinConstrain<double, psi::DEVICE_CPU>::cal_mw_from_lambda(int i_step)
{
}

template <>
ModuleBase::matrix SpinConstrain<double, psi::DEVICE_CPU>::cal_MW_k(
    LCAO_Matrix* LM,
    const std::vector<std::vector<std::complex<double>>>& dm)
{
    ModuleBase::matrix orbMulP;
    return orbMulP;
}

template <>
void SpinConstrain<double, psi::DEVICE_CPU>::cal_MW(const int& step, LCAO_Matrix* LM, bool print)
{
}

template <>
void SpinConstrain<double, psi::DEVICE_CPU>::calculate_MW(const std::vector<std::vector<std::vector<double>>>& AorbMulP)
{
}

template<>
std::vector<std::vector<std::vector<double>>> SpinConstrain<double, psi::DEVICE_CPU>::convert(const ModuleBase::matrix &orbMulP)
{
    std::vector<std::vector<std::vector<double>>> AorbMulP;
    return AorbMulP;
}

template <>
void SpinConstrain<double, psi::DEVICE_CPU>::run_lambda_loop(int outer_step)
{
}

template <>
bool SpinConstrain<double, psi::DEVICE_CPU>::check_rms_stop(int outer_step, int i_step, double rms_error)
{
    return false;
}

template <>
void SpinConstrain<double, psi::DEVICE_CPU>::print_termination()
{
}

template <>
void SpinConstrain<double, psi::DEVICE_CPU>::collect_MW(ModuleBase::matrix& MecMulP,
                                                        const ModuleBase::ComplexMatrix& mud,
                                                        int nw)
{
}