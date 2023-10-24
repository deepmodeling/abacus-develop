#include "spin_constrain.h"

template <>
void SpinConstrain<double, psi::DEVICE_CPU>::cal_h_lambda(std::complex<double>* h_lambda,
                                                          const std::vector<std::complex<double>>& Sloc2)
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
void SpinConstrain<double, psi::DEVICE_CPU>::cal_MW(const int& step, LCAO_Matrix* LM, const UnitCell& ucell, bool print)
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
void SpinConstrain<double, psi::DEVICE_CPU>::init_sc(const UnitCell& ucell,
                                            int NPOL,
                                            std::string sc_file,
                                            Parallel_Orbitals* ParaV_in,
                                            int nspin_in,
                                            double sc_thr_in,
                                            int nsc_in,
                                            int nsc_min_in,
                                            double alpha_trial_in,
                                            double sccut_in,
                                            bool decay_grad_switch_in,
                                            K_Vectors kv_in,
                                            std::string KS_SOLVER_in,
                                            LCAO_Matrix* LM_in,
                                            hsolver::HSolver<double, psi::DEVICE_CPU>* phsol_in,
                                            hamilt::Hamilt<double, psi::DEVICE_CPU>* p_hamilt_in,
                                            psi::Psi<double>* psi_in,
                                            elecstate::ElecState* pelec_in)
{
}

template <>
void SpinConstrain<double, psi::DEVICE_CPU>::set_input_parameters(double sc_thr_in,
                                                                                int nsc_in,
                                                                                int nsc_min_in,
                                                                                double alpha_trial_in,
                                                                                double sccut_in,
                                                                                bool decay_grad_switch_in)
{
}

template <>
void SpinConstrain<double, psi::DEVICE_CPU>::set_ParaV(Parallel_Orbitals* ParaV_in)
{
}

template <>
void SpinConstrain<double, psi::DEVICE_CPU>::set_solver_parameters(
    int nspin_in,
    K_Vectors kv_in,
    hsolver::HSolver<double, psi::DEVICE_CPU>* phsol_in,
    hamilt::Hamilt<double, psi::DEVICE_CPU>* p_hamilt_in,
    psi::Psi<double>* psi_in,
    elecstate::ElecState* pelec_in,
    std::string KS_SOLVER_in,
    LCAO_Matrix* LM_in)
{
}

template <>
void SpinConstrain<double, psi::DEVICE_CPU>::bcast_ScData(const UnitCell& ucell, int NPOL, std::string sc_file)
{
}