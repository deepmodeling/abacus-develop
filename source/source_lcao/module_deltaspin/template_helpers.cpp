#include "lambda_loop_helper.h"
#include "spin_constrain.h"

/**
 * @file template_helpers.cpp
 * @brief Stub implementations for the TK=double (nspin=2) template specialization.
 *
 * @par Why stubs?
 * Even for nspin=2 (collinear spin), ABACUS uses complex arithmetic internally
 * for the Hamiltonian and wavefunctions. The TK=double specialization exists
 * only to satisfy the linker when the code is compiled with nspin=2.
 * All actual computation is done by the TK=std::complex<double> specialization.
 *
 * @par Design rationale
 * - The SpinConstrain template is instantiated for both TK=double and TK=complex<double>
 * - For TK=double, all methods that perform actual computation are no-ops
 * - Simple getters/setters (nspin, npol, atomCounts, etc.) still work correctly
 *   because they are not template-specialized (they use the base template)
 *
 * @par History
 * The following routines were originally member-function stubs of
 * SpinConstrain<double>. They have been lifted to free-function stubs in
 * the spinconstrain namespace alongside their complex<double> counterparts
 * (see lambda_loop_helper.cpp). Each returns the safe no-op value.
 */

/// @brief cal_mw_from_lambda stub (TK=double): no-op
template <>
void spinconstrain::SpinConstrain<double>::cal_mw_from_lambda(int i_step,
        const ModuleBase::Vector3<double>* delta_lambda)
{
}

/// @brief cal_mi_lcao stub (TK=double): no-op
template <>
void spinconstrain::SpinConstrain<double>::cal_mi_lcao(const int& step, bool print)
{
}

/// @brief run_lambda_loop stub (TK=double): no-op
template <>
void spinconstrain::SpinConstrain<double>::run_lambda_loop(int outer_step,
        bool rerun,
        std::ostream& ofs_running)
{
}

/// @brief check_rms_stop stub (TK=double): always return false (continue)
template <>
bool spinconstrain::check_rms_stop<double>(const spinconstrain::SpinConstrain<double>& sc,
                                           int outer_step,
                                           int i_step,
                                           double rms_error,
                                           double duration,
                                           double total_duration,
                                           std::ostream& ofs_running)
{
    return false;
}

/// @brief check_restriction stub (TK=double): no-op
template <>
void spinconstrain::check_restriction<double>(const spinconstrain::SpinConstrain<double>& sc,
                                              const std::vector<ModuleBase::Vector3<double>>& search,
                                              double& alpha_trial,
                                              std::ostream& ofs_running)
{
}

/// @brief cal_alpha_opt stub (TK=double): return 0.0
template <>
double spinconstrain::cal_alpha_opt<double>(const spinconstrain::SpinConstrain<double>& sc,
                                            std::vector<ModuleBase::Vector3<double>> spin,
                                            std::vector<ModuleBase::Vector3<double>> spin_plus,
                                            const double alpha_trial)
{
    return 0.0;
}

/// @brief print_termination stub (TK=double): no-op
template <>
void spinconstrain::print_termination<double>(const spinconstrain::SpinConstrain<double>& sc,
                                               std::ostream& ofs_running)
{
}

/// @brief print_header stub (TK=double): no-op
template <>
void spinconstrain::print_header<double>(const spinconstrain::SpinConstrain<double>& sc,
                                         std::ostream& ofs_running)
{
}

/// @brief check_gradient_decay stub (TK=double): always return false (no early termination)
template <>
bool spinconstrain::check_gradient_decay<double>(const spinconstrain::SpinConstrain<double>& sc,
                                                 std::vector<ModuleBase::Vector3<double>> new_spin,
                                                 std::vector<ModuleBase::Vector3<double>> old_spin,
                                                 std::vector<ModuleBase::Vector3<double>> new_delta_lambda,
                                                 std::vector<ModuleBase::Vector3<double>> old_delta_lambda,
                                                 bool print,
                                                 std::ostream& ofs_running)
{
    return false;
}

/// @brief run_lambda_linear_scan stub (TK=double): no-op
template <>
void spinconstrain::SpinConstrain<double>::run_lambda_linear_scan(int outer_step, std::ostream& ofs_running)
{
}

/// @brief reset_dspin_operator stub (TK=double): no-op
template <>
void spinconstrain::SpinConstrain<double>::reset_dspin_operator()
{
}
