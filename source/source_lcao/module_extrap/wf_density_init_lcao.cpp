#include "source_lcao/module_extrap/wf_history_lcao.h"

#include "source_base/global_function.h"
#include "source_base/global_variable.h"
#include "source_base/timer.h"
#include "source_base/tool_quit.h"
#include "source_estate/module_dm/cal_dm_psi.h"
#include "source_lcao/hamilt_lcao.h"
#include "source_lcao/module_operator_lcao/operator_lcao.h"
#include "source_lcao/rho_tau_lcao.h"

#include <complex>
#include <iomanip>
#include <sstream>

namespace ModuleExtrap
{

namespace
{

std::string wfc_extrapolation_failure_message(const WfExtrapApplyResult& result)
{
    std::ostringstream oss;
    oss << std::scientific << std::setprecision(16);
    oss << "WFN extrapolation failed with status: " << to_string(result.status) << "\n\n"
        << " snapshot_istep=" << result.snapshot_istep << "\n"
        << " failed_state=" << result.failed_state << "\n"
        << " failed_pivot_index=" << result.failed_pivot_index << "\n"
        << " failed_pivot=" << result.failed_pivot << "\n"
        << " nstate=" << result.nstate << "\n"
        << " nbands=" << result.nbands << "\n"
        << " nbasis=" << result.nbasis << "\n"
        << " min_metric_diag=" << result.min_metric_diag << "\n"
        << " max_metric_diag=" << result.max_metric_diag << "\n"
        << " max_metric_abs=" << result.max_metric_abs << "\n"
        << " max_metric_asymmetry=" << result.max_metric_asymmetry << "\n"
        << " max_orthonormality_deviation=" << result.max_orthonormality_deviation << "\n";
    return oss.str();
}

} // namespace

template <typename TK, typename TR>
bool initialize_gamma_density_from_history(WfHistoryLCAO<TK>&,
                                           hamilt::HamiltLCAO<TK, TR>&,
                                           const Parallel_Orbitals&,
                                           psi::Psi<TK>&,
                                           const ModuleBase::matrix&,
                                           elecstate::DensityMatrix<TK, double>&,
                                           Charge&,
                                           const int,
                                           const std::string&)
{
    ModuleBase::WARNING_QUIT("initialize_gamma_density_from_history",
                             "WFN extrapolation is currently supported only for the real Gamma-only NAO path.");
    return false;
}

template <>
bool initialize_gamma_density_from_history<double, double>(WfHistoryLCAO<double>& history,
                                                            hamilt::HamiltLCAO<double, double>& hamiltonian,
                                                            const Parallel_Orbitals& pv,
                                                            psi::Psi<double>& psi,
                                                            const ModuleBase::matrix& wg_now,
                                                            elecstate::DensityMatrix<double, double>& dmat,
                                                            Charge& charge,
                                                            const int nspin,
                                                            const std::string& ks_solver)
{
    auto* lcao_op = dynamic_cast<hamilt::OperatorLCAO<double, double>*>(hamiltonian.getOperator());
    if (lcao_op == nullptr)
    {
        ModuleBase::WARNING_QUIT("initialize_gamma_density_from_history",
                                 "Failed to access the LCAO operator chain for WFN extrapolation.");
    }

    ModuleBase::timer::start("WFN_Extrap", "prepare_overlap");
    lcao_op->contributeHR();
    const int sk_layout = ModuleBase::GlobalFunc::IS_COLUMN_MAJOR_KS_SOLVER(ks_solver) ? 1 : 0;
    hamiltonian.updateSk(0, sk_layout);
    ModuleBase::timer::end("WFN_Extrap", "prepare_overlap");

    ModuleBase::timer::start("WFN_Extrap", "apply");
    const WfExtrapApplyResult result = history.try_use_prev_wf_gamma(hamiltonian.getSk(), pv, psi, wg_now);
    ModuleBase::timer::end("WFN_Extrap", "apply");
    if (!result.ok())
    {
        ModuleBase::WARNING_QUIT("initialize_gamma_density_from_history",
                                 wfc_extrapolation_failure_message(result));
    }

    ModuleBase::timer::start("WFN_Extrap", "rebuild_density");
    elecstate::cal_dm_psi(dmat.get_paraV_pointer(), wg_now, psi, dmat);
    dmat.cal_DMR();
    LCAO_domain::dm2rho(dmat.get_DMR_vector(), nspin, &charge);
    ModuleBase::timer::end("WFN_Extrap", "rebuild_density");

    GlobalV::ofs_running << " WFN extrapolation: use_prev_wf from ionic step " << result.snapshot_istep
                         << ", active bands = " << result.nactive_bands
                         << ", max |C^T S C - I| = " << result.max_orthonormality_deviation << std::endl;
    return true;
}

template bool initialize_gamma_density_from_history<std::complex<double>, double>(
    WfHistoryLCAO<std::complex<double>>&, hamilt::HamiltLCAO<std::complex<double>, double>&,
    const Parallel_Orbitals&, psi::Psi<std::complex<double>>&, const ModuleBase::matrix&,
    elecstate::DensityMatrix<std::complex<double>, double>&, Charge&, int, const std::string&);
template bool initialize_gamma_density_from_history<std::complex<double>, std::complex<double>>(
    WfHistoryLCAO<std::complex<double>>&, hamilt::HamiltLCAO<std::complex<double>, std::complex<double>>&,
    const Parallel_Orbitals&, psi::Psi<std::complex<double>>&, const ModuleBase::matrix&,
    elecstate::DensityMatrix<std::complex<double>, double>&, Charge&, int, const std::string&);

} // namespace ModuleExtrap
