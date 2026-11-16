#ifndef ABACUS_LR_NACV_WORKFLOW_H
#define ABACUS_LR_NACV_WORKFLOW_H

#include <iosfwd>
#include <string>
#include <vector>

class Parallel_Blkcyclic;
class pseudopot_cell_vl;
class Structure_Factor;
class TwoCenterBundle;
namespace ModuleBase
{
class matrix;
}
namespace ModulePW
{
class PW_Basis;
}
namespace container
{
class Tensor;
}
namespace module_dm
{
template <typename TK, typename TR> class DensityMatrix;
}

namespace LR
{
template <typename T> struct GradientInputs;

/**
 * @brief Immutable borrowed controls for one NAC request batch.
 * @details References remain owned by resolved INPUT for the complete ionic step.
 * No configuration is read globally; empty pairs disable NAC-specific validation.
 */
struct NacvRequest
{
    const std::vector<int> &pairs;
    const std::string &spin_name;
    const std::string &response_solver;
    bool local_potential;
    bool hartree_potential;
    double response_tolerance;
    int response_max_iterations;
};

/**
 * @brief Borrowed electronic and force data for one geometry's NAC batch.
 * @details Amplitudes use the original LR distribution; gradient inputs describe
 * the complete response window. Energies and residuals are in Ry. All references,
 * including the read-only ground density, outlive run_nac_requests.
 */
template <typename T> struct NacvInputs
{
    const GradientInputs<T> &gradient;
    const std::vector<container::Tensor> &amplitudes;
    const std::vector<Parallel_Blkcyclic> &amplitude_layouts;
    const std::vector<int> &virtuals;
    int local_amplitude_size;
    const ModuleBase::matrix &energies;
    const std::vector<std::vector<double>> &residuals;
    const module_dm::DensityMatrix<T, double> &ground_density;
    const ModuleBase::matrix &occupations;
    const ModulePW::PW_Basis &density_basis;
    const ModulePW::PW_Basis &smooth_basis;
    const pseudopot_cell_vl &local_potential;
    const Structure_Factor &structure_factor;
    const TwoCenterBundle &two_center_bundle;
    bool open_shell;
    const NacvRequest &request;
};

/**
 * @brief Collect and report residuals while the local LR operator is alive.
 * @param[in] hamiltonian Current distributed LR operator.
 * @param[in] amplitudes State-major local amplitudes.
 * @param[in] energies Excitation energies in Ry.
 * @param[in] local_size Local amplitude count per root.
 * @param[in] state_count Number of solved roots.
 * @param[in] requested Whether to perform the diagnostic; false returns empty.
 * @param[in] ionic_step One-based geometry index.
 * @param[in] spin_name Existing spin-adaptation output label.
 * @param[in] rank Writer rank is zero.
 * @param[in,out] output Borrowed running log.
 * @return Owned global normalized residual norms in Ry.
 */
template <typename T, typename THamilt>
std::vector<double> collect_nac_eigenpairs(
    const THamilt &hamiltonian,
    const T *amplitudes,
    const double *energies,
    int local_size,
    int state_count,
    bool requested,
    int ionic_step,
    const std::string &spin_name,
    int rank,
    std::ostream &output);

/**
 * @brief Evaluate and log a complete real-Gamma NAC batch in input pair order.
 * @param[in] data Borrowed geometry-local inputs; numerical contractions are unchanged.
 * @param[in] ionic_step One-based index written to the existing running log.
 * @throw std::logic_error For orbitals other than real Gamma.
 * @details Each pair owns its force evaluator. Ground density is prepared once
 * per geometry and remains read-only throughout the batch.
 */
template <typename T> void run_nac_requests(const NacvInputs<T> &data, int ionic_step);
} // namespace LR

#endif
