#ifndef ABACUS_LR_NACV_RESPONSE_H
#define ABACUS_LR_NACV_RESPONSE_H
#include <memory>
#include <string>
namespace LR
{
template <typename T> struct GradientInputs;
/** @brief INPUT-validated CG controls, with an absolute residual tolerance in Ry. */
struct NacvResponseControls
{
    double tolerance;
    int max_iterations;
};
/**
 * @brief Reuse orbital-response operators and direct factors within one NAC geometry.
 * @details The object is owned by a single run_nac_requests call and borrows
 * its complete orbital window and kernels. GS response and LR RHS kernels keep
 * their distinct XC purposes; mutable operator density buffers are used sequentially.
 */
class NacvResponse
{
  public:
    /**
     * @param[in] inputs Borrowed converged reference; must outlive this object.
     * @param[in] channel Selected singlet/triplet or unrestricted channel.
     * @param[in] open_shell Whether amplitudes concatenate alpha and beta blocks.
     * @param[in] solver Existing NAC response solver name.
     * @param[in] controls Existing accuracy settings, unchanged.
     */
    NacvResponse(
        const GradientInputs<double> &inputs,
        int channel,
        bool open_shell,
        const std::string &solver,
        const NacvResponseControls &controls);
    ~NacvResponse();
    /**
     * @param[in,out] rhs Packed local RHS in Ry, preserved by the solver.
     * @param[out] solution Dimensionless local Z, with count rows of inputs.nloc.
     * @param[in] count Number of RHS rows.
     */
    void solve(double *rhs, double *solution, int count);
    /**
     * @param[in] amplitudes Local state-major polarized TDA vectors.
     * @param[out] solution Corresponding Z rows in the same layout.
     * @param[in] count Number of polarized vectors.
     */
    void solve_polarized(const double *amplitudes, double *solution, int count);

  private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};
}
#endif
