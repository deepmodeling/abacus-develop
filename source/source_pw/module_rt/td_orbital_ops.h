#ifndef HSOLVER_TD_ORBITAL_OPS_H
#define HSOLVER_TD_ORBITAL_OPS_H
#include "source_base/matrix.h"
#include "source_hsolver/hs_operator.h"
#include "source_hsolver/module_linalg/linear_workspace.h"
#include "source_hsolver/module_linalg/orthonormal.h"
#include "source_psi/psi.h"

#include <iosfwd>

namespace hsolver
{
/** @brief Independent orthogonality maxima over accepted trial iterations and local k points. */
struct TDOrthStats
{
    double before = 0.0;
    double after = 0.0;
};

template <typename T, typename Device>
class TDOrbOps
{
  private:
    const diag_comm_info comm_;
    const int npol_;
    const OrthMethod orth_method_;
    const bool out_stat_;
    const std::vector<int> global_k_;
    ct::Tensor clean_;
    ct::Tensor hpsi_;
    LinearAlgebra<T, Device> algebra_;
    Orthonormal<T, Device> orthonormal_;
    TDOrthStats orth_stats_;
    std::vector<std::vector<double>> orth_norms_;
    std::ostream& log_;
    unsigned int warned_events_ = 0;
    LinearWorkspace<T, Device> band_products_;
    void prepare_buffers(int nbands, int nbasis);
    void correct_orbitals(T* current, int ld, int dim, int bands, int ik, int istep, int iter, bool detailed_output);
    void record_orth(const OrthResult& result, int ik, int istep, int iter);
    void report_orth_warning(const OrthResult& result, int ik, int istep, int iter);

  public:
    /** @brief Configure common orbital processing independently of a propagator. */
    TDOrbOps(int npol,
             OrthMethod orth_method,
             bool out_stat,
             const std::vector<int>& global_k,
             const diag_comm_info& comm,
             std::ostream& log);
    /** @brief Apply common postprocessing to a completed predictor or corrector. */
    void correct(psi::Psi<T, Device>* current, int step, int iter, bool detailed_output);
    /** @brief Reset diagnostics without invalidating preconditioner history. */
    void reset_orth_stats();
    /** @brief Publish diagnostics accumulated over accepted trial intervals only. */
    void set_orth_stats(const TDOrthStats& stats)
    {
        orth_stats_ = stats;
    }
    /** @brief Clear spinor padding without changing physical coefficients or the k-point cursor; no communication. */
    void zero_padding(psi::Psi<T, Device>* current);
    /** @brief Orthogonalize initial orbitals with clean padding; return whether this pool changed any physical orbital. */
    bool orth_initial(psi::Psi<T, Device>* current, int iter);
    /** @brief Check the converged initial state once; finite orthogonality errors are diagnostic only. */
    void check_initial(const psi::Psi<T, Device>& current, int iter);
    /** @brief Occupation-weighted electron count for the local k-point pool. */
    double wave_electrons(const ModuleBase::matrix& occupations) const;
    const TDOrthStats& orth_stats() const
    {
        return orth_stats_;
    }

    /** @brief Evaluate band expectations after the caller restores the endpoint Hamiltonian. */
    void cal_band_energy(HSOperator<T, Device>& op, const psi::Psi<T, Device>& current, ModuleBase::matrix* energies);

    /** @brief Return the pool-reduced, column-major Psi-adjoint H Psi for one local k point. */
    std::vector<std::complex<double>> cal_eij(HSOperator<T, Device>& op, const psi::Psi<T, Device>& current, const int ik);
};
} // namespace hsolver
#endif
