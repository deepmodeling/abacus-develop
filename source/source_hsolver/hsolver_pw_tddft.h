#ifndef HSOLVER_PW_TDDFT_H
#define HSOLVER_PW_TDDFT_H

#include "source_base/matrix.h"
#include "source_basis/module_pw/pw_basis_k.h"
#include "source_hsolver/hs_operator.h"
#include "source_hsolver/hsolver_linear.h"
#include "source_psi/psi.h"

#include <iosfwd>
#include <memory>
#include <string>

namespace hsolver
{

/** @brief PW Crank-Nicolson solves and endpoint Hamiltonian expectations. */
template <typename T, typename Device>
class HSolverPWTDDFT
{
  public:
    HSolverPWTDDFT(const ModulePW::PW_Basis_K& basis,
                   const std::string& method,
                   const std::string& preconditioner,
                   const double tolerance,
                   const int max_iterations,
                   const bool kinetic_enabled,
                   const diag_comm_info& comm,
                   std::ostream& log);

    /** @brief Propagate from the fixed previous step, retaining the current iterate as the initial guess.
     *  @param dt Electronic time step in Hartree atomic units.
     *  @param momentum_shift Propagation vector potential in Hartree atomic units (inverse Bohr).
     */
    void solve(HSOperator<T, Device>& op,
               const psi::Psi<T, Device>& previous,
               psi::Psi<T, Device>* current,
               const double dt,
               const ModuleBase::Vector3<double>& momentum_shift,
               const int istep,
               const int iter,
               const bool detailed_output,
               std::ostream& log);

    /** @brief Evaluate band expectations after the caller restores the endpoint Hamiltonian. */
    void cal_band_energy(HSOperator<T, Device>& op, const psi::Psi<T, Device>& current, ModuleBase::matrix* energies);

  private:
    using Real = typename GetTypeReal<T>::type;
    const ModulePW::PW_Basis_K& basis_;
    const bool kinetic_enabled_;
    const bool kinetic_preconditioner_;
    std::unique_ptr<HSolverLinear<T, Device>> linear_solver_;
    LinearWorkspace<T, Device> band_products_;
    ct::Tensor rhs_;
    ct::Tensor hpsi_;
    ct::Tensor inverse_kinetic_;

    void prepare_buffers(const int nbands, const int nbasis);
    void update_precond(const int ik, const int dim, const T coefficient, const ModuleBase::Vector3<double>& momentum_shift);
};

} // namespace hsolver
#endif
