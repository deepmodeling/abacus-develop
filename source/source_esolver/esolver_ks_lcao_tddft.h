#ifndef ESOLVER_KS_LCAO_TDDFT_H
#define ESOLVER_KS_LCAO_TDDFT_H
#include "esolver_ks.h"
#include "esolver_ks_lcao.h"
#include "source_base/module_container/ATen/core/tensor.h" // ct::Tensor
#include "source_lcao/module_rt/td_info.h"
#include "source_lcao/module_rt/td_moving_gauge.h"
#include "source_lcao/module_rt/velocity_op.h"

#include <memory>

namespace elecstate
{
//! Shared time-dependent field state for propagation and output.
class TDFieldManager;
} // namespace elecstate

namespace ModuleESolver
{

//! Real-time TDDFT energy solver using the LCAO basis.
template <typename TR, typename Device = base_device::DEVICE_CPU>
class ESolver_KS_LCAO_TDDFT : public ESolver_KS_LCAO<std::complex<double>, TR>
{
  public:
    //! Initialize the solver identity and basis name.
    ESolver_KS_LCAO_TDDFT();

    //! Release the saved wave functions, Hamiltonian, and solver-owned TD contexts.
    ~ESolver_KS_LCAO_TDDFT();

    //! Initialize the LCAO solver, shared field state, and optional restart wave functions.
    void before_all_runners(BaseCell& basecell, const Input_para& inp) override;

  protected:
    //! Return the solver-owned RT state used by LCAO operators and output.
    TD_info* get_td_info() const override { return td_p; }

    //! Run the electronic steps belonging to one MD call, including initialization or restart replay.
    virtual void runner(BaseCell& basecell, const int istep) override;

    //! Obtain the density from initial-state diagonalization, loaded wave functions, or time propagation.
    virtual void hamilt2rho_single(UnitCell& ucell, const int istep, const int iter, const double ethr) override;

    //! Save H, S, and wave functions when requested at convergence or the maximum iteration count.
    void store_h_s_psi(UnitCell& ucell, const int istep, const int iter, const bool conv_esolver);

    //! Finish an electronic iteration and update the saved state and requested energy-density matrix.
    void iter_finish(UnitCell& ucell, const int istep, const int estep, const int estep_max, int& iter, bool& conv_esolver);

    //! Write the electronic-step energy and TD observables after the SCF loop.
    virtual void after_scf(UnitCell& ucell, const int istep, const bool conv_esolver) override;

    //! Print the electronic-step header using a 1-based step number.
    void print_step();

    //! Wave functions for all local k-points of the last electronic step.
    psi::Psi<std::complex<double>>* psi_laststep = nullptr;

    //! Hamiltonian for all local k-points of the last electronic step.
    ct::Tensor Hk_laststep = ct::Tensor(ct::DataType::DT_COMPLEX_DOUBLE);

    //! Overlap matrix for all local k-points of the last electronic step.
    ct::Tensor Sk_laststep = ct::Tensor(ct::DataType::DT_COMPLEX_DOUBLE);

    //! Control tensor-based heterogeneous computing of the TDDFT solver.
    bool use_tensor = false;
    //! Use gathered full matrices rather than distributed matrices in tensor propagation.
    bool use_lapack = false;

    //! Store Hk_laststep and Sk_laststep on the CPU, including for GPU propagation.
    static constexpr ct::DeviceType ct_device_type_hs = ct::DeviceType::CpuDevice;

    //! Current 0-based electronic-step index, with -1 before the initial state.
    int totstep = -1;

    //! Velocity matrix for calculating current
    Velocity_op<TR>* velocity_mat = nullptr;

    //! Solver-owned TD context for field-dependent LCAO operators and observables.
    TD_info* td_p = nullptr;

    //! Moving spatial gauge for Ehrenfest dynamics, to calculate the correction term arising from the movement of basis
    bool use_td_moving_gauge = false;
    //! Moving-basis correction matrices, allocated when the moving spatial gauge is enabled.
    module_rt::TD_MovingGauge* td_mg_ = nullptr;

    //! Whether the saved restart endpoint has been replayed to rebuild the electronic state.
    bool restart_done = false;

    //! Field manager shared by propagation, potentials, and restart output.
    std::shared_ptr<elecstate::TDFieldManager> td_field_manager_;

  private:
    //! Fixed electronic-step origin loaded at restart, independent of field advancement.
    int restart_step_ = 0;

    //! Apply specified occupations and rebuild band energy, density matrices, and charge density.
    void weight_dm_rho(const UnitCell& ucell);
};

} // namespace ModuleESolver
#endif // ESOLVER_KS_LCAO_TDDFT_H
