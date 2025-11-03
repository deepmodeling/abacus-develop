#ifndef ESOLVER_KS_LCAO_TDDFT_H
#define ESOLVER_KS_LCAO_TDDFT_H
#include "esolver_ks.h"
#include "esolver_ks_lcao.h"
#include "source_base/module_container/ATen/core/tensor.h" // ct::Tensor
#include "source_lcao/module_rt/gather_mat.h"              // MPI gathering and distributing functions
#include "source_lcao/module_rt/td_info.h"
#include "source_lcao/module_rt/velocity_op.h"

namespace ModuleESolver
{

template <typename TR, typename Device = base_device::DEVICE_CPU>
class ESolver_KS_LCAO_TDDFT : public ESolver_KS_LCAO<std::complex<double>, TR>
{
  public:
    ESolver_KS_LCAO_TDDFT();

    ~ESolver_KS_LCAO_TDDFT();

    void before_all_runners(UnitCell& ucell, const Input_para& inp) override;

  protected:
    virtual void runner(UnitCell& cell, const int istep) override;

    virtual void hamilt2rho_single(UnitCell& ucell, const int istep, const int iter, const double ethr) override;

    void store_h_s_psi(UnitCell& ucell, const int istep, const int iter, const bool conv_esolver);

    void iter_finish(UnitCell& ucell,
                     const int istep,
                     const int estep,
                     const int estep_max,
                     int& iter,
                     bool& conv_esolver);

    virtual void after_scf(UnitCell& ucell, const int istep, const bool conv_esolver) override;

    void print_step();

    //! Wave function for all k-points of last time step
    psi::Psi<std::complex<double>>* psi_laststep = nullptr;

    //! Hamiltonian for all k-points of last time step
    ct::Tensor Hk_laststep = ct::Tensor(ct::DataType::DT_COMPLEX_DOUBLE);

    //! Overlap matrix for all k-points of last time step
    ct::Tensor Sk_laststep = ct::Tensor(ct::DataType::DT_COMPLEX_DOUBLE);

    //! Control heterogeneous computing of the TDDFT solver
    bool use_tensor = false;
    bool use_lapack = false;

    // Control the device type for Hk_laststep and Sk_laststep
    // Set to CPU temporarily, should wait for further GPU development
    static constexpr ct::DeviceType ct_device_type_hs = ct::DeviceType::CpuDevice;

    //! Total steps for evolving the wave function
    int totstep = -1;

    //! Velocity matrix for calculating current
    Velocity_op<TR>* velocity_mat = nullptr;

    TD_info* td_p = nullptr;

    //! Restart flag
    bool restart_done = false;

  private:
    void weight_dm_rho(const UnitCell& ucell);
};

} // namespace ModuleESolver
#endif // ESOLVER_KS_LCAO_TDDFT_H
