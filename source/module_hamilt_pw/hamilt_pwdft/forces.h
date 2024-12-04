#ifndef FORCES_H
#define FORCES_H

#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_base/matrix.h"
#include "module_base/module_device/memory_op.h"
#include "module_basis/module_pw/pw_basis_k.h"
#include "module_cell/klist.h"
#include "module_cell/module_symmetry/symmetry.h"
#include "module_elecstate/elecstate.h"
#include "module_hamilt_pw/hamilt_pwdft/kernels/force_op.h"
#include "module_hsolver/kernels/math_kernel_op.h"
#include "module_psi/psi.h"
#include "structure_factor.h"

template <typename FPTYPE, typename Device = base_device::DEVICE_CPU>
class Forces
{
  public:
    template <typename T>
    friend class Force_Stress_LCAO;
    /* This routine is a driver routine which compute the forces
     * acting on the atoms, the complete forces in plane waves
     * is computed from 4 main parts
     * (1) cal_force_loc: contribution due to local potential.
     * (2) cal_foce_ew: contribution due to ewald potential.
     * (3) cal_force_cc: contributino due to NLCC.
     * (4) cal_nl: contribution due to the non-local pseudopotential.
     * (5) cal_force_us: contribution due to US pseudopotential.
     * (6) cal_scc: contributino due to incomplete SCF calculation.
     */
    Forces(const int nat_in) : nat(nat_in){};
    ~Forces(){};

    void cal_force(const UnitCell& ucell,
                   ModuleBase::matrix& force,
                   const elecstate::ElecState& elec,
                   ModulePW::PW_Basis* rho_basis,
                   ModuleSymmetry::Symmetry* p_symm,
                   Structure_Factor* p_sf,
                   const pseudopot_cell_vl* locpp,
                   const pseudopot_cell_vnl* nlpp = nullptr,
                   K_Vectors* pkv = nullptr,
                   ModulePW::PW_Basis_K* psi_basis = nullptr,
                   const psi::Psi<std::complex<FPTYPE>, Device>* psi_in = nullptr);

  protected:
    int nat = 0;
    int npwx = 0;

    void cal_force_loc(const UnitCell& ucell,
                       ModuleBase::matrix& forcelc,
                       ModulePW::PW_Basis* rho_basis,
                       const ModuleBase::matrix& vloc,
                       const Charge* const chr);
    void cal_force_ew(const UnitCell& ucell,
                      ModuleBase::matrix& forceion, 
                      ModulePW::PW_Basis* rho_basis, 
                      const Structure_Factor* p_sf);
    void cal_force_cc(ModuleBase::matrix& forcecc,
                      ModulePW::PW_Basis* rho_basis,
                      const Charge* const chr,
                      const bool* numeric,
                      UnitCell& ucell_in);
    /**
     * @brief This routine computes the atomic force of non-local pseudopotential
     *    F^{NL}_i = \sum_{n,k}f_{nk}\sum_I \sum_{lm,l'm'}D_{l,l'}^{I} [
     *               \sum_G \langle c_{nk}(\mathbf{G+K})|\beta_{lm}^I(\mathbf{G+K})\rangle *
     *               \sum_{G'}\langle \beta_{lm}^I(\mathbf{G+K})*(-j)^l(\mathbf{G+K})_i |c_{nk}(\mathbf{G+K})\rangle ]
     *    there would be three parts in the above equation:
     *    (1) sum over becp and dbecp with D_{l,l'}^{I} ----- first line in the above equation
     *    (2) calculate becp = <psi | beta> ----- second line in the above equation
     *    (3) calculate dbecp = <psi | \nabla beta> ----- third line in the above equation
     */
    void cal_force_nl(ModuleBase::matrix& forcenl,
                      const ModuleBase::matrix& wg,
                      const ModuleBase::matrix& ekb,
                      const K_Vectors* p_kv,
                      const ModulePW::PW_Basis_K* psi_basis,
                      const Structure_Factor* p_sf,
                      const pseudopot_cell_vnl& nlpp_in,
                      const UnitCell& ucell_in,
                      const psi::Psi<std::complex<FPTYPE>, Device>* psi_in = nullptr);
    void cal_force_scc(ModuleBase::matrix& forcescc,
                       ModulePW::PW_Basis* rho_basis,
                       const ModuleBase::matrix& v_current,
                       const bool vnew_exist,
                       const bool* numeric,
                       const UnitCell& ucell_in);
    void cal_force_us(ModuleBase::matrix& forcenl,
                      ModulePW::PW_Basis* rho_basis,
                      const pseudopot_cell_vnl& ppcell_in,
                      const elecstate::ElecState& elec,
                      const UnitCell& ucell);
    void cal_ylm(int lmax, int npw, const FPTYPE* gk_in, FPTYPE* ylm);
    void deriv_drhoc(const bool& numeric,
                     const int mesh,
                     const FPTYPE* r,
                     const FPTYPE* rab,
                     const FPTYPE* rhoc,
                     FPTYPE* drhocg,
                     ModulePW::PW_Basis* rho_basis,
                     int type,
                     const UnitCell& ucell_in); // used in nonlinear core correction stress
    void deriv_drhoc_scc(const bool& numeric,
                     const int mesh,
                     const FPTYPE* r,
                     const FPTYPE* rab,
                     const FPTYPE* rhoc,
                     FPTYPE* drhocg,
                     ModulePW::PW_Basis* rho_basis,
                     const UnitCell& ucell_in); // used in nonlinear core correction stress
  protected:
    Device* ctx = {};
    base_device::DEVICE_CPU* cpu_ctx = {};
    base_device::AbacusDevice_t device = {};
  private:
    using gemm_op = hsolver::gemm_op<std::complex<FPTYPE>, Device>;

    using resmem_complex_op = base_device::memory::resize_memory_op<std::complex<FPTYPE>, Device>;
    using resmem_complex_h_op = base_device::memory::resize_memory_op<std::complex<FPTYPE>, base_device::DEVICE_CPU>;
    using delmem_complex_op = base_device::memory::delete_memory_op<std::complex<FPTYPE>, Device>;
    using delmem_complex_h_op = base_device::memory::delete_memory_op<std::complex<FPTYPE>, base_device::DEVICE_CPU>;
    using syncmem_complex_h2d_op
        = base_device::memory::synchronize_memory_op<std::complex<FPTYPE>, Device, base_device::DEVICE_CPU>;
    using syncmem_complex_d2h_op
        = base_device::memory::synchronize_memory_op<std::complex<FPTYPE>, base_device::DEVICE_CPU, Device>;

    using resmem_var_op = base_device::memory::resize_memory_op<FPTYPE, Device>;
    using delmem_var_op = base_device::memory::delete_memory_op<FPTYPE, Device>;
    using syncmem_var_h2d_op = base_device::memory::synchronize_memory_op<FPTYPE, Device, base_device::DEVICE_CPU>;
    using syncmem_var_d2h_op = base_device::memory::synchronize_memory_op<FPTYPE, base_device::DEVICE_CPU, Device>;

    using resmem_int_op = base_device::memory::resize_memory_op<int, Device>;
    using delmem_int_op = base_device::memory::delete_memory_op<int, Device>;
    using syncmem_int_h2d_op = base_device::memory::synchronize_memory_op<int, Device, base_device::DEVICE_CPU>;
};

#endif
