#ifndef FS_KIN_TOOLS_H
#define FS_KIN_TOOLS_H
#include "module_base/module_device/device.h"
#include "module_basis/module_pw/pw_basis_k.h"
#include "module_cell/klist.h"
#include "module_cell/unitcell.h"
#include "module_cell/module_symmetry/symmetry.h"
#include "module_hamilt_pw/hamilt_pwdft/kernels/stress_op.h"

#include <complex>
namespace hamilt
{
template <typename FPTYPE, typename Device>
class FS_Kin_tools
{
  public:
    FS_Kin_tools(const UnitCell& ucell_in,
                 const K_Vectors* kv_in,
                 const ModulePW::PW_Basis_K* wfc_basis_in,
                 const ModuleBase::matrix& wg);
    ~FS_Kin_tools();

    /**
     * @brief calculate G+k and store it in gk and also calculate kfac
     */
    void cal_gk(const int& ik);

    /**
     * @brief calculate stress tensor for kinetic energy
     *        stress = \sum_{G,k,i}  wk(k) * gk_l(G) * gk_m(G) * d_kfac(G) * occ_i*|ppsi_i(G)|^2
     * 
     * @param ik k-point index
     * @param npm number of bands
     * @param occ if use the occupation of the bands
     * @param psi wavefunctions
     */
    void cal_stress_kin(const int& ik, const int& npm, const bool& occ, const std::complex<FPTYPE>* psi);

    /**
     * @brief symmetrize the stress tensor
     */
    void symmetrize_stress(ModuleSymmetry::Symmetry* p_symm, ModuleBase::matrix& sigma);

  protected:
    Device* ctx = {};
    base_device::DEVICE_CPU* cpu_ctx = {};
    base_device::AbacusDevice_t device = {};
    std::vector<FPTYPE> gk3_;
    std::vector<FPTYPE*> gk;
    std::vector<FPTYPE> kfac;
    std::vector<FPTYPE> s_kin;
    FPTYPE* d_gk = nullptr;
    FPTYPE* d_kfac = nullptr;
    const FPTYPE* wg = nullptr;
    const FPTYPE* wk = nullptr;
    const ModulePW::PW_Basis_K* wfc_basis_ = nullptr;
    const UnitCell& ucell_;
    const int nksbands_;


  private:
    using resmem_var_op = base_device::memory::resize_memory_op<FPTYPE, Device>;
    using setmem_var_op = base_device::memory::set_memory_op<FPTYPE, Device>;
    using delmem_var_op = base_device::memory::delete_memory_op<FPTYPE, Device>;
    using syncmem_var_h2d_op = base_device::memory::synchronize_memory_op<FPTYPE, Device, base_device::DEVICE_CPU>;
    using syncmem_var_d2h_op = base_device::memory::synchronize_memory_op<FPTYPE, base_device::DEVICE_CPU, Device>;
    using cal_multi_dot_op = hamilt::cal_multi_dot_op<FPTYPE, Device>;
};

} // namespace hamilt

#endif