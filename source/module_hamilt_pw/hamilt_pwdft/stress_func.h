#ifndef STRESS_FUNC_H
#define STRESS_FUNC_H

#include "module_base/complexmatrix.h"
#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_base/matrix.h"
#include "module_base/parallel_reduce.h"
#include "module_base/realarray.h"
#include "module_base/vector3.h"
#include "module_basis/module_pw/pw_basis_k.h"
#include "module_cell/klist.h"
#include "module_elecstate/module_charge/charge.h"
#include "module_hamilt_pw/hamilt_pwdft/kernels/stress_op.h"
#include "module_hamilt_pw/hamilt_pwdft/structure_factor.h"
#include "module_hsolver/kernels/math_kernel_op.h"
#include "module_psi/psi.h"

//-------------------------------------------------------------------
// mohan reconstruction note: 2021-02-07
// the stress code needs reconstructions (by Daye Zheng)
// 1) add explanations for each function, each variable, for
// main procedures, make the code readable
// 2) divide the stress class into several files, each file
// deals with only one part of the stress, it is convenient for
// the next-step reconstruction, for example, we want to make
// pw as an external variable instead of a global variable
// 3) for PW and LCAO, keeps only one copy of the code, for example
// the ewald term needs only one copy, it will reduce the
// time to maintain both copies of stress codes.
// 4) remain openning interfaces for others to contribute stress
// codes, for example, molecular dynamics will have an ionic stress
// term, +U? exx? may have other stress terms.
// 5) delete useless comments and tests, if you have a useless code,
// please explicitly explain why you want to keep the test
// 6) format should be beautiful! code should be readable like a
// note (let readers be comfortable)
//-------------------------------------------------------------------

//----------------------------------------------------------------
// compute the stress terms in terms of the plane wave basis set
// the stress terms include:
// 1) the stress from the electron kinetic energy
// 2) the stress from the local pseudopotentials
// 3) the stress from the non-local pseudopotentials
// 4) the stress from the Hartree term
// 5) the stress from the non-linear core correction (if any)
// 6) the strees from the exchange-correlation functional term
// 7) the stress from the ewald term (ion-ion intraction under
//		periodic boundary conditions).
// 8) the stress from ionic contributions (for molecular dynamics)
//----------------------------------------------------------------

template <typename FPTYPE, typename Device = psi::DEVICE_CPU>
class Stress_Func
{
  public:
    Stress_Func(){};
    ~Stress_Func(){};

    // stress functions
    //  1) the stress from the electron kinetic energy
    void stress_kin(ModuleBase::matrix& sigma,
                    const ModuleBase::matrix& wg,
                    ModuleSymmetry::Symmetry* p_symm,
                    K_Vectors* p_kv,
                    ModulePW::PW_Basis_K* wfc_basis,
                    const psi::Psi<complex<FPTYPE>>* psi_in = nullptr); // electron kinetic part in PW basis

    // 2) the stress from the Hartree term
    void stress_har(ModuleBase::matrix& sigma,
                    ModulePW::PW_Basis* rho_basis,
                    const bool is_pw,
                    const Charge* const chr); // hartree part in PW or LCAO basis

    // 3) the stress from the ewald term (ion-ion intraction under
    //		periodic boundary conditions).
    void stress_ewa(ModuleBase::matrix& sigma,
                    ModulePW::PW_Basis* rho_basis,
                    const bool is_pw); // ewald part in PW or LCAO basis

    // 4) the stress from the local pseudopotentials
    void stress_loc(ModuleBase::matrix& sigma,
                    ModulePW::PW_Basis* rho_basis,
                    const Structure_Factor* p_sf,
                    const bool is_pw,
                    const Charge* const chr); // local pseudopotential part in PW or LCAO

    void dvloc_of_g(const int& msh,
                    const FPTYPE* rab,
                    const FPTYPE* r,
                    const FPTYPE* vloc_at,
                    const FPTYPE& zp,
                    FPTYPE* dvloc,
                    ModulePW::PW_Basis* rho_basis); // used in local pseudopotential stress

    void dvloc_coul(const FPTYPE& zp,
                    FPTYPE* dvloc,
                    ModulePW::PW_Basis* rho_basis); // used in local pseudopotential stress

    // 5) the stress from the non-linear core correction (if any)
    void stress_cc(ModuleBase::matrix& sigma,
                   ModulePW::PW_Basis* rho_basis,
                   const Structure_Factor* p_sf,
                   const bool is_pw,
                   const Charge* const chr); // nonlinear core correction stress in PW or LCAO basis

    void deriv_drhoc(const bool& numeric,
                     const int mesh,
                     const FPTYPE* r,
                     const FPTYPE* rab,
                     const FPTYPE* rhoc,
                     FPTYPE* drhocg,
                     ModulePW::PW_Basis* rho_basis); // used in nonlinear core correction stress

    // 6) the stress from the exchange-correlation functional term
    void stress_gga(ModuleBase::matrix& sigma,
                    ModulePW::PW_Basis* rho_basis,
                    const Charge* const chr); // gga part in both PW and LCAO basis
    void stress_mgga(ModuleBase::matrix& sigma,
                     const ModuleBase::matrix& wg,
                     const ModuleBase::matrix& v_ofk,
                     const Charge* const chr,
                     K_Vectors* p_kv,
                     ModulePW::PW_Basis_K* wfc_basis,
                     const psi::Psi<complex<FPTYPE>>* psi_in = nullptr); // gga part in PW basis

    // 7) the stress from the non-local pseudopotentials
    void stress_nl(ModuleBase::matrix& sigma,
                   const ModuleBase::matrix& wg,
                   Structure_Factor* p_sf,
                   K_Vectors* p_kv,
                   ModuleSymmetry::Symmetry* p_symm,
                   ModulePW::PW_Basis_K* wfc_basis,
                   const psi::Psi<complex<FPTYPE>, Device>* psi_in); // nonlocal part in PW basis

    void get_dvnl1(ModuleBase::ComplexMatrix& vkb,
                   const int ik,
                   const int ipol,
                   Structure_Factor* p_sf,
                   ModulePW::PW_Basis_K* wfc_basis); // used in nonlocal part in PW basis
    void dylmr2(const int nylm,
                const int ngy,
                ModuleBase::Vector3<FPTYPE>* gk,
                ModuleBase::matrix& dylm,
                const int ipol); // used in get_dvnl1()
    void get_dvnl2(ModuleBase::ComplexMatrix& vkb,
                   const int ik,
                   Structure_Factor* p_sf,
                   ModulePW::PW_Basis_K* wfc_basis); // used in nonlocal part in PW basis
    FPTYPE Polynomial_Interpolation_nl(const ModuleBase::realArray& table,
                                       const int& dim1,
                                       const int& dim2,
                                       const FPTYPE& table_interval,
                                       const FPTYPE& x); // used in get_dvnl2()

    // functions for stress print
    void print_stress(const std::string& name, const ModuleBase::matrix& f, const bool screen, bool ry) const;

    void printstress_total(const ModuleBase::matrix& scs, bool ry);

    static FPTYPE stress_invalid_threshold_ev;
    static FPTYPE output_acc;

  private:
    Device* ctx = {};
    psi::DEVICE_CPU* cpu_ctx = {};
    psi::AbacusDevice_t device = {};
    using gemm_op = hsolver::gemm_op<FPTYPE, Device>;
    using cal_stress_nl_op = hamilt::cal_stress_nl_op<FPTYPE, Device>;
    using cal_dbecp_noevc_nl_op = hamilt::cal_dbecp_noevc_nl_op<FPTYPE, Device>;

    using resmem_complex_op = psi::memory::resize_memory_op<std::complex<FPTYPE>, Device>;
    using resmem_complex_h_op = psi::memory::resize_memory_op<std::complex<FPTYPE>, psi::DEVICE_CPU>;
    using setmem_complex_op = psi::memory::set_memory_op<std::complex<FPTYPE>, Device>;
    using delmem_complex_op = psi::memory::delete_memory_op<std::complex<FPTYPE>, Device>;
    using delmem_complex_h_op = psi::memory::delete_memory_op<std::complex<FPTYPE>, psi::DEVICE_CPU>;
    using syncmem_complex_h2d_op = psi::memory::synchronize_memory_op<std::complex<FPTYPE>, Device, psi::DEVICE_CPU>;
    using syncmem_complex_d2h_op = psi::memory::synchronize_memory_op<std::complex<FPTYPE>, psi::DEVICE_CPU, Device>;

    using resmem_var_op = psi::memory::resize_memory_op<FPTYPE, Device>;
    using resmem_var_h_op = psi::memory::resize_memory_op<FPTYPE, psi::DEVICE_CPU>;
    using setmem_var_op = psi::memory::set_memory_op<FPTYPE, Device>;
    using delmem_var_op = psi::memory::delete_memory_op<FPTYPE, Device>;
    using delmem_var_h_op = psi::memory::delete_memory_op<FPTYPE, psi::DEVICE_CPU>;
    using syncmem_var_h2d_op = psi::memory::synchronize_memory_op<FPTYPE, Device, psi::DEVICE_CPU>;
    using syncmem_var_d2h_op = psi::memory::synchronize_memory_op<FPTYPE, psi::DEVICE_CPU, Device>;

    using resmem_int_op = psi::memory::resize_memory_op<int, Device>;
    using delmem_int_op = psi::memory::delete_memory_op<int, Device>;
    using syncmem_int_h2d_op = psi::memory::synchronize_memory_op<int, Device, psi::DEVICE_CPU>;
};

#endif
