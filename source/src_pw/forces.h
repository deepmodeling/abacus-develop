#ifndef FORCES_H
#define FORCES_H

#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"
#include "../module_base/matrix.h"
#include "../module_pw/pw_basis.h"
#include "module_hsolver/include/math_kernel.h"
#include "module_psi/psi.h"
#include "module_psi/include/memory.h"
#include "charge.h"
#include "src_pw/include/force_multi_device.h"

template<typename FPTYPE, typename Device = psi::DEVICE_CPU>
class Forces
{
public:
	friend class Force_Stress_LCAO;
    /* This routine is a driver routine which compute the forces
     * acting on the atoms, the complete forces in plane waves
     * is computed from 4 main parts
     * (1) cal_force_loc: contribution due to local potential.
     * (2) cal_foce_ew: contribution due to ewald potential.
     * (3) cal_force_cc: contributino due to NLCC.
     * (4) cal_nl: contribution due to the non-local pseudopotential.
     * (4) cal_scc: contributino due to incomplete SCF calculation.
     */
    Forces() {};
    ~Forces() {};

    void init(ModuleBase::matrix& force, const ModuleBase::matrix& wg, const Charge* const chr, const psi::Psi<std::complex<FPTYPE>, Device>* psi_in=nullptr);

protected:

    int nat = 0;
	static FPTYPE output_acc;

    void cal_force_loc(ModuleBase::matrix& forcelc, ModulePW::PW_Basis* rho_basis, const Charge* const chr);
    void cal_force_ew(ModuleBase::matrix& forceion, ModulePW::PW_Basis* rho_basis);
    void cal_force_cc(ModuleBase::matrix& forcecc, ModulePW::PW_Basis* rho_basis, const Charge* const chr);
    void cal_force_nl(ModuleBase::matrix& forcenl, const ModuleBase::matrix& wg, const psi::Psi<std::complex<FPTYPE>, Device>* psi_in=nullptr);
    void cal_force_scc(ModuleBase::matrix& forcescc, ModulePW::PW_Basis* rho_basis);

    static void print( const std::string &name, const ModuleBase::matrix &f, bool rv=true );
    static void print_to_files( std::ofstream &ofs, const std::string &name, const ModuleBase::matrix &f );

private:
    Device * ctx = {};
    psi::DEVICE_CPU * cpu_ctx = {};
    psi::AbacusDevice_t device = {};
    using gemm_op = hsolver::gemm_op<FPTYPE, Device>;
    using cal_vkb1_nl_op = src_pw::cal_vkb1_nl_op<FPTYPE, Device>;
    using cal_force_nl_op = src_pw::cal_force_nl_op<FPTYPE, Device>;

    using resmem_complex_op = psi::memory::resize_memory_op<std::complex<FPTYPE>, Device>;
    using resmem_complex_h_op = psi::memory::resize_memory_op<std::complex<FPTYPE>, psi::DEVICE_CPU>;
    using delmem_complex_op = psi::memory::delete_memory_op<std::complex<FPTYPE>, Device>;
    using delmem_complex_h_op = psi::memory::delete_memory_op<std::complex<FPTYPE>, psi::DEVICE_CPU>;
    using syncmem_complex_h2d_op = psi::memory::synchronize_memory_op<std::complex<FPTYPE>, Device, psi::DEVICE_CPU>;
    using syncmem_complex_d2h_op = psi::memory::synchronize_memory_op<std::complex<FPTYPE>, psi::DEVICE_CPU, Device>;

    using resmem_var_op = psi::memory::resize_memory_op<FPTYPE, Device>;
    using delmem_var_op = psi::memory::delete_memory_op<FPTYPE, Device>;
    using syncmem_var_h2d_op = psi::memory::synchronize_memory_op<FPTYPE, Device, psi::DEVICE_CPU>;
    using syncmem_var_d2h_op = psi::memory::synchronize_memory_op<FPTYPE, psi::DEVICE_CPU, Device>;

    using resmem_int_op = psi::memory::resize_memory_op<int, Device>;
    using delmem_int_op = psi::memory::delete_memory_op<int, Device>;
    using syncmem_int_h2d_op = psi::memory::synchronize_memory_op<int, Device, psi::DEVICE_CPU>;

};

#endif
