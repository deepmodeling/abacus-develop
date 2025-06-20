#ifndef VELOCITY_PW_H
#define VELOCITY_PW_H
#include "operator_pw.h"
#include "source_cell/unitcell.h"
#include "module_hamilt_pw/hamilt_pwdft/VNL_in_pw.h"
#include "source_basis/module_pw/pw_basis_k.h"
namespace hamilt
{

//velocity operator mv = im/\hbar * [H,r] =  p + im/\hbar [V_NL, r] 
template <typename FPTYPE, typename Device = base_device::DEVICE_CPU>
class Velocity
{
    public:
    Velocity(
        const ModulePW::PW_Basis_K* wfcpw_in,
        const int* isk_in,
        pseudopot_cell_vnl* ppcell_in,
        const UnitCell* ucell_in,
        const bool nonlocal_in = true
    );

    ~Velocity();

    void init(const int ik_in);

    /**
     * @brief calculate \hat{v}|\psi>
     * 
     * @param psi_in Psi class which contains some information
     * @param n_npwx nbands * NPOL
     * @param tmpsi_in |\psi_i>    size: n_npwx*npwx 
     * @param tmvpsi \hat{v}|\psi> size: 3*n_npwx*npwx
     * @param add true : tmvpsi = tmvpsi + v|\psi>  false: tmvpsi = v|\psi>
     * 
     */
    void act(const psi::Psi<std::complex<FPTYPE>, Device>* psi_in,
             const int n_npwx,
             const std::complex<FPTYPE>* tmpsi_in,
             std::complex<FPTYPE>* tmvpsi,
             const bool add = false) const;
    
    bool nonlocal = true;

    private:
    const ModulePW::PW_Basis_K* wfcpw = nullptr;

    const int* isk = nullptr;

    pseudopot_cell_vnl* ppcell = nullptr;

    const UnitCell* ucell = nullptr;

    int ik=0;

    double tpiba=0.0;
  
  private:
    FPTYPE* gx_ = nullptr; ///<[Device, npwx] x component of G+K
    FPTYPE* gy_ = nullptr; ///<[Device, npwx] y component of G+K
    FPTYPE* gz_ = nullptr; ///<[Device, npwx] z component of G+K
    std::complex<FPTYPE>* vkb_ = nullptr;     ///<[Device, nkb * npwk_max] nonlocal pseudopotential vkb
    std::complex<FPTYPE>* gradvkb_ = nullptr; ///<[Device, 3*nkb * npwk_max] gradient of nonlocal pseudopotential gradvkb
    FPTYPE* deeq_ = nullptr;                  ///<[Device] D matrix for nonlocal pseudopotential
    
    using Complex = std::complex<FPTYPE>;
    using resmem_var_op = base_device::memory::resize_memory_op<FPTYPE, Device>;
    using delmem_var_op = base_device::memory::delete_memory_op<FPTYPE, Device>;
    using syncmem_var_h2d_op = base_device::memory::synchronize_memory_op<FPTYPE, Device, base_device::DEVICE_CPU>;
    using castmem_var_h2d_op = base_device::memory::cast_memory_op<FPTYPE, double, Device, base_device::DEVICE_CPU>;
    using resmem_complex_op = base_device::memory::resize_memory_op<std::complex<FPTYPE>, Device>;
    using setmem_complex_op = base_device::memory::set_memory_op<std::complex<FPTYPE>, Device>;
    using delmem_complex_op = base_device::memory::delete_memory_op<std::complex<FPTYPE>, Device>;
    using castmem_complex_h2d_op = base_device::memory::cast_memory_op<std::complex<FPTYPE>, std::complex<double>, Device, base_device::DEVICE_CPU>;
    using syncmem_complex_d2h_op = base_device::memory::synchronize_memory_op<std::complex<FPTYPE>, base_device::DEVICE_CPU, Device>;
    using syncmem_complex_h2d_op = base_device::memory::synchronize_memory_op<std::complex<FPTYPE>, Device, base_device::DEVICE_CPU>;
};
}
#endif