#ifndef SETUP_PSI_PW_H
#define SETUP_PSI_PW_H

#include "source_psi/psi_prepare.h"
#include "source_cell/unitcell.h"
#include "source_cell/klist.h"
#include "source_pw/module_pwdft/structure_factor.h"
#include "source_basis/module_pw/pw_basis_k.h"
#include "source_pw/module_pwdft/vnl_pw.h"
#include "source_io/module_parameter/input_parameter.h"
#include "source_base/module_device/device.h"
#include "source_hamilt/hamilt.h"

class Setup_Psi_pw
{
    public:

    Setup_Psi_pw();
    ~Setup_Psi_pw();

    //------------
    // public types
    //------------
    
    // Precision type: 0 = float, 1 = double, 2 = complex<float>, 3 = complex<double>
    enum class PrecisionType {
        Float = 0,
        Double = 1,
        ComplexFloat = 2,
        ComplexDouble = 3
    };

    //------------
    // variables
    // psi_cpu, complex<double> on cpu
    //------------

    // originally, this term is psi
    // for PW, we have psi_cpu
    psi::Psi<std::complex<double>, base_device::DEVICE_CPU>* psi_cpu = nullptr;

    // psi_initializer controller
    psi::PSIPrepareBase* p_psi_init = nullptr;

    //------------
    // functions
    //------------

    void before_runner(
        const UnitCell &ucell,
        const K_Vectors &kv,
        const Structure_Factor &sf,
        const ModulePW::PW_Basis_K &pw_wfc, 
        const pseudopot_cell_vnl &ppcell,
        const Input_para &inp);

    void init(hamilt::HamiltBase* p_hamilt);

    void update_psi_d();

    // Transfer data from device to host in pw basis
    void copy_d2h();

    void clean();

    //------------
    // accessor functions
    //------------
    
    // Get basic information (no type conversion needed, use psi_cpu)
    int get_nbands() const { return this->psi_cpu->get_nbands(); }
    int get_nk() const { return this->psi_cpu->get_nk(); }
    int get_nbasis() const { return this->psi_cpu->get_nbasis(); }
    size_t size() const { return this->psi_cpu->size(); }
    
    // Get runtime type information
    base_device::AbacusDevice_t get_device_type() const { return device_type_; }
    PrecisionType get_precision_type() const { return precision_type_; }
    
    // Get psi_t pointer (template version, for backward compatibility)
    template <typename T, typename Device>
    psi::Psi<T, Device>* get_psi_t() { return static_cast<psi::Psi<T, Device>*>(psi_t); }
    
    template <typename T, typename Device>
    const psi::Psi<T, Device>* get_psi_t() const { return static_cast<const psi::Psi<T, Device>*>(psi_t); }
    
    // Get psi_d pointer (template version, for backward compatibility)
    template <typename T, typename Device>
    psi::Psi<std::complex<double>, Device>* get_psi_d() { 
        return static_cast<psi::Psi<std::complex<double>, Device>*>(psi_d); 
    }
    
    template <typename T, typename Device>
    const psi::Psi<std::complex<double>, Device>* get_psi_d() const { 
        return static_cast<const psi::Psi<std::complex<double>, Device>*>(psi_d); 
    }

    private:

    //------------
    // private variables
    //------------
    
    // originally, this term is kspw_psi
    // if CPU, kspw_psi = psi, otherwise, kspw_psi has a new copy
    void* psi_t = nullptr;  // Use void* to store pointer, runtime type information records actual type 

    // originally, this term is __kspw_psi
    void* psi_d = nullptr;  // Use void* to store pointer, runtime type information records actual type

    bool already_initpsi = false;

    //------------
    // runtime type information
    //------------
    base_device::AbacusDevice_t device_type_ = base_device::CpuDevice;
    PrecisionType precision_type_ = PrecisionType::ComplexDouble;

    //------------
    // private functions
    //------------

    template <typename T, typename Device>
    void before_runner_impl(
        const UnitCell &ucell,
        const K_Vectors &kv,
        const Structure_Factor &sf,
        const ModulePW::PW_Basis_K &pw_wfc, 
        const pseudopot_cell_vnl &ppcell,
        const Input_para &inp);

    template <typename T, typename Device>
    void init_impl(hamilt::Hamilt<T, Device>* p_hamilt);

    template <typename T, typename Device>
    void update_psi_d_impl();

    template <typename T, typename Device>
    void clean_impl();

    template <typename T, typename Device>
    void copy_d2h_impl();

    template <typename T, typename Device>
    void castmem_d2h_impl(std::complex<double>* dst, const std::complex<double>* src, const size_t size);
    
    template <typename T, typename Device>
    void castmem_d2h_impl(std::complex<double>* dst, const std::complex<float>* src, const size_t size);

};


#endif
