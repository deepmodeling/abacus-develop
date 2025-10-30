#ifndef SETUP_PSI_PW_H
#define SETUP_PSI_PW_H

#include "source_psi/psi_init.h"
#include "source_cell/unitcell.h"
#include "source_cell/klist.h"
#include "source_pw/module_pwdft/structure_factor.h"
#include "source_basis/module_pw/pw_basis_k.h"
#include "source_pw/module_pwdft/VNL_in_pw.h"
#include "source_io/module_parameter/input_parameter.h"
#include "source_base/module_device/device.h"
#include "source_hamilt/hamilt.h"

template <typename T, typename Device = base_device::DEVICE_CPU>
class Setup_Psi_pw
{
    public:

    Setup_Psi_pw();
    ~Setup_Psi_pw();

    //------------
    // variables
    // psi_cpu, complex<double> on cpu
    // psi_t, complex<T> on cpu/gpu
    // psi_d, complex<double> on cpu/gpu 
    //------------

    // originally, this term is psi
    // for PW, we have psi_cpu
    psi::Psi<std::complex<double>, base_device::DEVICE_CPU>* psi_cpu = nullptr;

    // originally, this term is kspw_psi
    // if CPU, kspw_psi = psi, otherwise, kspw_psi has a new copy
    psi::Psi<T, Device>* psi_t = nullptr; 

    // originally, this term is __kspw_psi
    psi::Psi<std::complex<double>, Device>* psi_d = nullptr;

    // psi_initializer controller
    psi::PSIInit<T, Device>* p_psi_init = nullptr;

    bool already_initpsi = false;

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

    void init(hamilt::Hamilt<T, Device>* p_hamilt);

    void update_psi_d();

    // Transfer data from device to host in pw basis
    void copy_d2h(const base_device::AbacusDevice_t &device);

    void clean();

    private:

    using castmem_2d_d2h_op
        = base_device::memory::cast_memory_op<std::complex<double>, T, base_device::DEVICE_CPU, Device>;

};


#endif
