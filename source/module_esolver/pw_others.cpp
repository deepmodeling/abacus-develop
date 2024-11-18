#include "esolver_ks_pw.h"

#include "module_base/global_variable.h"
#include "module_hamilt_pw/hamilt_pwdft/elecond.h"
#include "module_io/input_conv.h"
#include "module_io/nscf_band.h"
#include "module_io/output_log.h"
#include "module_io/write_dos_pw.h"
#include "module_io/write_istate_info.h"
#include "module_io/write_wfc_pw.h"

#include <iostream>

//--------------temporary----------------------------
#include "module_elecstate/module_charge/symmetry_rho.h"
#include "module_elecstate/occupy.h"
#include "module_hamilt_general/module_ewald/H_Ewald_pw.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_io/print_info.h"
//-----force-------------------
#include "module_hamilt_pw/hamilt_pwdft/forces.h"
//-----stress------------------
#include "module_hamilt_pw/hamilt_pwdft/stress_pw.h"
//---------------------------------------------------
#include "module_base/memory.h"
#include "module_base/module_device/device.h"
#include "module_elecstate/elecstate_pw.h"
#include "module_hamilt_general/module_vdw/vdw.h"
#include "module_hamilt_pw/hamilt_pwdft/hamilt_pw.h"
#include "module_hsolver/diago_iter_assist.h"
#include "module_hsolver/hsolver_pw.h"
#include "module_hsolver/kernels/dngvd_op.h"
#include "module_hsolver/kernels/math_kernel_op.h"
#include "module_io/berryphase.h"
#include "module_io/numerical_basis.h"
#include "module_io/numerical_descriptor.h"
#include "module_io/to_wannier90_pw.h"
#include "module_io/winput.h"
#include "module_io/write_elecstat_pot.h"
#include "module_io/write_wfc_r.h"
#include "module_parameter/parameter.h"
#ifdef USE_PAW
#include "module_cell/module_paw/paw_cell.h"
#endif
#include <ATen/kernels/blas.h>
#include <ATen/kernels/lapack.h>
#include "module_base/formatter.h"

namespace ModuleESolver {


template <typename T, typename Device>
void ESolver_KS_PW<T, Device>::others(const int istep) {
    ModuleBase::TITLE("ESolver_KS_PW", "others");

    const std::string cal_type = PARAM.inp.calculation;

    if (cal_type == "test_memory") {
        Cal_Test::test_memory(this->pw_rho,
                              this->pw_wfc,
                              this->p_chgmix->get_mixing_mode(),
                              this->p_chgmix->get_mixing_ndim());
    } else if (cal_type == "gen_bessel") {
        Numerical_Descriptor nc;
        nc.output_descriptor(this->psi[0],
                             PARAM.inp.bessel_descriptor_lmax,
                             PARAM.inp.bessel_descriptor_rcut,
                             PARAM.inp.bessel_descriptor_tolerence,
                             this->kv.get_nks());
        ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "GENERATE DESCRIPTOR FOR DEEPKS");
    } else {
        ModuleBase::WARNING_QUIT("ESolver_KS_PW::others",
                                 "CALCULATION type not supported");
    }

    return;
}

template class ESolver_KS_PW<std::complex<float>, base_device::DEVICE_CPU>;
template class ESolver_KS_PW<std::complex<double>, base_device::DEVICE_CPU>;
#if ((defined __CUDA) || (defined __ROCM))
template class ESolver_KS_PW<std::complex<float>, base_device::DEVICE_GPU>;
template class ESolver_KS_PW<std::complex<double>, base_device::DEVICE_GPU>;
#endif
} // namespace ModuleESolver
