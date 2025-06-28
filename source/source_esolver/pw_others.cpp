#include "esolver_ks_pw.h"

#include "source_base/global_variable.h"
#include "source_pw/hamilt_pwdft/elecond.h"
#include "source_io/input_conv.h"
#include "source_io/output_log.h"

#include <iostream>

//--------------temporary----------------------------
#include "source_estate/module_charge/symmetry_rho.h"
#include "source_estate/occupy.h"
#include "source_hamilt/module_ewald/H_Ewald_pw.h"
#include "source_pw/hamilt_pwdft/global.h"
#include "source_io/print_info.h"
//-----force-------------------
#include "source_pw/hamilt_pwdft/forces.h"
//-----stress------------------
#include "source_pw/hamilt_pwdft/stress_pw.h"
//---------------------------------------------------
#include "source_base/memory.h"
#include "source_base/module_device/device.h"
#include "source_estate/elecstate_pw.h"
#include "source_hamilt/module_vdw/vdw.h"
#include "source_pw/hamilt_pwdft/hamilt_pw.h"
#include "source_hsolver/diago_iter_assist.h"
#include "source_hsolver/hsolver_pw.h"
#include "source_hsolver/kernels/dngvd_op.h"
#include "source_base/kernels/math_kernel_op.h"
#include "source_io/berryphase.h"
#include "source_io/numerical_basis.h"
#include "source_io/numerical_descriptor.h"
#include "source_io/to_wannier90_pw.h"
#include "source_io/winput.h"
#include "source_io/write_elecstat_pot.h"
#include "source_io/write_wfc_r.h"
#include "module_parameter/parameter.h"

#include <ATen/kernels/blas.h>
#include <ATen/kernels/lapack.h>
#include "source_base/formatter.h"

// mohan add 2025-03-06
#include "source_io/cal_test.h"

namespace ModuleESolver {

template <typename T, typename Device>
void ESolver_KS_PW<T, Device>::others(UnitCell& ucell, const int istep)
{
    ModuleBase::TITLE("ESolver_KS_PW", "others");

    const std::string cal_type = PARAM.inp.calculation;

    if (cal_type == "test_memory") 
    {
        Cal_Test::test_memory(ucell.nat,
                              ucell.ntype,
                              ucell.GGT,
                              this->pw_rho,
                              this->pw_wfc,
                              this->p_chgmix->get_mixing_mode(),
                              this->p_chgmix->get_mixing_ndim());
    } 
    else if (cal_type == "gen_bessel") 
    {
        Numerical_Descriptor nc;
        nc.output_descriptor(ucell,
                             this->psi[0],
                             PARAM.inp.bessel_descriptor_lmax,
                             PARAM.inp.bessel_descriptor_rcut,
                             PARAM.inp.bessel_descriptor_tolerence,
                             this->kv.get_nks());
        ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "GENERATE DESCRIPTOR FOR DEEPKS");
    } 
    else 
    {
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
