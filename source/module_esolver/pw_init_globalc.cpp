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
#include "module_io/rho_io.h"
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
void ESolver_KS_PW<T, Device>::Init_GlobalC(const Input_para& inp,
                                            UnitCell& ucell,
                                            pseudopot_cell_vnl& ppcell) {
    // GlobalC is a historically left-over namespace, it is used to store global
    // classes, including: pseudopot_cell_vnl: pseudopotential in cell, V
    // non-local UnitCell: cell information with atomic properties Grid_Driver:
    // Parallel_Grid:
    // Parallel_Kpoints:
    // Restart:
    // Exx_Info:
    // Exx_Lip:

    // GlobalC would be refactored out in the future. If there is better idea
    // about how to organize information stored in classes above, please feel
    // free to discuss with issue or pull request.

    if (this->psi != nullptr) {
        delete this->psi;
    }

    //! init pseudopotential
    ppcell.init(ucell.ntype, &this->sf, this->pw_wfc);

    //! initalize local pseudopotential
    ppcell.init_vloc(ppcell.vloc, this->pw_rhod);
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "LOCAL POTENTIAL");

    //! Initalize non-local pseudopotential
    ppcell.init_vnl(ucell, this->pw_rhod);
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "NON-LOCAL POTENTIAL");

    //! Allocate psi
    this->p_wf_init->allocate_psi(this->psi,
                                  this->kv.get_nkstot(),
                                  this->kv.get_nks(),
                                  this->kv.ngk.data(),
                                  this->pw_wfc->npwk_max,
                                  &(this->sf));

    this->kspw_psi
        = PARAM.inp.device == "gpu" || PARAM.inp.precision == "single"
              ? new psi::Psi<T, Device>(this->psi[0])
              : reinterpret_cast<psi::Psi<T, Device>*>(this->psi);

    if (PARAM.inp.precision == "single") {
        ModuleBase::Memory::record("Psi_single",
                                   sizeof(T) * this->psi[0].size());
    }

    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "INIT BASIS");
}



template class ESolver_KS_PW<std::complex<float>, base_device::DEVICE_CPU>;
template class ESolver_KS_PW<std::complex<double>, base_device::DEVICE_CPU>;
#if ((defined __CUDA) || (defined __ROCM))
template class ESolver_KS_PW<std::complex<float>, base_device::DEVICE_GPU>;
template class ESolver_KS_PW<std::complex<double>, base_device::DEVICE_GPU>;
#endif
} // namespace ModuleESolver
