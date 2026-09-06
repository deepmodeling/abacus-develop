#include "esolver_ks_pw.h"
#include "source_base/formatter.h"
#include "source_base/module_device/device.h"
#include "source_io/module_bessel/numerical_descriptor.h"

// mohan add 2025-03-06
#include "source_io/module_output/cal_test.h"

namespace ModuleESolver
{

template <typename T, typename Device>
void ESolver_KS_PW<T, Device>::others(BaseCell& basecell, const int istep)
{
    basecell.require_kind(BaseCell::Kind::unitcell, __FUNCTION__);
    UnitCell& ucell = static_cast<UnitCell&>(basecell);

    ModuleBase::TITLE("ESolver_KS_PW", "others");

    const std::string cal_type = this->inp_->calculation;

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
                             *(this->stp.psi_cpu),
                             this->inp_->bessel_descriptor_lmax,
                             this->inp_->bessel_descriptor_rcut,
                             this->inp_->bessel_descriptor_tolerence,
                             this->kv.get_nks());
        ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "GENERATE DESCRIPTOR FOR DEEPKS");
    }
    else
    {
        ModuleBase::WARNING_QUIT("ESolver_KS_PW::others", "CALCULATION type not supported");
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
