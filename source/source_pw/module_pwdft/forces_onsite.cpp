#include "forces.h"
#include "source_base/parallel_reduce.h"
#include "source_base/timer.h"
#include "source_base/tool_title.h"
#include "source_pw/module_pwdft/onsite_proj.h"
#include "source_pw/module_pwdft/kernels/force_op.h"
#include "source_io/module_parameter/parameter.h"
#include "source_lcao/module_dftu/dftu.h"
#include "source_lcao/module_deltaspin/spin_constrain.h"

template <typename FPTYPE, typename Device>
void Forces<FPTYPE, Device>::cal_force_onsite(ModuleBase::matrix& force_onsite,
                                          const ModuleBase::matrix& wg,
                                          const ModulePW::PW_Basis_K* wfc_basis,
										  const UnitCell& ucell_in,
										  const Plus_U &dftu,
										  const psi::Psi <std::complex<FPTYPE>, Device>* psi_in)
{
    ModuleBase::TITLE("Forces", "cal_force_onsite");
    if(psi_in == nullptr || wfc_basis == nullptr)
    {
        return;
    }
    ModuleBase::timer::start("Forces", "cal_force_onsite");

    FPTYPE* force = nullptr;
    resmem_var_op()(force, ucell_in.nat * 3);
    base_device::memory::set_memory_op<FPTYPE, Device>()(force, 0.0, ucell_in.nat * 3);

    auto* onsite_p = projectors::OnsiteProjector<FPTYPE, Device>::get_instance();

    const int nks = wfc_basis->nks;
    for (int ik = 0; ik < nks; ik++)
    {
        int nbands_occ = wg.nc;
        while (wg(ik, nbands_occ - 1) == 0.0)
        {
            nbands_occ--;
            if (nbands_occ == 0)
            {
                break;
            }
        }
        const int npm = nbands_occ;
        onsite_p->get_fs_tools()->cal_becp(ik, npm);
        for (int ipol = 0; ipol < 3; ipol++)
        {
            onsite_p->get_fs_tools()->cal_dbecp_f(ik, npm, ipol);
        }
        if(PARAM.inp.dft_plus_u)
        {
            onsite_p->cal_force_onsite_dftu(ik, npm, force, dftu, nks, wg.c);
        }
        if(PARAM.inp.sc_mag_switch)
        {
            spinconstrain::SpinConstrain<std::complex<double>>& sc = 
              spinconstrain::SpinConstrain<std::complex<double>>::getScInstance();
            onsite_p->cal_force_onsite_dspin(ik, npm, force, sc.get_sc_lambda().data(), wg.c);
        }
        
    }

    syncmem_var_d2h_op()(force_onsite.c, force, force_onsite.nr * force_onsite.nc);
    delmem_var_op()(force);
    Parallel_Reduce::reduce_all(force_onsite.c, force_onsite.nr * force_onsite.nc);

    ModuleBase::timer::end("Forces", "cal_force_onsite");
}

template class Forces<double, base_device::DEVICE_CPU>;
#if ((defined __CUDA) || (defined __ROCM))
template class Forces<double, base_device::DEVICE_GPU>;
#endif
