#include "forces.h"
#include "module_base/timer.h"
#include "module_base/tool_title.h"
#include "module_hamilt_pw/hamilt_pwdft/onsite_projector.h"
#include "module_hamilt_pw/hamilt_pwdft/kernels/force_op.h"
#include "module_parameter/parameter.h"
#include "module_hamilt_lcao/module_dftu/dftu.h"
#include "module_hamilt_lcao/module_deltaspin/spin_constrain.h"

template <typename FPTYPE, typename Device>
void Forces<FPTYPE, Device>::cal_force_onsite(ModuleBase::matrix& force_onsite,
                                          const ModuleBase::matrix& wg,
                                          const ModulePW::PW_Basis_K* wfc_basis,
                                          const UnitCell& ucell_in,
                                          const psi::Psi<complex<FPTYPE>, Device>* psi_in)
{
    ModuleBase::TITLE("Forces", "cal_force_onsite");
    if(psi_in == nullptr || wfc_basis == nullptr)
    {
        return;
    }
    ModuleBase::timer::tick("Forces", "cal_force_onsite");

    // allocate memory for the force
    FPTYPE* force = nullptr;
    resmem_var_op()(force, ucell_in.nat * 3);
    base_device::memory::set_memory_op<FPTYPE, Device>()(force, 0.0, ucell_in.nat * 3);

    auto* onsite_p = projectors::OnsiteProjector<FPTYPE, Device>::get_instance();

    const int nks = wfc_basis->nks;
    for (int ik = 0; ik < nks; ik++) // loop k points
    {
        // skip zero weights to speed up
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
        // calculate becp = <psi|beta> for all beta functions
        for (int ipol = 0; ipol < 3; ipol++)
        {
            // calculate dbecp = <psi|\nabla beta> for all beta functions
            onsite_p->get_fs_tools()->cal_dbecp_f(ik, npm, ipol);
        }
        // calculate the force_i = \sum_{n,k}f_{nk}\sum_I \sum_{lm,l'm'}D_{l,l'}^{I} becp * dbecp_i
        // force for DFT+U
        if(PARAM.inp.dft_plus_u)
        {
            auto* dftu = ModuleDFTU::DFTU::get_instance();
            onsite_p->get_fs_tools()->cal_force_dftu(ik, npm, force, dftu->orbital_corr.data(), dftu->get_eff_pot_pw(0), dftu->get_size_eff_pot_pw(), wg.c);
        }
        if(PARAM.inp.sc_mag_switch)
        {
            spinconstrain::SpinConstrain<std::complex<double>>& sc = spinconstrain::SpinConstrain<std::complex<double>>::getScInstance();
            const std::vector<ModuleBase::Vector3<double>>& lambda = sc.get_sc_lambda();
            onsite_p->get_fs_tools()->cal_force_dspin(ik, npm, force, lambda.data(), wg.c);
        }
        
    } // end ik

    syncmem_var_d2h_op()(force_onsite.c, force, force_onsite.nr * force_onsite.nc);
    delmem_var_op()(force);
    // sum up force_onsite from all processors
    Parallel_Reduce::reduce_all(force_onsite.c, force_onsite.nr * force_onsite.nc);

    ModuleBase::timer::tick("Forces", "cal_force_onsite");
}

template class Forces<double, base_device::DEVICE_CPU>;
#if ((defined __CUDA) || (defined __ROCM))
template class Forces<double, base_device::DEVICE_GPU>;
#endif