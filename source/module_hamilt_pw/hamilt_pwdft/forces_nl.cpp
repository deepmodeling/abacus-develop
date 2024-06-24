#include "forces.h"
#include "module_base/timer.h"
#include "module_base/tool_title.h"
#include "module_hamilt_pw/hamilt_pwdft/fs_nonlocal_tools.h"

#ifdef _OPENMP
#include <omp.h>
#endif

template <typename FPTYPE, typename Device>
void Forces<FPTYPE, Device>::cal_force_nl(ModuleBase::matrix& forcenl,
                                          const ModuleBase::matrix& wg,
                                          const ModuleBase::matrix& ekb,
                                          const K_Vectors* p_kv,
                                          const ModulePW::PW_Basis_K* wfc_basis,
                                          const Structure_Factor* p_sf,
                                          pseudopot_cell_vnl* nlpp_in,
                                          const UnitCell& ucell_in,
                                          const psi::Psi<complex<FPTYPE>, Device>* psi_in)
{
    ModuleBase::TITLE("Forces", "cal_force_nl");
    if (nlpp_in->nkb == 0 || psi_in == nullptr || wfc_basis == nullptr)
    {
        return;
    }
    ModuleBase::timer::tick("Forces", "cal_force_nl");

    // allocate memory for the force
    FPTYPE* force = nullptr;
    resmem_var_op()(this->ctx, force, ucell_in.nat * 3);
    base_device::memory::set_memory_op<FPTYPE, Device>()(this->ctx, force, 0.0, ucell_in.nat * 3);

    hamilt::FS_Nonlocal_tools<FPTYPE, Device> nl_tools(nlpp_in, &ucell_in, psi_in, p_kv, wfc_basis, p_sf, wg, ekb);

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
        const int npm = ucell_in.get_npol() * nbands_occ;
        // calculate becp = <psi|beta> for all beta functions
        nl_tools.cal_becp(ik, npm);
        for (int ipol = 0; ipol < 3; ipol++)
        {
            nl_tools.cal_dbecp_f(ik, npm, ipol);
        }
        nl_tools.cal_force(ik, npm, force);
    } // end ik

    syncmem_var_d2h_op()(this->cpu_ctx, this->ctx, forcenl.c, force, forcenl.nr * forcenl.nc);
    delmem_var_op()(this->ctx, force);
    // sum up forcenl from all processors
    Parallel_Reduce::reduce_all(forcenl.c, forcenl.nr * forcenl.nc);

    ModuleBase::timer::tick("Forces", "cal_force_nl");
}

template class Forces<double, base_device::DEVICE_CPU>;
#if ((defined __CUDA) || (defined __ROCM))
template class Forces<double, base_device::DEVICE_GPU>;
#endif