#include "sto_forces.h"

#include "module_base/mathzone.h"
#include "module_cell/module_symmetry/symmetry.h"
#include "module_elecstate/elecstate.h"
#include "module_elecstate/potentials/efield.h"
#include "module_elecstate/potentials/gatefield.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_io/output_log.h"
#include "module_parameter/parameter.h"
#include "module_hamilt_pw/hamilt_pwdft/fs_nonlocal_tools.h"

// new
#include "module_base/math_integral.h"
#include "module_base/parallel_reduce.h"
#include "module_base/timer.h"
#include "module_hamilt_general/module_xc/xc_functional.h"

template <typename FPTYPE, typename Device>
void Sto_Forces<FPTYPE, Device>::cal_stoforce(ModuleBase::matrix& force,
                                              const elecstate::ElecState& elec,
                                              ModulePW::PW_Basis* rho_basis,
                                              ModuleSymmetry::Symmetry* p_symm,
                                              const Structure_Factor* p_sf,
                                              K_Vectors* pkv,
                                              ModulePW::PW_Basis_K* wfc_basis,
                                              const pseudopot_cell_vnl& nlpp,
                                              UnitCell& ucell,
                                              const psi::Psi<std::complex<FPTYPE>, Device>& psi,
                                              const Stochastic_WF<std::complex<FPTYPE>, Device>& stowf)
{
    ModuleBase::timer::tick("Sto_Forces", "cal_force");
    ModuleBase::TITLE("Sto_Forces", "init");
    this->device = base_device::get_device_type<Device>(this->ctx);
    const ModuleBase::matrix& wg = elec.wg;
    const Charge* chr = elec.charge;
    force.create(this->nat, 3);

    ModuleBase::matrix forcelc(this->nat, 3);
    ModuleBase::matrix forceion(this->nat, 3);
    ModuleBase::matrix forcecc(this->nat, 3);
    ModuleBase::matrix forcenl(this->nat, 3);
    ModuleBase::matrix forcescc(this->nat, 3);
    this->cal_force_loc(ucell,forcelc, rho_basis, nlpp.vloc, chr);
    this->cal_force_ew(ucell,forceion, rho_basis, p_sf);
    this->cal_sto_force_nl(forcenl, wg, pkv, wfc_basis, p_sf, nlpp, ucell, psi, stowf);
    this->cal_force_cc(forcecc, rho_basis, chr, nlpp.numeric, ucell);
    this->cal_force_scc(forcescc, rho_basis, elec.vnew, elec.vnew_exist, nlpp.numeric, ucell);

    // impose total force = 0
    ModuleBase::matrix force_e;
    if (PARAM.inp.efield_flag)
    {
        force_e.create(this->nat, 3);
        elecstate::Efield::compute_force(ucell, force_e);
    }

    ModuleBase::matrix force_gate;
    if (PARAM.inp.gate_flag)
    {
        force_gate.create(this->nat, 3);
        elecstate::Gatefield::compute_force(ucell, force_gate);
    }

    int iat = 0;
    for (int ipol = 0; ipol < 3; ipol++)
    {
        double sum = 0.0;
        iat = 0;

        for (int it = 0; it < ucell.ntype; it++)
        {
            for (int ia = 0; ia < ucell.atoms[it].na; ia++)
            {
                force(iat, ipol) = forcelc(iat, ipol) + forceion(iat, ipol) + forcenl(iat, ipol) + forcecc(iat, ipol)
                                   + forcescc(iat, ipol);

                if (PARAM.inp.efield_flag)
                {
                    force(iat, ipol) = force(iat, ipol) + force_e(iat, ipol);
                }

                if (PARAM.inp.gate_flag)
                {
                    force(iat, ipol) = force(iat, ipol) + force_gate(iat, ipol);
                }

                sum += force(iat, ipol);

                iat++;
            }
        }

        if (!(PARAM.inp.gate_flag || PARAM.inp.efield_flag))
        {
            double compen = sum / ucell.nat;
            for (int iat = 0; iat < ucell.nat; ++iat)
            {
                force(iat, ipol) = force(iat, ipol) - compen;
            }
        }
    }

    if (PARAM.inp.gate_flag || PARAM.inp.efield_flag)
    {
        GlobalV::ofs_running << "Atomic forces are not shifted if gate_flag or efield_flag == true!" << std::endl;
    }

    if (ModuleSymmetry::Symmetry::symm_flag == 1)
    {
        double d1, d2, d3;
        for (int iat = 0; iat < ucell.nat; iat++)
        {
            ModuleBase::Mathzone::Cartesian_to_Direct(force(iat, 0),
                                                      force(iat, 1),
                                                      force(iat, 2),
                                                      ucell.a1.x,
                                                      ucell.a1.y,
                                                      ucell.a1.z,
                                                      ucell.a2.x,
                                                      ucell.a2.y,
                                                      ucell.a2.z,
                                                      ucell.a3.x,
                                                      ucell.a3.y,
                                                      ucell.a3.z,
                                                      d1,
                                                      d2,
                                                      d3);

            force(iat, 0) = d1;
            force(iat, 1) = d2;
            force(iat, 2) = d3;
        }
        p_symm->symmetrize_vec3_nat(force.c);
        for (int iat = 0; iat < ucell.nat; iat++)
        {
            ModuleBase::Mathzone::Direct_to_Cartesian(force(iat, 0),
                                                      force(iat, 1),
                                                      force(iat, 2),
                                                      ucell.a1.x,
                                                      ucell.a1.y,
                                                      ucell.a1.z,
                                                      ucell.a2.x,
                                                      ucell.a2.y,
                                                      ucell.a2.z,
                                                      ucell.a3.x,
                                                      ucell.a3.y,
                                                      ucell.a3.z,
                                                      d1,
                                                      d2,
                                                      d3);
            force(iat, 0) = d1;
            force(iat, 1) = d2;
            force(iat, 2) = d3;
        }
    }

    GlobalV::ofs_running << setiosflags(std::ios::fixed) << std::setprecision(6) << std::endl;

    // output force in unit eV/Angstrom
    GlobalV::ofs_running << std::endl;

    if (PARAM.inp.test_force)
    {
        ModuleIO::print_force(GlobalV::ofs_running, ucell, "LOCAL    FORCE (eV/Angstrom)", forcelc, false);
        ModuleIO::print_force(GlobalV::ofs_running, ucell, "NONLOCAL FORCE (eV/Angstrom)", forcenl, false);
        ModuleIO::print_force(GlobalV::ofs_running, ucell, "NLCC     FORCE (eV/Angstrom)", forcecc, false);
        ModuleIO::print_force(GlobalV::ofs_running, ucell, "ION      FORCE (eV/Angstrom)", forceion, false);
        ModuleIO::print_force(GlobalV::ofs_running, ucell, "SCC      FORCE (eV/Angstrom)", forcescc, false);
        if (PARAM.inp.efield_flag)
        {
            ModuleIO::print_force(GlobalV::ofs_running, ucell, "EFIELD   FORCE (eV/Angstrom)", force_e, false);
        }
        if (PARAM.inp.gate_flag)
        {
            ModuleIO::print_force(GlobalV::ofs_running,
                                  ucell,
                                  "GATEFIELD   FORCE (eV/Angstrom)",
                                  force_gate,
                                  false);
        }
    }
    ModuleIO::print_force(GlobalV::ofs_running, ucell, "TOTAL-FORCE (eV/Angstrom)", force, false);
    ModuleBase::timer::tick("Sto_Forces", "cal_force");
    return;
}

template <typename FPTYPE, typename Device>
void Sto_Forces<FPTYPE, Device>::cal_sto_force_nl(
    ModuleBase::matrix& forcenl,
    const ModuleBase::matrix& wg,
    K_Vectors* p_kv,
    ModulePW::PW_Basis_K* wfc_basis,
    const Structure_Factor* p_sf,
    const pseudopot_cell_vnl& nlpp,
    const UnitCell& ucell,
    const psi::Psi<std::complex<FPTYPE>, Device>& psi_in,
    const Stochastic_WF<std::complex<FPTYPE>, Device>& stowf)
{
    ModuleBase::TITLE("Sto_Forces", "cal_force_nl");
    const int nkb = nlpp.nkb;
    if (nkb == 0)
    {
        return;
    }

    ModuleBase::timer::tick("Sto_Forces", "cal_force_nl");

    const int* nchip = stowf.nchip;
    const int npwx = wfc_basis->npwk_max;
    int nksbands = psi_in.get_nbands();
    if (GlobalV::MY_STOGROUP != 0)
    {
        nksbands = 0;
    }

   // allocate memory for the force
    FPTYPE* force = nullptr;
    resmem_var_op()(this->ctx, force, ucell.nat * 3);
    base_device::memory::set_memory_op<FPTYPE, Device>()(this->ctx, force, 0.0, ucell.nat * 3);

    hamilt::FS_Nonlocal_tools<FPTYPE, Device> nl_tools(&nlpp, &ucell, p_kv, wfc_basis, p_sf, wg, nullptr);


    for (int ik = 0; ik < wfc_basis->nks; ik++)
    {
        const int nstobands = nchip[ik];
        const int max_nbands = stowf.shchi->get_nbands() + nksbands;
        const int npw = wfc_basis->npwk[ik];
        psi_in.fix_k(ik);
        stowf.shchi->fix_k(ik);

        nl_tools.cal_vkb(ik, max_nbands); // vkb has dimension of nkb * max_nbands * npol

        // calculate becp = <psi|beta> for all beta functions
        nl_tools.cal_becp(ik, nksbands, psi_in.get_pointer(), 0);
        nl_tools.cal_becp(ik, nstobands, stowf.shchi->get_pointer(), nksbands);
        nl_tools.reduce_pool_becp(max_nbands);

        for (int ipol = 0; ipol < 3; ipol++)
        {
            nl_tools.cal_vkb_deri_f(ik, max_nbands, ipol); // vkb_deri has dimension of nkb * max_nbands * npol
            // calculate dbecp = <psi|\nabla beta> for all beta functions
            nl_tools.cal_dbecp_f(ik, max_nbands, nksbands, ipol, psi_in.get_pointer(), 0);
            nl_tools.cal_dbecp_f(ik, max_nbands, nstobands, ipol, stowf.shchi->get_pointer(), nksbands);
            nl_tools.revert_vkb(ik, ipol);
        }
        nl_tools.cal_force(ik, max_nbands, nksbands, true, force, 0);
        nl_tools.cal_force(ik, max_nbands, nstobands, false, force, nksbands);
    } // end ik

    syncmem_var_d2h_op()(this->cpu_ctx, this->ctx, forcenl.c, force, forcenl.nr * forcenl.nc);
    delmem_var_op()(this->ctx, force);
    // sum up forcenl from all processors
    Parallel_Reduce::reduce_all(forcenl.c, forcenl.nr * forcenl.nc);

    
    ModuleBase::timer::tick("Sto_Forces", "cal_force_nl");
    return;
}

template class Sto_Forces<double, base_device::DEVICE_CPU>;
#if ((defined __CUDA) || (defined __ROCM))
template class Sto_Forces<double, base_device::DEVICE_GPU>;
#endif
