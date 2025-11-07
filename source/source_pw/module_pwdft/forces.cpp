#include "forces.h"

#include "source_io/module_parameter/parameter.h"
#include "source_pw/module_pwdft/global.h"
#include "source_io/output_log.h"
// new
#include "source_base/complexmatrix.h"
#include "source_base/libm/libm.h"
#include "source_base/math_integral.h"
#include "source_base/mathzone.h"
#include "source_base/timer.h"
#include "source_base/tool_threading.h"
#include "source_estate/module_pot/efield.h"
#include "source_estate/module_pot/gatefield.h"
#include "source_hamilt/module_ewald/H_Ewald_pw.h"
#include "source_hamilt/module_surchem/surchem.h"
#include "source_hamilt/module_vdw/vdw.h"

#ifdef _OPENMP
#include <omp.h>
#endif


template <typename FPTYPE, typename Device>
void Forces<FPTYPE, Device>::cal_force(UnitCell& ucell,
                                       ModuleBase::matrix& force,
                                       const elecstate::ElecState& elec,
                                       const ModulePW::PW_Basis* const rho_basis,
                                       ModuleSymmetry::Symmetry* p_symm,
                                       Structure_Factor* p_sf,
                                       surchem& solvent,
                                       const Plus_U *p_dftu, //mohan add 2025-11-06
                                       const pseudopot_cell_vl* locpp,
                                       const pseudopot_cell_vnl* p_nlpp,
                                       K_Vectors* pkv,
                                       ModulePW::PW_Basis_K* wfc_basis,
                                       const psi::Psi<std::complex<FPTYPE>, Device>* psi_in)
{
    ModuleBase::timer::tick("Forces", "cal_force");
    ModuleBase::TITLE("Forces", "init");
    this->device = base_device::get_device_type<Device>(this->ctx);
    const ModuleBase::matrix& wg = elec.wg;
    const ModuleBase::matrix& ekb = elec.ekb;
    const Charge* const chr = elec.charge;
    force.create(nat, 3);

    ModuleBase::matrix forcelc(nat, 3);
    ModuleBase::matrix forceion(nat, 3);
    ModuleBase::matrix forcecc(nat, 3);
    ModuleBase::matrix forcenl(nat, 3);
    ModuleBase::matrix forcescc(nat, 3);
    ModuleBase::matrix forcepaw(nat, 3);
    ModuleBase::matrix forceonsite(nat, 3);

    // Force due to local ionic potential
    this->cal_force_loc(ucell,forcelc, rho_basis, locpp->vloc, chr);

    // Ewald
    this->cal_force_ew(ucell,forceion, rho_basis, p_sf);

    // Force due to nonlocal part of pseudopotential
    if (wfc_basis != nullptr)
    {

        this->npwx = wfc_basis->npwk_max;
        Forces::cal_force_nl(forcenl, wg, ekb, pkv, wfc_basis, p_sf, *p_nlpp, ucell, psi_in);

        if (PARAM.globalv.use_uspp)
        {
            this->cal_force_us(forcenl, rho_basis, *p_nlpp, elec, ucell);
        }

        // DFT+U and DeltaSpin
        // here maybe a bug when OFDFT calls +U, mohan add 20251107
        if(PARAM.inp.dft_plus_u || PARAM.inp.sc_mag_switch)
        {
            this->cal_force_onsite(forceonsite, wg, wfc_basis, ucell, *p_dftu, psi_in);
        }
    }

    // non-linear core correction
    Forces::cal_force_cc(forcecc, rho_basis, chr, locpp->numeric, ucell);

    // force due to core charge
    this->cal_force_scc(forcescc, rho_basis, elec.vnew, elec.vnew_exist, locpp->numeric, ucell);

    ModuleBase::matrix stress_vdw_pw; //.create(3,3);
    ModuleBase::matrix force_vdw;
    force_vdw.create(nat, 3);
    auto vdw_solver = vdw::make_vdw(ucell, PARAM.inp);
    if (vdw_solver != nullptr)
    {
        const std::vector<ModuleBase::Vector3<double>>& force_vdw_temp = vdw_solver->get_force();
        for (int iat = 0; iat < this->nat; ++iat)
        {
            force_vdw(iat, 0) = force_vdw_temp[iat].x;
            force_vdw(iat, 1) = force_vdw_temp[iat].y;
            force_vdw(iat, 2) = force_vdw_temp[iat].z;
        }
        if (PARAM.inp.test_force)
        {
            ModuleIO::print_force(GlobalV::ofs_running, ucell, "VDW      FORCE (Ry/Bohr)", force_vdw);
        }
    }

    ModuleBase::matrix force_e;
    if (PARAM.inp.efield_flag)
    {
        force_e.create(this->nat, 3);
        elecstate::Efield::compute_force(ucell, force_e);
        if (PARAM.inp.test_force)
        {
            ModuleIO::print_force(GlobalV::ofs_running, ucell, "EFIELD      FORCE (Ry/Bohr)", force_e);
        }
    }

    ModuleBase::matrix force_gate;
    if (PARAM.inp.gate_flag)
    {
        force_gate.create(this->nat, 3);
        elecstate::Gatefield::compute_force(ucell, force_gate);
        if (PARAM.inp.test_force)
        {
            ModuleIO::print_force(GlobalV::ofs_running, ucell, "GATEFIELD      FORCE (Ry/Bohr)", force_gate);
        }
    }

    ModuleBase::matrix forcesol;
    if (PARAM.inp.imp_sol)
    {
        forcesol.create(this->nat, 3);
        solvent.cal_force_sol(ucell, rho_basis, locpp->vloc, forcesol);
        if (PARAM.inp.test_force)
        {
            ModuleIO::print_force(GlobalV::ofs_running, ucell, "IMP_SOL      FORCE (Ry/Bohr)", forcesol);
        }
    }

    // impose total force = 0
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

                if (vdw_solver != nullptr) // linpz and jiyy added vdw force, modified by zhengdy
                {
                    force(iat, ipol) += force_vdw(iat, ipol);
                }

                if (PARAM.inp.efield_flag)
                {
                    force(iat, ipol) = force(iat, ipol) + force_e(iat, ipol);
                }

                if (PARAM.inp.gate_flag)
                {
                    force(iat, ipol) = force(iat, ipol) + force_gate(iat, ipol);
                }

                if (PARAM.inp.imp_sol)
                {
                    force(iat, ipol) = force(iat, ipol) + forcesol(iat, ipol);
                }

                if(PARAM.inp.dft_plus_u || PARAM.inp.sc_mag_switch)
                {
                    force(iat, ipol) += forceonsite(iat, ipol);
                }

                sum += force(iat, ipol);

                iat++;
            }
        }

        if (!(PARAM.inp.gate_flag || PARAM.inp.efield_flag))
        {
            double compen = sum / this->nat;
            for (int iat = 0; iat < this->nat; ++iat)
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
        for (int iat = 0; iat < this->nat; iat++)
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
        for (int iat = 0; iat < this->nat; iat++)
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

    GlobalV::ofs_running << std::setiosflags(std::ios::fixed) << std::setprecision(6) << std::endl;
    /*if(PARAM.inp.test_force)
    {
        ModuleIO::print_force(GlobalV::ofs_running, ucell,"LOCAL    FORCE (Ry/Bohr)", forcelc);
        ModuleIO::print_force(GlobalV::ofs_running, ucell,"NONLOCAL FORCE (Ry/Bohr)", forcenl);
        ModuleIO::print_force(GlobalV::ofs_running, ucell,"NLCC     FORCE (Ry/Bohr)", forcecc);
        ModuleIO::print_force(GlobalV::ofs_running, ucell,"ION      FORCE (Ry/Bohr)", forceion);
        ModuleIO::print_force(GlobalV::ofs_running, ucell,"SCC      FORCE (Ry/Bohr)", forcescc);
        if(GlobalV::EFIELD) ModuleIO::print_force(GlobalV::ofs_running, ucell,"EFIELD   FORCE (Ry/Bohr)",
    force_e);
    }*/

    /*
        ModuleIO::print_force(GlobalV::ofs_running, ucell,"   TOTAL-FORCE (Ry/Bohr)", force);

        if(INPUT.out_force)                                                   // pengfei 2016-12-20
        {
            std::ofstream ofs("FORCE.dat");
            if(!ofs)
            {
                std::cout << "open FORCE.dat error !" <<std::endl;
            }
            for(int iat=0; iat<this->nat; iat++)
            {
                ofs << "   " << force(iat,0)*ModuleBase::Ry_to_eV / 0.529177
                    << "   " << force(iat,1)*ModuleBase::Ry_to_eV / 0.529177
                    << "   " << force(iat,2)*ModuleBase::Ry_to_eV / 0.529177 << std::endl;
            }
            ofs.close();
        }
    */

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
        if (PARAM.inp.imp_sol)
        {
            ModuleIO::print_force(GlobalV::ofs_running,
                                  ucell,
                                  "IMP_SOL   FORCE (eV/Angstrom)",
                                  forcesol,
                                  false);
        }
        if (PARAM.inp.dft_plus_u || PARAM.inp.sc_mag_switch)
        {
            ModuleIO::print_force(GlobalV::ofs_running,
                                  ucell,
                                  "ONSITE_PROJ    FORCE (eV/Angstrom)",
                                  forceonsite,
                                  false);
        }
    }
    ModuleIO::print_force(GlobalV::ofs_running, ucell, "TOTAL-FORCE (eV/Angstrom)", force, false);
    ModuleBase::timer::tick("Forces", "cal_force");
    return;
}

template <typename FPTYPE, typename Device>
void Forces<FPTYPE, Device>::cal_force_loc(const UnitCell& ucell,
                                           ModuleBase::matrix& forcelc,
                                           const ModulePW::PW_Basis* const rho_basis,
                                           const ModuleBase::matrix& vloc,
                                           const Charge* const chr)
{
    ModuleBase::TITLE("Forces", "cal_force_loc");
    ModuleBase::timer::tick("Forces", "cal_force_loc");
    this->device = base_device::get_device_type<Device>(this->ctx);
    std::complex<double>* aux = new std::complex<double>[rho_basis->nmaxgr];
    // now, in all pools , the charge are the same,
    // so, the force calculated by each pool is equal.

    /*
        blocking rho_basis->nrxx for data locality.

        By blocking aux with block size 1024,
        we can keep the blocked aux in L1 cache when iterating PARAM.inp.nspin loop
        performance will be better when number of atom is quite huge
    */
    const int block_ir = 1024;
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int irb = 0; irb < rho_basis->nrxx; irb += block_ir)
    {
        // calculate the actual task length of this block
        int ir_end = std::min(irb + block_ir, rho_basis->nrxx);

        { // is = 0
            for (int ir = irb; ir < ir_end; ++ir)
            { // initialize aux
                aux[ir] = std::complex<double>(chr->rho[0][ir], 0.0);
            }
        }
        if (PARAM.inp.nspin == 2)
        {
            for (int ir = irb; ir < ir_end; ++ir)
            { // accumulate aux
                aux[ir] += std::complex<double>(chr->rho[1][ir], 0.0);
            }
        }
    }

    // to G space. maybe need fftw with OpenMP
    rho_basis->real2recip(aux, aux);

    if(this->device == base_device::GpuDevice)
    {
        std::vector<double> tau_h;
        std::vector<double> gcar_h;
        tau_h.resize(this->nat * 3);
        for(int iat = 0; iat < this->nat; ++iat)
        {
            int it = ucell.iat2it[iat];
            int ia = ucell.iat2ia[iat];
            tau_h[iat * 3] = ucell.atoms[it].tau[ia].x;
            tau_h[iat * 3 + 1] = ucell.atoms[it].tau[ia].y;
            tau_h[iat * 3 + 2] = ucell.atoms[it].tau[ia].z;
        }

        gcar_h.resize(rho_basis->npw * 3);
        for(int ig = 0; ig < rho_basis->npw; ++ig)
        {
            gcar_h[ig * 3] = rho_basis->gcar[ig].x;
            gcar_h[ig * 3 + 1] = rho_basis->gcar[ig].y;
            gcar_h[ig * 3 + 2] = rho_basis->gcar[ig].z;
        }

        int* iat2it_d = nullptr;
        int* ig2gg_d = nullptr;
        double* gcar_d = nullptr;
        double* tau_d = nullptr;
        std::complex<double>* aux_d = nullptr;
        double* forcelc_d  = nullptr;
        double* vloc_d = nullptr;

        resmem_int_op()(iat2it_d, this->nat);
        resmem_int_op()(ig2gg_d, rho_basis->npw);
        resmem_var_op()(gcar_d, rho_basis->npw * 3);
        resmem_var_op()(tau_d, this->nat * 3);
        resmem_complex_op()(aux_d, rho_basis->npw);
        resmem_var_op()(forcelc_d, this->nat * 3);
        resmem_var_op()(vloc_d, vloc.nr * vloc.nc);

        syncmem_int_h2d_op()(iat2it_d, ucell.iat2it, this->nat);
        syncmem_int_h2d_op()(ig2gg_d, rho_basis->ig2igg, rho_basis->npw);
        syncmem_var_h2d_op()(gcar_d, gcar_h.data(), rho_basis->npw * 3);
        syncmem_var_h2d_op()(tau_d, tau_h.data(), this->nat * 3);
        syncmem_complex_h2d_op()(aux_d, aux, rho_basis->npw);
        syncmem_var_h2d_op()(forcelc_d, forcelc.c, this->nat * 3);
        syncmem_var_h2d_op()(vloc_d, vloc.c, vloc.nr * vloc.nc);

        hamilt::cal_force_loc_op<FPTYPE, Device>()(
            this->nat,
            rho_basis->npw,
            ucell.tpiba * ucell.omega,
            iat2it_d,
            ig2gg_d,
            gcar_d,
            tau_d,
            aux_d,
            vloc_d,
            vloc.nc,
            forcelc_d);
        syncmem_var_d2h_op()(forcelc.c, forcelc_d, this->nat * 3);

        delmem_int_op()(iat2it_d);
        delmem_int_op()(ig2gg_d);
        delmem_var_op()(gcar_d);
        delmem_var_op()(tau_d);
        delmem_complex_op()(aux_d);
        delmem_var_op()(forcelc_d);
        delmem_var_op()(vloc_d);
    }
    else{  // calculate forces on CPU
        #ifdef _OPENMP
        #pragma omp parallel for
        #endif
        for (int iat = 0; iat < this->nat; ++iat)
        {
            // read `it` `ia` from the table
            int it = ucell.iat2it[iat];
            int ia = ucell.iat2ia[iat];
            for (int ig = 0; ig < rho_basis->npw; ig++)
            {
                const double phase = ModuleBase::TWO_PI * (rho_basis->gcar[ig] * ucell.atoms[it].tau[ia]);
                double sinp, cosp;
                ModuleBase::libm::sincos(phase, &sinp, &cosp);
                const double factor
                    = vloc(it, rho_basis->ig2igg[ig]) * (cosp * aux[ig].imag() + sinp * aux[ig].real());
                forcelc(iat, 0) += rho_basis->gcar[ig][0] * factor;
                forcelc(iat, 1) += rho_basis->gcar[ig][1] * factor;
                forcelc(iat, 2) += rho_basis->gcar[ig][2] * factor;
            }
            forcelc(iat, 0) *= (ucell.tpiba * ucell.omega);
            forcelc(iat, 1) *= (ucell.tpiba * ucell.omega);
            forcelc(iat, 2) *= (ucell.tpiba * ucell.omega);
        }
    }
    // this->print(GlobalV::ofs_running, "local forces", forcelc);
    Parallel_Reduce::reduce_pool(forcelc.c, forcelc.nr * forcelc.nc);
    delete[] aux;
    ModuleBase::timer::tick("Forces", "cal_force_loc");
    return;
}

template <typename FPTYPE, typename Device>
void Forces<FPTYPE, Device>::cal_force_ew(const UnitCell& ucell,
                                          ModuleBase::matrix& forceion,
                                          const ModulePW::PW_Basis* const rho_basis,
                                          const Structure_Factor* p_sf)
{
    ModuleBase::TITLE("Forces", "cal_force_ew");
    ModuleBase::timer::tick("Forces", "cal_force_ew");
    this->device = base_device::get_device_type<Device>(this->ctx);
    double fact = 2.0;
    std::vector<std::complex<double>> aux(rho_basis->npw);

    /*
        blocking rho_basis->nrxnpwx for data locality.

        By blocking aux with block size 1024,
        we can keep the blocked aux in L1 cache when iterating ucell.ntype loop
        performance will be better when number of atom is quite huge
    */
    const int block_ig = 1024;
#pragma omp parallel for
    for (int igb = 0; igb < rho_basis->npw; igb += block_ig)
    {
        // calculate the actual task length of this block
        int ig_end = std::min(igb + block_ig, rho_basis->npw);

        for (int ig = igb; ig < ig_end; ++ig)
        {
            aux[ig] = 0.0;
        }
        for (int it = 0; it < ucell.ntype; it++)
        {
            if (ucell.atoms[it].na != 0)
            {
                double dzv;
                {
                    dzv = ucell.atoms[it].ncpp.zv;
                }

                for (int ig = igb; ig < ig_end; ++ig)
                { // accumulate aux
                    aux[ig] += dzv * conj(p_sf->strucFac(it, ig));
                }
            }
        }
    }

    // calculate total ionic charge
    double charge = 0.0;
    for (int it = 0; it < ucell.ntype; it++)
    {
        {
            charge += ucell.atoms[it].na * ucell.atoms[it].ncpp.zv; // mohan modify 2007-11-7
        }
    }

    double alpha = 1.1;
    double upperbound;
    do
    {
        alpha -= 0.10;
        // choose alpha in order to have convergence in the sum over G
        // upperbound is a safe upper bound for the error in the sum over G

        if (alpha <= 0.0)
        {
            ModuleBase::WARNING_QUIT("ewald", "Can't find optimal alpha.");
        }
        upperbound = 2.0 * charge * charge * sqrt(2.0 * alpha / ModuleBase::TWO_PI)
                     * erfc(sqrt(ucell.tpiba2 * rho_basis->ggecut / 4.0 / alpha));
    } while (upperbound > 1.0e-6);
    const int ig0 = rho_basis->ig_gge0;
#pragma omp parallel for
    for (int ig = 0; ig < rho_basis->npw; ig++)
    {
        if (ig== ig0)
        {
            continue; // skip G=0
        }
        aux[ig] *= ModuleBase::libm::exp(-1.0 * rho_basis->gg[ig] * ucell.tpiba2 / alpha / 4.0)
                / (rho_basis->gg[ig] * ucell.tpiba2);
    }

    // set pos rho_basis->ig_gge0 to zero
    if (rho_basis->ig_gge0 >= 0 && rho_basis->ig_gge0 < rho_basis->npw)
    {
        aux[rho_basis->ig_gge0] = std::complex<double>(0.0, 0.0);
    }
    if(this->device == base_device::GpuDevice)
    {
        std::vector<double> tau_h(this->nat * 3);
        std::vector<double> gcar_h(rho_basis->npw * 3);
        for(int iat = 0; iat < this->nat; ++iat)
        {
            int it = ucell.iat2it[iat];
            int ia = ucell.iat2ia[iat];
            tau_h[iat * 3] = ucell.atoms[it].tau[ia].x;
            tau_h[iat * 3 + 1] = ucell.atoms[it].tau[ia].y;
            tau_h[iat * 3 + 2] = ucell.atoms[it].tau[ia].z;
        }
        for(int ig = 0; ig < rho_basis->npw; ++ig)
        {
            gcar_h[ig * 3] = rho_basis->gcar[ig].x;
            gcar_h[ig * 3 + 1] = rho_basis->gcar[ig].y;
            gcar_h[ig * 3 + 2] = rho_basis->gcar[ig].z;
        }
        std::vector<double> it_fact_h(ucell.ntype);
        for(int it = 0; it < ucell.ntype; ++it)
        {
            it_fact_h[it] = ucell.atoms[it].ncpp.zv * ModuleBase::e2 * ucell.tpiba * ModuleBase::TWO_PI / ucell.omega * fact;
        }

        int* iat2it_d = nullptr;
        double* gcar_d = nullptr;
        double* tau_d = nullptr;
        double* it_fact_d = nullptr;
        std::complex<double>* aux_d = nullptr;
        double* forceion_d  = nullptr;
        resmem_int_op()(iat2it_d, this->nat);
        resmem_var_op()(gcar_d, rho_basis->npw * 3);
        resmem_var_op()(tau_d, this->nat * 3);
        resmem_var_op()(it_fact_d, ucell.ntype);
        resmem_complex_op()(aux_d, rho_basis->npw);
        resmem_var_op()(forceion_d, this->nat * 3);

        syncmem_int_h2d_op()(iat2it_d, ucell.iat2it, this->nat);
        syncmem_var_h2d_op()(gcar_d, gcar_h.data(), rho_basis->npw * 3);
        syncmem_var_h2d_op()(tau_d, tau_h.data(), this->nat * 3);
        syncmem_var_h2d_op()(it_fact_d, it_fact_h.data(), ucell.ntype);
        syncmem_complex_h2d_op()(aux_d, aux.data(), rho_basis->npw);
        syncmem_var_h2d_op()(forceion_d, forceion.c, this->nat * 3);

        hamilt::cal_force_ew_op<FPTYPE, Device>()(
            this->nat,
            rho_basis->npw,
            rho_basis->ig_gge0,
            iat2it_d,
            gcar_d,
            tau_d,
            it_fact_d,
            aux_d,
            forceion_d);
        
        syncmem_var_d2h_op()(forceion.c, forceion_d, this->nat * 3);
        delmem_int_op()(iat2it_d);
        delmem_var_op()(gcar_d);
        delmem_var_op()(tau_d);
        delmem_var_op()(it_fact_d);
        delmem_complex_op()(aux_d);
        delmem_var_op()(forceion_d);
    } else // calculate forces on CPU
    {
    #pragma omp parallel for
        for(int iat = 0; iat < this->nat; ++iat)
        {
            const int it = ucell.iat2it[iat];
            const int ia = ucell.iat2ia[iat];
            double it_fact = ucell.atoms[it].ncpp.zv * ModuleBase::e2 * ucell.tpiba * ModuleBase::TWO_PI / ucell.omega * fact;

            for(int ig = 0; ig < rho_basis->npw; ++ig)
            {
                if(ig != rho_basis->ig_gge0) // skip G=0
                {
                    const ModuleBase::Vector3<double> gcar = rho_basis->gcar[ig];
                    const double arg = ModuleBase::TWO_PI * (gcar * ucell.atoms[it].tau[ia]);
                    double sinp, cosp;
                    ModuleBase::libm::sincos(arg, &sinp, &cosp);
                    double sumnb = -cosp * aux[ig].imag() + sinp * aux[ig].real();
                    forceion(iat, 0) += gcar[0] * sumnb;
                    forceion(iat, 1) += gcar[1] * sumnb;
                    forceion(iat, 2) += gcar[2] * sumnb;
                }
            }
            forceion(iat, 0) *= it_fact;
            forceion(iat, 1) *= it_fact;
            forceion(iat, 2) *= it_fact;
        }
    }
    // means that the processor contains G=0 term.
    #pragma omp parallel
    {
        if (rho_basis->ig_gge0 >= 0)
        {
            double rmax = 5.0 / (sqrt(alpha) * ucell.lat0);
            int nrm = 0;

            // output of rgen: the number of vectors in the sphere
            const int mxr = 200;
            // the maximum number of R vectors included in r
            std::vector<ModuleBase::Vector3<double>> r(mxr);
            std::vector<double> r2(mxr);
            std::vector<int> irr(mxr);
            // the square modulus of R_j-tau_s-tau_s'

            const double sqa = sqrt(alpha);
            const double sq8a_2pi = sqrt(8.0 * alpha / ModuleBase::TWO_PI);

            // iterating atoms.
            #pragma omp for
            for(int iat1 = 0; iat1 < this->nat; iat1++)
            {
                int T1 = ucell.iat2it[iat1];
                int I1 = ucell.iat2ia[iat1];
                for(int iat2 = 0; iat2 < this->nat; iat2++)
                {
                    int T2 = ucell.iat2it[iat2];
                    int I2 = ucell.iat2ia[iat2];
                    if (iat1 != iat2)
                    {
                        ModuleBase::Vector3<double> d_tau
                            = ucell.atoms[T1].tau[I1] - ucell.atoms[T2].tau[I2];
                        H_Ewald_pw::rgen(d_tau, rmax, irr.data(), ucell.latvec, ucell.G, r.data(), r2.data(), nrm);

                        for (int n = 0; n < nrm; n++)
                        {
                            const double rr = sqrt(r2[n]) * ucell.lat0;

                            double factor;
                            {
                                factor = ucell.atoms[T1].ncpp.zv * ucell.atoms[T2].ncpp.zv
                                            * ModuleBase::e2 / (rr * rr)
                                            * (erfc(sqa * rr) / rr + sq8a_2pi * ModuleBase::libm::exp(-alpha * rr * rr))
                                            * ucell.lat0;
                            }
                            forceion(iat1, 0) -= factor * r[n].x;
                            forceion(iat1, 1) -= factor * r[n].y;
                            forceion(iat1, 2) -= factor * r[n].z;
                        }
                    }
                } // atom b
            } // atom a
        }
    }
    Parallel_Reduce::reduce_pool(forceion.c, forceion.nr * forceion.nc);
    // this->print(GlobalV::ofs_running, "ewald forces", forceion);

    ModuleBase::timer::tick("Forces", "cal_force_ew");

    return;
}



template class Forces<double, base_device::DEVICE_CPU>;
#if ((defined __CUDA) || (defined __ROCM))
template class Forces<double, base_device::DEVICE_GPU>;
#endif
