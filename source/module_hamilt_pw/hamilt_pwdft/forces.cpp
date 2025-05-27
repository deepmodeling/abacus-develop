#include "forces.h"

#include "module_parameter/parameter.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_io/output_log.h"
// new
#include "module_base/complexmatrix.h"
#include "module_base/libm/libm.h"
#include "module_base/math_integral.h"
#include "module_base/mathzone.h"
#include "module_base/timer.h"
#include "module_base/tool_threading.h"
#include "module_elecstate/potentials/efield.h"
#include "module_elecstate/potentials/gatefield.h"
#include "module_hamilt_general/module_ewald/H_Ewald_pw.h"
#include "module_hamilt_general/module_surchem/surchem.h"
#include "module_hamilt_general/module_vdw/vdw.h"
#include "module_hamilt_pw/hamilt_pwdft/kernels/force_op.h"
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif
#ifdef USE_PAW
#include "module_cell/module_paw/paw_cell.h"
#endif

// Data preparation function for Ewald parallel operator
template <typename FPTYPE>
void prepare_ew_parallel_data(const UnitCell& ucell,
                              ModulePW::PW_Basis* rho_basis,
                              const std::complex<double>* aux,
                              const double alpha,
                              const double fact,
                              std::vector<FPTYPE>& gcar_flat,
                              std::vector<FPTYPE>& tau_flat,
                              std::vector<int>& iat2it_array,
                              std::vector<FPTYPE>& it_facts,
                              std::vector<int>& atom_na,
                              std::vector<FPTYPE>& zv_values,
                              std::vector<FPTYPE>& latvec_flat,
                              std::vector<FPTYPE>& G_flat)
{
    const int nat = ucell.nat;
    const int npw = rho_basis->npw;
    const int ntype = ucell.ntype;
    
    // Prepare gcar_flat
    gcar_flat.resize(npw * 3);
    for (int ig = 0; ig < npw; ig++)
    {
        gcar_flat[ig * 3 + 0] = static_cast<FPTYPE>(rho_basis->gcar[ig][0]);
        gcar_flat[ig * 3 + 1] = static_cast<FPTYPE>(rho_basis->gcar[ig][1]);
        gcar_flat[ig * 3 + 2] = static_cast<FPTYPE>(rho_basis->gcar[ig][2]);
    }
    
    // Prepare tau_flat and iat2it_array
    tau_flat.resize(nat * 3);
    iat2it_array.resize(nat);
    int iat = 0;
    for (int it = 0; it < ntype; it++)
    {
        for (int ia = 0; ia < ucell.atoms[it].na; ia++)
        {
            tau_flat[iat * 3 + 0] = static_cast<FPTYPE>(ucell.atoms[it].tau[ia].x);
            tau_flat[iat * 3 + 1] = static_cast<FPTYPE>(ucell.atoms[it].tau[ia].y);
            tau_flat[iat * 3 + 2] = static_cast<FPTYPE>(ucell.atoms[it].tau[ia].z);
            iat2it_array[iat] = it;
            iat++;
        }
    }
    
    // Prepare it_facts
    it_facts.resize(nat);
    iat = 0;
    for (int it = 0; it < ntype; it++)
    {
        double zv;
        if (PARAM.inp.use_paw)
        {
#ifdef USE_PAW
            zv = GlobalC::paw_cell.get_val(it);
#else
            zv = ucell.atoms[it].ncpp.zv;
#endif
        }
        else
        {
            zv = ucell.atoms[it].ncpp.zv;
        }
        
        const FPTYPE it_fact = static_cast<FPTYPE>(zv * ModuleBase::e2 * ucell.tpiba * ModuleBase::TWO_PI / ucell.omega * fact);
        
        for (int ia = 0; ia < ucell.atoms[it].na; ia++)
        {
            it_facts[iat] = it_fact;
            iat++;
        }
    }
    
    // Prepare atom_na
    atom_na.resize(ntype);
    for (int it = 0; it < ntype; it++)
    {
        atom_na[it] = ucell.atoms[it].na;
    }
    
    // Prepare zv_values
    zv_values.resize(ntype);
    for (int it = 0; it < ntype; it++)
    {
        if (PARAM.inp.use_paw)
        {
#ifdef USE_PAW
            zv_values[it] = static_cast<FPTYPE>(GlobalC::paw_cell.get_val(it));
#else
            zv_values[it] = static_cast<FPTYPE>(ucell.atoms[it].ncpp.zv);
#endif
        }
        else
        {
            zv_values[it] = static_cast<FPTYPE>(ucell.atoms[it].ncpp.zv);
        }
    }
    
    // Prepare latvec_flat (row-major order)
    latvec_flat.resize(9);
    latvec_flat[0] = static_cast<FPTYPE>(ucell.latvec.e11);
    latvec_flat[1] = static_cast<FPTYPE>(ucell.latvec.e12);
    latvec_flat[2] = static_cast<FPTYPE>(ucell.latvec.e13);
    latvec_flat[3] = static_cast<FPTYPE>(ucell.latvec.e21);
    latvec_flat[4] = static_cast<FPTYPE>(ucell.latvec.e22);
    latvec_flat[5] = static_cast<FPTYPE>(ucell.latvec.e23);
    latvec_flat[6] = static_cast<FPTYPE>(ucell.latvec.e31);
    latvec_flat[7] = static_cast<FPTYPE>(ucell.latvec.e32);
    latvec_flat[8] = static_cast<FPTYPE>(ucell.latvec.e33);
    
    // Prepare G_flat (row-major order)
    G_flat.resize(9);
    G_flat[0] = static_cast<FPTYPE>(ucell.G.e11);
    G_flat[1] = static_cast<FPTYPE>(ucell.G.e12);
    G_flat[2] = static_cast<FPTYPE>(ucell.G.e13);
    G_flat[3] = static_cast<FPTYPE>(ucell.G.e21);
    G_flat[4] = static_cast<FPTYPE>(ucell.G.e22);
    G_flat[5] = static_cast<FPTYPE>(ucell.G.e23);
    G_flat[6] = static_cast<FPTYPE>(ucell.G.e31);
    G_flat[7] = static_cast<FPTYPE>(ucell.G.e32);
    G_flat[8] = static_cast<FPTYPE>(ucell.G.e33);
}

// Data preparation function for local force sincos operator
template <typename FPTYPE>
void prepare_loc_sincos_data(const UnitCell& ucell,
                             ModulePW::PW_Basis* rho_basis,
                             const ModuleBase::matrix& vloc,
                             const std::complex<double>* aux,
                             std::vector<FPTYPE>& gcar_flat,
                             std::vector<FPTYPE>& tau_flat,
                             std::vector<int>& iat2it_array,
                             std::vector<FPTYPE>& vloc_factors)
{
    const int nat = ucell.nat;
    const int npw = rho_basis->npw;
    
    // Prepare gcar_flat
    gcar_flat.resize(npw * 3);
    for (int ig = 0; ig < npw; ig++)
    {
        gcar_flat[ig * 3 + 0] = static_cast<FPTYPE>(rho_basis->gcar[ig][0]);
        gcar_flat[ig * 3 + 1] = static_cast<FPTYPE>(rho_basis->gcar[ig][1]);
        gcar_flat[ig * 3 + 2] = static_cast<FPTYPE>(rho_basis->gcar[ig][2]);
    }
    
    // Prepare tau_flat and iat2it_array
    tau_flat.resize(nat * 3);
    iat2it_array.resize(nat);
    int iat = 0;
    for (int it = 0; it < ucell.ntype; it++)
    {
        for (int ia = 0; ia < ucell.atoms[it].na; ia++)
        {
            tau_flat[iat * 3 + 0] = static_cast<FPTYPE>(ucell.atoms[it].tau[ia].x);
            tau_flat[iat * 3 + 1] = static_cast<FPTYPE>(ucell.atoms[it].tau[ia].y);
            tau_flat[iat * 3 + 2] = static_cast<FPTYPE>(ucell.atoms[it].tau[ia].z);
            iat2it_array[iat] = it;
            iat++;
        }
    }
    
    // Prepare vloc_factors: vloc(it, ig2igg[ig]) for each G-vector and atom type
    vloc_factors.resize(npw * ucell.ntype);
    for (int it = 0; it < ucell.ntype; it++)
    {
        for (int ig = 0; ig < npw; ig++)
        {
            vloc_factors[it * npw + ig] = static_cast<FPTYPE>(vloc(it, rho_basis->ig2igg[ig]));
        }
    }
}

// Data copy back function for local force sincos operator
template <typename FPTYPE>
void copy_back_loc_force_data(const FPTYPE* force,
                              const int nat,
                              const FPTYPE& scale_factor,
                              ModuleBase::matrix& forcelc)
{
    for (int iat = 0; iat < nat; iat++)
    {
        forcelc(iat, 0) = static_cast<double>(force[iat * 3 + 0] * scale_factor);
        forcelc(iat, 1) = static_cast<double>(force[iat * 3 + 1] * scale_factor);
        forcelc(iat, 2) = static_cast<double>(force[iat * 3 + 2] * scale_factor);
    }
}

// Data copy back function for Ewald parallel operator
template <typename FPTYPE>
void copy_back_ew_force_data(const FPTYPE* force,
                             const int nat,
                             ModuleBase::matrix& forceion)
{
    for (int iat = 0; iat < nat; iat++)
    {
        forceion(iat, 0) = static_cast<double>(force[iat * 3 + 0]);
        forceion(iat, 1) = static_cast<double>(force[iat * 3 + 1]);
        forceion(iat, 2) = static_cast<double>(force[iat * 3 + 2]);
    }
}

template <typename FPTYPE, typename Device>
void Forces<FPTYPE, Device>::cal_force(UnitCell& ucell,
                                       ModuleBase::matrix& force,
                                       const elecstate::ElecState& elec,
                                       ModulePW::PW_Basis* rho_basis,
                                       ModuleSymmetry::Symmetry* p_symm,
                                       Structure_Factor* p_sf,
                                       surchem& solvent,
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
    // For PAW, calculated together in paw_cell.calculate_force
    if (!PARAM.inp.use_paw)
    {
        this->cal_force_loc(ucell,forcelc, rho_basis, locpp->vloc, chr);
    }
    else
    {
        forcelc.zero_out();
    }

    // Ewald
    this->cal_force_ew(ucell,forceion, rho_basis, p_sf);

    // Force due to nonlocal part of pseudopotential
    if (wfc_basis != nullptr)
    {
        if (!PARAM.inp.use_paw)
        {
            this->npwx = wfc_basis->npwk_max;
            Forces::cal_force_nl(forcenl, wg, ekb, pkv, wfc_basis, p_sf, *p_nlpp, ucell, psi_in);

            if (PARAM.globalv.use_uspp)
            {
                this->cal_force_us(forcenl, rho_basis, *p_nlpp, elec, ucell);
            }
        }
        else
        {
#ifdef USE_PAW
            for (int ik = 0; ik < wfc_basis->nks; ik++)
            {
                const int npw = wfc_basis->npwk[ik];
                ModuleBase::Vector3<double>* _gk = new ModuleBase::Vector3<double>[npw];
                for (int ig = 0; ig < npw; ig++)
                {
                    _gk[ig] = wfc_basis->getgpluskcar(ik, ig);
                }

                double* kpt;
                kpt = new double[3];
                kpt[0] = wfc_basis->kvec_c[ik].x;
                kpt[1] = wfc_basis->kvec_c[ik].y;
                kpt[2] = wfc_basis->kvec_c[ik].z;

                double** kpg;
                double** gcar;
                kpg = new double*[npw];
                gcar = new double*[npw];
                for (int ipw = 0; ipw < npw; ipw++)
                {
                    kpg[ipw] = new double[3];
                    kpg[ipw][0] = _gk[ipw].x;
                    kpg[ipw][1] = _gk[ipw].y;
                    kpg[ipw][2] = _gk[ipw].z;

                    gcar[ipw] = new double[3];
                    gcar[ipw][0] = wfc_basis->getgcar(ik, ipw).x;
                    gcar[ipw][1] = wfc_basis->getgcar(ik, ipw).y;
                    gcar[ipw][2] = wfc_basis->getgcar(ik, ipw).z;
                }

                GlobalC::paw_cell.set_paw_k(npw,
                                            wfc_basis->npwk_max,
                                            kpt,
                                            wfc_basis->get_ig2ix(ik).data(),
                                            wfc_basis->get_ig2iy(ik).data(),
                                            wfc_basis->get_ig2iz(ik).data(),
                                            (const double**)kpg,
                                            ucell.tpiba,
                                            (const double**)gcar);

                delete[] kpt;
                for (int ipw = 0; ipw < npw; ipw++)
                {
                    delete[] kpg[ipw];
                    delete[] gcar[ipw];
                }
                delete[] kpg;
                delete[] gcar;

                GlobalC::paw_cell.get_vkb();

                GlobalC::paw_cell.set_currentk(ik);

                psi_in[0].fix_k(ik);
                double *weight, *epsilon;
                weight = new double[PARAM.inp.nbands];
                epsilon = new double[PARAM.inp.nbands];
                for (int ib = 0; ib < PARAM.inp.nbands; ib++)
                {
                    weight[ib] = wg(ik, ib);
                    epsilon[ib] = ekb(ik, ib);
                }
                GlobalC::paw_cell.paw_nl_force(reinterpret_cast<std::complex<double>*>(psi_in[0].get_pointer()),
                                               epsilon,
                                               weight,
                                               PARAM.inp.nbands,
                                               forcenl.c);

                delete[] weight;
                delete[] epsilon;
            }
#endif
        }
        // DFT+U and DeltaSpin
        if(PARAM.inp.dft_plus_u || PARAM.inp.sc_mag_switch)
        {
            this->cal_force_onsite(forceonsite, wg, wfc_basis, ucell, psi_in);
        }
    }

    // non-linear core correction
    // not relevant for PAW
    if (!PARAM.inp.use_paw)
    {
        Forces::cal_force_cc(forcecc, rho_basis, chr, locpp->numeric, ucell);
    }
    else
    {
        forcecc.zero_out();
    }

    // force due to core charge
    // For PAW, calculated together in paw_cell.calculate_force
    if (!PARAM.inp.use_paw)
    {
        this->cal_force_scc(forcescc, rho_basis, elec.vnew, elec.vnew_exist, locpp->numeric, ucell);
    }
    else
    {
        forcescc.zero_out();
    }

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

#ifdef USE_PAW
    if (PARAM.inp.use_paw)
    {
        double* force_paw;
        double* rhor;
        rhor = new double[rho_basis->nrxx];
        for (int ir = 0; ir < rho_basis->nrxx; ir++)
        {
            rhor[ir] = 0.0;
        }
        for (int is = 0; is < PARAM.inp.nspin; is++)
        {
            for (int ir = 0; ir < rho_basis->nrxx; ir++)
            {
                rhor[ir] += chr->rho[is][ir] + chr->nhat[is][ir];
            }
        }

        force_paw = new double[3 * this->nat];
        ModuleBase::matrix v_xc, v_effective;
        v_effective.create(PARAM.inp.nspin, rho_basis->nrxx);
        v_effective.zero_out();
        elec.pot->update_from_charge(elec.charge, &ucell);
        v_effective = elec.pot->get_effective_v();

        v_xc.create(PARAM.inp.nspin, rho_basis->nrxx);
        v_xc.zero_out();
        const std::tuple<double, double, ModuleBase::matrix> etxc_vtxc_v
            = XC_Functional::v_xc(rho_basis->nrxx, elec.charge, &ucell);
        v_xc = std::get<2>(etxc_vtxc_v);

        GlobalC::paw_cell.calculate_force(v_effective.c, v_xc.c, rhor, force_paw);

        for (int iat = 0; iat < this->nat; iat++)
        {
            // Ha to Ry
            forcepaw(iat, 0) = force_paw[3 * iat] * 2.0;
            forcepaw(iat, 1) = force_paw[3 * iat + 1] * 2.0;
            forcepaw(iat, 2) = force_paw[3 * iat + 2] * 2.0;
        }

        delete[] force_paw;
        delete[] rhor;
    }
#endif

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

                if (PARAM.inp.use_paw)
                {
                    force(iat, ipol) += forcepaw(iat, ipol);
                }

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
        if (PARAM.inp.use_paw)
        {
            ModuleIO::print_force(GlobalV::ofs_running,
                                  ucell,
                                  "PAW      FORCE (eV/Angstrom)",
                                  forcepaw,
                                  false);
        }
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
                                           ModulePW::PW_Basis* rho_basis,
                                           const ModuleBase::matrix& vloc,
                                           const Charge* const chr)
{
    ModuleBase::TITLE("Forces", "cal_force_loc");
    ModuleBase::timer::tick("Forces", "cal_force_loc");

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

    // Prepare data for the sincos operator
    std::vector<FPTYPE> gcar_flat, tau_flat, vloc_factors;
    std::vector<int> iat2it_array;
    prepare_loc_sincos_data<FPTYPE>(ucell, rho_basis, vloc, aux,
                                    gcar_flat, tau_flat, iat2it_array, vloc_factors);

    // Prepare force array for the operator
    std::vector<FPTYPE> force_flat(this->nat * 3, 0.0);

    // Convert aux to FPTYPE complex array
    std::vector<std::complex<FPTYPE>> aux_fptype(rho_basis->npw);
    for (int ig = 0; ig < rho_basis->npw; ig++)
    {
        aux_fptype[ig] = std::complex<FPTYPE>(static_cast<FPTYPE>(aux[ig].real()),
                                              static_cast<FPTYPE>(aux[ig].imag()));
    }

    // Call the sincos operator
    hamilt::cal_force_loc_sincos_op<FPTYPE, base_device::DEVICE_CPU> sincos_op;
    sincos_op(this->ctx,
              this->nat,
              rho_basis->npw,
              ucell.ntype,
              gcar_flat.data(),
              tau_flat.data(),
              iat2it_array.data(),
              vloc_factors.data(),
              aux_fptype.data(),
              force_flat.data());

    // Copy back the results with scaling
    const FPTYPE scale_factor = static_cast<FPTYPE>(ucell.tpiba * ucell.omega);
    copy_back_loc_force_data<FPTYPE>(force_flat.data(), this->nat, scale_factor, forcelc);

    // this->print(GlobalV::ofs_running, "local forces", forcelc);
    Parallel_Reduce::reduce_pool(forcelc.c, forcelc.nr * forcelc.nc);
    delete[] aux;
    ModuleBase::timer::tick("Forces", "cal_force_loc");
    return;
}

template <typename FPTYPE, typename Device>
void Forces<FPTYPE, Device>::cal_force_ew(const UnitCell& ucell,
                                          ModuleBase::matrix& forceion,
                                          ModulePW::PW_Basis* rho_basis,
                                          const Structure_Factor* p_sf)
{
    ModuleBase::TITLE("Forces", "cal_force_ew");
    ModuleBase::timer::tick("Forces", "cal_force_ew");

    double fact = 2.0;
    std::complex<double>* aux = new std::complex<double>[rho_basis->npw];

    /*
        blocking rho_basis->nrxnpwx for data locality.

        By blocking aux with block size 1024,
        we can keep the blocked aux in L1 cache when iterating ucell.ntype loop
        performance will be better when number of atom is quite huge
    */
    const int block_ig = 1024;
#ifdef _OPENMP
#pragma omp parallel for
#endif
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
                if (PARAM.inp.use_paw)
                {
    #ifdef USE_PAW
                    dzv = GlobalC::paw_cell.get_val(it);
    #endif
                }
                else
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
        if (PARAM.inp.use_paw)
        {
#ifdef USE_PAW
            charge += ucell.atoms[it].na * GlobalC::paw_cell.get_val(it);
#endif
        }
        else
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

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int ig = 0; ig < rho_basis->npw; ig++)
    {
        aux[ig] *= ModuleBase::libm::exp(-1.0 * rho_basis->gg[ig] * ucell.tpiba2 / alpha / 4.0)
                   / (rho_basis->gg[ig] * ucell.tpiba2);
    }

    // set pos rho_basis->ig_gge0 to zero
    if (rho_basis->ig_gge0 >= 0 && rho_basis->ig_gge0 < rho_basis->npw)
    {
        aux[rho_basis->ig_gge0] = std::complex<double>(0.0, 0.0);
    }

    // Prepare data for the parallel operator
    std::vector<FPTYPE> gcar_flat, tau_flat, it_facts, zv_values, latvec_flat, G_flat;
    std::vector<int> iat2it_array, atom_na;
    prepare_ew_parallel_data<FPTYPE>(ucell, rho_basis, aux, alpha, fact,
                                     gcar_flat, tau_flat, iat2it_array, it_facts,
                                     atom_na, zv_values, latvec_flat, G_flat);

    // Prepare force array for the operator
    std::vector<FPTYPE> force_flat(this->nat * 3, 0.0);

    // Convert aux to FPTYPE complex array
    std::vector<std::complex<FPTYPE>> aux_fptype(rho_basis->npw);
    for (int ig = 0; ig < rho_basis->npw; ig++)
    {
        aux_fptype[ig] = std::complex<FPTYPE>(static_cast<FPTYPE>(aux[ig].real()),
                                              static_cast<FPTYPE>(aux[ig].imag()));
    }

    // Call the parallel operator
    hamilt::cal_force_ew_parallel_op<FPTYPE, base_device::DEVICE_CPU> parallel_op;
    parallel_op(this->ctx,
                this->nat,
                rho_basis->npw,
                ucell.ntype,
                rho_basis->ig_gge0,
                gcar_flat.data(),
                tau_flat.data(),
                iat2it_array.data(),
                it_facts.data(),
                atom_na.data(),
                zv_values.data(),
                aux_fptype.data(),
                static_cast<FPTYPE>(alpha),
                static_cast<FPTYPE>(fact),
                static_cast<FPTYPE>(ucell.lat0),
                latvec_flat.data(),
                G_flat.data(),
                PARAM.inp.use_paw,
                force_flat.data());

    // Copy back the results
    copy_back_ew_force_data<FPTYPE>(force_flat.data(), this->nat, forceion);

    Parallel_Reduce::reduce_pool(forceion.c, forceion.nr * forceion.nc);

    // this->print(GlobalV::ofs_running, "ewald forces", forceion);

    ModuleBase::timer::tick("Forces", "cal_force_ew");

    delete[] aux;

    return;
}



template class Forces<double, base_device::DEVICE_CPU>;
#if ((defined __CUDA) || (defined __ROCM))
template class Forces<double, base_device::DEVICE_GPU>;
#endif
