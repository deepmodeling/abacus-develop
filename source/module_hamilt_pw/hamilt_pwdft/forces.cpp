#include "forces.h"

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
#include "module_psi/kernels/device.h"
#ifdef _OPENMP
#include <omp.h>
#endif

template <typename FPTYPE, typename Device>
void Forces<FPTYPE, Device>::cal_force(ModuleBase::matrix& force,
                                       const elecstate::ElecState& elec,
                                       ModulePW::PW_Basis* rho_basis,
                                       ModuleSymmetry::Symmetry* p_symm,
                                       Structure_Factor* p_sf,
                                       K_Vectors* pkv,
                                       ModulePW::PW_Basis_K* wfc_basis,
                                       const psi::Psi<std::complex<FPTYPE>, Device>* psi_in)
{
    ModuleBase::TITLE("Forces", "init");
    this->device = psi::device::get_device_type<Device>(this->ctx);
    const ModuleBase::matrix& wg = elec.wg;
    const ModuleBase::matrix& ekb = elec.ekb;
    const Charge* const chr = elec.charge;
    force.create(nat, 3);

    ModuleBase::matrix forcelc(nat, 3);
    ModuleBase::matrix forceion(nat, 3);
    ModuleBase::matrix forcecc(nat, 3);
    ModuleBase::matrix forcenl(nat, 3);
    ModuleBase::matrix forcescc(nat, 3);
    this->cal_force_loc(forcelc, rho_basis, chr);
    this->cal_force_ew(forceion, rho_basis, p_sf);
    if(wfc_basis != nullptr)
    {
        this->npwx = wfc_basis->npwk_max;
        this->cal_force_nl(forcenl, wg, ekb, pkv, wfc_basis, psi_in);
        if (GlobalV::use_uspp)
        {
            this->cal_force_us(forcenl, rho_basis, &GlobalC::ppcell, elec, GlobalC::ucell);
        }
    }
    this->cal_force_cc(forcecc, rho_basis, chr);
    this->cal_force_scc(forcescc, rho_basis, elec.vnew, elec.vnew_exist);

    ModuleBase::matrix stress_vdw_pw; //.create(3,3);
    ModuleBase::matrix force_vdw;
    force_vdw.create(nat, 3);
    auto vdw_solver = vdw::make_vdw(GlobalC::ucell, INPUT);
    if (vdw_solver != nullptr)
    {
        const std::vector<ModuleBase::Vector3<double>>& force_vdw_temp = vdw_solver->get_force();
        for (int iat = 0; iat < this->nat; ++iat)
        {
            force_vdw(iat, 0) = force_vdw_temp[iat].x;
            force_vdw(iat, 1) = force_vdw_temp[iat].y;
            force_vdw(iat, 2) = force_vdw_temp[iat].z;
        }
        if (GlobalV::TEST_FORCE)
        {
            ModuleIO::print_force(GlobalV::ofs_running, GlobalC::ucell, "VDW      FORCE (Ry/Bohr)", force_vdw);
        }
    }

    ModuleBase::matrix force_e;
    if (GlobalV::EFIELD_FLAG)
    {
        force_e.create(this->nat, 3);
        elecstate::Efield::compute_force(GlobalC::ucell, force_e);
        if (GlobalV::TEST_FORCE)
        {
            ModuleIO::print_force(GlobalV::ofs_running, GlobalC::ucell, "EFIELD      FORCE (Ry/Bohr)", force_e);
        }
    }

    ModuleBase::matrix force_gate;
    if (GlobalV::GATE_FLAG)
    {
        force_gate.create(this->nat, 3);
        elecstate::Gatefield::compute_force(GlobalC::ucell, force_gate);
        if (GlobalV::TEST_FORCE)
        {
            ModuleIO::print_force(GlobalV::ofs_running, GlobalC::ucell, "GATEFIELD      FORCE (Ry/Bohr)", force_gate);
        }
    }

    ModuleBase::matrix forcesol;
    if (GlobalV::imp_sol)
    {
        forcesol.create(this->nat, 3);
        GlobalC::solvent_model.cal_force_sol(GlobalC::ucell, rho_basis, forcesol);
        if (GlobalV::TEST_FORCE)
        {
            ModuleIO::print_force(GlobalV::ofs_running, GlobalC::ucell, "IMP_SOL      FORCE (Ry/Bohr)", forcesol);
        }
    }

    // impose total force = 0
    int iat = 0;
    for (int ipol = 0; ipol < 3; ipol++)
    {
        double sum = 0.0;
        iat = 0;

        for (int it = 0; it < GlobalC::ucell.ntype; it++)
        {
            for (int ia = 0; ia < GlobalC::ucell.atoms[it].na; ia++)
            {
                force(iat, ipol) = forcelc(iat, ipol) + forceion(iat, ipol) + forcenl(iat, ipol) + forcecc(iat, ipol)
                                   + forcescc(iat, ipol);

                if (vdw_solver != nullptr) // linpz and jiyy added vdw force, modified by zhengdy
                {
                    force(iat, ipol) += force_vdw(iat, ipol);
                }

                if (GlobalV::EFIELD_FLAG)
                {
                    force(iat, ipol) = force(iat, ipol) + force_e(iat, ipol);
                }

                if (GlobalV::GATE_FLAG)
                {
                    force(iat, ipol) = force(iat, ipol) + force_gate(iat, ipol);
                }

                if (GlobalV::imp_sol)
                {
                    force(iat, ipol) = force(iat, ipol) + forcesol(iat, ipol);
                }

                sum += force(iat, ipol);

                iat++;
            }
        }

        if (!(GlobalV::GATE_FLAG || GlobalV::EFIELD_FLAG))
        {
            double compen = sum / this->nat;
            for (int iat = 0; iat < this->nat; ++iat)
            {
                force(iat, ipol) = force(iat, ipol) - compen;
            }
        }
    }

    if (GlobalV::GATE_FLAG || GlobalV::EFIELD_FLAG)
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
                                                      GlobalC::ucell.a1.x,
                                                      GlobalC::ucell.a1.y,
                                                      GlobalC::ucell.a1.z,
                                                      GlobalC::ucell.a2.x,
                                                      GlobalC::ucell.a2.y,
                                                      GlobalC::ucell.a2.z,
                                                      GlobalC::ucell.a3.x,
                                                      GlobalC::ucell.a3.y,
                                                      GlobalC::ucell.a3.z,
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
                                                      GlobalC::ucell.a1.x,
                                                      GlobalC::ucell.a1.y,
                                                      GlobalC::ucell.a1.z,
                                                      GlobalC::ucell.a2.x,
                                                      GlobalC::ucell.a2.y,
                                                      GlobalC::ucell.a2.z,
                                                      GlobalC::ucell.a3.x,
                                                      GlobalC::ucell.a3.y,
                                                      GlobalC::ucell.a3.z,
                                                      d1,
                                                      d2,
                                                      d3);
            force(iat, 0) = d1;
            force(iat, 1) = d2;
            force(iat, 2) = d3;
        }
    }

    GlobalV::ofs_running << std::setiosflags(std::ios::fixed) << std::setprecision(6) << std::endl;
    /*if(GlobalV::TEST_FORCE)
    {
        ModuleIO::print_force(GlobalV::ofs_running, GlobalC::ucell,"LOCAL    FORCE (Ry/Bohr)", forcelc);
        ModuleIO::print_force(GlobalV::ofs_running, GlobalC::ucell,"NONLOCAL FORCE (Ry/Bohr)", forcenl);
        ModuleIO::print_force(GlobalV::ofs_running, GlobalC::ucell,"NLCC     FORCE (Ry/Bohr)", forcecc);
        ModuleIO::print_force(GlobalV::ofs_running, GlobalC::ucell,"ION      FORCE (Ry/Bohr)", forceion);
        ModuleIO::print_force(GlobalV::ofs_running, GlobalC::ucell,"SCC      FORCE (Ry/Bohr)", forcescc);
        if(GlobalV::EFIELD) ModuleIO::print_force(GlobalV::ofs_running, GlobalC::ucell,"EFIELD   FORCE (Ry/Bohr)",
    force_e);
    }*/

    /*
        ModuleIO::print_force(GlobalV::ofs_running, GlobalC::ucell,"   TOTAL-FORCE (Ry/Bohr)", force);

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

    if (GlobalV::TEST_FORCE)
    {
        ModuleIO::print_force(GlobalV::ofs_running, GlobalC::ucell, "LOCAL    FORCE (eV/Angstrom)", forcelc, 0);
        ModuleIO::print_force(GlobalV::ofs_running, GlobalC::ucell, "NONLOCAL FORCE (eV/Angstrom)", forcenl, 0);
        ModuleIO::print_force(GlobalV::ofs_running, GlobalC::ucell, "NLCC     FORCE (eV/Angstrom)", forcecc, 0);
        ModuleIO::print_force(GlobalV::ofs_running, GlobalC::ucell, "ION      FORCE (eV/Angstrom)", forceion, 0);
        ModuleIO::print_force(GlobalV::ofs_running, GlobalC::ucell, "SCC      FORCE (eV/Angstrom)", forcescc, 0);
        if (GlobalV::EFIELD_FLAG)
            ModuleIO::print_force(GlobalV::ofs_running, GlobalC::ucell, "EFIELD   FORCE (eV/Angstrom)", force_e, 0);
        if (GlobalV::GATE_FLAG)
            ModuleIO::print_force(GlobalV::ofs_running,
                                  GlobalC::ucell,
                                  "GATEFIELD   FORCE (eV/Angstrom)",
                                  force_gate,
                                  0);
        if (GlobalV::imp_sol)
            ModuleIO::print_force(GlobalV::ofs_running, GlobalC::ucell, "IMP_SOL   FORCE (eV/Angstrom)", forcesol, 0);
    }
    ModuleIO::print_force(GlobalV::ofs_running, GlobalC::ucell, "TOTAL-FORCE (eV/Angstrom)", force, 0);

    return;
}

template <typename FPTYPE, typename Device>
void Forces<FPTYPE, Device>::cal_force_loc(ModuleBase::matrix& forcelc,
                                           ModulePW::PW_Basis* rho_basis,
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
        we can keep the blocked aux in L1 cache when iterating GlobalV::NSPIN loop
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
        for (int is = 1; is < GlobalV::NSPIN; is++)
        {
            for (int ir = irb; ir < ir_end; ++ir)
            { // accumulate aux
                aux[ir] += std::complex<double>(chr->rho[is][ir], 0.0);
            }
        }
    }

    // to G space. maybe need fftw with OpenMP
    rho_basis->real2recip(aux, aux);

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int iat = 0; iat < this->nat; ++iat)
    {
        // read `it` `ia` from the table
        int it = GlobalC::ucell.iat2it[iat];
        int ia = GlobalC::ucell.iat2ia[iat];
        for (int ig = 0; ig < rho_basis->npw; ig++)
        {
            const double phase = ModuleBase::TWO_PI * (rho_basis->gcar[ig] * GlobalC::ucell.atoms[it].tau[ia]);
            double sinp, cosp;
            ModuleBase::libm::sincos(phase, &sinp, &cosp);
            const double factor
                = GlobalC::ppcell.vloc(it, rho_basis->ig2igg[ig]) * (cosp * aux[ig].imag() + sinp * aux[ig].real());
            forcelc(iat, 0) += rho_basis->gcar[ig][0] * factor;
            forcelc(iat, 1) += rho_basis->gcar[ig][1] * factor;
            forcelc(iat, 2) += rho_basis->gcar[ig][2] * factor;
        }
        forcelc(iat, 0) *= (GlobalC::ucell.tpiba * GlobalC::ucell.omega);
        forcelc(iat, 1) *= (GlobalC::ucell.tpiba * GlobalC::ucell.omega);
        forcelc(iat, 2) *= (GlobalC::ucell.tpiba * GlobalC::ucell.omega);
    }

    // this->print(GlobalV::ofs_running, "local forces", forcelc);
    Parallel_Reduce::reduce_pool(forcelc.c, forcelc.nr * forcelc.nc);
    delete[] aux;
    ModuleBase::timer::tick("Forces", "cal_force_loc");
    return;
}

template <typename FPTYPE, typename Device>
void Forces<FPTYPE, Device>::cal_force_ew(ModuleBase::matrix& forceion,
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
        we can keep the blocked aux in L1 cache when iterating GlobalC::ucell.ntype loop
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

        { // it = 0
            const double dzv = static_cast<double>(GlobalC::ucell.atoms[0].ncpp.zv);
            for (int ig = igb; ig < ig_end; ++ig)
            { // initialize aux
                aux[ig] = dzv * conj(p_sf->strucFac(0, ig));
            }
        }
        for (int it = 1; it < GlobalC::ucell.ntype; it++)
        {
            const double dzv = static_cast<double>(GlobalC::ucell.atoms[it].ncpp.zv);
            for (int ig = igb; ig < ig_end; ++ig)
            { // accumulate aux
                aux[ig] += dzv * conj(p_sf->strucFac(it, ig));
            }
        }
    }

    // calculate total ionic charge
    double charge = 0.0;
    for (int it = 0; it < GlobalC::ucell.ntype; it++)
    {
        charge += GlobalC::ucell.atoms[it].na * GlobalC::ucell.atoms[it].ncpp.zv; // mohan modify 2007-11-7
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
                     * erfc(sqrt(GlobalC::ucell.tpiba2 * rho_basis->ggecut / 4.0 / alpha));
    } while (upperbound > 1.0e-6);

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int ig = 0; ig < rho_basis->npw; ig++)
    {
        aux[ig] *= ModuleBase::libm::exp(-1.0 * rho_basis->gg[ig] * GlobalC::ucell.tpiba2 / alpha / 4.0)
                   / (rho_basis->gg[ig] * GlobalC::ucell.tpiba2);
    }

    // set pos rho_basis->ig_gge0 to zero
    if (rho_basis->ig_gge0 >= 0 && rho_basis->ig_gge0 < rho_basis->npw)
    {
        aux[rho_basis->ig_gge0] = std::complex<double>(0.0, 0.0);
    }

#ifdef _OPENMP
#pragma omp parallel
    {
        int num_threads = omp_get_num_threads();
        int thread_id = omp_get_thread_num();
#else
    int num_threads = 1;
    int thread_id = 0;
#endif

        /* Here is task distribution for multi-thread,
            0. atom will be iterated both in main nat loop and the loop in `if (rho_basis->ig_gge0 >= 0)`.
                To avoid syncing, we must calculate work range of each thread by our self
            1. Calculate the iat range [iat_beg, iat_end) by each thread
                a. when it is single thread stage, [iat_beg, iat_end) will be [0, nat)
            2. each thread iterate atoms form `iat_beg` to `iat_end-1`
        */
        int iat_beg, iat_end;
        int it_beg, ia_beg;
        ModuleBase::TASK_DIST_1D(num_threads, thread_id, this->nat, iat_beg, iat_end);
        iat_end = iat_beg + iat_end;
        GlobalC::ucell.iat2iait(iat_beg, &ia_beg, &it_beg);

        int iat = iat_beg;
        int it = it_beg;
        int ia = ia_beg;

        // preprocess ig_gap for skipping the ig point
        int ig_gap = (rho_basis->ig_gge0 >= 0 && rho_basis->ig_gge0 < rho_basis->npw) ? rho_basis->ig_gge0 : -1;

        double it_fact = 0.;
        int last_it = -1;

        // iterating atoms
        while (iat < iat_end)
        {
            if (it != last_it)
            { // calculate it_tact when it is changed
                it_fact = GlobalC::ucell.atoms[it].ncpp.zv * ModuleBase::e2 * GlobalC::ucell.tpiba * ModuleBase::TWO_PI
                          / GlobalC::ucell.omega * fact;
                last_it = it;
            }

            const auto ig_loop = [&](int ig_beg, int ig_end) {
                for (int ig = ig_beg; ig < ig_end; ig++)
                {
                    const ModuleBase::Vector3<double> gcar = rho_basis->gcar[ig];
                    const double arg = ModuleBase::TWO_PI * (gcar * GlobalC::ucell.atoms[it].tau[ia]);
                    double sinp, cosp;
                    ModuleBase::libm::sincos(arg, &sinp, &cosp);
                    double sumnb = -cosp * aux[ig].imag() + sinp * aux[ig].real();
                    forceion(iat, 0) += gcar[0] * sumnb;
                    forceion(iat, 1) += gcar[1] * sumnb;
                    forceion(iat, 2) += gcar[2] * sumnb;
                }
            };

            // skip ig_gge0 point by separating ig loop into two part
            ig_loop(0, ig_gap);
            ig_loop(ig_gap + 1, rho_basis->npw);

            forceion(iat, 0) *= it_fact;
            forceion(iat, 1) *= it_fact;
            forceion(iat, 2) *= it_fact;

            ++iat;
            GlobalC::ucell.step_iait(&ia, &it);
        }

        // means that the processor contains G=0 term.
        if (rho_basis->ig_gge0 >= 0)
        {
            double rmax = 5.0 / (sqrt(alpha) * GlobalC::ucell.lat0);
            int nrm = 0;

            // output of rgen: the number of vectors in the sphere
            const int mxr = 200;
            // the maximum number of R vectors included in r
            ModuleBase::Vector3<double>* r = new ModuleBase::Vector3<double>[mxr];
            double* r2 = new double[mxr];
            ModuleBase::GlobalFunc::ZEROS(r2, mxr);
            int* irr = new int[mxr];
            ModuleBase::GlobalFunc::ZEROS(irr, mxr);
            // the square modulus of R_j-tau_s-tau_s'

            int iat1 = iat_beg;
            int T1 = it_beg;
            int I1 = ia_beg;
            const double sqa = sqrt(alpha);
            const double sq8a_2pi = sqrt(8.0 * alpha / ModuleBase::TWO_PI);

            // iterating atoms.
            // do not need to sync threads because task range of each thread is isolated
            while (iat1 < iat_end)
            {
                int iat2 = 0; // mohan fix bug 2011-06-07
                int I2 = 0;
                int T2 = 0;
                while (iat2 < this->nat)
                {
                    if (iat1 != iat2)
                    {
                        ModuleBase::Vector3<double> d_tau
                            = GlobalC::ucell.atoms[T1].tau[I1] - GlobalC::ucell.atoms[T2].tau[I2];
                        H_Ewald_pw::rgen(d_tau, rmax, irr, GlobalC::ucell.latvec, GlobalC::ucell.G, r, r2, nrm);

                        for (int n = 0; n < nrm; n++)
                        {
                            const double rr = sqrt(r2[n]) * GlobalC::ucell.lat0;

                            double factor = GlobalC::ucell.atoms[T1].ncpp.zv * GlobalC::ucell.atoms[T2].ncpp.zv
                                            * ModuleBase::e2 / (rr * rr)
                                            * (erfc(sqa * rr) / rr + sq8a_2pi * ModuleBase::libm::exp(-alpha * rr * rr))
                                            * GlobalC::ucell.lat0;

                            forceion(iat1, 0) -= factor * r[n].x;
                            forceion(iat1, 1) -= factor * r[n].y;
                            forceion(iat1, 2) -= factor * r[n].z;
                        }
                    }
                    ++iat2;
                    GlobalC::ucell.step_iait(&I2, &T2);
                } // atom b
                ++iat1;
                GlobalC::ucell.step_iait(&I1, &T1);
            } // atom a

            delete[] r;
            delete[] r2;
            delete[] irr;
        }
#ifdef _OPENMP
    }
#endif

Parallel_Reduce::reduce_pool(forceion.c, forceion.nr* forceion.nc);

    // this->print(GlobalV::ofs_running, "ewald forces", forceion);

    ModuleBase::timer::tick("Forces", "cal_force_ew");

    delete[] aux;

    return;
}

template <typename FPTYPE, typename Device>
void Forces<FPTYPE, Device>::cal_force_cc(ModuleBase::matrix& forcecc,
                                          ModulePW::PW_Basis* rho_basis,
                                          const Charge* const chr)
{
    ModuleBase::TITLE("Forces", "cal_force_cc");
    // recalculate the exchange-correlation potential.
    ModuleBase::timer::tick("Forces", "cal_force_cc");

    int total_works = 0;
    // cal total works for skipping preprocess
    for (int it = 0; it < GlobalC::ucell.ntype; ++it)
    {
        if (GlobalC::ucell.atoms[it].ncpp.nlcc)
        {
            total_works += GlobalC::ucell.atoms[it].na;
        }
    }
    if (total_works == 0)
    {
        ModuleBase::timer::tick("Forces", "cal_force_cc");
        return;
    }

    ModuleBase::matrix v(GlobalV::NSPIN, rho_basis->nrxx);

    if (XC_Functional::get_func_type() == 3 || XC_Functional::get_func_type() == 5)
    {
#ifdef USE_LIBXC
        const auto etxc_vtxc_v
            = XC_Functional::v_xc_meta(rho_basis->nrxx, GlobalC::ucell.omega, GlobalC::ucell.tpiba, chr);

        // etxc = std::get<0>(etxc_vtxc_v);
        // vtxc = std::get<1>(etxc_vtxc_v);
        v = std::get<2>(etxc_vtxc_v);
#else
        ModuleBase::WARNING_QUIT("cal_force_cc", "to use mGGA, compile with LIBXC");
#endif
    }
    else
    {
        if (GlobalV::NSPIN == 4)
            GlobalC::ucell.cal_ux();
        const auto etxc_vtxc_v = XC_Functional::v_xc(rho_basis->nrxx, chr, &GlobalC::ucell);

        // etxc = std::get<0>(etxc_vtxc_v);
        // vtxc = std::get<1>(etxc_vtxc_v);
        v = std::get<2>(etxc_vtxc_v);
    }

    const ModuleBase::matrix vxc = v;
    std::complex<double>* psiv = new std::complex<double>[rho_basis->nmaxgr];
    if (GlobalV::NSPIN == 1 || GlobalV::NSPIN == 4)
    {
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
        for (int ir = 0; ir < rho_basis->nrxx; ir++)
        {
            psiv[ir] = std::complex<double>(vxc(0, ir), 0.0);
        }
    }
    else
    {
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
        for (int ir = 0; ir < rho_basis->nrxx; ir++)
        {
            psiv[ir] = 0.5 * (vxc(0, ir) + vxc(1, ir));
        }
    }

    // to G space
    rho_basis->real2recip(psiv, psiv);

    // psiv contains now Vxc(G)
    double* rhocg = new double[rho_basis->ngg];
    ModuleBase::GlobalFunc::ZEROS(rhocg, rho_basis->ngg);

    for (int it = 0; it < GlobalC::ucell.ntype; ++it)
    {
        if (GlobalC::ucell.atoms[it].ncpp.nlcc)
        {
            chr->non_linear_core_correction(GlobalC::ppcell.numeric,
                                            GlobalC::ucell.atoms[it].ncpp.msh,
                                            GlobalC::ucell.atoms[it].ncpp.r,
                                            GlobalC::ucell.atoms[it].ncpp.rab,
                                            GlobalC::ucell.atoms[it].ncpp.rho_atc,
                                            rhocg);
#ifdef _OPENMP
#pragma omp parallel
            {
#endif
                for (int ia = 0; ia < GlobalC::ucell.atoms[it].na; ++ia)
                {
                    // get iat form table
                    int iat = GlobalC::ucell.itia2iat(it, ia);
                    double force[3] = {0, 0, 0};
#ifdef _OPENMP
#pragma omp for nowait
#endif
                    for (int ig = 0; ig < rho_basis->npw; ig++)
                    {
                        const ModuleBase::Vector3<double> gv = rho_basis->gcar[ig];
                        const ModuleBase::Vector3<double> pos = GlobalC::ucell.atoms[it].tau[ia];
                        const double rhocgigg = rhocg[rho_basis->ig2igg[ig]];
                        const std::complex<double> psiv_conj = conj(psiv[ig]);

                        const double arg = ModuleBase::TWO_PI * (gv.x * pos.x + gv.y * pos.y + gv.z * pos.z);
                        double sinp, cosp;
                        ModuleBase::libm::sincos(arg, &sinp, &cosp);
                        const std::complex<double> expiarg = std::complex<double>(sinp, cosp);

                        auto ipol0
                            = GlobalC::ucell.tpiba * GlobalC::ucell.omega * rhocgigg * gv.x * psiv_conj * expiarg;
                        force[0] += ipol0.real();

                        auto ipol1
                            = GlobalC::ucell.tpiba * GlobalC::ucell.omega * rhocgigg * gv.y * psiv_conj * expiarg;
                        force[1] += ipol1.real();

                        auto ipol2
                            = GlobalC::ucell.tpiba * GlobalC::ucell.omega * rhocgigg * gv.z * psiv_conj * expiarg;
                        force[2] += ipol2.real();
                    }
#ifdef _OPENMP
                    if (omp_get_num_threads() > 1)
                    {
#pragma omp atomic
                        forcecc(iat, 0) += force[0];
#pragma omp atomic
                        forcecc(iat, 1) += force[1];
#pragma omp atomic
                        forcecc(iat, 2) += force[2];
                    }
                    else
#endif
                    {
                        forcecc(iat, 0) += force[0];
                        forcecc(iat, 1) += force[1];
                        forcecc(iat, 2) += force[2];
                    }
                }
#ifdef _OPENMP
            } // omp parallel
#endif
        }
    }

    delete[] rhocg;

    delete[] psiv;                                                           // mohan fix bug 2012-03-22
    Parallel_Reduce::reduce_pool(forcecc.c, forcecc.nr * forcecc.nc); // qianrui fix a bug for kpar > 1
    ModuleBase::timer::tick("Forces", "cal_force_cc");
    return;
}

template <typename FPTYPE, typename Device>
void Forces<FPTYPE, Device>::cal_force_nl(ModuleBase::matrix& forcenl,
                                          const ModuleBase::matrix& wg,
                                          const ModuleBase::matrix& ekb,
                                          K_Vectors* p_kv,
                                          ModulePW::PW_Basis_K* wfc_basis,
                                          const psi::Psi<complex<FPTYPE>, Device>* psi_in)
{
    ModuleBase::TITLE("Forces", "cal_force_nl");
    ModuleBase::timer::tick("Forces", "cal_force_nl");

    const int nkb = GlobalC::ppcell.nkb;
    if (nkb == 0 || psi_in == nullptr || wfc_basis == nullptr)
    {
        return; // mohan add 2010-07-25
    }

    // There is a contribution for jh<>ih in US case or multi projectors case
    // Actually, the judge of nondiagonal should be done on every atom type
    const bool nondiagonal = (GlobalV::use_uspp || GlobalC::ppcell.multi_proj) ? true : false;

    // dbecp: conj( -iG * <Beta(nkb,npw)|psi(nbnd,npw)> )
    // ModuleBase::ComplexArray dbecp(3, GlobalV::NBANDS, nkb);
    // ModuleBase::ComplexMatrix becp(GlobalV::NBANDS, nkb);
    std::complex<FPTYPE>*dbecp = nullptr, *becp = nullptr, *vkb1 = nullptr, *vkb = nullptr;
    resmem_complex_op()(this->ctx, becp, GlobalV::NBANDS * nkb, "Force::becp");
    resmem_complex_op()(this->ctx, dbecp, 3 * GlobalV::NBANDS * nkb, "Force::dbecp");
    // vkb1: |Beta(nkb,npw)><Beta(nkb,npw)|psi(nbnd,npw)>
    // ModuleBase::ComplexMatrix vkb1(nkb, this->npwx);
    resmem_complex_op()(this->ctx, vkb1, this->npwx * nkb, "Force::vkb1");
    // init additional params
    FPTYPE *force = nullptr;
    FPTYPE *d_wg = nullptr;
    FPTYPE* d_ekb = nullptr;
    FPTYPE *gcar = nullptr;
    auto *deeq = GlobalC::ppcell.get_deeq_data<FPTYPE>();
    auto* qq_nt = GlobalC::ppcell.get_qq_nt_data<FPTYPE>();
    int wg_nc = wg.nc;
    int *atom_nh = nullptr, *atom_na = nullptr;
    int* h_atom_nh = new int[GlobalC::ucell.ntype];
    int* h_atom_na = new int[GlobalC::ucell.ntype];
    for (int ii = 0; ii < GlobalC::ucell.ntype; ii++)
    {
        h_atom_nh[ii] = GlobalC::ucell.atoms[ii].ncpp.nh;
        h_atom_na[ii] = GlobalC::ucell.atoms[ii].na;
    }
    if (this->device == psi::GpuDevice)
    {
        resmem_var_op()(this->ctx, d_wg, wg.nr * wg.nc);
        resmem_var_op()(this->ctx, d_ekb, ekb.nr * ekb.nc);
        resmem_var_op()(this->ctx, force, forcenl.nr * forcenl.nc);
        resmem_var_op()(this->ctx, gcar, 3 * wfc_basis->nks * wfc_basis->npwk_max);
        syncmem_var_h2d_op()(this->ctx, this->cpu_ctx, d_wg, wg.c, wg.nr * wg.nc);
        syncmem_var_h2d_op()(this->ctx, this->cpu_ctx, d_ekb, ekb.c, ekb.nr * ekb.nc);
        syncmem_var_h2d_op()(this->ctx, this->cpu_ctx, force, forcenl.c, forcenl.nr * forcenl.nc);
        syncmem_var_h2d_op()(this->ctx,
                             this->cpu_ctx,
                             gcar,
                             &wfc_basis->gcar[0][0],
                             3 * wfc_basis->nks * wfc_basis->npwk_max);

        resmem_int_op()(this->ctx, atom_nh, GlobalC::ucell.ntype);
        resmem_int_op()(this->ctx, atom_na, GlobalC::ucell.ntype);
        syncmem_int_h2d_op()(this->ctx, this->cpu_ctx, atom_nh, h_atom_nh, GlobalC::ucell.ntype);
        syncmem_int_h2d_op()(this->ctx, this->cpu_ctx, atom_na, h_atom_na, GlobalC::ucell.ntype);
    }
    else
    {
        d_wg = wg.c;
        d_ekb = ekb.c;
        force = forcenl.c;
        gcar = &wfc_basis->gcar[0][0];
        atom_nh = h_atom_nh;
        atom_na = h_atom_na;
    }

    for (int ik = 0; ik < wfc_basis->nks; ik++)
    {
        if (GlobalV::NSPIN == 2)
            GlobalV::CURRENT_SPIN = p_kv->isk[ik];
        const int nbasis = wfc_basis->npwk[ik];
        // generate vkb
        if (GlobalC::ppcell.nkb > 0)
        {
            vkb = GlobalC::ppcell.get_vkb_data<FPTYPE>();
            GlobalC::ppcell.getvnl(ctx, ik, vkb);
        }

        // get becp according to wave functions and vkb
        // important here ! becp must set zero!!
        // vkb: Beta(nkb,npw)
        // becp(nkb,nbnd): <Beta(nkb,npw)|psi(nbnd,npw)>
        // becp.zero_out();
        psi_in[0].fix_k(ik);
        char transa = 'C';
        char transb = 'N';
        ///
        /// only occupied band should be calculated.
        ///
        int nbands_occ = GlobalV::NBANDS;
        const double threshold = ModuleBase::threshold_wg * wg(ik, 0);
        while (std::fabs(wg(ik, nbands_occ - 1)) < threshold)
        {
            nbands_occ--;
            if (nbands_occ == 0)
            {
                break;
            }
        }
        int npm = GlobalV::NPOL * nbands_occ;
        gemm_op()(this->ctx,
                  transa,
                  transb,
                  nkb,
                  npm,
                  nbasis,
                  &ModuleBase::ONE,
                  vkb,
                  this->npwx,
                  psi_in[0].get_pointer(),
                  this->npwx,
                  &ModuleBase::ZERO,
                  becp,
                  nkb);
        if (this->device == psi::GpuDevice)
        {
            std::complex<FPTYPE>* h_becp = nullptr;
            resmem_complex_h_op()(this->cpu_ctx, h_becp, GlobalV::NBANDS * nkb);
            syncmem_complex_d2h_op()(this->cpu_ctx, this->ctx, h_becp, becp, GlobalV::NBANDS * nkb);
            Parallel_Reduce::reduce_pool(h_becp, GlobalV::NBANDS * nkb);
            syncmem_complex_h2d_op()(this->ctx, this->cpu_ctx, becp, h_becp, GlobalV::NBANDS * nkb);
            delmem_complex_h_op()(this->cpu_ctx, h_becp);
        }
        else
        {
            Parallel_Reduce::reduce_pool(becp, GlobalV::NBANDS * nkb);
        }
        // out.printcm_real("becp",becp,1.0e-4);
        //  Calculate the derivative of beta,
        //  |dbeta> =  -ig * |beta>
        // dbecp.zero_out();
        for (int ipol = 0; ipol < 3; ipol++)
        {
            cal_vkb1_nl_op()(this->ctx,
                             nkb,
                             this->npwx,
                             wfc_basis->npwk_max,
                             GlobalC::ppcell.vkb.nc,
                             nbasis,
                             ik,
                             ipol,
                             ModuleBase::NEG_IMAG_UNIT,
                             vkb,
                             gcar,
                             vkb1);
            std::complex<double>* pdbecp = dbecp + ipol * GlobalV::NBANDS * nkb;
            gemm_op()(this->ctx,
                      transa,
                      transb,
                      nkb,
                      npm,
                      nbasis,
                      &ModuleBase::ONE,
                      vkb1,
                      this->npwx,
                      psi_in[0].get_pointer(),
                      this->npwx,
                      &ModuleBase::ZERO,
                      pdbecp,
                      nkb);
        } // end ipol

        //		don't need to reduce here, keep dbecp different in each processor,
        //		and at last sum up all the forces.
        //		Parallel_Reduce::reduce_pool( dbecp.ptr, dbecp.ndata);
        cal_force_nl_op()(this->ctx,
                          nondiagonal,
                          nbands_occ,
                          wg_nc,
                          GlobalC::ucell.ntype,
                          GlobalV::CURRENT_SPIN,
                          GlobalC::ppcell.deeq.getBound2(),
                          GlobalC::ppcell.deeq.getBound3(),
                          GlobalC::ppcell.deeq.getBound4(),
                          forcenl.nc,
                          GlobalV::NBANDS,
                          ik,
                          nkb,
                          atom_nh,
                          atom_na,
                          GlobalC::ucell.tpiba,
                          d_wg,
                          d_ekb,
                          qq_nt,
                          deeq,
                          becp,
                          dbecp,
                          force);
    } // end ik

    if (this->device == psi::GpuDevice)
    {
        syncmem_var_d2h_op()(this->cpu_ctx, this->ctx, forcenl.c, force, forcenl.nr * forcenl.nc);
    }
    // sum up forcenl from all processors
    Parallel_Reduce::reduce_all(forcenl.c, forcenl.nr* forcenl.nc);

    delete[] h_atom_nh;
    delete[] h_atom_na;
    delmem_complex_op()(this->ctx, vkb1);
    delmem_complex_op()(this->ctx, becp);
    delmem_complex_op()(this->ctx, dbecp);
    if (this->device == psi::GpuDevice)
    {
        delmem_var_op()(this->ctx, d_wg);
        delmem_var_op()(this->ctx, gcar);
        delmem_var_op()(this->ctx, force);
        delmem_int_op()(this->ctx, atom_nh);
        delmem_int_op()(this->ctx, atom_na);
    }
    //  this->print(GlobalV::ofs_running, "nonlocal forces", forcenl);
    ModuleBase::timer::tick("Forces", "cal_force_nl");
}

template <typename FPTYPE, typename Device>
void Forces<FPTYPE, Device>::cal_force_scc(ModuleBase::matrix& forcescc,
                                           ModulePW::PW_Basis* rho_basis,
                                           const ModuleBase::matrix& vnew,
                                           const bool vnew_exist)
{
    ModuleBase::TITLE("Forces", "cal_force_scc");
    ModuleBase::timer::tick("Forces", "cal_force_scc");

    // for orbital free case
    if (!vnew_exist)
    {
        ModuleBase::timer::tick("Forces", "cal_force_scc");
        return;
    }

    std::complex<double>* psic = new std::complex<double>[rho_basis->nmaxgr];

    const int nrxx = vnew.nc;
    const int nspin = vnew.nr;

    if (nspin == 1 || nspin == 4)
    {
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
        for (int ir = 0; ir < nrxx; ir++)
        {
            psic[ir] = vnew(0, ir);
        }
    }
    else
    {
        int isup = 0;
        int isdw = 1;
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
        for (int ir = 0; ir < nrxx; ir++)
        {
            psic[ir] = (vnew(isup, ir) + vnew(isdw, ir)) * 0.5;
        }
    }

    int ndm = 0;

    for (int it = 0; it < GlobalC::ucell.ntype; it++)
    {
        if (ndm < GlobalC::ucell.atoms[it].ncpp.msh)
        {
            ndm = GlobalC::ucell.atoms[it].ncpp.msh;
        }
    }

    // work space
    double* rhocgnt = new double[rho_basis->ngg];
    ModuleBase::GlobalFunc::ZEROS(rhocgnt, rho_basis->ngg);

    rho_basis->real2recip(psic, psic);

    int igg0 = 0;
    const int ig0 = rho_basis->ig_gge0;
    if (rho_basis->gg_uniq[0] < 1.0e-8)
        igg0 = 1;

    double fact = 2.0;
    for (int nt = 0; nt < GlobalC::ucell.ntype; nt++)
    {
        //		Here we compute the G.ne.0 term
        const int mesh = GlobalC::ucell.atoms[nt].ncpp.msh;
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (int ig = igg0; ig < rho_basis->ngg; ++ig)
        {
            double* aux = new double[mesh];
            const double gx = sqrt(rho_basis->gg_uniq[ig]) * GlobalC::ucell.tpiba;
            for (int ir = 0; ir < mesh; ir++)
            {
                if (GlobalC::ucell.atoms[nt].ncpp.r[ir] < 1.0e-8)
                {
                    aux[ir] = GlobalC::ucell.atoms[nt].ncpp.rho_at[ir];
                }
                else
                {
                    const double gxx = gx * GlobalC::ucell.atoms[nt].ncpp.r[ir];
                    aux[ir] = GlobalC::ucell.atoms[nt].ncpp.rho_at[ir] * ModuleBase::libm::sin(gxx) / gxx;
                }
            }
            ModuleBase::Integral::Simpson_Integral(mesh, aux, GlobalC::ucell.atoms[nt].ncpp.rab, rhocgnt[ig]);
            delete[] aux; // mohan fix bug 2012-03-22
        }

        int iat = 0;
        for (int it = 0; it < GlobalC::ucell.ntype; it++)
        {
            for (int ia = 0; ia < GlobalC::ucell.atoms[it].na; ia++)
            {
                if (nt == it)
                {
                    const ModuleBase::Vector3<double> pos = GlobalC::ucell.atoms[it].tau[ia];
                    double &force0 = forcescc(iat, 0), &force1 = forcescc(iat, 1), &force2 = forcescc(iat, 2);
#ifdef _OPENMP
#pragma omp parallel for reduction(+ : force0) reduction(+ : force1) reduction(+ : force2)
#endif
                    for (int ig = 0; ig < rho_basis->npw; ++ig)
                    {
                        if (ig == ig0)
                            continue;
                        const ModuleBase::Vector3<double> gv = rho_basis->gcar[ig];
                        const double rhocgntigg = rhocgnt[rho_basis->ig2igg[ig]];
                        const double arg = ModuleBase::TWO_PI * (gv * pos);
                        double sinp, cosp;
                        ModuleBase::libm::sincos(arg, &sinp, &cosp);
                        const std::complex<double> cpm = std::complex<double>(sinp, cosp) * conj(psic[ig]);

                        force0 += fact * rhocgntigg * GlobalC::ucell.tpiba * gv.x * cpm.real();
                        force1 += fact * rhocgntigg * GlobalC::ucell.tpiba * gv.y * cpm.real();
                        force2 += fact * rhocgntigg * GlobalC::ucell.tpiba * gv.z * cpm.real();
                    }
                    // std::cout << " forcescc = " << forcescc(iat,0) << " " << forcescc(iat,1) << " " <<
                    // forcescc(iat,2) << std::endl;
                }
                iat++;
            }
        }
    }

    Parallel_Reduce::reduce_pool(forcescc.c, forcescc.nr* forcescc.nc);

    delete[] psic;    // mohan fix bug 2012-03-22
    delete[] rhocgnt; // mohan fix bug 2012-03-22

    ModuleBase::timer::tick("Forces", "cal_force_scc");
    return;
}

template class Forces<double, psi::DEVICE_CPU>;
#if ((defined __CUDA) || (defined __ROCM))
template class Forces<double, psi::DEVICE_GPU>;
#endif