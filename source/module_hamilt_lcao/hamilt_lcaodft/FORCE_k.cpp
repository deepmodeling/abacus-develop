#include "FORCE.h"
#include "source_base/memory.h"
#include "source_base/parallel_reduce.h"
#include "source_base/timer.h"
#include "source_base/tool_threading.h"
#include "module_basis/module_ao/ORB_read.h"
#include "module_cell/module_neighbor/sltk_grid_driver.h"
#include "module_elecstate/cal_dm.h"
#include "module_elecstate/elecstate_lcao.h"
#include "module_elecstate/module_dm/cal_dm_psi.h"
#include "module_hamilt_lcao/hamilt_lcaodft/LCAO_domain.h"
#include "module_hamilt_lcao/hamilt_lcaodft/pulay_force_stress.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_io/write_HS.h"
#include "module_parameter/parameter.h"

#include <map>
#include <unordered_map>

#ifdef __MLALGO
#include "module_hamilt_lcao/module_deepks/LCAO_deepks.h"
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

template <>
void Force_LCAO<std::complex<double>>::allocate(const UnitCell& ucell,
                                                const Grid_Driver& gd,
                                                const Parallel_Orbitals& pv,
                                                ForceStressArrays& fsr, // mohan add 2024-06-15
                                                const TwoCenterBundle& two_center_bundle,
                                                const LCAO_Orbitals& orb,
                                                const int& nks,
                                                const std::vector<ModuleBase::Vector3<double>>& kvec_d)
{
    ModuleBase::TITLE("Forces", "allocate");
    ModuleBase::timer::tick("Forces", "allocate");

    const int nnr = pv.nnr;

    assert(nnr >= 0);

    //--------------------------------
    // (1) allocate for dSx dSy & dSz
    //--------------------------------
    fsr.DSloc_Rx = new double[nnr];
    fsr.DSloc_Ry = new double[nnr];
    fsr.DSloc_Rz = new double[nnr];

    const auto init_DSloc_Rxyz = [this, nnr, &fsr](int num_threads, int thread_id) {
        int beg = 0;
        int len = 0;
        ModuleBase::BLOCK_TASK_DIST_1D(num_threads, thread_id, nnr, 1024, beg, len);
        ModuleBase::GlobalFunc::ZEROS(fsr.DSloc_Rx + beg, len);
        ModuleBase::GlobalFunc::ZEROS(fsr.DSloc_Ry + beg, len);
        ModuleBase::GlobalFunc::ZEROS(fsr.DSloc_Rz + beg, len);
    };

    ModuleBase::OMP_PARALLEL(init_DSloc_Rxyz);
    ModuleBase::Memory::record("Force::dS_K", sizeof(double) * nnr * 3);

    if (PARAM.inp.cal_stress)
    {
        fsr.DH_r = new double[3 * nnr];
        fsr.stvnl11 = new double[nnr];
        fsr.stvnl12 = new double[nnr];
        fsr.stvnl13 = new double[nnr];
        fsr.stvnl22 = new double[nnr];
        fsr.stvnl23 = new double[nnr];
        fsr.stvnl33 = new double[nnr];
        const auto init_DH_r_stvnl = [this, nnr, &fsr](int num_threads, int thread_id) {
            int beg, len;
            ModuleBase::BLOCK_TASK_DIST_1D(num_threads, thread_id, nnr, 1024, beg, len);
            ModuleBase::GlobalFunc::ZEROS(fsr.DH_r + 3 * beg, 3 * len);
            ModuleBase::GlobalFunc::ZEROS(fsr.stvnl11 + beg, len);
            ModuleBase::GlobalFunc::ZEROS(fsr.stvnl12 + beg, len);
            ModuleBase::GlobalFunc::ZEROS(fsr.stvnl13 + beg, len);
            ModuleBase::GlobalFunc::ZEROS(fsr.stvnl22 + beg, len);
            ModuleBase::GlobalFunc::ZEROS(fsr.stvnl23 + beg, len);
            ModuleBase::GlobalFunc::ZEROS(fsr.stvnl33 + beg, len);
        };
        ModuleBase::OMP_PARALLEL(init_DH_r_stvnl);

        ModuleBase::Memory::record("Stress::dHr", sizeof(double) * nnr * 3);
        ModuleBase::Memory::record("Stress::dSR", sizeof(double) * nnr * 6);
    }

    //-----------------------------
    // calculate dS = <phi | dphi>
    //-----------------------------
    bool cal_deri = true;
    LCAO_domain::build_ST_new(fsr,
                              'S',
                              cal_deri,
                              PARAM.inp.cal_stress,
                              ucell,
                              orb,
                              pv,
                              two_center_bundle,
                              &gd,
                              nullptr); // delete lm.SlocR

    //-----------------------------------------
    // (2) allocate for <phi | T + Vnl | dphi>
    //-----------------------------------------
    fsr.DHloc_fixedR_x = new double[nnr];
    fsr.DHloc_fixedR_y = new double[nnr];
    fsr.DHloc_fixedR_z = new double[nnr];

    const auto init_DHloc_fixedR_xyz = [this, nnr, &fsr](int num_threads, int thread_id) {
        int beg = 0;
        int len = 0;
        ModuleBase::BLOCK_TASK_DIST_1D(num_threads, thread_id, nnr, 1024, beg, len);
        ModuleBase::GlobalFunc::ZEROS(fsr.DHloc_fixedR_x + beg, len);
        ModuleBase::GlobalFunc::ZEROS(fsr.DHloc_fixedR_y + beg, len);
        ModuleBase::GlobalFunc::ZEROS(fsr.DHloc_fixedR_z + beg, len);
    };
    ModuleBase::OMP_PARALLEL(init_DHloc_fixedR_xyz);
    ModuleBase::Memory::record("Force::dTVNL", sizeof(double) * nnr * 3);

    // calculate dT=<phi|kin|dphi> in LCAO
    // calculate T + VNL(P1) in LCAO basis
    LCAO_domain::build_ST_new(fsr,
                              'T',
                              cal_deri,
                              PARAM.inp.cal_stress,
                              ucell,
                              orb,
                              pv,
                              two_center_bundle,
                              &gd,
                              nullptr); // delete lm.Hloc_fixedR

    // calculate asynchronous S matrix to output for Hefei-NAMD
    if (PARAM.inp.cal_syns)
    {
        cal_deri = false;

        ModuleBase::WARNING_QUIT("cal_syns", "This function has been broken and will be fixed later.");

        LCAO_domain::build_ST_new(fsr,
                                  'S',
                                  cal_deri,
                                  PARAM.inp.cal_stress,
                                  ucell,
                                  orb,
                                  pv,
                                  two_center_bundle,
                                  &(gd),
                                  nullptr, // delete lm.SlocR
                                  PARAM.inp.cal_syns,
                                  PARAM.inp.dmax);

        for (int ik = 0; ik < nks; ik++)
        {

            bool bit = false; // LiuXh, 2017-03-21
        }
    }

    ModuleBase::timer::tick("Forces", "allocate");
    return;
}

template <>
void Force_LCAO<std::complex<double>>::finish_ftable(ForceStressArrays& fsr)
{
    delete[] fsr.DSloc_Rx;
    delete[] fsr.DSloc_Ry;
    delete[] fsr.DSloc_Rz;
    delete[] fsr.DHloc_fixedR_x;
    delete[] fsr.DHloc_fixedR_y;
    delete[] fsr.DHloc_fixedR_z;

    if (PARAM.inp.cal_stress)
    {
        delete[] fsr.DH_r;
        delete[] fsr.stvnl11;
        delete[] fsr.stvnl12;
        delete[] fsr.stvnl13;
        delete[] fsr.stvnl22;
        delete[] fsr.stvnl23;
        delete[] fsr.stvnl33;
    }
    return;
}

// be called in Force_LCAO::start_force_calculation
template <>
void Force_LCAO<std::complex<double>>::ftable(const bool isforce,
                                              const bool isstress,
                                              ForceStressArrays& fsr, // mohan add 2024-06-15
                                              const UnitCell& ucell,
                                              const Grid_Driver& gd,
                                              const psi::Psi<std::complex<double>>* psi,
                                              const elecstate::ElecState* pelec,
                                              ModuleBase::matrix& foverlap,
                                              ModuleBase::matrix& ftvnl_dphi,
                                              ModuleBase::matrix& fvnl_dbeta,
                                              ModuleBase::matrix& fvl_dphi,
                                              ModuleBase::matrix& soverlap,
                                              ModuleBase::matrix& stvnl_dphi,
                                              ModuleBase::matrix& svnl_dbeta,
                                              ModuleBase::matrix& svl_dphi,
#ifdef __MLALGO
                                              ModuleBase::matrix& fvnl_dalpha,
                                              ModuleBase::matrix& svnl_dalpha,
                                              LCAO_Deepks<std::complex<double>>& ld,
#endif
                                              TGint<std::complex<double>>::type& gint,
                                              const TwoCenterBundle& two_center_bundle,
                                              const LCAO_Orbitals& orb,
                                              const Parallel_Orbitals& pv,
                                              const K_Vectors* kv,
                                              Record_adj* ra)
{
    ModuleBase::TITLE("Forces", "ftable");
    ModuleBase::timer::tick("Forces", "ftable");

    elecstate::DensityMatrix<complex<double>, double>* dm
        = dynamic_cast<const elecstate::ElecStateLCAO<std::complex<double>>*>(pelec)->get_DM();

    this->allocate(ucell,
                   gd,
                   pv,
                   fsr, // mohan add 2024-06-16
                   two_center_bundle,
                   orb,
                   kv->get_nks(),
                   kv->kvec_d);

    const double* dSx[3] = {fsr.DSloc_Rx, fsr.DSloc_Ry, fsr.DSloc_Rz};
    // calculate the energy density matrix
    // and the force related to overlap matrix and energy density matrix.
    PulayForceStress::cal_pulay_fs(
        foverlap,
        soverlap,
        this->cal_edm(pelec, *psi, *dm, *kv, pv, PARAM.inp.nspin, PARAM.inp.nbands, ucell, *ra),
        ucell,
        pv,
        dSx,
        fsr.DH_r,
        isforce,
        isstress,
        ra,
        -1.0,
        1.0);

    const double* dHx[3] = {fsr.DHloc_fixedR_x, fsr.DHloc_fixedR_y, fsr.DHloc_fixedR_z};                    // T+Vnl
    const double* dHxy[6] = {fsr.stvnl11, fsr.stvnl12, fsr.stvnl13, fsr.stvnl22, fsr.stvnl23, fsr.stvnl33}; // T
    // tvnl_dphi
    PulayForceStress::cal_pulay_fs(ftvnl_dphi, stvnl_dphi, *dm, ucell, pv, dHx, dHxy, isforce, isstress, ra, 1.0, -1.0);

    // doing on the real space grid.
    // vl_dphi
    PulayForceStress::cal_pulay_fs(fvl_dphi,
                                   svl_dphi,
                                   *dm,
                                   ucell,
                                   pelec->pot,
                                   gint,
                                   isforce,
                                   isstress,
                                   false /*reset dm to gint*/);

#ifdef __MLALGO
    if (PARAM.inp.deepks_scf)
    {
        // No need to update E_delta since it have been done in LCAO_Deepks_Interface in after_scf
        DeePKS_domain::cal_f_delta<std::complex<double>>(ld.dm_r,
                                                         ucell,
                                                         orb,
                                                         gd,
                                                         pv,
                                                         kv->get_nks(),
                                                         kv->kvec_d,
                                                         ld.phialpha,
                                                         ld.gedm,
                                                         ld.inl_index,
                                                         fvnl_dalpha,
                                                         isstress,
                                                         svnl_dalpha);
    }
#endif

    //----------------------------------------------------------------
    // reduce the force according to 2D distribution of H & S matrix.
    //----------------------------------------------------------------
    if (isforce)
    {
        Parallel_Reduce::reduce_pool(foverlap.c, foverlap.nr * foverlap.nc);
        Parallel_Reduce::reduce_pool(ftvnl_dphi.c, ftvnl_dphi.nr * ftvnl_dphi.nc);
        Parallel_Reduce::reduce_pool(fvnl_dbeta.c, fvnl_dbeta.nr * fvnl_dbeta.nc);
        Parallel_Reduce::reduce_pool(fvl_dphi.c, fvl_dphi.nr * fvl_dphi.nc);
#ifdef __MLALGO
        Parallel_Reduce::reduce_pool(fvnl_dalpha.c, fvnl_dalpha.nr * fvnl_dalpha.nc);
#endif
    }
    if (isstress)
    {
        Parallel_Reduce::reduce_pool(soverlap.c, soverlap.nr * soverlap.nc);
        Parallel_Reduce::reduce_pool(stvnl_dphi.c, stvnl_dphi.nr * stvnl_dphi.nc);
        Parallel_Reduce::reduce_pool(svnl_dbeta.c, svnl_dbeta.nr * svnl_dbeta.nc);
        Parallel_Reduce::reduce_pool(svl_dphi.c, svl_dphi.nr * svl_dphi.nc);
#ifdef __MLALGO
        Parallel_Reduce::reduce_pool(svnl_dalpha.c, svnl_dalpha.nr * svnl_dalpha.nc);
#endif
    }

#ifdef __MLALGO
    if (PARAM.inp.deepks_scf && PARAM.inp.deepks_out_unittest)
    {
        std::ofstream ofs_f("F_delta.dat");
        std::ofstream ofs_s("stress_delta.dat");
        ofs_f << std::setprecision(10);
        ofs_s << std::setprecision(10);
        fvnl_dalpha.print(ofs_f);
        ofs_f.close();
        svnl_dalpha.print(ofs_s);
        ofs_s.close();
    }
#endif

    ModuleBase::timer::tick("Forces", "ftable");
    return;
}
