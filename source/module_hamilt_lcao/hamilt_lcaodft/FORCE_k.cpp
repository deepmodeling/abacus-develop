#include "FORCE_k.h"

#include "module_base/memory.h"
#include "module_base/timer.h"
#include "module_base/parallel_reduce.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_cell/module_neighbor/sltk_grid_driver.h"
#include "module_elecstate/cal_dm.h"
#include "module_base/tool_threading.h"
#include "module_io/write_HS.h"

#include <map>
#include <unordered_map>

#ifdef __DEEPKS
#include "module_hamilt_lcao/module_deepks/LCAO_deepks.h"
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

Force_LCAO_k::Force_LCAO_k()
{
}

Force_LCAO_k::~Force_LCAO_k()
{
}

// be called in Force_LCAO::start_force_calculation
void Force_LCAO_k::ftable_k(const bool isforce,
                            const bool isstress,
                            Record_adj& ra,
                            const psi::Psi<std::complex<double>>* psi,
                            Local_Orbital_Charge& loc,
                            const elecstate::ElecState* pelec,
                            ModuleBase::matrix& foverlap,
                            ModuleBase::matrix& ftvnl_dphi,
                            ModuleBase::matrix& fvnl_dbeta,
                            ModuleBase::matrix& fvl_dphi,
                            ModuleBase::matrix& soverlap,
                            ModuleBase::matrix& stvnl_dphi,
                            ModuleBase::matrix& svnl_dbeta,
#ifdef __DEEPKS
                            ModuleBase::matrix& svl_dphi,
                            ModuleBase::matrix& svnl_dalpha,
#else
                            ModuleBase::matrix& svl_dphi,
#endif
                            LCAO_Hamilt& uhm,
                            const K_Vectors& kv)
{
    ModuleBase::TITLE("Force_LCAO_k", "ftable_k");
    ModuleBase::timer::tick("Force_LCAO_k", "ftable_k");

    this->UHM = &uhm;

    const Parallel_Orbitals* pv = loc.ParaV;
    this->allocate_k(*pv, kv.nks, kv.kvec_d);

    // calculate the energy density matrix
    // and the force related to overlap matrix and energy density matrix.
    this->cal_foverlap_k(isforce, isstress, ra, psi, loc, foverlap, soverlap, pelec, kv.nks, kv);

    // calculate the density matrix
    double** dm2d = new double*[GlobalV::NSPIN];
    for (int is = 0; is < GlobalV::NSPIN; is++)
    {
        dm2d[is] = new double[pv->nnr];
    }
    ModuleBase::Memory::record("Force::dm2d_K", sizeof(double) * GlobalV::NSPIN * pv->nnr);
    auto init_dm2d = [dm2d, pv](int num_threads, int thread_id)
    {
        int beg, len;
        ModuleBase::BLOCK_TASK_DIST_1D(num_threads, thread_id, pv->nnr, 1024, beg, len);
        for (int is = 0; is < GlobalV::NSPIN; is++)
        {
            ModuleBase::GlobalFunc::ZEROS(dm2d[is] + beg, len);
        }
    };
    ModuleBase::OMP_PARALLEL(init_dm2d);

    loc.cal_dm_R(loc.dm_k, ra, dm2d, kv);

    this->cal_ftvnl_dphi_k(dm2d, isforce, isstress, ra, ftvnl_dphi, stvnl_dphi);

    // ---------------------------------------
    // doing on the real space grid.
    // ---------------------------------------
    this->cal_fvl_dphi_k(isforce, isstress, pelec->pot, fvl_dphi, svl_dphi, loc.DM_R);

    this->cal_fvnl_dbeta_k_new(dm2d, isforce, isstress, fvnl_dbeta, svnl_dbeta);

#ifdef __DEEPKS
    if (GlobalV::deepks_scf)
    {
        GlobalC::ld.cal_projected_DM_k(loc.dm_k,
                                       GlobalC::ucell,
                                       GlobalC::ORB,
                                       GlobalC::GridD,
                                       pv->trace_loc_row,
                                       pv->trace_loc_col,
                                       kv.nks,
                                       kv.kvec_d);
        GlobalC::ld.cal_descriptor();
        GlobalC::ld.cal_gedm(GlobalC::ucell.nat);

        GlobalC::ld.cal_f_delta_k(loc.dm_k,
                                  GlobalC::ucell,
                                  GlobalC::ORB,
                                  GlobalC::GridD,
                                  pv->trace_loc_row,
                                  pv->trace_loc_col,
                                  kv.nks,
                                  kv.kvec_d,
                                  isstress,
                                  svnl_dalpha);
#ifdef __MPI
        Parallel_Reduce::reduce_double_all(GlobalC::ld.F_delta.c, GlobalC::ld.F_delta.nr * GlobalC::ld.F_delta.nc);
        if (isstress)
        {
            Parallel_Reduce::reduce_double_pool(svnl_dalpha.c, svnl_dalpha.nr * svnl_dalpha.nc);
        }
#endif
        if (GlobalV::deepks_out_unittest)
        {
            GlobalC::ld.print_dm_k(kv.nks, loc.dm_k);
            GlobalC::ld.check_projected_dm();
            GlobalC::ld.check_descriptor(GlobalC::ucell);
            GlobalC::ld.check_gedm();
            GlobalC::ld.add_v_delta_k(GlobalC::ucell,
                                      GlobalC::ORB,
                                      GlobalC::GridD,
                                      pv->trace_loc_row,
                                      pv->trace_loc_col,
                                      pv->nnr);
            GlobalC::ld.check_v_delta_k(pv->nnr);
            for (int ik = 0; ik < kv.nks; ik++)
            {
                uhm.LM->folding_fixedH(ik, kv.kvec_d);
            }
            GlobalC::ld.cal_e_delta_band_k(loc.dm_k,
                                           pv->trace_loc_row,
                                           pv->trace_loc_col,
                                           kv.nks,
                                           pv->nrow,
                                           pv->ncol);
            ofstream ofs("E_delta_bands.dat");
            ofs << std::setprecision(10) << GlobalC::ld.e_delta_band;
            ofstream ofs1("E_delta.dat");
            ofs1 << std::setprecision(10) << GlobalC::ld.E_delta;
            GlobalC::ld.check_f_delta(GlobalC::ucell.nat, svnl_dalpha);
        }
    }
#endif

    for (int is = 0; is < GlobalV::NSPIN; is++)
    {
        delete[] dm2d[is];
    }
    delete[] dm2d;

    //----------------------------------------------------------------
    // reduce the force according to 2D distribution of H & S matrix.
    //----------------------------------------------------------------
    if (isforce)
    {
        Parallel_Reduce::reduce_double_pool(foverlap.c, foverlap.nr * foverlap.nc);
        Parallel_Reduce::reduce_double_pool(ftvnl_dphi.c, ftvnl_dphi.nr * ftvnl_dphi.nc);
        Parallel_Reduce::reduce_double_pool(fvnl_dbeta.c, fvnl_dbeta.nr * fvnl_dbeta.nc);
        Parallel_Reduce::reduce_double_pool(fvl_dphi.c, fvl_dphi.nr * fvl_dphi.nc);
    }
    if (isstress)
    {
        Parallel_Reduce::reduce_double_pool(soverlap.c, soverlap.nr * soverlap.nc);
        Parallel_Reduce::reduce_double_pool(stvnl_dphi.c, stvnl_dphi.nr * stvnl_dphi.nc);
        Parallel_Reduce::reduce_double_pool(svnl_dbeta.c, svnl_dbeta.nr * svnl_dbeta.nc);
        Parallel_Reduce::reduce_double_pool(svl_dphi.c, svl_dphi.nr * svl_dphi.nc);
    }


    ModuleBase::timer::tick("Force_LCAO_k", "ftable_k");
    return;
}

void Force_LCAO_k::allocate_k(const Parallel_Orbitals& pv,
                            const int& nks,
                            const std::vector<ModuleBase::Vector3<double>>& kvec_d)
{
    ModuleBase::TITLE("Force_LCAO_k", "allocate_k");
    ModuleBase::timer::tick("Force_LCAO_k", "allocate_k");

    this->ParaV = &pv;
    const int nnr = pv.nnr;
    //--------------------------------
    // (1) allocate for dSx dSy & dSz
    //--------------------------------
    this->UHM->LM->DSloc_Rx = new double[nnr];
    this->UHM->LM->DSloc_Ry = new double[nnr];
    this->UHM->LM->DSloc_Rz = new double[nnr];
    const auto init_DSloc_Rxyz = [this, nnr](int num_threads, int thread_id) {
        int beg, len;
        ModuleBase::BLOCK_TASK_DIST_1D(num_threads, thread_id, nnr, 1024, beg, len);
        ModuleBase::GlobalFunc::ZEROS(this->UHM->LM->DSloc_Rx + beg, len);
        ModuleBase::GlobalFunc::ZEROS(this->UHM->LM->DSloc_Ry + beg, len);
        ModuleBase::GlobalFunc::ZEROS(this->UHM->LM->DSloc_Rz + beg, len);
    };
    ModuleBase::OMP_PARALLEL(init_DSloc_Rxyz);
    ModuleBase::Memory::record("Force::dS_K", sizeof(double) * nnr * 3);

    if (GlobalV::CAL_STRESS)
    {
        this->UHM->LM->DH_r = new double[3 * nnr];
        this->UHM->LM->stvnl11 = new double[nnr];
        this->UHM->LM->stvnl12 = new double[nnr];
        this->UHM->LM->stvnl13 = new double[nnr];
        this->UHM->LM->stvnl22 = new double[nnr];
        this->UHM->LM->stvnl23 = new double[nnr];
        this->UHM->LM->stvnl33 = new double[nnr];
        const auto init_DH_r_stvnl = [this, nnr](int num_threads, int thread_id) {
            int beg, len;
            ModuleBase::BLOCK_TASK_DIST_1D(num_threads, thread_id, nnr, 1024, beg, len);
            ModuleBase::GlobalFunc::ZEROS(this->UHM->LM->DH_r + 3 * beg, 3 * len);
            ModuleBase::GlobalFunc::ZEROS(this->UHM->LM->stvnl11 + beg, len);
            ModuleBase::GlobalFunc::ZEROS(this->UHM->LM->stvnl12 + beg, len);
            ModuleBase::GlobalFunc::ZEROS(this->UHM->LM->stvnl13 + beg, len);
            ModuleBase::GlobalFunc::ZEROS(this->UHM->LM->stvnl22 + beg, len);
            ModuleBase::GlobalFunc::ZEROS(this->UHM->LM->stvnl23 + beg, len);
            ModuleBase::GlobalFunc::ZEROS(this->UHM->LM->stvnl33 + beg, len);
        };
        ModuleBase::OMP_PARALLEL(init_DH_r_stvnl);
        
        ModuleBase::Memory::record("Stress::dHr", sizeof(double) * nnr * 3);
        ModuleBase::Memory::record("Stress::dSR", sizeof(double) * nnr * 6);
    }

    //-----------------------------
    // calculate dS = <phi | dphi>
    //-----------------------------
    bool cal_deri = true;
    this->UHM->genH.build_ST_new('S', cal_deri, GlobalC::ucell, this->UHM->genH.LM->SlocR.data());

    //-----------------------------------------
    // (2) allocate for <phi | T + Vnl | dphi>
    //-----------------------------------------
    this->UHM->LM->DHloc_fixedR_x = new double[nnr];
    this->UHM->LM->DHloc_fixedR_y = new double[nnr];
    this->UHM->LM->DHloc_fixedR_z = new double[nnr];
    const auto init_DHloc_fixedR_xyz = [this, nnr](int num_threads, int thread_id) {
        int beg, len;
        ModuleBase::BLOCK_TASK_DIST_1D(num_threads, thread_id, nnr, 1024, beg, len);
        ModuleBase::GlobalFunc::ZEROS(this->UHM->LM->DHloc_fixedR_x + beg, len);
        ModuleBase::GlobalFunc::ZEROS(this->UHM->LM->DHloc_fixedR_y + beg, len);
        ModuleBase::GlobalFunc::ZEROS(this->UHM->LM->DHloc_fixedR_z + beg, len);
    };
    ModuleBase::OMP_PARALLEL(init_DHloc_fixedR_xyz);
    ModuleBase::Memory::record("Force::dTVNL", sizeof(double) * nnr * 3);

    // calculate dT=<phi|kin|dphi> in LCAO
    // calculate T + VNL(P1) in LCAO basis
    this->UHM->genH.build_ST_new('T', cal_deri, GlobalC::ucell, this->UHM->genH.LM->Hloc_fixedR.data());
    // test(this->UHM->LM->DHloc_fixedR_x,"this->UHM->LM->DHloc_fixedR_x T part");

    // calculate dVnl=<phi|dVnl|dphi> in LCAO
    this->UHM->genH.build_Nonlocal_mu_new (this->UHM->genH.LM->Hloc_fixed.data(), cal_deri);
    // test(this->UHM->LM->DHloc_fixedR_x,"this->UHM->LM->DHloc_fixedR_x Vnl part");

    // calculate asynchronous S matrix to output for Hefei-NAMD
    if (INPUT.cal_syns)
    {
        cal_deri = false;
        //this->UHM->genH.build_ST_new('S', cal_deri, GlobalC::ucell, this->UHM->genH.LM->SlocR.data(), INPUT.cal_syns);
        this->UHM->genH.build_ST_new('S', cal_deri, GlobalC::ucell, this->UHM->genH.LM->SlocR.data(), INPUT.cal_syns, INPUT.dmax);
        for (int ik = 0; ik < nks; ik++)
        {
            this->UHM->genH.LM->zeros_HSk('S');
            this->UHM->genH.LM->folding_fixedH(ik, kvec_d, 1);
            bool bit = false; // LiuXh, 2017-03-21
            ModuleIO::saving_HS(0,
                                this->UHM->genH.LM->Hloc2.data(),
                                this->UHM->genH.LM->Sloc2.data(),
                                bit,
                                1,
                                "data-" + std::to_string(ik),
                                this->ParaV[0],
                                0); // LiuXh, 2017-03-21
        }
    }

    ModuleBase::timer::tick("Force_LCAO_k", "allocate_k");
    return;
}

void Force_LCAO_k::finish_k(void)
{
    delete[] this->UHM->LM->DSloc_Rx;
    delete[] this->UHM->LM->DSloc_Ry;
    delete[] this->UHM->LM->DSloc_Rz;
    delete[] this->UHM->LM->DHloc_fixedR_x;
    delete[] this->UHM->LM->DHloc_fixedR_y;
    delete[] this->UHM->LM->DHloc_fixedR_z;
    if (GlobalV::CAL_STRESS)
    {
        delete[] this->UHM->LM->DH_r;
        delete[] this->UHM->LM->stvnl11;
        delete[] this->UHM->LM->stvnl12;
        delete[] this->UHM->LM->stvnl13;
        delete[] this->UHM->LM->stvnl22;
        delete[] this->UHM->LM->stvnl23;
        delete[] this->UHM->LM->stvnl33;
    }
    return;
}

#include "module_hamilt_lcao/hamilt_lcaodft/record_adj.h"
void Force_LCAO_k::cal_foverlap_k(const bool isforce,
                                  const bool isstress,
                                  Record_adj& ra,
                                  const psi::Psi<std::complex<double>>* psi,
                                  Local_Orbital_Charge& loc,
                                  ModuleBase::matrix& foverlap,
                                  ModuleBase::matrix& soverlap,
                                  const elecstate::ElecState* pelec,
                                  const int& nks,
                                  const K_Vectors& kv
                                  )
{
    ModuleBase::TITLE("Force_LCAO_k", "cal_foverlap_k");
    ModuleBase::timer::tick("Force_LCAO_k", "cal_foverlap_k");

    const Parallel_Orbitals* pv = this->ParaV;
    //--------------------------------------------
    // (1) allocate energy density matrix (nnr)
    //--------------------------------------------
    double** edm2d = new double*[GlobalV::NSPIN];
    for (int is = 0; is < GlobalV::NSPIN; is++)
    {
        edm2d[is] = new double[pv->nnr];
    }
    ModuleBase::Memory::record("Force::edm2d", sizeof(double) * pv->nnr * GlobalV::NSPIN);
    auto init_edm2d = [edm2d, pv](int num_threads, int thread_id)
    {
        int beg, len;
        ModuleBase::BLOCK_TASK_DIST_1D(num_threads, thread_id, pv->nnr, 1024, beg, len);
        for (int is = 0; is < GlobalV::NSPIN; is++)
        {
            ModuleBase::GlobalFunc::ZEROS(edm2d[is] + beg, len);
        }
    };
    ModuleBase::OMP_PARALLEL(init_edm2d);

    //--------------------------------------------
    // calculate the energy density matrix here.
    //--------------------------------------------
    ModuleBase::timer::tick("Force_LCAO_k", "cal_edm_2d");

    ModuleBase::matrix wgEkb;
    wgEkb.create(nks, GlobalV::NBANDS);
    ModuleBase::Memory::record("Force::wgEkb", sizeof(double) * nks * GlobalV::NBANDS);
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static, 1024)
#endif
    for (int ik = 0; ik < nks; ik++)
    {
        for (int ib = 0; ib < GlobalV::NBANDS; ib++)
        {
            wgEkb(ik, ib) = pelec->wg(ik, ib) * pelec->ekb(ik, ib);
        }
    }
    std::vector<ModuleBase::ComplexMatrix> edm_k(nks);

    // use the original formula (Hamiltonian matrix) to calculate energy density matrix
    if (loc.edm_k_tddft.size())
    {
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
        for (int ik = 0; ik < nks; ++ik)
        {
            edm_k[ik] = loc.edm_k_tddft[ik];
        }
    }
    else
    {
        elecstate::cal_dm(loc.ParaV, wgEkb, psi[0], edm_k);
    }

    loc.cal_dm_R(edm_k, ra, edm2d, kv);
    ModuleBase::timer::tick("Force_LCAO_k", "cal_edm_2d");

    //--------------------------------------------
    // summation \sum_{i,j} E(i,j)*dS(i,j)
    // BEGIN CALCULATION OF FORCE OF EACH ATOM
    //--------------------------------------------
    int total_irr = 0;
#ifdef _OPENMP
#pragma omp parallel
{
    int num_threads = omp_get_num_threads();
    ModuleBase::matrix local_soverlap(3, 3);
    int local_total_irr = 0;
#else
    ModuleBase::matrix& local_soverlap = soverlap;
    int& local_total_irr = total_irr;
#endif

    ModuleBase::Vector3<double> tau1, dtau, tau2;

#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
    for (int iat = 0; iat < GlobalC::ucell.nat; iat++)
    {
        const int T1 = GlobalC::ucell.iat2it[iat];
        Atom* atom1 = &GlobalC::ucell.atoms[T1];
        const int I1 = GlobalC::ucell.iat2ia[iat];

        double *foverlap_iat;
        if(isforce) foverlap_iat = &foverlap(iat, 0);

#ifdef _OPENMP
        // using local stack to avoid false sharing in multi-threaded case
        double foverlap_temp[3] = {0.0, 0.0, 0.0};
        if (num_threads > 1)
        {
            foverlap_iat = foverlap_temp;
        }
#endif
        int irr = pv->nlocstart[iat];
        const int start1 = GlobalC::ucell.itiaiw2iwt(T1, I1, 0);
        for (int cb = 0; cb < ra.na_each[iat]; ++cb)
        {
            const int T2 = ra.info[iat][cb][3];
            const int I2 = ra.info[iat][cb][4];
            const int start2 = GlobalC::ucell.itiaiw2iwt(T2, I2, 0);

            Atom* atom2 = &GlobalC::ucell.atoms[T2];

            for (int jj = 0; jj < atom1->nw; jj++)
            {
                const int iw1_all = start1 + jj;

                // HPSEPS
                const int mu = pv->trace_loc_row[iw1_all];
                if (mu < 0)
                    continue;

                for (int kk = 0; kk < atom2->nw; kk++)
                {
                    const int iw2_all = start2 + kk;

                    // HPSEPS
                    const int nu = pv->trace_loc_col[iw2_all];
                    if (nu < 0)
                        continue;
                    //==============================================================
                    // here we use 'minus', but in GlobalV::GAMMA_ONLY_LOCAL we use 'plus',
                    // both are correct because the 'DSloc_Rx' is used in 'row' (-),
                    // however, the 'DSloc_x' in GAMMA is used in 'col' (+),
                    // mohan update 2011-06-16
                    //==============================================================
                    for (int is = 0; is < GlobalV::NSPIN; ++is)
                    {
                        double edm2d2 = 2.0 * edm2d[is][irr];
                        if (isforce)
                        {
                            foverlap_iat[0] -= edm2d2 * this->UHM->LM->DSloc_Rx[irr];
                            foverlap_iat[1] -= edm2d2 * this->UHM->LM->DSloc_Ry[irr];
                            foverlap_iat[2] -= edm2d2 * this->UHM->LM->DSloc_Rz[irr];
                        }
                        if (isstress)
                        {
                            for (int ipol = 0; ipol < 3; ipol++)
                            {
                                local_soverlap(0, ipol) += edm2d[is][irr] * this->UHM->LM->DSloc_Rx[irr]
                                                        * this->UHM->LM->DH_r[irr * 3 + ipol];
                                if (ipol < 1)
                                    continue;
                                local_soverlap(1, ipol) += edm2d[is][irr] * this->UHM->LM->DSloc_Ry[irr]
                                                        * this->UHM->LM->DH_r[irr * 3 + ipol];
                                if (ipol < 2)
                                    continue;
                                local_soverlap(2, ipol) += edm2d[is][irr] * this->UHM->LM->DSloc_Rz[irr]
                                                        * this->UHM->LM->DH_r[irr * 3 + ipol];
                            }
                        }
                    }
                    ++local_total_irr;
                    ++irr;
                } // end kk
            } // end jj
        } // end cb
#ifdef _OPENMP
        if (isforce && num_threads > 1)
        {
            foverlap(iat, 0) += foverlap_iat[0];
            foverlap(iat, 1) += foverlap_iat[1];
            foverlap(iat, 2) += foverlap_iat[2];
        }
#endif
    } // end iat
#ifdef _OPENMP
    #pragma omp critical(cal_foverlap_k_reduce)
    {
        total_irr += local_total_irr;
        if (isstress)
        {
            for (int ipol = 0; ipol < 3; ipol++)
            {
                soverlap(0, ipol) += local_soverlap(0, ipol);
                if (ipol < 1)
                    continue;
                soverlap(1, ipol) += local_soverlap(1, ipol);
                if (ipol < 2)
                    continue;
                soverlap(2, ipol) += local_soverlap(2, ipol);
            }
        }
    }
}
#endif

    if (isstress)
    {
        StressTools::stress_fill(GlobalC::ucell.lat0, GlobalC::ucell.omega, soverlap);
    }

    if (total_irr != pv->nnr)
    {
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "wrong irr", total_irr);
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "wrong LNNR.nnr", pv->nnr);
        ModuleBase::WARNING_QUIT("Force_LCAO_k::cal_foverlap_k", "irr!=LNNR.nnr");
    }

    for (int is = 0; is < GlobalV::NSPIN; is++)
    {
        delete[] edm2d[is];
    }
    delete[] edm2d;

    ModuleBase::timer::tick("Force_LCAO_k", "cal_foverlap_k");
    return;
}

void Force_LCAO_k::cal_ftvnl_dphi_k(double** dm2d,
                                    const bool isforce,
                                    const bool isstress,
                                    Record_adj& ra,
                                    ModuleBase::matrix& ftvnl_dphi,
                                    ModuleBase::matrix& stvnl_dphi)
{
    ModuleBase::TITLE("Force_LCAO_k", "cal_ftvnl_dphi_k");
    ModuleBase::timer::tick("Force_LCAO_k", "cal_ftvnl_dphi_k");

    int total_irr = 0;
    const Parallel_Orbitals* pv = this->ParaV;
    // get the adjacent atom's information.

    //	GlobalV::ofs_running << " calculate the ftvnl_dphi_k force" << std::endl;
#ifdef _OPENMP
#pragma omp parallel
{
    int num_threads = omp_get_num_threads();
    ModuleBase::matrix local_stvnl_dphi(3, 3);
    int local_total_irr = 0;
    #pragma omp for schedule(dynamic)
#else
    ModuleBase::matrix& local_stvnl_dphi = stvnl_dphi;
    int& local_total_irr = total_irr;
#endif
    for (int iat = 0; iat < GlobalC::ucell.nat; iat++)
    {
        const int T1 = GlobalC::ucell.iat2it[iat];
        Atom* atom1 = &GlobalC::ucell.atoms[T1];
        const int I1 = GlobalC::ucell.iat2ia[iat];

        int irr = pv->nlocstart[iat];
        const int start1 = GlobalC::ucell.itiaiw2iwt(T1, I1, 0);
        double *ftvnl_dphi_iat;
        if(isforce) ftvnl_dphi_iat = &ftvnl_dphi(iat, 0);
#ifdef _OPENMP
        // using local stack to avoid false sharing in multi-threaded case
        double ftvnl_dphi_temp[3] = {0.0, 0.0, 0.0};
        if (num_threads > 1)
        {
            ftvnl_dphi_iat = ftvnl_dphi_temp;
        }
#endif
        for (int cb = 0; cb < ra.na_each[iat]; ++cb)
        {
            const int T2 = ra.info[iat][cb][3];
            const int I2 = ra.info[iat][cb][4];
            const int start2 = GlobalC::ucell.itiaiw2iwt(T2, I2, 0);
            Atom* atom2 = &GlobalC::ucell.atoms[T2];

            for (int jj = 0; jj < atom1->nw; ++jj)
            {
                const int iw1_all = start1 + jj;
                const int mu = pv->trace_loc_row[iw1_all];
                if (mu < 0)
                    continue;
                for (int kk = 0; kk < atom2->nw; ++kk)
                {
                    const int iw2_all = start2 + kk;
                    const int nu = pv->trace_loc_col[iw2_all];
                    if (nu < 0)
                        continue;
                    //==============================================================
                    // here we use 'minus', but in GlobalV::GAMMA_ONLY_LOCAL we use 'plus',
                    // both are correct because the 'DSloc_Rx' is used in 'row' (-),
                    // however, the 'DSloc_x' is used in 'col' (+),
                    // mohan update 2011-06-16
                    //==============================================================
                    for (int is = 0; is < GlobalV::NSPIN; ++is)
                    {
                        double dm2d2 = 2.0 * dm2d[is][irr];
                        if (isforce)
                        {
                            ftvnl_dphi_iat[0] += dm2d2 * this->UHM->LM->DHloc_fixedR_x[irr];
                            ftvnl_dphi_iat[1] += dm2d2 * this->UHM->LM->DHloc_fixedR_y[irr];
                            ftvnl_dphi_iat[2] += dm2d2 * this->UHM->LM->DHloc_fixedR_z[irr];
                        }
                        if (isstress)
                        {
                            local_stvnl_dphi(0, 0) -= dm2d[is][irr] * this->UHM->LM->stvnl11[irr];
                            local_stvnl_dphi(0, 1) -= dm2d[is][irr] * this->UHM->LM->stvnl12[irr];
                            local_stvnl_dphi(0, 2) -= dm2d[is][irr] * this->UHM->LM->stvnl13[irr];
                            local_stvnl_dphi(1, 1) -= dm2d[is][irr] * this->UHM->LM->stvnl22[irr];
                            local_stvnl_dphi(1, 2) -= dm2d[is][irr] * this->UHM->LM->stvnl23[irr];
                            local_stvnl_dphi(2, 2) -= dm2d[is][irr] * this->UHM->LM->stvnl33[irr];
                        }
                    }
                    ++local_total_irr;
                    ++irr;
                } // end kk
            } // end jj
        } // end cb
#ifdef _OPENMP
        if (isforce && num_threads > 1)
        {
            ftvnl_dphi(iat, 0) += ftvnl_dphi_iat[0];
            ftvnl_dphi(iat, 1) += ftvnl_dphi_iat[1];
            ftvnl_dphi(iat, 2) += ftvnl_dphi_iat[2];
        }
#endif
    } // end iat
#ifdef _OPENMP
    #pragma omp critical(cal_ftvnl_dphi_k_reduce)
    {
        total_irr += local_total_irr;
        if (isstress)
        {
            stvnl_dphi(0, 0) += local_stvnl_dphi(0, 0);
            stvnl_dphi(0, 1) += local_stvnl_dphi(0, 1);
            stvnl_dphi(0, 2) += local_stvnl_dphi(0, 2);
            stvnl_dphi(1, 1) += local_stvnl_dphi(1, 1);
            stvnl_dphi(1, 2) += local_stvnl_dphi(1, 2);
            stvnl_dphi(2, 2) += local_stvnl_dphi(2, 2);
        }
    }
}
#endif
    assert(total_irr == pv->nnr);

    //	test(this->UHM->LM->DSloc_Rx);
    //	test(dm2d[0],"dm2d");

    if (isstress)
    {
        StressTools::stress_fill(GlobalC::ucell.lat0, GlobalC::ucell.omega, stvnl_dphi);
    }

    ModuleBase::timer::tick("Force_LCAO_k", "cal_ftvnl_dphi_k");
    return;
}

void Force_LCAO_k::test(double* mmm, const std::string& name)
{
    const Parallel_Orbitals* pv = this->ParaV;
    if (GlobalV::NPROC != 1)
        return;
    std::cout << "test!" << std::endl;

    int irr = 0;
    int ca = 0;

    GlobalV::ofs_running << " Calculate the test in Force_LCAO_k" << std::endl;
    Record_adj RA;
    RA.for_2d(*this->UHM->LM->ParaV, GlobalV::GAMMA_ONLY_LOCAL);

    double* test;
    test = new double[GlobalV::NLOCAL * GlobalV::NLOCAL];
    ModuleBase::GlobalFunc::ZEROS(test, GlobalV::NLOCAL * GlobalV::NLOCAL);

    for (int T1 = 0; T1 < GlobalC::ucell.ntype; T1++)
    {
        Atom* atom1 = &GlobalC::ucell.atoms[T1];
        for (int I1 = 0; I1 < atom1->na; I1++)
        {
            // const int iat = GlobalC::ucell.itia2iat(T1,I1);
            const int start1 = GlobalC::ucell.itiaiw2iwt(T1, I1, 0);
            for (int cb = 0; cb < RA.na_each[ca]; cb++)
            {
                const int T2 = RA.info[ca][cb][3];
                const int I2 = RA.info[ca][cb][4];
                Atom* atom2 = &GlobalC::ucell.atoms[T2];
                const int start2 = GlobalC::ucell.itiaiw2iwt(T2, I2, 0);

                for (int jj = 0; jj < atom1->nw; jj++)
                {
                    const int iw1_all = start1 + jj;
                    for (int kk = 0; kk < atom2->nw; kk++)
                    {
                        const int iw2_all = start2 + kk;
                        assert(irr < pv->nnr);
                        // test[iw1_all*GlobalV::NLOCAL+iw2_all] += this->UHM->LM->DHloc_fixedR_x[irr];
                        test[iw1_all * GlobalV::NLOCAL + iw2_all] += mmm[irr];
                        ++irr;
                    }
                }
            }
            ++ca;
        }
    }

    std::cout << "\n " << name << std::endl;
    std::cout << std::setprecision(4);
    for (int i = 0; i < GlobalV::NLOCAL; i++)
    {
        for (int j = 0; j < GlobalV::NLOCAL; j++)
        {
            if (std::abs(test[i * GlobalV::NLOCAL + j]) > 1.0e-5)
                std::cout << std::setw(12) << test[i * GlobalV::NLOCAL + j];
            else
                std::cout << std::setw(12) << "0";
        }
        std::cout << std::endl;
    }
    delete[] test;

    RA.delete_grid(); // xiaohui add 2015-02-04
    return;
}

typedef std::tuple<int, int, int, int> key_tuple;

// must consider three-center H matrix.
void Force_LCAO_k::cal_fvnl_dbeta_k_new(double** dm2d,
                                        const bool isforce,
                                        const bool isstress,
                                        ModuleBase::matrix& fvnl_dbeta,
                                        ModuleBase::matrix& svnl_dbeta)
{
    ModuleBase::TITLE("Force_LCAO_k", "cal_fvnl_dbeta_k_new");
    ModuleBase::timer::tick("Force_LCAO_k", "cal_fvnl_dbeta_k_new");
    const Parallel_Orbitals* pv = this->ParaV;

    // Data structure for storing <psi|beta>, for a detailed description
    // check out the same data structure in build_Nonlocal_mu_new
    std::vector<std::map<key_tuple, std::unordered_map<int, std::vector<std::vector<double>>>>> nlm_tot;

    nlm_tot.resize(GlobalC::ucell.nat);

#ifdef _OPENMP
// use schedule(dynamic) for load balancing because adj_num is various
#pragma omp parallel for schedule(dynamic)
#endif
    for (int iat = 0; iat < GlobalC::ucell.nat; iat++)
    {

        const int it = GlobalC::ucell.iat2it[iat];
        const int ia = GlobalC::ucell.iat2ia[iat];

        // Step 1 : generate <psi|beta>
        // type of atom; distance; atomic basis; projectors

        const double Rcut_Beta = GlobalC::ucell.infoNL.Beta[it].get_rcut_max();
        const ModuleBase::Vector3<double> tau = GlobalC::ucell.atoms[it].tau[ia];
        AdjacentAtomInfo adjs;
        GlobalC::GridD.Find_atom(GlobalC::ucell, tau, it, ia, &adjs);

        nlm_tot[iat].clear();

        for (int ad = 0; ad < adjs.adj_num + 1; ++ad)
        {
            const int T1 = adjs.ntype[ad];
            const int I1 = adjs.natom[ad];
            const int start1 = GlobalC::ucell.itiaiw2iwt(T1, I1, 0);
            const double Rcut_AO1 = GlobalC::ORB.Phi[T1].getRcut();

            const ModuleBase::Vector3<double> &tau1 = adjs.adjacent_tau[ad];
            const Atom* atom1 = &GlobalC::ucell.atoms[T1];
            const int nw1_tot = atom1->nw * GlobalV::NPOL;

            const ModuleBase::Vector3<double> dtau = tau1 - tau;
            const double dist1 = dtau.norm2() * pow(GlobalC::ucell.lat0, 2);
            if (dist1 > pow(Rcut_Beta + Rcut_AO1, 2))
            {
                continue;
            }

            std::unordered_map<int, std::vector<std::vector<double>>> nlm_cur;
            nlm_cur.clear();

            for (int iw1 = 0; iw1 < nw1_tot; ++iw1)
            {
                const int iw1_all = start1 + iw1;
                const int iw1_local = pv->trace_loc_row[iw1_all];
                const int iw2_local = pv->trace_loc_col[iw1_all];
                if (iw1_local < 0 && iw2_local < 0)
                    continue;
                const int iw1_0 = iw1 / GlobalV::NPOL;
                std::vector<std::vector<double>> nlm;
                GlobalC::UOT.snap_psibeta_half(GlobalC::ORB,
                                               GlobalC::ucell.infoNL,
                                               nlm,
                                               tau1,
                                               T1,
                                               atom1->iw2l[iw1_0], // L1
                                               atom1->iw2m[iw1_0], // m1
                                               atom1->iw2n[iw1_0], // N1
                                               tau,
                                               it,
                                               1); // R0,T0

                nlm_cur.insert({iw1_all, nlm});
            } // end iw
            const int iat1 = GlobalC::ucell.itia2iat(T1, I1);
            const int rx1 = adjs.box[ad].x;
            const int ry1 = adjs.box[ad].y;
            const int rz1 = adjs.box[ad].z;
            key_tuple key_1(iat1, rx1, ry1, rz1);
            nlm_tot[iat][key_1] = nlm_cur;
        } // end ad
    }

    //=======================================================
    // Step2:
    // calculate sum_(L0,M0) beta<psi_i|beta><beta|psi_j>
    // and accumulate the value to Hloc_fixedR(i,j)
    //=======================================================
    int total_nnr = 0;
#ifdef _OPENMP
#pragma omp parallel reduction(+:total_nnr)
{
    ModuleBase::matrix local_svnl_dbeta(3, 3);
    const int num_threads = omp_get_num_threads();
#else
    ModuleBase::matrix &local_svnl_dbeta = svnl_dbeta;
#endif

    ModuleBase::Vector3<double> tau1;
    ModuleBase::Vector3<double> tau2;
    ModuleBase::Vector3<double> dtau;
    ModuleBase::Vector3<double> tau0;
    ModuleBase::Vector3<double> dtau1;
    ModuleBase::Vector3<double> dtau2;

    double rcut;
    double distance;

    double rcut1;
    double rcut2;
    double distance1;
    double distance2;

#ifdef _OPENMP
// use schedule(dynamic) for load balancing because adj_num is various
#pragma omp for schedule(dynamic)
#endif
    for (int iat1 = 0; iat1 < GlobalC::ucell.nat; iat1++)
    {
        const int T1 = GlobalC::ucell.iat2it[iat1];
        const Atom* atom1 = &GlobalC::ucell.atoms[T1];

        {
            const int I1 = GlobalC::ucell.iat2ia[iat1];
            tau1 = atom1->tau[I1];
            AdjacentAtomInfo adjs;
            GlobalC::GridD.Find_atom(GlobalC::ucell, tau1, T1, I1, &adjs);
            const int start1 = GlobalC::ucell.itiaiw2iwt(T1, I1, 0);
            int nnr = this->ParaV->nlocstart[iat1];

            /*
                !!!!!!!!!!!!!!!!
                This optimization is also improving the performance of single thread.
                Making memory access more linearly in the core loop
            */
            bool iat_recorded = false;
            bool force_updated = false;
            // record iat of adjs
            std::vector<int> adj_iat;
            // record fvnl_dbeta diff of adjs
            std::vector<double> adj_fvnl_dbeta;
            if (isforce)
            {
                adj_iat.resize(adjs.adj_num + 1);
                adj_fvnl_dbeta.resize((adjs.adj_num + 1) * 3, 0.0);
            }

            for (int ad2 = 0; ad2 < adjs.adj_num + 1; ++ad2)
            {
                const int T2 = adjs.ntype[ad2];
                const Atom* atom2 = &GlobalC::ucell.atoms[T2];

                const int I2 = adjs.natom[ad2];
                const int iat2 = GlobalC::ucell.itia2iat(T2, I2);
                const int start2 = GlobalC::ucell.itiaiw2iwt(T2, I2, 0);
                tau2 = adjs.adjacent_tau[ad2];

                const int rx2 = adjs.box[ad2].x;
                const int ry2 = adjs.box[ad2].y;
                const int rz2 = adjs.box[ad2].z;

                dtau = tau2 - tau1;
                distance = dtau.norm2() * pow(GlobalC::ucell.lat0, 2);
                // this rcut is in order to make nnr consistent
                // with other matrix.
                rcut = pow(GlobalC::ORB.Phi[T1].getRcut() + GlobalC::ORB.Phi[T2].getRcut(), 2);

                // check if this a adjacent atoms.
                bool is_adj = false;
                if (distance < rcut)
                    is_adj = true;
                else if (distance >= rcut)
                {
                    for (int ad0 = 0; ad0 < adjs.adj_num + 1; ++ad0)
                    {
                        const int T0 = adjs.ntype[ad0];
                        if (GlobalC::ucell.infoNL.nproj[T0] == 0)
                            continue;
                        const int I0 = adjs.natom[ad0];
                        // const int iat0 = GlobalC::ucell.itia2iat(T0, I0);
                        // const int start0 = GlobalC::ucell.itiaiw2iwt(T0, I0, 0);

                        tau0 = adjs.adjacent_tau[ad0];
                        dtau1 = tau0 - tau1;
                        distance1 = dtau1.norm() * GlobalC::ucell.lat0;
                        rcut1 = GlobalC::ORB.Phi[T1].getRcut() + GlobalC::ucell.infoNL.Beta[T0].get_rcut_max();

                        dtau2 = tau0 - tau2;
                        distance2 = dtau2.norm() * GlobalC::ucell.lat0;
                        rcut2 = GlobalC::ORB.Phi[T2].getRcut() + GlobalC::ucell.infoNL.Beta[T0].get_rcut_max();

                        if (distance1 < rcut1 && distance2 < rcut2)
                        {
                            is_adj = true;
                            break;
                        }
                    }
                }

                if (is_adj)
                {
                    for (int ad0 = 0; ad0 < adjs.adj_num + 1; ++ad0)
                    {
                        const int T0 = adjs.ntype[ad0];
                        const int I0 = adjs.natom[ad0];
                        const int iat = GlobalC::ucell.itia2iat(T0, I0);
                        if (!iat_recorded && isforce) adj_iat[ad0] = iat;

                        // mohan add 2010-12-19
                        if (GlobalC::ucell.infoNL.nproj[T0] == 0)
                            continue;

                        // const int I0 = GlobalC::GridD.getNatom(ad0);
                        // const int start0 = GlobalC::ucell.itiaiw2iwt(T0, I0, 0);
                        tau0 = adjs.adjacent_tau[ad0];

                        dtau1 = tau0 - tau1;
                        dtau2 = tau0 - tau2;
                        const double distance1 = dtau1.norm2() * pow(GlobalC::ucell.lat0, 2);
                        const double distance2 = dtau2.norm2() * pow(GlobalC::ucell.lat0, 2);

                        // seems a bug here!! mohan 2011-06-17
                        rcut1 = pow(GlobalC::ORB.Phi[T1].getRcut() + GlobalC::ucell.infoNL.Beta[T0].get_rcut_max(), 2);
                        rcut2 = pow(GlobalC::ORB.Phi[T2].getRcut() + GlobalC::ucell.infoNL.Beta[T0].get_rcut_max(), 2);

                        double r0[3];
                        double r1[3];
                        r1[0] = (tau1.x - tau0.x);
                        r1[1] = (tau1.y - tau0.y);
                        r1[2] = (tau1.z - tau0.z);
                        r0[0] = (tau2.x - tau0.x);
                        r0[1] = (tau2.y - tau0.y);
                        r0[2] = (tau2.z - tau0.z);

                        if (distance1 >= rcut1 || distance2 >= rcut2)
                        {
                            continue;
                        }

                        const int rx0 = adjs.box[ad0].x;
                        const int ry0 = adjs.box[ad0].y;
                        const int rz0 = adjs.box[ad0].z;
                        key_tuple key1(iat1, -rx0, -ry0, -rz0);
                        key_tuple key2(iat2, rx2 - rx0, ry2 - ry0, rz2 - rz0);

                        int nnr_inner = 0;
                        for (int j = 0; j < atom1->nw * GlobalV::NPOL; j++)
                        {
                            const int j0 = j / GlobalV::NPOL; // added by zhengdy-soc
                            const int iw1_all = start1 + j;
                            const int mu = pv->trace_loc_row[iw1_all];
                            if (mu < 0)
                                continue;

                            for (int k = 0; k < atom2->nw * GlobalV::NPOL; k++)
                            {
                                const int k0 = k / GlobalV::NPOL;
                                const int iw2_all = start2 + k;
                                const int nu = pv->trace_loc_col[iw2_all];
                                if (nu < 0)
                                    continue;

                                // const Atom* atom0 = &GlobalC::ucell.atoms[T0];
                                double nlm[3] = {0, 0, 0};
                                std::vector<double> nlm_1 = nlm_tot[iat][key2][iw2_all][0];
                                std::vector<std::vector<double>> nlm_2;
                                nlm_2.resize(3);
                                for (int i = 0; i < 3; i++)
                                {
                                    nlm_2[i] = nlm_tot[iat][key1][iw1_all][i + 1];
                                }

                                assert(nlm_1.size() == nlm_2[0].size());

                                const int nproj = GlobalC::ucell.infoNL.nproj[T0];
                                int ib = 0;
                                for (int nb = 0; nb < nproj; nb++)
                                {
                                    const int L0 = GlobalC::ucell.infoNL.Beta[T0].Proj[nb].getL();
                                    for (int m = 0; m < 2 * L0 + 1; m++)
                                    {
                                        for (int ir = 0; ir < 3; ir++)
                                        {
                                            nlm[ir]
                                                += nlm_2[ir][ib] * nlm_1[ib] * GlobalC::ucell.atoms[T0].ncpp.dion(nb, nb);
                                        }
                                        ib += 1;
                                    }
                                }
                                assert(ib == nlm_1.size());

                                double nlm1[3] = {0, 0, 0};
                                if (isstress)
                                {
                                    std::vector<double> nlm_1 = nlm_tot[iat][key1][iw1_all][0];
                                    std::vector<std::vector<double>> nlm_2;
                                    nlm_2.resize(3);
                                    for (int i = 0; i < 3; i++)
                                    {
                                        nlm_2[i] = nlm_tot[iat][key2][iw2_all][i + 1];
                                    }

                                    assert(nlm_1.size() == nlm_2[0].size());

                                    const int nproj = GlobalC::ucell.infoNL.nproj[T0];
                                    int ib = 0;
                                    for (int nb = 0; nb < nproj; nb++)
                                    {
                                        const int L0 = GlobalC::ucell.infoNL.Beta[T0].Proj[nb].getL();
                                        for (int m = 0; m < 2 * L0 + 1; m++)
                                        {
                                            for (int ir = 0; ir < 3; ir++)
                                            {
                                                nlm1[ir] += nlm_2[ir][ib] * nlm_1[ib]
                                                            * GlobalC::ucell.atoms[T0].ncpp.dion(nb, nb);
                                            }
                                            ib += 1;
                                        }
                                    }
                                    assert(ib == nlm_1.size());
                                }

                                /// only one projector for each atom force, but another projector for stress
                                force_updated = true;
                                for (int is = 0; is < GlobalV::NSPIN; ++is)
                                {
                                    double dm2d2 = 2.0 * dm2d[is][nnr + nnr_inner];
                                    for (int jpol = 0; jpol < 3; jpol++)
                                    {
                                        if (isforce)
                                        {
                                            adj_fvnl_dbeta[ad0 * 3 + jpol] -= dm2d2 * nlm[jpol];
                                        }
                                        if (isstress)
                                        {
                                            for (int ipol = jpol; ipol < 3; ipol++)
                                            {
                                                local_svnl_dbeta(jpol, ipol)
                                                    += dm2d[is][nnr + nnr_inner]
                                                       * (nlm[jpol] * r1[ipol] + nlm1[jpol] * r0[ipol]);
                                            }
                                        }
                                    }
                                }
                                nnr_inner++;
                            } // k
                        } // j
                    } // ad0

                    // outer circle : accumulate nnr
                    for (int j = 0; j < atom1->nw * GlobalV::NPOL; j++)
                    {
                        const int j0 = j / GlobalV::NPOL; // added by zhengdy-soc
                        const int iw1_all = start1 + j;
                        const int mu = pv->trace_loc_row[iw1_all];
                        if (mu < 0)
                            continue;

                        // fix a serious bug: atom2[T2] -> atom2
                        // mohan 2010-12-20
                        for (int k = 0; k < atom2->nw * GlobalV::NPOL; k++)
                        {
                            const int k0 = k / GlobalV::NPOL;
                            const int iw2_all = start2 + k;
                            const int nu = pv->trace_loc_col[iw2_all];
                            if (nu < 0)
                                continue;
                            total_nnr++;
                            nnr++;
                        }
                    }
                    iat_recorded = true;
                } // is_adj
            } // ad2

            // sum the diff to fvnl_dbeta
            if (force_updated && isforce)
            {
#ifdef _OPENMP
                if (num_threads > 1)
                {
                    for (int ad0 = 0; ad0 < adjs.adj_num + 1; ++ad0)
                    {
                        #pragma omp atomic
                        fvnl_dbeta(adj_iat[ad0], 0) += adj_fvnl_dbeta[ad0 * 3 + 0];
                        #pragma omp atomic
                        fvnl_dbeta(adj_iat[ad0], 1) += adj_fvnl_dbeta[ad0 * 3 + 1];
                        #pragma omp atomic
                        fvnl_dbeta(adj_iat[ad0], 2) += adj_fvnl_dbeta[ad0 * 3 + 2];
                    }
                }
                else
#endif
                {
                    for (int ad0 = 0; ad0 < adjs.adj_num + 1; ++ad0)
                    {
                        fvnl_dbeta(adj_iat[ad0], 0) += adj_fvnl_dbeta[ad0 * 3 + 0];
                        fvnl_dbeta(adj_iat[ad0], 1) += adj_fvnl_dbeta[ad0 * 3 + 1];
                        fvnl_dbeta(adj_iat[ad0], 2) += adj_fvnl_dbeta[ad0 * 3 + 2];
                    }
                }
            }
        } // I1
    } // T1

#ifdef _OPENMP
    if (isstress)
    {
        #pragma omp critical(cal_fvnl_dbeta_k_new_reduce)
        {
            for(int l=0;l<3;l++)
            {
                for(int m=0;m<3;m++)
                {
                    svnl_dbeta(l,m) += local_svnl_dbeta(l, m);
                }
            }
        }
    }
}
#endif

    assert(total_nnr == pv->nnr);

    if (isstress)
    {
        StressTools::stress_fill(GlobalC::ucell.lat0, GlobalC::ucell.omega, svnl_dbeta);
    }

    ModuleBase::timer::tick("Force_LCAO_k", "cal_fvnl_dbeta_k_new");
    return;
}

// calculate the force due to < phi | Vlocal | dphi >
void Force_LCAO_k::cal_fvl_dphi_k(const bool isforce,
                                  const bool isstress,
                                  const elecstate::Potential* pot_in,
                                  ModuleBase::matrix& fvl_dphi,
                                  ModuleBase::matrix& svl_dphi,
                                  double** DM_R)
{
    ModuleBase::TITLE("Force_LCAO_k", "cal_fvl_dphi_k");
    ModuleBase::timer::tick("Force_LCAO_k", "cal_fvl_dphi_k");

    if (!isforce && !isstress)
    {
        ModuleBase::timer::tick("Force_LCAO_k", "cal_fvl_dphi_k");
        return;
    }
    assert(this->UHM->LM->DHloc_fixedR_x != NULL);
    assert(this->UHM->LM->DHloc_fixedR_y != NULL);
    assert(this->UHM->LM->DHloc_fixedR_z != NULL);

    int istep = 1;

    for (int is = 0; is < GlobalV::NSPIN; ++is)
    {
        GlobalV::CURRENT_SPIN = is;

        const double* vr_eff1 = pot_in->get_effective_v(GlobalV::CURRENT_SPIN);
        const double* vofk_eff1 = nullptr;
        if(XC_Functional::get_func_type()==3 || XC_Functional::get_func_type()==5)
        {
            vofk_eff1 = pot_in->get_effective_vofk(GlobalV::CURRENT_SPIN);
        }

        //--------------------------------
        // Grid integration here.
        //--------------------------------
        // fvl_dphi can not be set to zero here if Vna is used
        if (isstress || isforce)
        {
            if (XC_Functional::get_func_type() == 3 || XC_Functional::get_func_type() == 5)
            {
                Gint_inout inout(DM_R,
                                 vr_eff1,
                                 vofk_eff1,
                                 isforce,
                                 isstress,
                                 &fvl_dphi,
                                 &svl_dphi,
                                 Gint_Tools::job_type::force_meta);
                this->UHM->GK.cal_gint(&inout);
            }
            else
            {
                Gint_inout inout(DM_R,
                                 vr_eff1,
                                 isforce,
                                 isstress,
                                 &fvl_dphi,
                                 &svl_dphi,
                                 Gint_Tools::job_type::force);
                this->UHM->GK.cal_gint(&inout);
            }
        }
    }

    if (isstress)
    {
        StressTools::stress_fill(-1.0, GlobalC::ucell.omega, svl_dphi);
    }

    ModuleBase::timer::tick("Force_LCAO_k", "cal_fvl_dphi_k");
    return;
}
