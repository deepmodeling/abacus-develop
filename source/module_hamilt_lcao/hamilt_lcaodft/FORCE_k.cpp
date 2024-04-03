#include "FORCE_k.h"

#include <map>
#include <unordered_map>

#include "module_base/memory.h"
#include "module_base/parallel_reduce.h"
#include "module_base/timer.h"
#include "module_base/tool_threading.h"
#include "module_cell/module_neighbor/sltk_grid_driver.h"
#include "module_elecstate/cal_dm.h"
#include "module_elecstate/module_dm/cal_dm_psi.h"
#include "module_elecstate/elecstate_lcao.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_io/write_HS.h"

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
                            LCAO_Hamilt &uhm,
                            Gint_k &gint_k,
                            Parallel_Orbitals &pv,
                            LCAO_Matrix &lm,
                            const K_Vectors& kv)
{
    ModuleBase::TITLE("Force_LCAO_k", "ftable_k");
    ModuleBase::timer::tick("Force_LCAO_k", "ftable_k");

    elecstate::DensityMatrix<complex<double>,double>* DM
    = dynamic_cast<const elecstate::ElecStateLCAO<std::complex<double>>*>(pelec)->get_DM();

	this->allocate_k(
			pv, 
			lm, 
            uhm.genH,
			kv.nks, 
			kv.kvec_d);

    // calculate the energy density matrix
    // and the force related to overlap matrix and energy density matrix.
	this->cal_foverlap_k(
			isforce, 
			isstress, 
			ra, 
			psi, 
			loc, 
			pv, 
			lm, 
			DM, 
			foverlap, 
			soverlap, 
			pelec, 
			kv.nks, 
			kv);

	this->cal_ftvnl_dphi_k(
			DM, 
			pv, 
			lm, 
			isforce, 
			isstress, 
			ra, 
			ftvnl_dphi, 
			stvnl_dphi);

	// doing on the real space grid.
	this->cal_fvl_dphi_k(
			isforce, 
			isstress, 
			lm, 
			gint_k, 
			pelec->pot, 
			fvl_dphi, 
			svl_dphi, 
			loc.DM_R);

	this->cal_fvnl_dbeta_k(
			DM, 
			isforce, 
			isstress, 
			pv, 
			fvnl_dbeta, 
			svnl_dbeta);

#ifdef __DEEPKS
    if (GlobalV::deepks_scf)
    {
        const std::vector<std::vector<std::complex<double>>>& dm_k = DM->get_DMK_vector();
        GlobalC::ld.cal_projected_DM_k(DM, GlobalC::ucell, GlobalC::ORB, GlobalC::GridD, kv.nks, kv.kvec_d);
        GlobalC::ld.cal_descriptor();
        GlobalC::ld.cal_gedm(GlobalC::ucell.nat);

        GlobalC::ld.cal_f_delta_k(dm_k,
                                  GlobalC::ucell,
                                  GlobalC::ORB,
                                  GlobalC::GridD,
                                  kv.nks,
                                  kv.kvec_d,
                                  isstress,
                                  svnl_dalpha);
#ifdef __MPI
        Parallel_Reduce::reduce_all(GlobalC::ld.F_delta.c, GlobalC::ld.F_delta.nr * GlobalC::ld.F_delta.nc);
        if (isstress)
        {
            Parallel_Reduce::reduce_pool(svnl_dalpha.c, svnl_dalpha.nr * svnl_dalpha.nc);
        }
#endif
        /*if (GlobalV::deepks_out_unittest)
        {
            GlobalC::ld.print_dm_k(kv.nks, dm_k);
            GlobalC::ld.check_projected_dm();
            GlobalC::ld.check_descriptor(GlobalC::ucell);
            GlobalC::ld.check_gedm();
            GlobalC::ld.add_v_delta_k(GlobalC::ucell, GlobalC::ORB, GlobalC::GridD, pv->nnr);
            GlobalC::ld.check_v_delta_k(pv->nnr);
            for (int ik = 0; ik < kv.nks; ik++)
            {
                uhm.LM->folding_fixedH(ik, kv.kvec_d);
            }
            GlobalC::ld.cal_e_delta_band_k(dm_k, kv.nks);
            std::ofstream ofs("E_delta_bands.dat");
            ofs << std::setprecision(10) << GlobalC::ld.e_delta_band;
            std::ofstream ofs1("E_delta.dat");
            ofs1 << std::setprecision(10) << GlobalC::ld.E_delta;
            GlobalC::ld.check_f_delta(GlobalC::ucell.nat, svnl_dalpha);
        }*/
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
    }
    if (isstress)
    {
        Parallel_Reduce::reduce_pool(soverlap.c, soverlap.nr * soverlap.nc);
        Parallel_Reduce::reduce_pool(stvnl_dphi.c, stvnl_dphi.nr * stvnl_dphi.nc);
        Parallel_Reduce::reduce_pool(svnl_dbeta.c, svnl_dbeta.nr * svnl_dbeta.nc);
        Parallel_Reduce::reduce_pool(svl_dphi.c, svl_dphi.nr * svl_dphi.nc);
    }

    ModuleBase::timer::tick("Force_LCAO_k", "ftable_k");
    return;
}

void Force_LCAO_k::allocate_k(const Parallel_Orbitals& pv,
                              LCAO_Matrix &lm,
                              LCAO_gen_fixedH &gen_h,
                              const int& nks,
                              const std::vector<ModuleBase::Vector3<double>>& kvec_d)
{
    ModuleBase::TITLE("Force_LCAO_k", "allocate_k");
    ModuleBase::timer::tick("Force_LCAO_k", "allocate_k");

    const int nnr = pv.nnr;

    assert(nnr>=0);

    //--------------------------------
    // (1) allocate for dSx dSy & dSz
    //--------------------------------
    lm.DSloc_Rx = new double[nnr];
    lm.DSloc_Ry = new double[nnr];
    lm.DSloc_Rz = new double[nnr];

    // mohan add lm on 2024-03-31
	const auto init_DSloc_Rxyz = [this, nnr, &lm](int num_threads, int thread_id) 
	{
		int beg=0;
		int len=0;
		ModuleBase::BLOCK_TASK_DIST_1D(num_threads, thread_id, nnr, 1024, beg, len);
		ModuleBase::GlobalFunc::ZEROS(lm.DSloc_Rx + beg, len);
		ModuleBase::GlobalFunc::ZEROS(lm.DSloc_Ry + beg, len);
		ModuleBase::GlobalFunc::ZEROS(lm.DSloc_Rz + beg, len);
	};

    ModuleBase::OMP_PARALLEL(init_DSloc_Rxyz);
    ModuleBase::Memory::record("Force::dS_K", sizeof(double) * nnr * 3);

    if (GlobalV::CAL_STRESS)
    {
        lm.DH_r = new double[3 * nnr];
        lm.stvnl11 = new double[nnr];
        lm.stvnl12 = new double[nnr];
        lm.stvnl13 = new double[nnr];
        lm.stvnl22 = new double[nnr];
        lm.stvnl23 = new double[nnr];
        lm.stvnl33 = new double[nnr];
        // mohan add lm on 2024-03-31
        const auto init_DH_r_stvnl = [this, nnr, &lm](int num_threads, int thread_id) {
            int beg, len;
            ModuleBase::BLOCK_TASK_DIST_1D(num_threads, thread_id, nnr, 1024, beg, len);
            ModuleBase::GlobalFunc::ZEROS(lm.DH_r + 3 * beg, 3 * len);
            ModuleBase::GlobalFunc::ZEROS(lm.stvnl11 + beg, len);
            ModuleBase::GlobalFunc::ZEROS(lm.stvnl12 + beg, len);
            ModuleBase::GlobalFunc::ZEROS(lm.stvnl13 + beg, len);
            ModuleBase::GlobalFunc::ZEROS(lm.stvnl22 + beg, len);
            ModuleBase::GlobalFunc::ZEROS(lm.stvnl23 + beg, len);
            ModuleBase::GlobalFunc::ZEROS(lm.stvnl33 + beg, len);
        };
        ModuleBase::OMP_PARALLEL(init_DH_r_stvnl);

        ModuleBase::Memory::record("Stress::dHr", sizeof(double) * nnr * 3);
        ModuleBase::Memory::record("Stress::dSR", sizeof(double) * nnr * 6);
    }

    //-----------------------------
    // calculate dS = <phi | dphi>
    //-----------------------------
    bool cal_deri = true;
    gen_h.build_ST_new('S', cal_deri, GlobalC::ucell, gen_h.LM->SlocR.data());

    //-----------------------------------------
    // (2) allocate for <phi | T + Vnl | dphi>
    //-----------------------------------------
    lm.DHloc_fixedR_x = new double[nnr];
    lm.DHloc_fixedR_y = new double[nnr];
    lm.DHloc_fixedR_z = new double[nnr];

    // mohan add lm on 2024-03-31
    const auto init_DHloc_fixedR_xyz = [this, nnr, &lm](int num_threads, int thread_id) {
        int beg=0;
        int len=0;
        ModuleBase::BLOCK_TASK_DIST_1D(num_threads, thread_id, nnr, 1024, beg, len);
        ModuleBase::GlobalFunc::ZEROS(lm.DHloc_fixedR_x + beg, len);
        ModuleBase::GlobalFunc::ZEROS(lm.DHloc_fixedR_y + beg, len);
        ModuleBase::GlobalFunc::ZEROS(lm.DHloc_fixedR_z + beg, len);
    };
    ModuleBase::OMP_PARALLEL(init_DHloc_fixedR_xyz);
    ModuleBase::Memory::record("Force::dTVNL", sizeof(double) * nnr * 3);

    // calculate dT=<phi|kin|dphi> in LCAO
    // calculate T + VNL(P1) in LCAO basis
    gen_h.build_ST_new('T', cal_deri, GlobalC::ucell, gen_h.LM->Hloc_fixedR.data());

    // calculate dVnl=<phi|dVnl|dphi> in LCAO
    gen_h.build_Nonlocal_mu_new(gen_h.LM->Hloc_fixed.data(), cal_deri);

    // calculate asynchronous S matrix to output for Hefei-NAMD
    if (INPUT.cal_syns)
    {
        cal_deri = false;
        // gen_h.build_ST_new('S', cal_deri, GlobalC::ucell, gen_h.LM->SlocR.data(),
		// INPUT.cal_syns);
		gen_h.build_ST_new('S', 
				cal_deri, 
				GlobalC::ucell, 
				lm.SlocR.data(), 
				INPUT.cal_syns, 
				INPUT.dmax);

        for (int ik = 0; ik < nks; ik++)
        {
            lm.zeros_HSk('S');
            lm.folding_fixedH(ik, kvec_d, 1);
            bool bit = false; // LiuXh, 2017-03-21
			ModuleIO::save_mat(0, 
					lm.Hloc2.data(), 
					GlobalV::NLOCAL, 
					bit, 
					GlobalV::out_ndigits, 
					0, 
					GlobalV::out_app_flag, 
					"H", 
					"data-" + std::to_string(ik), 
					pv, 
					GlobalV::DRANK);

            ModuleIO::save_mat(0, 
					lm.Sloc2.data(), 
					GlobalV::NLOCAL, 
					bit, 
					GlobalV::out_ndigits, 
					0, 
					GlobalV::out_app_flag, 
					"S", 
					"data-" + std::to_string(ik), 
					pv, 
					GlobalV::DRANK);
        }
    }

    ModuleBase::timer::tick("Force_LCAO_k", "allocate_k");
    return;
}


void Force_LCAO_k::finish_k(LCAO_Matrix &lm)
{
    delete[] lm.DSloc_Rx;
    delete[] lm.DSloc_Ry;
    delete[] lm.DSloc_Rz;
    delete[] lm.DHloc_fixedR_x;
    delete[] lm.DHloc_fixedR_y;
    delete[] lm.DHloc_fixedR_z;

    if (GlobalV::CAL_STRESS)
    {
        delete[] lm.DH_r;
        delete[] lm.stvnl11;
        delete[] lm.stvnl12;
        delete[] lm.stvnl13;
        delete[] lm.stvnl22;
        delete[] lm.stvnl23;
        delete[] lm.stvnl33;
    }
    return;
}


void Force_LCAO_k::test(
		Parallel_Orbitals &pv,
		double* mmm, 
		const std::string& name)
{
    // mohan remove 'const' for pv, 2024-03-31
    if (GlobalV::NPROC != 1)
    {
        return;
    }

    std::cout << "test!" << std::endl;

    int irr = 0;
    int ca = 0;

    GlobalV::ofs_running << " Calculate the test in Force_LCAO_k" << std::endl;
    Record_adj RA;
    
    // mohan update 2024-03-31
    RA.for_2d(pv, GlobalV::GAMMA_ONLY_LOCAL);

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
                        assert(irr < pv.nnr);
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
			{
				std::cout << std::setw(12) << test[i * GlobalV::NLOCAL + j];
			}
			else
			{
				std::cout << std::setw(12) << "0";
			}
		}
        std::cout << std::endl;
    }
    delete[] test;

    RA.delete_grid(); // xiaohui add 2015-02-04
    return;
}



