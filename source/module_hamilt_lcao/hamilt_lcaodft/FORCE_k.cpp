#include "FORCE.h"

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
#include "module_hamilt_lcao/hamilt_lcaodft/LCAO_domain.h"

#ifdef __DEEPKS
#include "module_hamilt_lcao/module_deepks/LCAO_deepks.h"
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

template<>
void Force_LCAO<std::complex<double>>::allocate(
    const Parallel_Orbitals& pv,
    LCAO_Matrix& lm,
    ForceStressArrays& fsr, // mohan add 2024-06-15
    const ORB_gen_tables* uot,
    const int& nks,
    const std::vector<ModuleBase::Vector3<double>>& kvec_d)
{
    ModuleBase::TITLE("Force_LCAO", "allocate");
    ModuleBase::timer::tick("Force_LCAO", "allocate");

    const int nnr = pv.nnr;

    assert(nnr>=0);

    //--------------------------------
    // (1) allocate for dSx dSy & dSz
    //--------------------------------
    fsr.DSloc_Rx = new double[nnr];
    fsr.DSloc_Ry = new double[nnr];
    fsr.DSloc_Rz = new double[nnr];

	const auto init_DSloc_Rxyz = [this, nnr, &fsr](int num_threads, int thread_id) 
	{
		int beg=0;
		int len=0;
		ModuleBase::BLOCK_TASK_DIST_1D(num_threads, thread_id, nnr, 1024, beg, len);
		ModuleBase::GlobalFunc::ZEROS(fsr.DSloc_Rx + beg, len);
		ModuleBase::GlobalFunc::ZEROS(fsr.DSloc_Ry + beg, len);
		ModuleBase::GlobalFunc::ZEROS(fsr.DSloc_Rz + beg, len);
	};

    ModuleBase::OMP_PARALLEL(init_DSloc_Rxyz);
    ModuleBase::Memory::record("Force::dS_K", sizeof(double) * nnr * 3);

    if (GlobalV::CAL_STRESS)
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
    LCAO_domain::build_ST_new(
        lm,
        fsr,
        'S', 
        cal_deri, 
        GlobalC::ucell, 
        GlobalC::ORB, 
        pv,
        *uot, 
        &GlobalC::GridD, 
        nullptr); // delete lm.SlocR

    //-----------------------------------------
    // (2) allocate for <phi | T + Vnl | dphi>
    //-----------------------------------------
    fsr.DHloc_fixedR_x = new double[nnr];
    fsr.DHloc_fixedR_y = new double[nnr];
    fsr.DHloc_fixedR_z = new double[nnr];

    const auto init_DHloc_fixedR_xyz = [this, nnr, &fsr](int num_threads, int thread_id) {
        int beg=0;
        int len=0;
        ModuleBase::BLOCK_TASK_DIST_1D(num_threads, thread_id, nnr, 1024, beg, len);
        ModuleBase::GlobalFunc::ZEROS(fsr.DHloc_fixedR_x + beg, len);
        ModuleBase::GlobalFunc::ZEROS(fsr.DHloc_fixedR_y + beg, len);
        ModuleBase::GlobalFunc::ZEROS(fsr.DHloc_fixedR_z + beg, len);
    };
    ModuleBase::OMP_PARALLEL(init_DHloc_fixedR_xyz);
    ModuleBase::Memory::record("Force::dTVNL", sizeof(double) * nnr * 3);

    // calculate dT=<phi|kin|dphi> in LCAO
    // calculate T + VNL(P1) in LCAO basis
    LCAO_domain::build_ST_new(
            lm,
            fsr,
			'T', 
			cal_deri, 
			GlobalC::ucell, 
			GlobalC::ORB, 
			pv,
			*uot, 
			&GlobalC::GridD, 
			lm.Hloc_fixedR.data());

    // calculate dVnl=<phi|dVnl|dphi> in LCAO
	LCAO_domain::build_Nonlocal_mu_new(
			lm,
            fsr,
			lm.Hloc_fixed.data(), 
			cal_deri, 
			GlobalC::ucell, 
			GlobalC::ORB, 
			*uot, 
			&GlobalC::GridD);

	// calculate asynchronous S matrix to output for Hefei-NAMD
	if (INPUT.cal_syns)
	{
		cal_deri = false;
        
		LCAO_domain::build_ST_new(
				lm,
                fsr,
				'S', 
				cal_deri, 
				GlobalC::ucell,
				GlobalC::ORB,
				pv,
				*uot,
                &(GlobalC::GridD),
				nullptr, // delete lm.SlocR
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

template<>
void Force_LCAO<std::complex<double>>::finish_ftable(ForceStressArrays &fsr)
{
    delete[] fsr.DSloc_Rx;
    delete[] fsr.DSloc_Ry;
    delete[] fsr.DSloc_Rz;
    delete[] fsr.DHloc_fixedR_x;
    delete[] fsr.DHloc_fixedR_y;
    delete[] fsr.DHloc_fixedR_z;

    if (GlobalV::CAL_STRESS)
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

template<>
void Force_LCAO<std::complex<double>>::test(
    Parallel_Orbitals& pv,
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

    // be called in Force_LCAO::start_force_calculation
template<>
void Force_LCAO<std::complex<double>>::ftable(
	const bool isforce,
	const bool isstress,
	ForceStressArrays &fsr, // mohan add 2024-06-15
	const UnitCell& ucell,
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
#ifdef __DEEPKS
	ModuleBase::matrix& svnl_dalpha,
#endif
	TGint<std::complex<double>>::type& gint,
	const ORB_gen_tables* uot,
	const Parallel_Orbitals& pv,
	LCAO_Matrix& lm,
	const K_Vectors* kv,
	Record_adj* ra)
{
	ModuleBase::TITLE("Force_LCAO", "ftable");
	ModuleBase::timer::tick("Force_LCAO", "ftable");

	elecstate::DensityMatrix<complex<double>, double>* dm
		= dynamic_cast<const elecstate::ElecStateLCAO<std::complex<double>>*>(pelec)->get_DM();

	this->allocate(
			pv,
			lm,
            fsr, // mohan add 2024-06-16
			uot,
			kv->get_nks(),
			kv->kvec_d);

	// calculate the energy density matrix
	// and the force related to overlap matrix and energy density matrix.
	this->cal_fedm(
			isforce,
			isstress,
            fsr,
			ucell,
			dm,
			psi,
			pv,
			pelec,
			lm,
			foverlap,
			soverlap,
			kv,
			ra);

	this->cal_ftvnl_dphi(
			dm,
			pv,
			ucell,
			fsr,
			isforce,
			isstress,
			ftvnl_dphi,
			stvnl_dphi,
			ra);

	// doing on the real space grid.
	this->cal_fvl_dphi(
			isforce,
			isstress,
			pelec->pot,
			gint,
			fvl_dphi,
			svl_dphi);

	this->cal_fvnl_dbeta(
			dm,
			pv,
			ucell,
			GlobalC::ORB,
			*uot,
			GlobalC::GridD,
			isforce,
			isstress,
			fvnl_dbeta,
			svnl_dbeta);

#ifdef __DEEPKS
	if (GlobalV::deepks_scf)
	{
		const std::vector<std::vector<std::complex<double>>>& dm_k = dm->get_DMK_vector();

		GlobalC::ld.cal_projected_DM_k(dm, ucell, GlobalC::ORB, GlobalC::GridD);

		GlobalC::ld.cal_descriptor(ucell.nat);

		GlobalC::ld.cal_gedm(ucell.nat);

		GlobalC::ld.cal_f_delta_k(dm_k,
				ucell,
				GlobalC::ORB,
				GlobalC::GridD,
				kv->get_nks(),
				kv->kvec_d,
				isstress,
				svnl_dalpha);
#ifdef __MPI
		Parallel_Reduce::reduce_all(GlobalC::ld.F_delta.c, GlobalC::ld.F_delta.nr * GlobalC::ld.F_delta.nc);
		if (isstress)
		{
			Parallel_Reduce::reduce_pool(svnl_dalpha.c, svnl_dalpha.nr * svnl_dalpha.nc);
		}
#endif
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

    ModuleBase::timer::tick("Force_LCAO", "ftable");
    return;
}
