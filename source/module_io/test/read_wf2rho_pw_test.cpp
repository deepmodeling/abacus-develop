#include "gmock/gmock.h"
#include "gtest/gtest.h"

#undef __LCAO

#define private public
#include "module_cell/klist.h"
#include "module_cell/unitcell.h"
#include "module_elecstate/module_charge/charge.h"
#include "module_elecstate/module_charge/symmetry_rho.h"
#include "module_hamilt_general/module_xc/xc_functional.h"
#include "module_hamilt_pw/hamilt_pwdft/parallel_grid.h"
#include "module_io/read_wf2rho_pw.h"
#include "module_io/write_wfc_pw.h"
#include "module_io/filename.h" // mohan add 2025-05-17
#include "module_parameter/parameter.h"
#include "module_psi/psi.h"

#ifdef __MPI
#include "module_base/parallel_global.h"
#include "module_basis/module_pw/test/test_tool.h"
#include "mpi.h"
#endif

Parallel_Grid::Parallel_Grid()
{
}
Parallel_Grid::~Parallel_Grid()
{
}
Charge::Charge()
{
}
Charge::~Charge()
{
}
UnitCell::UnitCell()
{
}
UnitCell::~UnitCell()
{
}
Magnetism::Magnetism()
{
}
Magnetism::~Magnetism()
{
}
int XC_Functional::func_type = 0;
bool XC_Functional::ked_flag = false;

Symmetry_rho::Symmetry_rho()
{
}
Symmetry_rho::~Symmetry_rho()
{
}
void Symmetry_rho::begin(const int& spin_now,
                         const Charge& CHR,
                         const ModulePW::PW_Basis* rho_basis,
                         ModuleSymmetry::Symmetry& symm) const
{
    return;
}

void cal_ik2iktot(std::vector<int>& ik2iktot, const int& nks, const int& nkstot)
{
    if(PARAM.inp.kpar==1)
	{
		for(int ik = 0; ik < nks; ++ik)
		{
			ik2iktot[ik] = ik;
		}
    }
    else if(PARAM.inp.kpar==2)
    {
        if(GlobalV::MY_POOL==0)
		{
			for(int ik = 0; ik < nks; ++ik)
			{
				ik2iktot[ik] = ik;
			}
		}
		else if(GlobalV::MY_POOL==1)
		{
			for(int ik = 0; ik < nks; ++ik)
			{
				ik2iktot[ik] = ik+2; // only works for this test
			}
		}
    }
}

namespace GlobalC
{
	Parallel_Grid Pgrid;
} // namespace GlobalC

/**
 * - Tested Functions:
 *  - write_wfc_pw()
 *  - read_wf2rho_pw()
 */

class ReadWfcRhoTest : public ::testing::Test
{
  protected:
    ModulePW::PW_Basis_K* wfcpw = nullptr;
    ModulePW::PW_Basis* rhopw = nullptr;
    K_Vectors* kv = nullptr;
    psi::Psi<std::complex<double>>* psi = nullptr;
    Charge chg;
    ModuleSymmetry::Symmetry symm;
    virtual void SetUp()
    {
        wfcpw = new ModulePW::PW_Basis_K;
        rhopw = new ModulePW::PW_Basis;
        kv = new K_Vectors;
        // output .dat file
        PARAM.input.out_wfc_pw = 2;
    }
    virtual void TearDown()
    {
        delete wfcpw;
        delete rhopw;
        delete kv;
    }
};

TEST_F(ReadWfcRhoTest, ReadWfcRho)
{
    // kpar=1, if nproc=1
    // kpar=2, if nproc>1
    const int kpar = GlobalV::KPAR;
    const int nspin = 1;
    const int nbands = 4;
    const int my_pool = GlobalV::MY_POOL;
    const int my_rank = GlobalV::MY_RANK;
    const int nks = 2;
    const int nkstot = GlobalV::KPAR * nks;

    //-------------------------
    // Initialize the k-points
    //-------------------------
    kv->set_nkstot(nkstot);
    kv->set_nks(nks);
    kv->isk = {0, 0};
    const double shift = my_pool * 0.1;
    kv->kvec_d = {ModuleBase::Vector3<double>(shift, shift, shift),
                  ModuleBase::Vector3<double>(0.5 + shift, 0.5 + shift, 0.5 + shift)};
    kv->ik2iktot.resize(nks);
    cal_ik2iktot(kv->ik2iktot, nks, nkstot);

    //-------------------------
    // Initialize the pw basis
    //-------------------------
#ifdef __MPI
    wfcpw->initmpi(GlobalV::NPROC_IN_POOL, GlobalV::RANK_IN_POOL, POOL_WORLD);
    rhopw->initmpi(GlobalV::NPROC_IN_POOL, GlobalV::RANK_IN_POOL, POOL_WORLD);
#endif
    rhopw->initgrids(5.3233, ModuleBase::Matrix3(-0.5, 0.0, 0.5, 0.0, 0.5, 0.5, -0.5, 0.5, 0.0), 80);
    rhopw->initparameters(false, 80);
    rhopw->setuptransform();
    rhopw->collect_local_pw();

    wfcpw->initgrids(5.3233, ModuleBase::Matrix3(-0.5, 0.0, 0.5, 0.0, 0.5, 0.5, -0.5, 0.5, 0.0), 80);
    wfcpw->initparameters(false, 20, nks, kv->kvec_d.data());
    wfcpw->setuptransform();
    wfcpw->collect_local_pw();

    //-------------------------
    // Initialize k points
    //-------------------------
    kv->kvec_c.clear();
    for (int ik = 0; ik < nks; ++ik)
    {
        kv->kvec_c.push_back(wfcpw->kvec_c[ik]);
    }
    kv->ngk = {wfcpw->npwk[0], wfcpw->npwk[1]};
    kv->wk = {1.0, 1.0};

    //----------------------------------------
    // Initialize weights for wave functions
    //----------------------------------------
    ModuleBase::matrix wg(nkstot, nbands);
    wg.fill_out(1.0);
    if (GlobalV::MY_RANK == 0)
    {
        std::ofstream ofs("istate.info");
        for (int ik = 0; ik < nkstot; ++ik)
        {
            ofs << "BAND               Energy(ev)               Occupation                Kpoint" << std::endl;
            for (int ib = 0; ib < nbands; ++ib)
            {
                ofs << "  " << ib + 1 << "                  0.0000000                " << 1.0 << std::endl;
            }
            ofs << std::endl;
        }
        ofs.close();
    }

    //----------------------------------------
    // Initialize wave functions Psi
    //----------------------------------------
    psi = new psi::Psi<std::complex<double>>(nks, nbands, wfcpw->npwk_max, kv->ngk, true);
    std::complex<double>* ptr = psi->get_pointer();
    for (int i = 0; i < nks * nbands * wfcpw->npwk_max; i++)
    {
        ptr[i] = std::complex<double>((i + my_pool * 100) / 100.0, (i + my_pool) / 100.0);
    }

    //----------------------------------------
    // Initialize charge density
    //----------------------------------------
    chg.rho = new double*[nspin];
    chg._space_rho = new double[rhopw->nrxx];
    chg.rho[0] = chg._space_rho;
    ModuleBase::GlobalFunc::ZEROS(chg.rho[0], rhopw->nrxx);
    chg.rhopw = rhopw;
    chg.nrxx = rhopw->nrxx;

    //----------------------------------------
    // set charge_ref
    //----------------------------------------
    Charge chg_ref;
    chg_ref.rho = new double*[nspin];
    chg_ref._space_rho = new double[rhopw->nrxx];
    chg_ref.rho[0] = chg_ref._space_rho;
    ModuleBase::GlobalFunc::ZEROS(chg_ref.rho[0], rhopw->nrxx);
    std::vector<std::complex<double>> rho_tmp(rhopw->nrxx);
    chg_ref.nrxx = rhopw->nrxx;

    for (int ik = 0; ik < nks; ++ik)
    {
        for (int ib = 0; ib < nbands; ++ib)
        {
            const std::complex<double>* wfc_ib = ptr + ik * nbands * wfcpw->npwk_max + ib * wfcpw->npwk_max;
            wfcpw->recip2real(wfc_ib, rho_tmp.data(), ik);

            const double w1 = wg(ik, ib) / wfcpw->omega;

            for (int ir = 0; ir < rhopw->nrxx; ir++)
            {
                chg_ref.rho[0][ir] += w1 * std::norm(rho_tmp[ir]);
            }
        }
    }

#ifdef __MPI
    chg_ref.init_chgmpi();
    chg_ref.reduce_diff_pools(chg_ref.rho[0]);
#endif

    // for spin=1 or 2, npol=1
    const int npol=1;

    // Write the wave functions to file
	const std::string out_dir = "./";

    // Read the wave functions to charge density
    std::stringstream ss;
    ss << "running_log" << GlobalV::MY_RANK << ".txt";
	std::ofstream running_log(ss.str().c_str());

    running_log << " rank=" << GlobalV::MY_RANK << std::endl;


    const double ecutwfc = 20; // this is a fake number

	ModuleIO::write_wfc_pw(
			kpar, my_pool, my_rank, nbands, nspin, npol,
			GlobalV::RANK_IN_POOL, GlobalV::NPROC_IN_POOL, 
			PARAM.input.out_wfc_pw, ecutwfc, out_dir, *psi, *kv, wfcpw,
			running_log);

	ModuleIO::read_wf2rho_pw(wfcpw, symm, chg, 
			out_dir, kpar, my_pool, my_rank, 
			GlobalV::NPROC_IN_POOL, GlobalV::RANK_IN_POOL,
			nbands, nspin, npol, 
			nkstot, kv->ik2iktot, kv->isk, running_log);

    // compare the charge density
    for (int ir = 0; ir < rhopw->nrxx; ++ir)
    {
        EXPECT_NEAR(chg.rho[0][ir], chg_ref.rho[0][ir], 1e-8);
    }

	if (GlobalV::NPROC == 1) 
	{
		EXPECT_NEAR(chg.rho[0][0], 8617.076357957576, 1e-8);
	} 
	else if (GlobalV::NPROC == 4)
	{
		const std::vector<double> ref = {8207.849135313403, 35.34776105132742, 8207.849135313403, 35.34776105132742};
		EXPECT_NEAR(chg.rho[0][0], ref[GlobalV::MY_RANK], 1e-8);
		// for (int ip = 0; ip < GlobalV::NPROC; ++ip)
		// {
		//     if (GlobalV::MY_RANK == ip)
		//     {
		//         std::cout.precision(16);
		//         std::cout << GlobalV::MY_RANK << " " << chg.rho[0][0] << std::endl;
		//     }
		//     MPI_Barrier(MPI_COMM_WORLD);
		// }
	}

    delete[] chg.rho;
    delete[] chg._space_rho;
    delete[] chg_ref.rho;
    delete[] chg_ref._space_rho;
    delete psi;

    if (GlobalV::MY_RANK == 0)
    {
        remove("running_log0.txt");
        remove("istate.info");
        remove("wfs1k1_pw.dat");
        remove("wfs1k2_pw.dat");
        if (GlobalV::KPAR == 2)
        {
            remove("wfs1k3_pw.dat");
			remove("wfs1k4_pw.dat");
			remove("running_log1.txt");
			remove("running_log2.txt");
			remove("running_log3.txt");
		}
    }
}

int main(int argc, char** argv)
{
#ifdef __MPI
    setupmpi(argc, argv, GlobalV::NPROC, GlobalV::MY_RANK);

    // when kpar == 2, nspin == 2
    PARAM.input.kpar = (GlobalV::NPROC > 1) ? 2 : 1;
    GlobalV::KPAR = PARAM.input.kpar;
    PARAM.input.bndpar = 1;
    Parallel_Global::divide_pools(GlobalV::NPROC,
                                  GlobalV::MY_RANK,
                                  PARAM.inp.bndpar,
                                  GlobalV::KPAR,
                                  GlobalV::NPROC_IN_BNDGROUP,
                                  GlobalV::RANK_IN_BPGROUP,
                                  GlobalV::MY_BNDGROUP,
                                  GlobalV::NPROC_IN_POOL,
                                  GlobalV::RANK_IN_POOL,
                                  GlobalV::MY_POOL);
#endif

    testing::InitGoogleTest(&argc, argv);
    int result = RUN_ALL_TESTS();

#ifdef __MPI
    finishmpi();
#endif
    return result;
}
