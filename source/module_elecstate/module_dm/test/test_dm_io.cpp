#include <fstream>

#include "module_hamilt_lcao/module_hcontainer/hcontainer.h"
#include "module_hamilt_lcao/module_hcontainer/output_hcontainer.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "module_io/csr_reader.h"
#include "prepare_unitcell.h"
#include "module_elecstate/module_dm/density_matrix.h"
#include "module_base/libm/libm.h"
#include <chrono>

// mock functions
#ifdef __LCAO
InfoNonlocal::InfoNonlocal()
{
}
InfoNonlocal::~InfoNonlocal()
{
}
LCAO_Orbitals::LCAO_Orbitals()
{
}
LCAO_Orbitals::~LCAO_Orbitals()
{
}
#endif
Magnetism::Magnetism()
{
    this->tot_magnetization = 0.0;
    this->abs_magnetization = 0.0;
    this->start_magnetization = nullptr;
}
Magnetism::~Magnetism()
{
    delete[] this->start_magnetization;
}
// mocke functions

/************************************************
 *  unit test of read and output DMK 
 ***********************************************/

/**
 * This unit test read sparse matrices of DMK
 * from nks files *.dmk, 
 * and output the matrices of DMR to a file output*.dmk.
 */

class DMIOTest : public testing::Test
{
  protected:
    UnitCell* ucell;
    UcellTestPrepare utp = UcellTestLib["Si"];
	std::vector<ModuleBase::ComplexMatrix> DMK;
	K_Vectors* kv = nullptr;
    // nw is the number of orbitals of each atom
    // it should container ucell.nat elements
    std::vector<int> nw = {13};
    int nks = 2;
	int nlocal = 0;
    void SetUp() override
    {
        // initalize a unitcell
        ucell = utp.SetUcellInfo(nw, nlocal);
        // initalize a DMK
		DMK.reserve(nks);
		ModuleBase::ComplexMatrix zero_dmk(nlocal,nlocal);
		for (int ik = 0;ik < nks;ik++)
		{
			DMK.push_back(zero_dmk);
		}
		// initalize a kvectors
		kv = new K_Vectors;
		kv->nks = nks;
		kv->kvec_d.resize(nks);
		kv->kvec_d[1].x = 0.5;
    }
    void TearDown()
	{
		DMK.clear();
		delete kv;
	}
};

TEST_F(DMIOTest, DMKIO)
{
    // construct paraV
    Parallel_Orbitals paraV;
    std::ofstream ofs("test.log");
    paraV.set_global2local(nlocal, nlocal, false, ofs);
    ofs.close();
    remove("test.log");
    paraV.set_atomic_trace(ucell->get_iat2iwt(), ucell->nat, nlocal);
    // construct DM
    elecstate::DensityMatrix<double> DM(nlocal,kv,&paraV);
    // read DMK
    DM.read_all_dmk("./support/");
    // write DMK
    DM.write_all_dmk("./support/output");
    // construct a new DM
    elecstate::DensityMatrix<double> DM1(nlocal,kv,&paraV);
    DM1.read_all_dmk("./support/output");
    // compare DMK1 with DMK
    EXPECT_NEAR(DM.get_dmK(0,0,0),DM1.get_dmK(0,0,0),1e-6);
	EXPECT_NEAR(DM.get_dmK(1,25,25),DM1.get_dmK(1,25,25),1e-6);
}
