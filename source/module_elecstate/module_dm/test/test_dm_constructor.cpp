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
 *  unit test of DensityMatrix constructor
 ***********************************************/

/**
 * This unit test construct a DensityMatrix object
 */

class DMTest : public testing::Test
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

TEST_F(DMTest, DMConstructor1)
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
    // compare
    EXPECT_NEAR(DM.get_dmK(0,0,0),0,1e-6);
    EXPECT_NEAR(DM.get_dmK(1,25,25),0,1e-6);
}
