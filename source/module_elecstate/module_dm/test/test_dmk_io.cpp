#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include "module_elecstate/module_dm/dmk_io.h"

/************************************************
 *  unit test of read_dmk and write_dmk
 ***********************************************/

/**
 * - Tested Functions:
 *   - read_dmk()
 *     - the function to read density matrix from file
 *     - the serial version without MPI
 *   - write_dmk()
 *     - the function to write density matrix to file
 *     - the serial version without MPI
 */

class DMKIOTest : public ::testing::Test
{
protected:
	int nks = 2;
	int nlocal = 26;
	std::vector<ModuleBase::ComplexMatrix> DMK;
	K_Vectors* kv = nullptr;
	void SetUp()
	{
		// initalize a DMK
		DMK.resize(nks);
		ModuleBase::ComplexMatrix zero_dmk(nlocal,nlocal);
		for (int ik = 0;ik < nks;ik++)
		{
			DMK[ik] = zero_dmk;
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

TEST_F(DMKIOTest,IO)
{
	// read
	std::string ssdk;
	for (int ik = 0;ik < kv->nks;++ik){
		ssdk = "./support/" + std::to_string(ik) + ".dmk";
		elecstate::read_dmk(*kv,ik,nlocal,ssdk,DMK);
	}

	// write
	int precision = 3;
	for (int ik = 0;ik < kv->nks;++ik){
		ssdk = "./support/" + std::to_string(ik) + ".dmk1";
		elecstate::write_dmk(*kv,ik,nlocal,ssdk,DMK);
	}

	// read again
	auto dmk = DMK;
	for (int ik = 0;ik < kv->nks;++ik){
		ssdk = "./support/" + std::to_string(ik) + ".dmk";
		elecstate::read_dmk(*kv,ik,nlocal,ssdk,dmk);
	}

	// compare DMK and dmk
	EXPECT_NEAR(DMK[0](0,0).real(),dmk[0](0,0).real(),1e-6);
	EXPECT_NEAR(DMK[1](25,25).real(),dmk[1](25,25).real(),1e-6);
	//for (int ik = 0;ik < kv->nks;++ik){
		//ssdk = "./support/" + std::to_string(ik) + ".dmk1";
		//remove(ssdk);
	//}
}
