#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include "module_elecstate/module_dm/dmk_io.h"
#include "prepare_unitcell.h"
#include "module_cell/parallel_kpoints.h"
#include "module_io/berryphase.h"

#ifdef __LCAO
InfoNonlocal::InfoNonlocal(){}
InfoNonlocal::~InfoNonlocal(){}
LCAO_Orbitals::LCAO_Orbitals(){}
LCAO_Orbitals::~LCAO_Orbitals(){}
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

namespace GlobalC
{
	Parallel_Kpoints Pkpoints;
	UnitCell ucell;
}
bool berryphase::berry_phase_flag=0;


/************************************************
 *  unit test of read_dmk and write_dmk
 ***********************************************/

/**
 * - Tested Functions:
 *   - read_dm()
 *     - the function to read density matrix from file
 *     - the serial version without MPI
 *   - write_dm()
 *     - the function to write density matrix to file
 *     - the serial version without MPI
 */

class DMKIOTest : public ::testing::Test
{
protected:
	int nspin = 1;
	int lgd = 26;
	int nnrg = 26*26;
	int nks = 2;
	std::vector<int> nw = {13};
	int nlocal = 0;
	UnitCell* ucell;
	std::vector<ModuleBase::ComplexMatrix> DMK;
	K_Vectors* kv = nullptr;
	void SetUp()
	{
		// initalize a DMK
		DMK.resize(nks);
		ModuleBase::ComplexMatrix zero_dmk(lgd,lgd);
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
		//delete ucell;
	}
};

TEST_F(DMKIOTest,IO)
{
	GlobalV::NLOCAL = lgd;
	UcellTestPrepare utp = UcellTestLib["Si"];
	ucell = utp.SetUcellInfo(nw, nlocal);
	//std::cout << kv->nks << std::endl;
	//for(int i=0; i<kv->nks; ++i)
    //{
    //    std::cout << "k" << i+1 << " " << kv->kvec_d[i].x << " " << kv->kvec_d[i].y<< " " <<kv->kvec_d[i].z<< std::endl;
    //}

	// read
	std::string ssdk;
	for (int ik = 0;ik < kv->nks;++ik){
        ssdk = "./support/" + std::to_string(ik) + ".dmk";
		//std::cout << ssdk << std::endl;
        elecstate::read_dmk(*kv,ik,ssdk,DMK);
    }

	// write
    int precision = 3;
    for (int ik = 0;ik < kv->nks;++ik){
        ssdk = "./support/" + std::to_string(ik) + ".dmk1";
        elecstate::write_dmk(*kv,ik,ssdk,precision,DMK);
    }

	// read again
	auto dmk = DMK;
	for (int ik = 0;ik < kv->nks;++ik){
        ssdk = "./support/" + std::to_string(ik) + ".dmk";
        elecstate::read_dmk(*kv,ik,ssdk,dmk);
    }

	// compare DMK and dmk
	EXPECT_NEAR(DMK[0](0,0).real(),dmk[0](0,0).real(),1e-6);
	EXPECT_NEAR(DMK[1](25,25).real(),dmk[1](25,25).real(),1e-6);
	//for (int ik = 0;ik < kv->nks;++ik){
        //ssdk = "./support/" + std::to_string(ik) + ".dmk1";
        //remove(ssdk);
    //}
}
