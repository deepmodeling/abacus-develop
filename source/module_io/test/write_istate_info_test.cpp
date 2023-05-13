#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include <streambuf>
#include "module_base/global_variable.h"
#ifdef __MPI
#include "mpi.h"
#include "module_base/parallel_global.h"
#include "module_cell/parallel_kpoints.h"
#endif
#include "../write_istate_info.h"
#include "for_testing_klist.h"

/************************************************
 *  unit test of write_istate_info
 ***********************************************/

/**
 * - Tested Functions:
 *   - write_istate_info()
 *     - print out electronic eigen energies and
 *     - occupation
 */

class IstateInfoTest : public ::testing::Test
{
protected:
	K_Vectors* kv = nullptr;
	Parallel_Kpoints* Pkpoints = nullptr;
	ModuleBase::matrix ekb;
	ModuleBase::matrix wg;
	void SetUp()
	{
		kv = new K_Vectors;
		Pkpoints = new Parallel_Kpoints;
	}
	void TearDown()
	{
		delete kv;
		delete Pkpoints;
	}
};

TEST_F(IstateInfoTest,OutIstateInfoS1)
{
	//preconditions
	GlobalV::KPAR = 1;
	GlobalV::NBANDS = 4;
	GlobalV::NSPIN = 1;
	GlobalV::global_out_dir = "./";
	//mpi setting
	Parallel_Global::init_pools();
	kv->nkstot = 100;
	Pkpoints->kinfo(kv->nkstot);
	//std::cout<<"my_rank "<<GlobalV::MY_RANK<<" pool rank/size: "
	//	<<GlobalV::RANK_IN_POOL<<"/"<<GlobalV::NPROC_IN_POOL<<std::endl;
	//std::cout<<"MY_POOL "<<GlobalV::MY_POOL<<std::endl;
	kv->nks = Pkpoints->nks_pool[GlobalV::MY_POOL];
	//std::cout<<"nks "<<kv->nks<<std::endl;
	ekb.create(kv->nks,GlobalV::NBANDS);
	wg.create(kv->nks,GlobalV::NBANDS);
	ekb.fill_out(0.15);
	wg.fill_out(0.0);
	kv->kvec_d.resize(kv->nkstot);
	int i=0;
	for(auto& kd : kv->kvec_d){
		kd.set(0.01*i,0.01*i,0.01*i);
		++i;
	}
	ModuleIO::write_istate_info(ekb,wg,*kv,Pkpoints);
	std::ifstream ifs;
	ifs.open("istate.info");
	std::string str((std::istreambuf_iterator<char>(ifs)),std::istreambuf_iterator<char>());
	EXPECT_THAT(str, testing::HasSubstr("BAND               Energy(ev)               Occupation                Kpoint = 100"));
	EXPECT_THAT(str, testing::HasSubstr("(0.99 0.99 0.99)"));
	EXPECT_THAT(str, testing::HasSubstr("4                  2.04085                        0"));
	ifs.close();
	remove("istate.info");
}

TEST_F(IstateInfoTest,OutIstateInfoS2)
{
	//preconditions
	GlobalV::KPAR = 1;
	GlobalV::NBANDS = 4;
	GlobalV::NSPIN = 2;
	GlobalV::global_out_dir = "./";
	//mpi setting
	Parallel_Global::init_pools();
	kv->nkstot = 100;
	Pkpoints->kinfo(kv->nkstot);
	//std::cout<<"my_rank "<<GlobalV::MY_RANK<<" pool rank/size: "
	//	<<GlobalV::RANK_IN_POOL<<"/"<<GlobalV::NPROC_IN_POOL<<std::endl;
	//std::cout<<"MY_POOL "<<GlobalV::MY_POOL<<std::endl;
	kv->nks = Pkpoints->nks_pool[GlobalV::MY_POOL];
	//std::cout<<"nks "<<kv->nks<<std::endl;
	ekb.create(kv->nks,GlobalV::NBANDS);
	wg.create(kv->nks,GlobalV::NBANDS);
	ekb.fill_out(0.15);
	wg.fill_out(0.0);
	kv->kvec_d.resize(kv->nkstot);
	int i=0;
	for(auto& kd : kv->kvec_d){
		kd.set(0.01*i,0.01*i,0.01*i);
		++i;
	}
	ModuleIO::write_istate_info(ekb,wg,*kv,Pkpoints);
	std::ifstream ifs;
	ifs.open("istate.info");
	std::string str((std::istreambuf_iterator<char>(ifs)),std::istreambuf_iterator<char>());
	EXPECT_THAT(str, testing::HasSubstr("BAND       Spin up Energy(ev)               Occupation     Spin down Energy(ev)               Occupation                Kpoint = 50"));
	EXPECT_THAT(str, testing::HasSubstr("4                  2.04085                        0                  2.04085                        0"));
	EXPECT_THAT(str, testing::HasSubstr("(0.49 0.49 0.49)"));
	ifs.close();
	remove("istate.info");
}

#ifdef __MPI
int main(int argc, char **argv)
{
	MPI_Init(&argc,&argv);

	testing::InitGoogleTest(&argc,argv);
	MPI_Comm_size(MPI_COMM_WORLD,&GlobalV::NPROC);
	MPI_Comm_rank(MPI_COMM_WORLD,&GlobalV::MY_RANK);
	int result = RUN_ALL_TESTS();

	MPI_Finalize();

	return result;
}
#endif
