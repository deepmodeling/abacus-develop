#include "gtest/gtest.h"
#include "module_base/global_function.h"
#include "module_base/constants.h"
#include "module_base/matrix3.h"

/************************************************
 *  serial unit test of functions in pw_basis.cpp
 ***********************************************/

/**
 * - Tested Functions:
 *   - Constructor
 *     - PW_Basis() and ~PW_Basis()
 *   - Initgrids1
 *     - initgrids() from gridecut
 *   - Initgrids2
 *     - initgrids() from nx,ny,nz
 *   - Initparameters
 *     - initparameters() for fft
 *   - Setfullpw
 *     - setfullpw(): set controlling parameters to get full pw in orbital free calculations
 *   - DistributeR
 *     - distribute_r(): distribute real space grids in z direction
 *   - DistributeMethod1
 *     - distribute_g() and distribution_method1(): set and distribute sticks
 *   - DistributeMethod2
 *     - distribute_g() and distribution_method2(): set and distribute sticks
 *   - GetStartGR
 *     - getstartgr(): get nmaxgr, numg, numr, startg, startr
 *   - SetupTransform
 *     - setuptransform(): for fft transform
 *   - CollectLocalPW
 *     - collect_local_pw: get gg, gdirect, gcar for local npw plane waves
 *   - CollectUniqgg
 *     - collect_uniqgg: get uniq gg without duplication in length
 */

#define protected public
#define private public
#include "../pw_basis.h"
#include "../fft.h"
class PWBasisTEST: public testing::Test
{
public:
	std::string precision_flag = "double";
	std::string device_flag = "cpu";
	ModulePW::PW_Basis pwb;
	ModulePW::PW_Basis pwb1;
};

TEST_F(PWBasisTEST,Constructor)
{
	ModulePW::PW_Basis pwb2(device_flag, precision_flag);
	EXPECT_EQ(pwb1.classname,"PW_Basis");
	EXPECT_EQ(pwb2.classname,"PW_Basis");
	EXPECT_EQ(pwb2.device,"cpu");
	EXPECT_EQ(pwb2.precision,"double");
	EXPECT_EQ(pwb2.ft.device,"cpu");
	EXPECT_EQ(pwb2.ft.precision,"double");
}

TEST_F(PWBasisTEST,Initgrids1)
{
	double lat0 = 1.8897261254578281;
	ModuleBase::Matrix3 latvec(10.0,0.0,0.0,
				0.0,10.0,0.0,
				0.0,0.0,10.0);
	double gridecut=10.0;
	pwb.initgrids(lat0,latvec,gridecut);
	EXPECT_DOUBLE_EQ(pwb.lat0,lat0);
	EXPECT_DOUBLE_EQ(pwb.tpiba,ModuleBase::TWO_PI/lat0);
	EXPECT_DOUBLE_EQ(pwb.tpiba2,pwb.tpiba*pwb.tpiba);
	EXPECT_DOUBLE_EQ(pwb.latvec.e11,latvec.e11);
	EXPECT_DOUBLE_EQ(pwb.GT.e11,latvec.Inverse().e11);
	EXPECT_DOUBLE_EQ(pwb.G.e11,pwb.GT.Transpose().e11);
	EXPECT_DOUBLE_EQ(pwb.GGT.e11,(pwb.G*pwb.GT).e11);
	EXPECT_DOUBLE_EQ(pwb.gridecut_lat,gridecut/pwb.tpiba2);
	EXPECT_NEAR(pwb.gridecut_lat,0.904561,1e-4);
	EXPECT_EQ(pwb.nx,20);
	EXPECT_EQ(pwb.ny,20);
	EXPECT_EQ(pwb.nz,20);
	EXPECT_TRUE(pwb.nx%2==0 || pwb.nx%3==0 || pwb.nx%5==0);
	EXPECT_TRUE(pwb.ny%2==0 || pwb.ny%3==0 || pwb.ny%5==0);
	EXPECT_TRUE(pwb.nz%2==0 || pwb.nz%3==0 || pwb.nz%5==0);
}

TEST_F(PWBasisTEST,Initgrids2)
{
	double lat0 = 1.8897261254578281;
	ModuleBase::Matrix3 latvec(10.0,0.0,0.0,
				0.0,10.0,0.0,
				0.0,0.0,10.0);
	int nx_in = 20;
	int ny_in = 20;
	int nz_in = 20;
	pwb.initgrids(lat0,latvec,nx_in,ny_in,nz_in);
	EXPECT_DOUBLE_EQ(pwb.lat0,lat0);
	EXPECT_DOUBLE_EQ(pwb.tpiba,ModuleBase::TWO_PI/lat0);
	EXPECT_DOUBLE_EQ(pwb.tpiba2,pwb.tpiba*pwb.tpiba);
	EXPECT_DOUBLE_EQ(pwb.latvec.e11,latvec.e11);
	EXPECT_DOUBLE_EQ(pwb.GT.e11,latvec.Inverse().e11);
	EXPECT_DOUBLE_EQ(pwb.G.e11,pwb.GT.Transpose().e11);
	EXPECT_DOUBLE_EQ(pwb.GGT.e11,(pwb.G*pwb.GT).e11);
	EXPECT_EQ(pwb.nx,nx_in);
	EXPECT_EQ(pwb.ny,ny_in);
	EXPECT_EQ(pwb.nz,nz_in);
	EXPECT_NEAR(pwb.gridecut_lat,0.999999,1e-4);
	EXPECT_NEAR(pwb.gridecut_lat*pwb.tpiba2,11.0551,1e-4);
}

TEST_F(PWBasisTEST,Initparameters)
{
	double lat0 = 1.8897261254578281;
	ModuleBase::Matrix3 latvec(10.0,0.0,0.0,
				0.0,10.0,0.0,
				0.0,0.0,10.0);
	double gridecut=10.0;
	//initparameters is always called after initgrids
	//because of nx,ny,nz, and tpiba2
	pwb.initgrids(lat0,latvec,gridecut);
	bool gamma_only_in = true;
	double pwecut_in = 11.0;
	int distribution_type_in = 1;
	bool xprime_in = true;
	pwb.initparameters(gamma_only_in,pwecut_in,distribution_type_in,xprime_in);
	EXPECT_EQ(pwb.xprime,xprime_in);
	EXPECT_EQ(pwb.gamma_only,gamma_only_in);
	EXPECT_EQ(pwb.xprime,xprime_in);
	EXPECT_TRUE(pwb.gamma_only);
	EXPECT_TRUE(pwb.xprime);
	EXPECT_EQ(pwb.fftnx,int(pwb.nx/2)+1);
	EXPECT_EQ(pwb.fftny,pwb.ny);
	EXPECT_EQ(pwb.fftnz,pwb.nz);
	EXPECT_EQ(pwb.ggecut,pwb.gridecut_lat);
	EXPECT_EQ(pwb.distribution_type,distribution_type_in);
}

TEST_F(PWBasisTEST,Setfullpw)
{
	bool inpt_full_pw = false;
	int inpt_full_pw_dim = 2;
	pwb.setfullpw(inpt_full_pw,inpt_full_pw_dim);
	EXPECT_FALSE(pwb.full_pw);
	EXPECT_EQ(pwb.full_pw_dim,0);
}

TEST_F(PWBasisTEST,DistributeR)
{
	//distribute_r depends on initgrids
	//because of nz
	double lat0 = 1.8897261254578281;
	ModuleBase::Matrix3 latvec(10.0,0.0,0.0,
				0.0,10.0,0.0,
				0.0,0.0,10.0);
	double gridecut=10.0;
	//initparameters is always called after initgrids
	//because of nx,ny,nz, and tpiba2
	pwb.initgrids(lat0,latvec,gridecut);
	//this is serial test, so that
	EXPECT_EQ(pwb.poolrank,0);
	EXPECT_EQ(pwb.poolnproc,1);
	pwb.distribute_r();
	EXPECT_EQ(pwb.startz[0],0);
	EXPECT_EQ(pwb.numz[0],pwb.nz);
	EXPECT_EQ(pwb.nplane,pwb.nz);
	EXPECT_EQ(pwb.nplane,20);
	EXPECT_EQ(pwb.nxy,400);
	EXPECT_EQ(pwb.nrxx,pwb.numz[0]*pwb.nxy);
}

TEST_F(PWBasisTEST,DistributeMethod1)
{
	double lat0 = 1.8897261254578281;
	ModuleBase::Matrix3 latvec(10.0,0.0,0.0,
				0.0,10.0,0.0,
				0.0,0.0,10.0);
	double gridecut=10.0;
	//initparameters is always called after initgrids
	//because of nx,ny,nz, and tpiba2
	//call initgrids
	pwb.initgrids(lat0,latvec,gridecut);
	bool gamma_only_in = true;
	double pwecut_in = 11.0;
	int distribution_type_in = 1;
	bool xprime_in = true;
	//call initparameters
	pwb.initparameters(gamma_only_in,pwecut_in,distribution_type_in,xprime_in);
	EXPECT_TRUE(pwb.gamma_only);
	EXPECT_TRUE(pwb.xprime);
	EXPECT_EQ(pwb.fftnx,int(pwb.nx/2)+1);
	EXPECT_EQ(pwb.fftny,pwb.ny);
	EXPECT_EQ(pwb.fftnz,pwb.nz);
	EXPECT_EQ(pwb.distribution_type,distribution_type_in);
	EXPECT_EQ(pwb.fftnxy,pwb.fftnx*pwb.fftny);
	EXPECT_EQ(pwb.fftnx,11);
	EXPECT_EQ(pwb.fftny,20);
	EXPECT_EQ(pwb.fftnz,20);
	//distribute_method1 depends on initparamters
	//because of fftnxy 
	EXPECT_EQ(pwb.fftnxy,220);
	EXPECT_EQ(pwb.distribution_type,1);
	//call distribute_g
	pwb.distribute_g();
	EXPECT_EQ(pwb.npwtot,1994);
	EXPECT_EQ(pwb.nstot,156);
}

TEST_F(PWBasisTEST,DistributeMethod2)
{
	double lat0 = 1.8897261254578281;
	ModuleBase::Matrix3 latvec(10.0,0.0,0.0,
				0.0,10.0,0.0,
				0.0,0.0,10.0);
	double gridecut=10.0;
	//initparameters is always called after initgrids
	//because of nx,ny,nz, and tpiba2
	//call initgrids
	pwb.initgrids(lat0,latvec,gridecut);
	bool gamma_only_in = true;
	double pwecut_in = 11.0;
	int distribution_type_in = 2;
	bool xprime_in = true;
	//call initparameters
	pwb.initparameters(gamma_only_in,pwecut_in,distribution_type_in,xprime_in);
	EXPECT_TRUE(pwb.gamma_only);
	EXPECT_TRUE(pwb.xprime);
	EXPECT_EQ(pwb.fftnx,int(pwb.nx/2)+1);
	EXPECT_EQ(pwb.fftny,pwb.ny);
	EXPECT_EQ(pwb.fftnz,pwb.nz);
	EXPECT_EQ(pwb.distribution_type,distribution_type_in);
	EXPECT_EQ(pwb.fftnxy,pwb.fftnx*pwb.fftny);
	EXPECT_EQ(pwb.fftnx,11);
	EXPECT_EQ(pwb.fftny,20);
	EXPECT_EQ(pwb.fftnz,20);
	//distribute_method1 depends on initparamters
	//because of fftnxy 
	EXPECT_EQ(pwb.fftnxy,220);
	EXPECT_EQ(pwb.distribution_type,2);
	//call distribute_g
	pwb.distribute_g();
	EXPECT_EQ(pwb.npwtot,1994);
	EXPECT_EQ(pwb.nstot,156);
	EXPECT_EQ(pwb.npw,1994);
	EXPECT_EQ(pwb.nst,156);
	EXPECT_EQ(pwb.nstnz,3120);
}

TEST_F(PWBasisTEST,GetStartGR)
{
	//getstartgr is called after distribute_r and distribute_g in setuptransform
	double lat0 = 1.8897261254578281;
	ModuleBase::Matrix3 latvec(10.0,0.0,0.0,
				0.0,10.0,0.0,
				0.0,0.0,10.0);
	double gridecut=10.0;
	//initparameters is always called after initgrids
	//because of nx,ny,nz, and tpiba2
	//call initgrids
	pwb.initgrids(lat0,latvec,gridecut);
	//call distribute_r
	pwb.distribute_r();
	bool gamma_only_in = true;
	double pwecut_in = 11.0;
	int distribution_type_in = 2;
	bool xprime_in = true;
	//call initparameters
	pwb.initparameters(gamma_only_in,pwecut_in,distribution_type_in,xprime_in);
	//call distribute_g
	pwb.distribute_g();
	//call getstartgr
	pwb.getstartgr();
	EXPECT_TRUE(pwb.gamma_only);
	EXPECT_EQ(pwb.npw,1994);
	EXPECT_EQ(pwb.nz,20);
	EXPECT_EQ(pwb.nst,156);
	EXPECT_EQ(pwb.nrxx,8000);
	EXPECT_EQ(pwb.nxy,400);
	EXPECT_EQ(pwb.nplane,20);
	EXPECT_EQ(pwb.nmaxgr,8000);
	EXPECT_EQ(pwb.numg[0],3120);
	EXPECT_EQ(pwb.numr[0],3120);
	EXPECT_EQ(pwb.startg[0],0);
	EXPECT_EQ(pwb.startr[0],0);
}

TEST_F(PWBasisTEST,SetupTransform)
{
	//getstartgr is called after distribute_r and distribute_g in setuptransform
	double lat0 = 1.8897261254578281;
	ModuleBase::Matrix3 latvec(10.0,0.0,0.0,
				0.0,10.0,0.0,
				0.0,0.0,10.0);
	double gridecut=10.0;
	//initparameters is always called after initgrids
	//because of nx,ny,nz, and tpiba2
	//call initgrids
	pwb.initgrids(lat0,latvec,gridecut);
	bool gamma_only_in = true;
	double pwecut_in = 11.0;
	int distribution_type_in = 2;
	bool xprime_in = true;
	pwb.initparameters(gamma_only_in,pwecut_in,distribution_type_in,xprime_in);
	//setuptransform for FFT
	//which calls fft planning functions
	//currently this is just a trivial test to see its successfull calling
	EXPECT_NO_THROW(pwb.setuptransform());
	EXPECT_EQ(pwb.npw,1994);
}

TEST_F(PWBasisTEST,CollectLocalPW)
{
	//getstartgr is called after distribute_r and distribute_g in setuptransform
	double lat0 = 1.8897261254578281;
	ModuleBase::Matrix3 latvec(10.0,0.0,0.0,
				0.0,10.0,0.0,
				0.0,0.0,10.0);
	double gridecut=10.0;
	//initparameters is always called after initgrids
	//because of nx,ny,nz, and tpiba2
	//call initgrids
	pwb.initgrids(lat0,latvec,gridecut);
	bool gamma_only_in = true;
	double pwecut_in = 11.0;
	int distribution_type_in = 2;
	bool xprime_in = true;
	pwb.initparameters(gamma_only_in,pwecut_in,distribution_type_in,xprime_in);
	//setuptransform for FFT
	//which calls fft planning functions
	//currently this is just a trivial test to see its successfull calling
	EXPECT_NO_THROW(pwb.setuptransform());
	EXPECT_EQ(pwb.npw,1994);
	pwb.collect_local_pw();
	EXPECT_EQ(pwb.ig_gge0,9);
}

TEST_F(PWBasisTEST,CollectUniqgg)
{
	//getstartgr is called after distribute_r and distribute_g in setuptransform
	double lat0 = 1.8897261254578281;
	ModuleBase::Matrix3 latvec(10.0,0.0,0.0,
				0.0,10.0,0.0,
				0.0,0.0,10.0);
	double gridecut=10.0;
	//initparameters is always called after initgrids
	//because of nx,ny,nz, and tpiba2
	//call initgrids
	pwb.initgrids(lat0,latvec,gridecut);
	bool gamma_only_in = true;
	double pwecut_in = 11.0;
	int distribution_type_in = 2;
	bool xprime_in = true;
	pwb.initparameters(gamma_only_in,pwecut_in,distribution_type_in,xprime_in);
	//setuptransform for FFT
	//which calls fft planning functions
	//currently this is just a trivial test to see its successfull calling
	EXPECT_NO_THROW(pwb.setuptransform());
	EXPECT_EQ(pwb.npw,1994);
	pwb.collect_local_pw();
	EXPECT_EQ(pwb.ig_gge0,9);
	pwb.collect_uniqgg();
	EXPECT_EQ(pwb.ngg,78);
}
#undef private
#undef protected

// It is weird to init MPI again in this serial test, but it is necessary
// to avoid MPICH not initialized error in intel compiler.
// The error happens somewhere in the DistributeMethod1 test. 
#include "mpi.h"
int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    testing::InitGoogleTest(&argc, argv);

    int nproc;
    int my_rank;

    MPI_Comm_size(MPI_COMM_WORLD,&nproc);
    MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

    int result = RUN_ALL_TESTS();

    MPI_Finalize();

    return result;
}