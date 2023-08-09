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
 *     - PW_Basis_K() and ~PW_Basis_K()
 *   - Initgrids1
 *     - initgrids() from gridecut, derived from class PW_Basis
 *   - Initgrids2
 *     - initgrids() from nx,ny,nz, derived from class PW_Basis
 *   - Initparameters
 *     - initparameters(), including k coordinates
 *   - SetupTransform
 *     - setuptransform(): for fft transform
 *   - CollectLocalPW
 *     - collect_local_pw: get gk2, gcar for local npw plane waves
 */

#define protected public
#define private public
#include "../pw_basis_k.h"
#include "../pw_basis.h"
#include "../fft.h"
class PWBasisKTEST: public testing::Test
{
public:
	std::string precision_double = "double";
	std::string precision_single = "single";
	std::string device_flag = "cpu";
};

TEST_F(PWBasisKTEST,Constructor)
{
	ModulePW::PW_Basis_K basis_k1;
	ModulePW::PW_Basis_K basis_k2(device_flag, precision_double);
	EXPECT_EQ(basis_k1.classname,"PW_Basis_K");
	EXPECT_EQ(basis_k2.classname,"PW_Basis_K");
	EXPECT_EQ(basis_k2.device,"cpu");
	EXPECT_EQ(basis_k2.ft.device,"cpu");
	EXPECT_EQ(basis_k2.precision,"double");
	EXPECT_EQ(basis_k2.ft.precision,"double");
	ModulePW::PW_Basis_K basis_k3(device_flag, precision_single);
	EXPECT_EQ(basis_k3.ft.precision,"single");
}

TEST_F(PWBasisKTEST,Initgrids1)
{
	ModulePW::PW_Basis_K basis_k;
	double lat0 = 1.8897261254578281;
	ModuleBase::Matrix3 latvec(10.0,0.0,0.0,
				0.0,10.0,0.0,
				0.0,0.0,10.0);
	double gridecut=10.0;
	basis_k.initgrids(lat0,latvec,gridecut);
	EXPECT_DOUBLE_EQ(basis_k.lat0,lat0);
	EXPECT_DOUBLE_EQ(basis_k.tpiba,ModuleBase::TWO_PI/lat0);
	EXPECT_DOUBLE_EQ(basis_k.tpiba2,basis_k.tpiba*basis_k.tpiba);
	EXPECT_DOUBLE_EQ(basis_k.latvec.e11,latvec.e11);
	EXPECT_DOUBLE_EQ(basis_k.GT.e11,latvec.Inverse().e11);
	EXPECT_DOUBLE_EQ(basis_k.G.e11,basis_k.GT.Transpose().e11);
	EXPECT_DOUBLE_EQ(basis_k.GGT.e11,(basis_k.G*basis_k.GT).e11);
	EXPECT_DOUBLE_EQ(basis_k.gridecut_lat,gridecut/basis_k.tpiba2);
	EXPECT_NEAR(basis_k.gridecut_lat,0.904561,1e-4);
	EXPECT_EQ(basis_k.nx,20);
	EXPECT_EQ(basis_k.ny,20);
	EXPECT_EQ(basis_k.nz,20);
	EXPECT_TRUE(basis_k.nx%2==0 || basis_k.nx%3==0 || basis_k.nx%5==0);
	EXPECT_TRUE(basis_k.ny%2==0 || basis_k.ny%3==0 || basis_k.ny%5==0);
	EXPECT_TRUE(basis_k.nz%2==0 || basis_k.nz%3==0 || basis_k.nz%5==0);
}

TEST_F(PWBasisKTEST,Initgrids2)
{
	ModulePW::PW_Basis_K basis_k;
	double lat0 = 1.8897261254578281;
	ModuleBase::Matrix3 latvec(10.0,0.0,0.0,
				0.0,10.0,0.0,
				0.0,0.0,10.0);
	int nx_in = 20;
	int ny_in = 20;
	int nz_in = 20;
	basis_k.initgrids(lat0,latvec,nx_in,ny_in,nz_in);
	EXPECT_DOUBLE_EQ(basis_k.lat0,lat0);
	EXPECT_DOUBLE_EQ(basis_k.tpiba,ModuleBase::TWO_PI/lat0);
	EXPECT_DOUBLE_EQ(basis_k.tpiba2,basis_k.tpiba*basis_k.tpiba);
	EXPECT_DOUBLE_EQ(basis_k.latvec.e11,latvec.e11);
	EXPECT_DOUBLE_EQ(basis_k.GT.e11,latvec.Inverse().e11);
	EXPECT_DOUBLE_EQ(basis_k.G.e11,basis_k.GT.Transpose().e11);
	EXPECT_DOUBLE_EQ(basis_k.GGT.e11,(basis_k.G*basis_k.GT).e11);
	EXPECT_EQ(basis_k.nx,nx_in);
	EXPECT_EQ(basis_k.ny,ny_in);
	EXPECT_EQ(basis_k.nz,nz_in);
	EXPECT_NEAR(basis_k.gridecut_lat,0.999999,1e-4);
	EXPECT_NEAR(basis_k.gridecut_lat*basis_k.tpiba2,11.0551,1e-4);
}

TEST_F(PWBasisKTEST, Initparameters) 
{
	ModulePW::PW_Basis_K basis_k(device_flag, precision_single);
	double lat0 = 1.8897261254578281;
	ModuleBase::Matrix3 latvec(10.0,0.0,0.0,
				0.0,10.0,0.0,
				0.0,0.0,10.0);
	int nx_in = 20;
	int ny_in = 20;
	int nz_in = 20;
	basis_k.initgrids(lat0,latvec,nx_in,ny_in,nz_in);
	const bool gamma_only_in = true;
	const double gk_ecut_in = 2.0;
	const int nks_in = 3;
	const ModuleBase::Vector3<double> kvec_d_in[3] = { {0.0, 0.0, 0.0}, {0.1, 0.2, 0.3}, {0.4, 0.5, 0.6} };
	const int distribution_type_in = 1;
	const bool xprime_in = true;	
	basis_k.initparameters(gamma_only_in, gk_ecut_in, nks_in,kvec_d_in, distribution_type_in, xprime_in);	
	EXPECT_EQ(basis_k.nks, nks_in);	
	EXPECT_NE(basis_k.kvec_d, nullptr);
	for(int i=0; i<nks_in; i++) {
	    EXPECT_EQ(basis_k.kvec_d[i], kvec_d_in[i]);
	}	
	EXPECT_NE(basis_k.kvec_c, nullptr);
	for(int i=0; i<nks_in; i++) {
	    EXPECT_EQ(basis_k.kvec_c[i], kvec_d_in[i] * basis_k.G);
	}	
	EXPECT_GT(basis_k.gk_ecut, 0.0);
	EXPECT_GT(basis_k.ggecut, 0.0);
	EXPECT_LE(basis_k.ggecut, basis_k.gridecut_lat);	
	EXPECT_FALSE(basis_k.gamma_only);
	EXPECT_EQ(basis_k.xprime, xprime_in);	
	if(basis_k.gamma_only) {
	    EXPECT_EQ(basis_k.fftny, basis_k.ny);
	    EXPECT_EQ(basis_k.fftnx, int(basis_k.nx / 2) + 1);
	} else {
	    EXPECT_EQ(basis_k.fftny, basis_k.ny);
	    EXPECT_EQ(basis_k.fftnx, basis_k.nx);
	}	
	EXPECT_EQ(basis_k.fftnz, basis_k.nz);
	EXPECT_EQ(basis_k.fftnxy, basis_k.fftnx * basis_k.fftny);
	EXPECT_EQ(basis_k.fftnxyz, basis_k.fftnxy * basis_k.fftnz);	
	EXPECT_EQ(basis_k.distribution_type, distribution_type_in);
}

TEST_F(PWBasisKTEST, SetupTransform) 
{
	ModulePW::PW_Basis_K basis_k(device_flag, precision_single);
	double lat0 = 1.8897261254578281;
	ModuleBase::Matrix3 latvec(10.0,0.0,0.0,
				0.0,10.0,0.0,
				0.0,0.0,10.0);
	double gridecut=10.0;
	basis_k.initgrids(lat0,latvec,gridecut);
	const bool gamma_only_in = true;
	const double gk_ecut_in = 10.0;
	const int nks_in = 3;
	const ModuleBase::Vector3<double> kvec_d_in[3] = { {0.0, 0.0, 0.0}, {0.1, 0.2, 0.3}, {0.4, 0.5, 0.6} };
	const int distribution_type_in = 1;
	const bool xprime_in = true;	
	basis_k.initparameters(gamma_only_in, gk_ecut_in, nks_in,kvec_d_in, distribution_type_in, xprime_in);	
	EXPECT_NO_THROW(basis_k.setuptransform());
	EXPECT_EQ(basis_k.npw,3695);
}

TEST_F(PWBasisKTEST, CollectLocalPW) 
{
	ModulePW::PW_Basis_K basis_k(device_flag, precision_single);
	double lat0 = 1.8897261254578281;
	ModuleBase::Matrix3 latvec(10.0,0.0,0.0,
				0.0,10.0,0.0,
				0.0,0.0,10.0);
	double gridecut=10.0;
	basis_k.initgrids(lat0,latvec,gridecut);
	const bool gamma_only_in = true;
	const double gk_ecut_in = 11.0;
	const int nks_in = 3;
	const ModuleBase::Vector3<double> kvec_d_in[3] = { {0.0, 0.0, 0.0}, {0.1, 0.2, 0.3}, {0.4, 0.5, 0.6} };
	const int distribution_type_in = 1;
	const bool xprime_in = true;	
	basis_k.initparameters(gamma_only_in, gk_ecut_in, nks_in,kvec_d_in, distribution_type_in, xprime_in);	
	EXPECT_NO_THROW(basis_k.setuptransform());
	EXPECT_NO_THROW(basis_k.collect_local_pw());
	EXPECT_EQ(basis_k.npw,3695);
	EXPECT_EQ(basis_k.npwk_max,2721);
}

#undef private
#undef protected

// It is weird to init MPI again in this serial test, but it is necessary
// to avoid "Attempting to use an MPI routine before initializing MPICH" in intel compiler.
// The error happens somewhere in the SetupTransform test. 
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