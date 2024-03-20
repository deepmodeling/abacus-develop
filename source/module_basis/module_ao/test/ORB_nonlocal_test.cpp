#include "gtest/gtest.h"
#include "module_base/global_variable.h"

#define private public
#include "module_basis/module_ao/ORB_nonlocal.h"
#undef private

#ifdef __MPI
#include <mpi.h>
#endif

/***********************************************************
 *      unit test of class "Numerical_Nonlocal"
 ***********************************************************/

/**
 * Tested functions
 *
 * - set_type_info
 *   copies the input arguments to class members;
 *   set rcut_max to the maximum of rcut of all Numerical_Nonlocal_Lm
 *
 * - all "getters"
 *   get access to class members
 *
 */


class NumericalNonlocalTest : public ::testing::Test
{
protected:

	void SetUp();
	void TearDown();

	// object under unit test
	Numerical_Nonlocal nn;

	// parameters used to initialize the Numerical_Nonlocal object
	std::string elem_label_;
	int ielem_;
	int lmax_;
	double rcut_max_;
	std::string type_ps_;
	int nproj_;
	std::vector<Numerical_Nonlocal_Lm> nnl;
};


void NumericalNonlocalTest::SetUp() {
	elem_label_ = "O";
	ielem_ = 1;
	lmax_ = 2;
	type_ps_ = "NC";
	nproj_ = 4;

	nnl.resize(nproj_);
	nnl[0].rcut = 1.0;
	nnl[1].rcut = 3.0;
	nnl[2].rcut = 4.0;
	nnl[3].rcut = 2.0;
	rcut_max_ = 4.0;
}


void NumericalNonlocalTest::TearDown() {

}


TEST_F(NumericalNonlocalTest, SetTypeInfo) {

	nn.set_type_info(ielem_, elem_label_, type_ps_, lmax_, nproj_, &nnl[0]);

	EXPECT_EQ(nn.label, elem_label_);
	EXPECT_EQ(nn.type, ielem_);
	EXPECT_EQ(nn.lmax, lmax_);
	EXPECT_DOUBLE_EQ(nn.rcut_max, rcut_max_);
	EXPECT_EQ(nn.nproj, nproj_);
}


TEST_F(NumericalNonlocalTest, Getters) {

	nn.set_type_info(ielem_, elem_label_, type_ps_, lmax_, nproj_, &nnl[0]);

	EXPECT_EQ(nn.getLmax(), nn.lmax);
	EXPECT_EQ(nn.getType(), nn.type);
	EXPECT_EQ(nn.getLabel(), nn.label);
	EXPECT_EQ(nn.getType_ps(), nn.type_ps);
	EXPECT_EQ(nn.get_rcut_max(), rcut_max_);
}


int main(int argc, char **argv)
{

#ifdef __MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD,&GlobalV::NPROC);
    MPI_Comm_rank(MPI_COMM_WORLD,&GlobalV::MY_RANK);
#endif
    testing::InitGoogleTest(&argc, argv);
    int result = RUN_ALL_TESTS();

#ifdef __MPI
    MPI_Finalize();
#endif

    return result;
}


