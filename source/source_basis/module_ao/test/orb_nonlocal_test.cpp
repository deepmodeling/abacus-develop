#include "gtest/gtest.h"
#include "source_base/global_variable.h"

#include "source_basis/module_ao/orb_nonlocal.h"

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

	// set_type_info takes rcut_max to be the largest rcut among the projectors.
	// Numerical_Nonlocal_Lm derives its rcut from the last point of the radial
	// mesh it is given, so build each projector through the public
	// set_NL_proj() instead of assigning the private member directly.
	const double rcut[4] = {1.0, 3.0, 4.0, 2.0};
	rcut_max_ = 4.0;

	// set_NL_proj asserts that both meshes are odd and longer than one point;
	// the beta values are irrelevant here, only rcut is read back.
	const int nr = 3;
	const int nk = 3;
	const double dk = 0.01;
	const double dr_uniform = 0.01;
	const double beta_r[nr] = {0.0, 0.0, 0.0};

	nnl.resize(nproj_);
	for (int i = 0; i < nproj_; ++i) {
		const double r_radial[nr] = {0.0, 0.5 * rcut[i], rcut[i]};
		const double rab[nr] = {0.5 * rcut[i], 0.5 * rcut[i], 0.5 * rcut[i]};
		nnl[i].set_NL_proj(elem_label_, ielem_, 0, nr, rab, r_radial, beta_r,
				nk, dk, dr_uniform);
	}
}


void NumericalNonlocalTest::TearDown() {

}


TEST_F(NumericalNonlocalTest, SetTypeInfo) {

	nn.set_type_info(ielem_, elem_label_, type_ps_, lmax_, nproj_, &nnl[0]);

	EXPECT_EQ(nn.getLabel(), elem_label_);
	EXPECT_EQ(nn.getType(), ielem_);
	EXPECT_EQ(nn.getLmax(), lmax_);
	EXPECT_DOUBLE_EQ(nn.get_rcut_max(), rcut_max_);
	EXPECT_EQ(nn.get_nproj(), nproj_);
}


TEST_F(NumericalNonlocalTest, Getters) {

	nn.set_type_info(ielem_, elem_label_, type_ps_, lmax_, nproj_, &nnl[0]);

	// Anchored to the values set_type_info was given, not to the members the
	// getters return -- comparing a getter against its own member can only ever
	// catch a getter wired to the wrong field, and cannot fail otherwise.
	EXPECT_EQ(nn.getLmax(), lmax_);
	EXPECT_EQ(nn.getType(), ielem_);
	EXPECT_EQ(nn.getLabel(), elem_label_);
	EXPECT_EQ(nn.getType_ps(), type_ps_);
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


