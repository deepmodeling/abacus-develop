#include "module_basis/module_nao/atomic_radials.h"

#include "gtest/gtest.h"
#include "module_base/spherical_bessel_transformer.h"

#ifdef __MPI
#include <mpi.h>
#endif

#include "module_base/constants.h"
#include "module_base/global_variable.h"
#include "module_basis/module_ao/ORB_atomic.h"

using ModuleBase::SphericalBesselTransformer;

/***********************************************************
 *      Unit test of class "AtomicRadials"
 ***********************************************************/
/*!
 *  Tested functions:
 *
 *  - build
 *      - parse an orbital file and initialize the NumericalRadial objects
 *
 *  - copy constructor, assignment operator & polymorphic clone
 *      - enabling deep copy
 *
 *  - cbegin & cend
 *      - pointers to the first and one-past-last read-only NumericalRadial objects
 *
 *  - all "getters"
 *      - Get access to private members.
 *
 *  - all "batch setters"
 *      - Set a property for all NumericalRadial objects at once
 *
 *  - to_numerical_orbital
 *      - Overwrites the content of a Numerical_Orbital object with the current object.
 *
 *                                                                      */
class AtomicRadialsTest : public ::testing::Test
{
  protected:
    void SetUp();
    void TearDown(){};

    AtomicRadials Ti_radials;                                         //!< object under test
    std::string file = "../../../../../tests/PP_ORB/Ti_gga_10au_100Ry_4s2p2d1f.orb"; //!< orbital file to read from
    std::string log_file = "./test_files/atomic_radials.log";         //!< file for logging

    double tol = 1e-12; //!< numerical tolerance for grid & values
};

void AtomicRadialsTest::SetUp()
{
#ifdef __MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &GlobalV::MY_RANK);
#endif
}

TEST_F(AtomicRadialsTest, ReadAndGet)
{
    Ti_radials.build(file, 0, nullptr, GlobalV::MY_RANK);

    EXPECT_EQ(Ti_radials.lmax(), 3);
    EXPECT_EQ(Ti_radials.nzeta(0), 4);
    EXPECT_EQ(Ti_radials.nzeta(1), 2);
    EXPECT_EQ(Ti_radials.nzeta(2), 2);
    EXPECT_EQ(Ti_radials.nzeta(3), 1);
    EXPECT_EQ(Ti_radials.nzeta_max(), 4);
    EXPECT_EQ(Ti_radials.nchi(), 9);
    EXPECT_DOUBLE_EQ(Ti_radials.rcut_max(), 10.0);
    EXPECT_DOUBLE_EQ(Ti_radials.orb_ecut(), 100.0);

    EXPECT_NEAR(Ti_radials.chi(0, 0).rvalue(0), -1.581711853170e-01, tol);
    EXPECT_NEAR(Ti_radials.chi(0, 0).rvalue(4), -1.583907030513e-01, tol);
    EXPECT_NEAR(Ti_radials.chi(0, 0).rvalue(996), -4.183526380009e-05, tol);
    EXPECT_NEAR(Ti_radials.chi(0, 0).rvalue(1000), 0, tol);

    EXPECT_NEAR(Ti_radials.chi(0, 3).rvalue(0), -1.166292682541e+00, tol);
    EXPECT_NEAR(Ti_radials.chi(0, 3).rvalue(4), -1.164223359672e+00, tol);
    EXPECT_NEAR(Ti_radials.chi(0, 3).rvalue(996), -3.183325576529e-04, tol);
    EXPECT_NEAR(Ti_radials.chi(0, 3).rvalue(1000), 0, tol);

    EXPECT_NEAR(Ti_radials.chi(3, 0).rvalue(0), 0, tol);
    EXPECT_NEAR(Ti_radials.chi(3, 0).rvalue(4), 3.744878535962e-05, tol);
    EXPECT_NEAR(Ti_radials.chi(3, 0).rvalue(996), 7.495357740660e-05, tol);
    EXPECT_NEAR(Ti_radials.chi(3, 0).rvalue(1000), 0, tol);
}

TEST_F(AtomicRadialsTest, BatchSet)
{
    int itype = 5;
    Ti_radials.build(file, itype, nullptr, GlobalV::MY_RANK);

    EXPECT_EQ(Ti_radials.itype(), 5);
    EXPECT_EQ(Ti_radials.chi(0, 0).itype(), 5);
    EXPECT_EQ(Ti_radials.chi(0, 3).itype(), 5);
    EXPECT_EQ(Ti_radials.chi(3, 0).itype(), 5);

    SphericalBesselTransformer sbt;
    Ti_radials.set_transformer(sbt);
    EXPECT_EQ(sbt, Ti_radials.chi(0, 0).sbt());
    EXPECT_EQ(sbt, Ti_radials.chi(0, 3).sbt());
    EXPECT_EQ(sbt, Ti_radials.chi(3, 0).sbt());

    Ti_radials.set_uniform_grid(true, 2001, 20.0);
    EXPECT_EQ(Ti_radials.chi(0, 0).nr(), 2001);
    EXPECT_EQ(Ti_radials.chi(0, 3).nr(), 2001);
    EXPECT_EQ(Ti_radials.chi(3, 0).nr(), 2001);
    EXPECT_NEAR(Ti_radials.chi(0, 0).rgrid(2000), 20, tol);
    EXPECT_NEAR(Ti_radials.chi(0, 3).rgrid(1500), 15, tol);
    EXPECT_NEAR(Ti_radials.chi(3, 0).rgrid(1200), 12, tol);
    EXPECT_NEAR(Ti_radials.chi(0, 0).rvalue(2000), 0, tol);
    EXPECT_NEAR(Ti_radials.chi(0, 3).rvalue(1500), 0, tol);
    EXPECT_NEAR(Ti_radials.chi(3, 0).rvalue(1200), 0, tol);

    double grid[5] = {0.0, 1.1, 2.2, 3.3, 4.4};
    Ti_radials.set_grid(true, 5, grid);
    EXPECT_EQ(Ti_radials.chi(0, 0).rgrid(1), 1.1);
    EXPECT_EQ(Ti_radials.chi(0, 3).rgrid(2), 2.2);
    EXPECT_EQ(Ti_radials.chi(3, 0).rgrid(3), 3.3);
}

TEST_F(AtomicRadialsTest, Copy)
{
    /*
     * This test checks whether
     *
     * 1. copy constructor
     * 2. assignment operator
     * 3. polymorphic clone
     *
     * work as expected.
     *                                                                  */
    int itype = 5;
    Ti_radials.build(file, itype, nullptr, GlobalV::MY_RANK);

    // copy constructor
    AtomicRadials Ti_copy(Ti_radials);

    EXPECT_EQ(Ti_copy.itype(), itype);
    EXPECT_EQ(Ti_copy.lmax(), 3);
    EXPECT_EQ(Ti_copy.nzeta(0), 4);
    EXPECT_EQ(Ti_copy.nzeta(1), 2);
    EXPECT_EQ(Ti_copy.nzeta(2), 2);
    EXPECT_EQ(Ti_copy.nzeta(3), 1);
    EXPECT_EQ(Ti_copy.nzeta_max(), 4);
    EXPECT_EQ(Ti_copy.nchi(), 9);
    EXPECT_DOUBLE_EQ(Ti_copy.rcut_max(), 10.0);
    EXPECT_DOUBLE_EQ(Ti_copy.orb_ecut(), 100.0);

    EXPECT_NEAR(Ti_copy.chi(0, 0).rvalue(0), -1.581711853170e-01, tol);
    EXPECT_NEAR(Ti_copy.chi(0, 0).rvalue(4), -1.583907030513e-01, tol);
    EXPECT_NEAR(Ti_copy.chi(0, 0).rvalue(996), -4.183526380009e-05, tol);
    EXPECT_NEAR(Ti_copy.chi(0, 0).rvalue(1000), 0, tol);

    EXPECT_NEAR(Ti_copy.chi(0, 3).rvalue(0), -1.166292682541e+00, tol);
    EXPECT_NEAR(Ti_copy.chi(0, 3).rvalue(4), -1.164223359672e+00, tol);
    EXPECT_NEAR(Ti_copy.chi(0, 3).rvalue(996), -3.183325576529e-04, tol);
    EXPECT_NEAR(Ti_copy.chi(0, 3).rvalue(1000), 0, tol);

    EXPECT_NEAR(Ti_copy.chi(3, 0).rvalue(0), 0, tol);
    EXPECT_NEAR(Ti_copy.chi(3, 0).rvalue(4), 3.744878535962e-05, tol);
    EXPECT_NEAR(Ti_copy.chi(3, 0).rvalue(996), 7.495357740660e-05, tol);
    EXPECT_NEAR(Ti_copy.chi(3, 0).rvalue(1000), 0, tol);

    // assignment operator
    AtomicRadials Ti_assign;
    Ti_assign = Ti_radials;

    EXPECT_EQ(Ti_assign.itype(), itype);
    EXPECT_EQ(Ti_assign.lmax(), 3);
    EXPECT_EQ(Ti_assign.nzeta(0), 4);
    EXPECT_EQ(Ti_assign.nzeta(1), 2);
    EXPECT_EQ(Ti_assign.nzeta(2), 2);
    EXPECT_EQ(Ti_assign.nzeta(3), 1);
    EXPECT_EQ(Ti_assign.nzeta_max(), 4);
    EXPECT_EQ(Ti_assign.nchi(), 9);
    EXPECT_DOUBLE_EQ(Ti_assign.rcut_max(), 10.0);
    EXPECT_DOUBLE_EQ(Ti_assign.orb_ecut(), 100.0);

    EXPECT_NEAR(Ti_assign.chi(0, 0).rvalue(0), -1.581711853170e-01, tol);
    EXPECT_NEAR(Ti_assign.chi(0, 0).rvalue(4), -1.583907030513e-01, tol);
    EXPECT_NEAR(Ti_assign.chi(0, 0).rvalue(996), -4.183526380009e-05, tol);
    EXPECT_NEAR(Ti_assign.chi(0, 0).rvalue(1000), 0, tol);

    EXPECT_NEAR(Ti_assign.chi(0, 3).rvalue(0), -1.166292682541e+00, tol);
    EXPECT_NEAR(Ti_assign.chi(0, 3).rvalue(4), -1.164223359672e+00, tol);
    EXPECT_NEAR(Ti_assign.chi(0, 3).rvalue(996), -3.183325576529e-04, tol);
    EXPECT_NEAR(Ti_assign.chi(0, 3).rvalue(1000), 0, tol);

    EXPECT_NEAR(Ti_assign.chi(3, 0).rvalue(0), 0, tol);
    EXPECT_NEAR(Ti_assign.chi(3, 0).rvalue(4), 3.744878535962e-05, tol);
    EXPECT_NEAR(Ti_assign.chi(3, 0).rvalue(996), 7.495357740660e-05, tol);
    EXPECT_NEAR(Ti_assign.chi(3, 0).rvalue(1000), 0, tol);

    // polymorphic clone
    RadialSet* ptr_Ti_polyclone = Ti_radials.clone();

    EXPECT_EQ(ptr_Ti_polyclone->itype(), itype);
    EXPECT_EQ(ptr_Ti_polyclone->lmax(), 3);
    EXPECT_EQ(ptr_Ti_polyclone->nzeta(0), 4);
    EXPECT_EQ(ptr_Ti_polyclone->nzeta(1), 2);
    EXPECT_EQ(ptr_Ti_polyclone->nzeta(2), 2);
    EXPECT_EQ(ptr_Ti_polyclone->nzeta(3), 1);
    EXPECT_EQ(ptr_Ti_polyclone->nzeta_max(), 4);
    EXPECT_EQ(ptr_Ti_polyclone->nchi(), 9);
    EXPECT_DOUBLE_EQ(ptr_Ti_polyclone->rcut_max(), 10.0);

    EXPECT_NEAR(ptr_Ti_polyclone->chi(0, 0).rvalue(0), -1.581711853170e-01, tol);
    EXPECT_NEAR(ptr_Ti_polyclone->chi(0, 0).rvalue(4), -1.583907030513e-01, tol);
    EXPECT_NEAR(ptr_Ti_polyclone->chi(0, 0).rvalue(996), -4.183526380009e-05, tol);
    EXPECT_NEAR(ptr_Ti_polyclone->chi(0, 0).rvalue(1000), 0, tol);

    EXPECT_NEAR(ptr_Ti_polyclone->chi(0, 3).rvalue(0), -1.166292682541e+00, tol);
    EXPECT_NEAR(ptr_Ti_polyclone->chi(0, 3).rvalue(4), -1.164223359672e+00, tol);
    EXPECT_NEAR(ptr_Ti_polyclone->chi(0, 3).rvalue(996), -3.183325576529e-04, tol);
    EXPECT_NEAR(ptr_Ti_polyclone->chi(0, 3).rvalue(1000), 0, tol);

    EXPECT_NEAR(ptr_Ti_polyclone->chi(3, 0).rvalue(0), 0, tol);
    EXPECT_NEAR(ptr_Ti_polyclone->chi(3, 0).rvalue(4), 3.744878535962e-05, tol);
    EXPECT_NEAR(ptr_Ti_polyclone->chi(3, 0).rvalue(996), 7.495357740660e-05, tol);
    EXPECT_NEAR(ptr_Ti_polyclone->chi(3, 0).rvalue(1000), 0, tol);

    delete ptr_Ti_polyclone;

    // normal clone
    AtomicRadials* ptr_Ti_clone = Ti_radials.clone();
    EXPECT_DOUBLE_EQ(ptr_Ti_clone->orb_ecut(), 100.0);

    delete ptr_Ti_clone;
}

TEST_F(AtomicRadialsTest, BeginAndEnd)
{
    int itype = 5;
    Ti_radials.build(file, itype, nullptr, GlobalV::MY_RANK);

    EXPECT_EQ(Ti_radials.cbegin(), &Ti_radials.chi(0, 0));
    EXPECT_EQ(Ti_radials.cend() - 1, &Ti_radials.chi(3, 0));
}

TEST_F(AtomicRadialsTest, ToNumericalOrbital)
{
    int itype = 5;
    Ti_radials.build(file, itype, nullptr, GlobalV::MY_RANK);

    ModuleBase::SphericalBesselTransformer sbt;
    Ti_radials.set_transformer(sbt);
    Ti_radials.set_uniform_grid(false, 1001, 30.0, 't');

    Numerical_Orbital no;
    Ti_radials.to_numerical_orbital(no);

    EXPECT_EQ(Ti_radials.lmax(), no.getLmax());
    EXPECT_EQ(Ti_radials.rcut_max(), no.getRcut());
    EXPECT_EQ(Ti_radials.itype(), no.getType());
    EXPECT_EQ(Ti_radials.nchi(), no.getTotal_nchi());
    EXPECT_EQ(Ti_radials.symbol(), no.getLabel());

    for (int l = 0; l <= Ti_radials.lmax(); l++)
    {
        EXPECT_EQ(Ti_radials.nzeta(l), no.getNchi(l));
    }
}

int main(int argc, char** argv)
{

#ifdef __MPI
    MPI_Init(&argc, &argv);
#endif

    testing::InitGoogleTest(&argc, argv);
    int result = RUN_ALL_TESTS();

#ifdef __MPI
    MPI_Finalize();
#endif

    return result;
}
