#include "module_basis/module_nao/radial_collection.h"
#include "gtest/gtest.h"
#include "module_base/spherical_bessel_transformer.h"

#ifdef __MPI
#include <mpi.h>
#endif

/***********************************************************
 *      Unit test of class "RadialCollection"
 ***********************************************************/
/*!
 *  Tested functions:
 *
 *  - build
 *      - parse a series of orbital or pseudopotential files and
 *        initialize AtomicRadials/BetaRadials objects accordingly.
 *
 *  - copy constructor and assignment operator
 *      - enable deep copy
 *
 *  - cbegin & cend
 *      - pointer-to-pointers that enables the iteration through
 *        non-contiguous read-only NumericalRadial objects
 *
 *  - all "getters"
 *      - Get access to private members.
 *
 *  - all "batch setters"
 *      - Set a property for all RadialSet objects at once
 *                                                                      */
class RadialCollectionTest : public ::testing::Test
{
  protected:
    void SetUp();
    void TearDown();

    RadialCollection orb;                                         //!< object under test
    int nfile = 0; // number of orbital/pseudopotential files
    std::string* file = nullptr; //!< orbitals file to read from
    std::string log_file = "./test_files/radial_collection.log";         //!< file for logging
};

void RadialCollectionTest::SetUp() {
    std::string dir = "../../../../../tests/PP_ORB/";
    nfile = 4;
    file = new std::string[nfile];
    file[0] = dir + "C_gga_8au_100Ry_2s2p1d.orb";
    file[1] = dir + "H_gga_8au_60Ry_2s1p.orb";
    file[2] = dir + "O_gga_10au_100Ry_2s2p1d.orb";
    file[3] = dir + "Fe_gga_9au_100Ry_4s2p2d1f.orb";
}

void RadialCollectionTest::TearDown() {
    delete[] file;
}

TEST_F(RadialCollectionTest, BuildAndGet) {
    orb.build(nfile, file, 'o');

    EXPECT_EQ(orb.symbol(0), "C");
    EXPECT_EQ(orb.symbol(1), "H");
    EXPECT_EQ(orb.symbol(2), "O");
    EXPECT_EQ(orb.symbol(3), "Fe");

    EXPECT_EQ(orb.ntype(), 4);
    EXPECT_EQ(orb.lmax(), 3);
    EXPECT_DOUBLE_EQ(orb.rcut_max(), 10.0);

    EXPECT_EQ(orb.nzeta(0,0), 2);
    EXPECT_EQ(orb.nzeta(0,1), 2);
    EXPECT_EQ(orb.nzeta(0,2), 1);

    EXPECT_EQ(orb.nzeta(1,0), 2);
    EXPECT_EQ(orb.nzeta(1,1), 1);

    EXPECT_EQ(orb.nzeta(2,0), 2);
    EXPECT_EQ(orb.nzeta(2,1), 2);
    EXPECT_EQ(orb.nzeta(2,2), 1);

    EXPECT_EQ(orb.nzeta(3,0), 4);
    EXPECT_EQ(orb.nzeta(3,1), 2);
    EXPECT_EQ(orb.nzeta(3,2), 2);
    EXPECT_EQ(orb.nzeta(3,3), 1);

    EXPECT_EQ(orb.nzeta_max(0), 2);
    EXPECT_EQ(orb.nzeta_max(1), 2);
    EXPECT_EQ(orb.nzeta_max(2), 2);
    EXPECT_EQ(orb.nzeta_max(3), 4);

    EXPECT_EQ(orb.nchi(0), 5);
    EXPECT_EQ(orb.nchi(1), 3);
    EXPECT_EQ(orb.nchi(2), 5);
    EXPECT_EQ(orb.nchi(3), 9);
    EXPECT_EQ(orb.nchi(), 22);

    for (int itype = 0; itype <= 3; ++itype) {
        EXPECT_EQ(orb(itype).itype(), itype);
    }

    for (int itype = 0; itype <= 3; ++itype) {
        for (int l = 0; l <= orb(itype).lmax(); ++l) {
            for (int izeta = 0; izeta != orb(itype).nzeta(l); ++izeta) {
                EXPECT_EQ(orb(itype, l, izeta).l(), l);
            }
        }
    }
}

TEST_F(RadialCollectionTest, BatchSet) {
    orb.build(nfile, file, 'o');

    ModuleBase::SphericalBesselTransformer sbt;
    orb.set_transformer(sbt);
    orb.set_uniform_grid(true, 2001, 20.0);

    // NOTE: cutoff radius is not necessarily the last rgrid point. This is
    // because the grid might have zero padding for the sake of FFT. rcut
    // keeps track of the "actual" cutoff radius.

    EXPECT_EQ(orb.rcut_max(), 10.0);
    std::array<int, 4> rcut = {8, 8, 10, 9};
    for (int itype = 0; itype != orb.ntype(); ++itype) {
        for (int l = 0; l <= orb(itype).lmax(); ++l) {
            for (int izeta = 0; izeta != orb.nzeta(itype, l); ++izeta) {
                EXPECT_EQ(sbt, orb(itype, l, izeta).sbt());
                EXPECT_DOUBLE_EQ(orb(itype, l, izeta).rcut(), rcut[itype]);
            }
        }
    }

    double* grid = new double[3];
    grid[0] = 0.0;
    grid[1] = 1.0;
    grid[2] = 3.14;

    orb.set_grid(true, 3, grid, 'i');
    for (int itype = 0; itype != orb.ntype(); ++itype) {
        for (int l = 0; l <= orb(itype).lmax(); ++l) {
            for (int izeta = 0; izeta != orb.nzeta(itype, l); ++izeta) {
                EXPECT_DOUBLE_EQ(orb(itype, l, izeta).rcut(), 3.14);
            }
        }
    }
    delete[] grid;
}

TEST_F(RadialCollectionTest, Copy)
{
    orb.build(nfile, file, 'o');

    // copy constructor
    RadialCollection orb2(orb);

    EXPECT_EQ(orb2.symbol(0), "C");
    EXPECT_EQ(orb2.symbol(1), "H");
    EXPECT_EQ(orb2.symbol(2), "O");
    EXPECT_EQ(orb2.symbol(3), "Fe");

    EXPECT_EQ(orb2.ntype(), 4);
    EXPECT_EQ(orb2.lmax(), 3);
    EXPECT_DOUBLE_EQ(orb2.rcut_max(), 10.0);

    EXPECT_EQ(orb2.nzeta(0, 0), 2);
    EXPECT_EQ(orb2.nzeta(0, 1), 2);
    EXPECT_EQ(orb2.nzeta(0, 2), 1);

    EXPECT_EQ(orb2.nzeta(1, 0), 2);
    EXPECT_EQ(orb2.nzeta(1, 1), 1);

    EXPECT_EQ(orb2.nzeta(2, 0), 2);
    EXPECT_EQ(orb2.nzeta(2, 1), 2);
    EXPECT_EQ(orb2.nzeta(2, 2), 1);

    EXPECT_EQ(orb2.nzeta(3, 0), 4);
    EXPECT_EQ(orb2.nzeta(3, 1), 2);
    EXPECT_EQ(orb2.nzeta(3, 2), 2);
    EXPECT_EQ(orb2.nzeta(3, 3), 1);

    EXPECT_EQ(orb2.nzeta_max(0), 2);
    EXPECT_EQ(orb2.nzeta_max(1), 2);
    EXPECT_EQ(orb2.nzeta_max(2), 2);
    EXPECT_EQ(orb2.nzeta_max(3), 4);

    EXPECT_EQ(orb2.nchi(0), 5);
    EXPECT_EQ(orb2.nchi(1), 3);
    EXPECT_EQ(orb2.nchi(2), 5);
    EXPECT_EQ(orb2.nchi(3), 9);
    EXPECT_EQ(orb2.nchi(), 22);

    for (int itype = 0; itype <= 3; ++itype)
    {
        EXPECT_EQ(orb2(itype).itype(), itype);
    }

    for (int itype = 0; itype <= 3; ++itype)
    {
        for (int l = 0; l <= orb2(itype).lmax(); ++l)
        {
            for (int izeta = 0; izeta != orb2(itype).nzeta(l); ++izeta)
            {
                EXPECT_EQ(orb2(itype, l, izeta).l(), l);
            }
        }
    }

    // assignment operator
    RadialCollection orb3;
    orb3 = orb;

    EXPECT_EQ(orb3.symbol(0), "C");
    EXPECT_EQ(orb3.symbol(1), "H");
    EXPECT_EQ(orb3.symbol(2), "O");
    EXPECT_EQ(orb3.symbol(3), "Fe");

    EXPECT_EQ(orb3.ntype(), 4);
    EXPECT_EQ(orb3.lmax(), 3);
    EXPECT_DOUBLE_EQ(orb3.rcut_max(), 10.0);

    EXPECT_EQ(orb3.nzeta(0, 0), 2);
    EXPECT_EQ(orb3.nzeta(0, 1), 2);
    EXPECT_EQ(orb3.nzeta(0, 2), 1);

    EXPECT_EQ(orb3.nzeta(1, 0), 2);
    EXPECT_EQ(orb3.nzeta(1, 1), 1);

    EXPECT_EQ(orb3.nzeta(2, 0), 2);
    EXPECT_EQ(orb3.nzeta(2, 1), 2);
    EXPECT_EQ(orb3.nzeta(2, 2), 1);

    EXPECT_EQ(orb3.nzeta(3, 0), 4);
    EXPECT_EQ(orb3.nzeta(3, 1), 2);
    EXPECT_EQ(orb3.nzeta(3, 2), 2);
    EXPECT_EQ(orb3.nzeta(3, 3), 1);

    EXPECT_EQ(orb3.nzeta_max(0), 2);
    EXPECT_EQ(orb3.nzeta_max(1), 2);
    EXPECT_EQ(orb3.nzeta_max(2), 2);
    EXPECT_EQ(orb3.nzeta_max(3), 4);

    EXPECT_EQ(orb3.nchi(0), 5);
    EXPECT_EQ(orb3.nchi(1), 3);
    EXPECT_EQ(orb3.nchi(2), 5);
    EXPECT_EQ(orb3.nchi(3), 9);
    EXPECT_EQ(orb3.nchi(), 22);

    for (int itype = 0; itype <= 3; ++itype)
    {
        EXPECT_EQ(orb3(itype).itype(), itype);
    }

    for (int itype = 0; itype <= 3; ++itype)
    {
        for (int l = 0; l <= orb3(itype).lmax(); ++l)
        {
            for (int izeta = 0; izeta != orb3(itype).nzeta(l); ++izeta)
            {
                EXPECT_EQ(orb3(itype, l, izeta).l(), l);
            }
        }
    }
}

TEST_F(RadialCollectionTest, Iteration)
{
    orb.build(nfile, file, 'o');
    EXPECT_EQ(*orb.cbegin(), &orb(0, 0, 0));
    EXPECT_EQ(*(orb.cbegin() + 2), &orb(1, 0, 0));
    EXPECT_EQ(*(orb.cbegin() + 9), &orb(3, 0, 3));
    EXPECT_EQ(*(orb.cbegin() + 10), &orb(0, 1, 0));
    EXPECT_EQ(*(orb.cbegin() + 17), &orb(0, 2, 0));
    EXPECT_EQ(*(orb.cbegin() + 21), &orb(3, 3, 0));
    EXPECT_EQ(*(orb.cend() - 1), &orb(3, 3, 0));
    //EXPECT_EQ(*(orb.cbegin() + 5), &orb(1, 0, 0));
    //EXPECT_EQ(*(orb.cbegin() + 8), &orb(2, 0, 0));
    //EXPECT_EQ(*(orb.cbegin() + 13), &orb(3, 0, 0));
    //EXPECT_EQ(*(orb.cend() - 1), &orb(3, 3, 0));
}

TEST_F(RadialCollectionTest, Build2) {
    // build a collection of truncated spherical Bessel functions
    int lmax = 3;
    int nbes = 10;
    double rcut = 10.0;
    double sigma = 0.0;
    double dr = 0.01;
    orb.build(lmax, nbes, rcut, sigma, dr);

    orb.lmax();

    EXPECT_EQ(orb.ntype(), 1);
    EXPECT_EQ(orb.lmax(), lmax);
    EXPECT_DOUBLE_EQ(orb.rcut_max(), rcut);

    for (int l = 0; l <= lmax; ++l) {
        EXPECT_EQ(orb.nzeta(0, l), nbes);
    }

    EXPECT_EQ(orb.nzeta_max(0), nbes);
    EXPECT_EQ(orb.nchi(0), nbes*(lmax+1));
    EXPECT_EQ(orb.nchi(), nbes*(lmax+1));
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
