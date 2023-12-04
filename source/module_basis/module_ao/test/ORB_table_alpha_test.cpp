#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include "module_base/sph_bessel_recursive.h"

#ifdef __MPI
#include <mpi.h>
#endif

#define private public
#include "module_basis/module_ao/ORB_table_alpha.h"
#undef private


/***********************************************************
 *      unit test of class "ORB_table_alpha"
 ***********************************************************/

/* Tested functions:
 *
 * - allocate
 *   copy the input parameters to class members & initialize k space mesh
 *
 * - init_DS_Opair
 *   make a 1-d index from the composite index (l, ld, ichi, ichid) of each
 *   element where l/ld is the angular momentum of atomic/descriptor orbital 
 *   and ichi/ichid is the index of basis function within that angular momentum;
 *   save the index map to member DS_Opair.
 *
 * - init_DS_2Lplus1
 *   calculate 2*max(lmax, lmaxd)+1 for each element, where lmax/lmaxd is 
 *   the maximum angular momentum of atomic/descriptor orbital.
 *
 * - init_Table_Alpha
 *   make a table for the radial overlap (and its derivative) between the atomic 
 *   and descriptor orbitals.
 *
 * - Destroy_Table_Alpha
 *   deallocate the radial overlap (and its derivative) table (Table_DSR)
 *
 * - get_rmesh
 *   calculate the number of real space mesh points given two radial cutoffs.
 *   the result is made odd to accomodate to Simpson_Integral.
 *
 * - cal_S_PhiAlpha_R
 *   core subroutine for calculating the radial overlap and its derivative
 *   (Eq. A. 3 of the ABACUS2016 paper)
 *
 * - print_Table_DSR (unit test incomplete)
 *   save S(R) table to file
 *
 ***********************************************************/

class OrbTableAlphaTest : public ::testing::Test
{
protected:

    // object under unit test
    ORB_table_alpha ota;

    // orbitals from which the table is calculated
    LCAO_Orbitals lcao_;

    // table for spherical bessel functions
    ModuleBase::Sph_Bessel_Recursive::D2 sbr_;

    void SetUp();
    void TearDown();

    // parameters to initialize lcao_
    std::ofstream ofs_log_;
    int ntype_;
    int lmax_;
    int out_mat_r_; // unused variable
    bool force_flag_;
    int my_rank_;
    bool deepks_setorb_;
    bool read_in_flag_;
    std::string descriptor_file_;
    std::vector<std::string> orbital_file_;
    double ecutwfc_;
    double dk_;
    double dR_;
    double Rmax_;

    // helper
    int calc_nr(double const& Rmax, double const& dR);
    void init_sph_bessel();
};


int OrbTableAlphaTest::calc_nr(double const& Rmax, double const& dR) {
    int nr = static_cast<int>(Rmax / dR) + 4;
    return (nr%2) ? nr : nr+1;
}


void OrbTableAlphaTest::SetUp() {

    // prepare data required to initialized ORB_table_alpha

    ofs_log_.open("ORB_table_alpha_test.log");
    ntype_ = 2;
    lmax_ = 2;
    out_mat_r_ = 0; // unused variable
    force_flag_ = true;
    my_rank_ = 0;
    deepks_setorb_ = true;

    read_in_flag_ = true;
    descriptor_file_ = "./lcao_H2O/jle.orb";
    orbital_file_.push_back("./lcao_H2O/H_gga_8au_60Ry_2s1p.orb");
    orbital_file_.push_back("./lcao_H2O/O_gga_7au_60Ry_2s2p1d.orb");

    // below we mimic ORB_control::read_orb_first
    ecutwfc_ = 123.0;
    dk_ = 0.01;
    dR_ = 0.01;
    Rmax_ = 20;

    lcao_.read_in_flag = read_in_flag_;
    lcao_.descriptor_file = descriptor_file_;
    lcao_.orbital_file = orbital_file_;
#ifdef __MPI
    lcao_.bcast_files(ntype_, GlobalV::MY_RANK);
#endif

    // see ORB_control::read_orb_first
    lcao_.ecutwfc = ecutwfc_;
    lcao_.dk = dk_;
    lcao_.dR = dR_;
    lcao_.Rmax = Rmax_;

    lcao_.Read_Orbitals(ofs_log_, ntype_, lmax_, deepks_setorb_, out_mat_r_, 
            force_flag_, my_rank_);

}


void OrbTableAlphaTest::TearDown() {

}


void OrbTableAlphaTest::init_sph_bessel() {
    // initialize spherical bessel table
	sbr_.set_dx(dR_*dk_);
    // max(l+ld) = 4, but the derivative calculation of j_l relies on j_{l+1}
	sbr_.cal_jlx(5, calc_nr(Rmax_, dR_), lcao_.get_kmesh());
}


TEST_F(OrbTableAlphaTest, Allocate) {

    ota.allocate(ntype_, lmax_, lcao_.get_kmesh(), Rmax_, dR_, dk_);

    EXPECT_EQ(ntype_, ota.ntype);
    EXPECT_EQ(lmax_, ota.lmax);
    EXPECT_EQ(lcao_.get_kmesh(), ota.kmesh);
    EXPECT_DOUBLE_EQ(Rmax_, ota.Rmax);
    EXPECT_DOUBLE_EQ(dR_, ota.dr);
    EXPECT_DOUBLE_EQ(dk_, ota.dk);

    EXPECT_EQ((2*lmax_+1)*(2*lmax_+1), ota.nlm);

    int nr = calc_nr(Rmax_, dR_);
    EXPECT_EQ(nr, ota.Rmesh);

    for (int ik = 0; ik != lcao_.get_kmesh(); ++ik) {
        EXPECT_DOUBLE_EQ(ota.kpoint[ik], ik*dk_);
        EXPECT_DOUBLE_EQ(ota.kab[ik], dk_);
    }

//    for (int ir = 0; ir != nr; ++ir) {
//        EXPECT_DOUBLE_EQ(ota.r[ir], ir*dR_);
//        EXPECT_DOUBLE_EQ(ota.rab[ir], dR_);
//    }
}


TEST_F(OrbTableAlphaTest, AllocateSanityCheck) {

    // make sure allocate() would abort with unphysical input

    EXPECT_EXIT(ota.allocate(0, lmax_, lcao_.get_kmesh(), Rmax_, dR_, dk_),
            testing::KilledBySignal(SIGABRT), "");

    EXPECT_EXIT(ota.allocate(ntype_, -1, lcao_.get_kmesh(), Rmax_, dR_, dk_),
            testing::KilledBySignal(SIGABRT), "");

    EXPECT_EXIT(ota.allocate(ntype_, lmax_, lcao_.get_kmesh(), -1.0, dR_, dk_),
            testing::KilledBySignal(SIGABRT), "");

    EXPECT_EXIT(ota.allocate(ntype_, lmax_, lcao_.get_kmesh(), Rmax_, -0.1, dk_),
            testing::KilledBySignal(SIGABRT), "");

    EXPECT_EXIT(ota.allocate(ntype_, lmax_, lcao_.get_kmesh(), Rmax_, dR_, -0.1),
            testing::KilledBySignal(SIGABRT), "");
}


TEST_F(OrbTableAlphaTest, GetRmesh) {

    dR_ = 0.07;
    ota.allocate(ntype_, lmax_, lcao_.get_kmesh(), Rmax_, dR_, dk_);

    double R1 = 0.0;
    double R2 = 0.0;
    int rmesh = 0;

    // check whether rmesh is an odd number
    for (size_t i = 0; i != 8; ++i) {
        R2 = 0.11*i;
        rmesh = ota.get_rmesh(R1, R2);
        EXPECT_TRUE( (rmesh % 2 == 1) && (rmesh > 0) );
    }

}


TEST_F(OrbTableAlphaTest, GetRmeshAbnormal) {

    ota.allocate(ntype_, lmax_, lcao_.get_kmesh(), Rmax_, dR_, dk_);

	testing::internal::CaptureStdout();

    EXPECT_EXIT(ota.get_rmesh(0.5, -1.0), testing::ExitedWithCode(0), "");

	std::string output = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output, testing::HasSubstr("rmesh <= 0"));
}


TEST_F(OrbTableAlphaTest, DS2Lplus1) {
    ota.allocate(ntype_, lmax_, lcao_.get_kmesh(), Rmax_, dR_, dk_);
    ota.init_DS_2Lplus1(lcao_);
    EXPECT_EQ(ota.DS_2Lplus1[0], 5);
}


TEST_F(OrbTableAlphaTest, CheckOpair) {

    // this test checks whether the indices in ota.DS_Opair are expected

    ota.allocate(ntype_, lmax_, lcao_.get_kmesh(), Rmax_, dR_, dk_);
    ota.init_DS_Opair(lcao_);

    // DS_Opair(elem, l, ld, ichi, ichid)
    //
    // note that the index increase order is different:
    // elem <-- l <-- ichi <-- ld <-- ichid

    // size check
    EXPECT_EQ(ota.DS_Opair.bound1, 2);
    EXPECT_EQ(ota.DS_Opair.bound2, 3);
    EXPECT_EQ(ota.DS_Opair.bound3, 3);
    EXPECT_EQ(ota.DS_Opair.bound4, 2);
    EXPECT_EQ(ota.DS_Opair.bound5, 2);

    // value check
    EXPECT_EQ(ota.DS_Opair(0, 0, 0, 0, 0), 0);
    EXPECT_EQ(ota.DS_Opair(0, 0, 0, 0, 1), 1);
    EXPECT_EQ(ota.DS_Opair(0, 0, 1, 0, 0), 2);
    EXPECT_EQ(ota.DS_Opair(0, 0, 0, 1, 0), 6);
    EXPECT_EQ(ota.DS_Opair(0, 1, 0, 0, 0), 12);

    EXPECT_EQ(ota.DS_Opair(1, 0, 0, 0, 0), 0);
    EXPECT_EQ(ota.DS_Opair(1, 0, 0, 0, 1), 1);
    EXPECT_EQ(ota.DS_Opair(1, 0, 1, 0, 0), 2);
    EXPECT_EQ(ota.DS_Opair(1, 0, 0, 1, 0), 6);
    EXPECT_EQ(ota.DS_Opair(1, 1, 0, 0, 0), 12);
    EXPECT_EQ(ota.DS_Opair(1, 2, 0, 0, 0), 24);
}


TEST_F(OrbTableAlphaTest, FiniteDiffOverlap) {

    // this test checks whether dS(R)/dR agrees with 
    // the finite difference derivative of S(R)

    ota.allocate(ntype_, lmax_, lcao_.get_kmesh(), Rmax_, dR_, dk_);
    ota.init_DS_Opair(lcao_);

    double* S_R = new double[ota.Rmesh];
    double* dS_R = new double[ota.Rmesh];

    init_sph_bessel();

    double max_tol = 1e-4;

    // loop over element type
    for (int ielem = 0; ielem != ntype_; ++ielem) {

        int rmesh = ota.get_rmesh(lcao_.Phi[ielem].getRcut(), 
                lcao_.Alpha[0].getRcut());

        // loop over each numerical atomic orbital (radial part)
        for (int l = 0; l <= lcao_.Phi[ielem].getLmax(); ++l) {
            for (int ichi = 0; ichi != lcao_.Phi[ielem].getNchi(l); ++ichi) {

                // loop over descriptor basis (radial part)
                for (int ld = 0; ld <= lcao_.get_lmax_d(); ++ld) {
                    for (int ichid = 0; ichid != lcao_.Alpha[0].getNchi(ld); ++ichid) {
                        // loop over possible angular momentum
                        // L = |L1-L2|, |L1-L2|+2, ..., L1+L2
                        for (int L = std::abs(l-ld); L <= l+ld; L = L+2) {
                            ota.cal_S_PhiAlpha_R(&sbr_, L, 
                                    lcao_.Phi[ielem].PhiLN(l,ichi), 
                                    lcao_.Alpha[0].PhiLN(ld, ichid), 
                                    rmesh, S_R, dS_R);

                            // finite difference 
                            // note that S_R[0] is unique
                            for (int ir = 5; ir != rmesh-4; ++ir) {
                                double dS_R_fd = 
                                    ( +1.0/280*S_R[ir-4] - 4.0/105*S_R[ir-3] 
                                      +1.0/5*S_R[ir-2] - 4.0/5*S_R[ir-1] 
                                      -1.0/280*S_R[ir+4] + 4.0/105*S_R[ir+3] 
                                      -1.0/5*S_R[ir+2] + 4.0/5*S_R[ir+1]
                                    ) / ota.dr;
                                EXPECT_NEAR(dS_R_fd, dS_R[ir], max_tol);
                            }
                        }
                    }
                }
            }
        }
    }
	delete [] S_R;
	delete [] dS_R;
}


TEST_F(OrbTableAlphaTest, DestroyTable) {

    ota.allocate(ntype_, lmax_, lcao_.get_kmesh(), Rmax_, dR_, dk_);
    ota.init_DS_Opair(lcao_);

    init_sph_bessel();
    ota.init_Table_Alpha(&sbr_, lcao_);

    // should do nothing
    ota.table_allocated = false;
    ota.Destroy_Table_Alpha(lcao_);
    EXPECT_NE(ota.Table_DSR, nullptr);

    // should destroy the table
    ota.table_allocated = true;
    ota.Destroy_Table_Alpha(lcao_);
    EXPECT_EQ(ota.Table_DSR, nullptr);
}


TEST_F(OrbTableAlphaTest, PrintTable) {

    ota.allocate(ntype_, lmax_, lcao_.get_kmesh(), Rmax_, dR_, dk_);
    ota.init_DS_Opair(lcao_);
    ota.init_DS_2Lplus1(lcao_);

    init_sph_bessel();

    ota.init_Table_Alpha(&sbr_, lcao_);
    EXPECT_NO_THROW(ota.print_Table_DSR(lcao_));
    
    // TODO content/format check to be done after code refactoring

    remove("./S_I_mu_alpha.dat");
	ota.Destroy_Table_Alpha(lcao_);
}


TEST_F(OrbTableAlphaTest, AutoDestroyTable) {
    ota.allocate(ntype_, lmax_, lcao_.get_kmesh(), Rmax_, dR_, dk_);
    ota.init_DS_Opair(lcao_);

    init_sph_bessel();
    ota.init_Table_Alpha(&sbr_, lcao_);

    // should do nothing
    ota.table_allocated = false;
    ota._destroy_table();
    EXPECT_NE(ota.Table_DSR, nullptr);

    // should destroy the table
    ota.table_allocated = true;
    ota._destroy_table();
    EXPECT_EQ(ota.Table_DSR, nullptr);
}


int main(int argc, char **argv)
{

#ifdef __MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD,&GlobalV::NPROC);
    MPI_Comm_rank(MPI_COMM_WORLD,&GlobalV::MY_RANK);
#endif

    testing::InitGoogleTest(&argc, argv);
    //GTEST_FLAG_SET(death_test_style, "threadsafe");
    int result = RUN_ALL_TESTS();

#ifdef __MPI
    MPI_Finalize();
#endif

    return result;
} 


