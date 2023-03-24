#include "gtest/gtest.h"
#include "gmock/gmock.h"

#ifdef __MPI
#include <mpi.h>
#endif

#define private public
#include "module_basis/module_ao/ORB_table_phi.h"
#undef private


/***********************************************************
 *      unit test of class "ORB_table_phi"
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

class OrbTablePhiTest : public ::testing::Test
{
protected:

	// object under unit test
	ORB_table_phi otp;

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


int OrbTablePhiTest::calc_nr(double const& Rmax, double const& dR) {
    int nr = static_cast<int>(Rmax / dR) + 4;
    return (nr%2) ? nr : nr+1;
}


void OrbTablePhiTest::SetUp() {

    // prepare data required to initialized ORB_table_phi

    ofs_log_.open("ORB_table_phi_test.log");
    ntype_ = 2;
    lmax_ = 2;
    out_mat_r_ = 0; // unused variable
    force_flag_ = true;
    my_rank_ = 0;
    deepks_setorb_ = false;

    read_in_flag_ = true;
    orbital_file_.push_back("./lcao_H2O/H_gga_8au_60Ry_2s1p.orb");
    orbital_file_.push_back("./lcao_H2O/O_gga_7au_60Ry_2s2p1d.orb");

    // below we mimic ORB_control::read_orb_first
    ecutwfc_ = 123.0;
    dk_ = 0.01;
    dR_ = 0.01;
    Rmax_ = 20;

    lcao_.read_in_flag = read_in_flag_;
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


void OrbTablePhiTest::TearDown() {

}


void OrbTablePhiTest::init_sph_bessel() {
    // initialize spherical bessel table
	sbr_.set_dx(dR_*dk_);

    // max(l+l) = 4, but the derivative calculation of j_l relies on j_{l+1}
	sbr_.cal_jlx(5, calc_nr(Rmax_, dR_), lcao_.get_kmesh());

	otp.pSB = &sbr_;
}


TEST_F(OrbTablePhiTest, TwoCenterIntegral) {

	// job == 1: calculate overlap



	// job == 2: calculate kinetic integral

}


TEST_F(OrbTablePhiTest, TwoCenterIntegralSelected) {

}


TEST_F(OrbTablePhiTest, InitTable) {

}


TEST_F(OrbTablePhiTest, DestroyTable) {

}


TEST_F(OrbTablePhiTest, AutoDestroyTable) {

}


TEST_F(OrbTablePhiTest, PlotTable) {

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


