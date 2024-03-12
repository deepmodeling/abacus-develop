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
 * - init_Table
 *   make tables for the radial overlap and kinetic energy integral (as well as
 *   their derivatives) between the atomic orbitals.
 *
 * - Destroy_Table
 *   deallocate the table allocated in init_Table.
 *
 * - _destroy_table
 *   same as Destroy_Table, but does not need input arguments and is automatically
 *   called by the destructor.
 *
 * - init_OV_Tpair
 *   squeeze the element pair index (T1, T2) into a 1-d index and save to OV_Tpair;
 *   calculate 2*lmax+1 for each element pair and save to OV_L2plus1.
 *
 * - init_OV_Opair
 *   squeeze the composite index (l1,l2,ichi1,ichi2) of each element pair into a 1-d
 *   index and save to OV_Opair. The index starts from 0 for each element pair.
 *
 * - get_rmesh
 *   calculate the number of real space mesh points given two radial cutoffs.
 *   the result is made odd to accomodate to Simpson_Integral.
 *
 * - cal_ST_Phi12_R
 *   core subroutine for calculating the radial overlap, kinetic integral, and their
 *   derivatives.
 *   (Eq. A.3 & A.7 of the ABACUS2016 paper)
 *
 * - plot_Table
 *   save input array to file
 *   (only one 1-d array, this function may need some refactoring in the future)
 *
 * Fuctions to be refactored in future:
 * - init_Table_Spherical_Bessel
 * - init_Lmax
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


TEST_F(OrbTablePhiTest, Allocate) {

	otp.allocate(ntype_, lmax_, lcao_.get_kmesh(), Rmax_, dR_, dk_);

	EXPECT_EQ(ntype_, otp.ntype);
    EXPECT_EQ(lmax_, otp.lmax);
    EXPECT_EQ(lcao_.get_kmesh(), otp.kmesh);
    EXPECT_DOUBLE_EQ(Rmax_, otp.Rmax);
    EXPECT_DOUBLE_EQ(dR_, otp.dr);
    EXPECT_DOUBLE_EQ(dk_, otp.dk);

    EXPECT_EQ((2*lmax_+1)*(2*lmax_+1), otp.nlm);

    int nr = calc_nr(Rmax_, dR_);
    EXPECT_EQ(nr, otp.Rmesh);

    for (int ik = 0; ik != lcao_.get_kmesh(); ++ik) {
        EXPECT_DOUBLE_EQ(otp.kpoint[ik], ik*dk_);
        EXPECT_DOUBLE_EQ(otp.kab[ik], dk_);
    }

    for (int ir = 0; ir != nr; ++ir) {
        EXPECT_DOUBLE_EQ(otp.r[ir], ir*dR_);
        EXPECT_DOUBLE_EQ(otp.rab[ir], dR_);
    }

}


TEST_F(OrbTablePhiTest, TwoCenterIntegral) {

	otp.allocate(ntype_, lmax_, lcao_.get_kmesh(), Rmax_, dR_, dk_);
	init_sph_bessel();

	int rmesh = otp.get_rmesh(lcao_.Phi[0].getRcut(), lcao_.Phi[1].getRcut());

	double* rs1 = new double[rmesh];
	double* rs2 = new double[rmesh];

	double* drs1 = new double[rmesh];
	double* drs2 = new double[rmesh];

	std::set<size_t> idx_r = {0, 1, 4};

	// job == 1: calculate overlap
	for (int ir = 0; ir != rmesh; ++ir) {
		rs1[ir] = rs2[ir] = drs1[ir] = drs2[ir] = 0.0;
	}

	otp.cal_ST_Phi12_R(1, 0, lcao_.Phi[0].PhiLN(1,0), lcao_.Phi[1].PhiLN(2,1), rmesh, rs1, drs1);
	otp.cal_ST_Phi12_R(1, 0, lcao_.Phi[0].PhiLN(1,0), lcao_.Phi[1].PhiLN(2,1), idx_r, rs2, drs2);

	for (int ir = 0; ir != rmesh; ++ir) {
		// 1. check the two overloaded versions agree with each other
		if (idx_r.find(ir) != idx_r.end()) {
			EXPECT_DOUBLE_EQ(rs1[ir], rs2[ir]);
			EXPECT_DOUBLE_EQ(drs1[ir], drs2[ir]);
		} else {
			EXPECT_DOUBLE_EQ(rs2[ir], 0.0);
			EXPECT_DOUBLE_EQ(drs2[ir], 0.0);
		}

		// 2. check whether the derivative agrees with finite difference
		double max_tol = 1e-4;
        for (int ir = 5; ir != rmesh-4; ++ir) {
            double drs1_fd = 
                ( +1.0/280*rs1[ir-4] - 4.0/105*rs1[ir-3] 
                  +1.0/5*rs1[ir-2] - 4.0/5*rs1[ir-1] 
                  -1.0/280*rs1[ir+4] + 4.0/105*rs1[ir+3] 
                  -1.0/5*rs1[ir+2] + 4.0/5*rs1[ir+1]
                ) / otp.dr;
            EXPECT_NEAR(drs1_fd, drs1[ir], max_tol);
        }
	}

	// job == 2: calculate kinetic integral
	for (int ir = 0; ir != rmesh; ++ir) {
		rs1[ir] = rs2[ir] = drs1[ir] = drs2[ir] = 0.0;
	}

	otp.cal_ST_Phi12_R(1, 0, lcao_.Phi[0].PhiLN(1,0), lcao_.Phi[1].PhiLN(2,1), rmesh, rs1, drs1);
	otp.cal_ST_Phi12_R(1, 0, lcao_.Phi[0].PhiLN(1,0), lcao_.Phi[1].PhiLN(2,1), idx_r, rs2, drs2);

	for (int ir = 0; ir != rmesh; ++ir) {
		if (idx_r.find(ir) != idx_r.end()) {
			EXPECT_DOUBLE_EQ(rs1[ir], rs2[ir]);
			EXPECT_DOUBLE_EQ(drs1[ir], drs2[ir]);
		} else {
			EXPECT_DOUBLE_EQ(rs2[ir], 0.0);
			EXPECT_DOUBLE_EQ(drs2[ir], 0.0);
		}
	}

	delete[] rs1;
	delete[] rs2;
	delete[] drs1;
	delete[] drs2;
}

TEST_F(OrbTablePhiTest, InitAndDestroyTable) {
	otp.allocate(ntype_, lmax_, lcao_.get_kmesh(), Rmax_, dR_, dk_);
	init_sph_bessel();
	otp.init_OV_Tpair(lcao_);
	otp.init_OV_Opair(lcao_);

	EXPECT_NO_THROW(otp.init_Table(lcao_));
	EXPECT_NE(otp.Table_SR, nullptr);
	EXPECT_NE(otp.Table_TR, nullptr);

	otp.Destroy_Table(lcao_);
	EXPECT_EQ(otp.Table_SR, nullptr);
	EXPECT_EQ(otp.Table_TR, nullptr);
}

TEST_F(OrbTablePhiTest, AutoDestroyTable) {
	otp.allocate(ntype_, lmax_, lcao_.get_kmesh(), Rmax_, dR_, dk_);
	init_sph_bessel();
	otp.init_OV_Tpair(lcao_);
	otp.init_OV_Opair(lcao_);

	otp.init_Table(lcao_);
	EXPECT_NE(otp.Table_SR, nullptr);
	EXPECT_NE(otp.Table_TR, nullptr);

	otp._destroy_table();
	EXPECT_EQ(otp.Table_SR, nullptr);
	EXPECT_EQ(otp.Table_TR, nullptr);
}


TEST_F(OrbTablePhiTest, GetRmesh) {

    dR_ = 0.07;
	otp.allocate(ntype_, lmax_, lcao_.get_kmesh(), Rmax_, dR_, dk_);

    double R1 = 0.0;
    double R2 = 0.0;
    int rmesh = 0;

    // check whether rmesh is an odd number
    for (size_t i = 0; i != 8; ++i) {
        R2 = 0.11*i;
        rmesh = otp.get_rmesh(R1, R2);
        EXPECT_TRUE( (rmesh % 2 == 1) && (rmesh > 0) );
    }

}


TEST_F(OrbTablePhiTest, GetRmeshAbnormal) {

	otp.allocate(ntype_, lmax_, lcao_.get_kmesh(), Rmax_, dR_, dk_);

	testing::internal::CaptureStdout();

    EXPECT_EXIT(otp.get_rmesh(0.5, -1.0), testing::ExitedWithCode(0), "");

	std::string output = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output, testing::HasSubstr("rmesh <= 0"));
}


TEST_F(OrbTablePhiTest, PlotTable) {

	otp.allocate(ntype_, lmax_, lcao_.get_kmesh(), Rmax_, dR_, dk_);
	init_sph_bessel();

	std::string fname = "ORB_table_phi_test.txt";
	int rmesh = otp.get_rmesh(lcao_.Phi[0].getRcut(), lcao_.Phi[1].getRcut());

	double* rs = new double[rmesh];
	double* drs = new double[rmesh];
	otp.cal_ST_Phi12_R(1, 0, lcao_.Phi[0].PhiLN(1,0), lcao_.Phi[1].PhiLN(2,1), rmesh, rs, drs);

	otp.plot_table(fname, rmesh, rs);

	std::ifstream ifs;
	ifs.open(fname.c_str());

	std::string str;
	std::getline(ifs, str); // skip the first line
	
	int ir;
	double val;
	double tol = 1e-5;
	while (std::getline(ifs, str)) {
		std::stringstream ss(str);
		ss >> ir >> val;
		EXPECT_NEAR(rs[ir], val, tol);
	}

	remove(fname.c_str());
}


TEST_F(OrbTablePhiTest, CheckTpair) {
	otp.allocate(ntype_, lmax_, lcao_.get_kmesh(), Rmax_, dR_, dk_);
	otp.init_OV_Tpair(lcao_);

	EXPECT_EQ(otp.OV_Tpair(0,0), 0);
	EXPECT_EQ(otp.OV_Tpair(0,1), 1);
	EXPECT_EQ(otp.OV_Tpair(1,1), 2);
}


TEST_F(OrbTablePhiTest, CheckOpair) {
	otp.allocate(ntype_, lmax_, lcao_.get_kmesh(), Rmax_, dR_, dk_);
	otp.init_OV_Tpair(lcao_);
	otp.init_OV_Opair(lcao_);

	// ***NOTE***
	// OV_Opair(L1,L2,N1,N2)
	// but internal index ordering is L1-N1-L2-N2

	// (H,H)
	EXPECT_EQ(otp.OV_Opair(0, 0, 0, 0, 0), 0);
	EXPECT_EQ(otp.OV_Opair(0, 0, 0, 0, 1), 1);
	EXPECT_EQ(otp.OV_Opair(0, 0, 1, 0, 0), 2);
	EXPECT_EQ(otp.OV_Opair(0, 0, 0, 1, 0), 3);
	EXPECT_EQ(otp.OV_Opair(0, 1, 0, 0, 0), 6);
	EXPECT_EQ(otp.OV_Opair(0, 1, 1, 0, 0), 8);

	// (H,O)
	EXPECT_EQ(otp.OV_Opair(1, 0, 0, 0, 0), 0);
	EXPECT_EQ(otp.OV_Opair(1, 0, 0, 0, 1), 1);
	EXPECT_EQ(otp.OV_Opair(1, 0, 1, 0, 0), 2);
	EXPECT_EQ(otp.OV_Opair(1, 0, 1, 0, 1), 3);
	EXPECT_EQ(otp.OV_Opair(1, 0, 2, 0, 0), 4);
	EXPECT_EQ(otp.OV_Opair(1, 0, 0, 1, 0), 5);
	EXPECT_EQ(otp.OV_Opair(1, 1, 0, 0, 0), 10);
	EXPECT_EQ(otp.OV_Opair(1, 1, 2, 0, 0), 14);

	// (O,O)
	EXPECT_EQ(otp.OV_Opair(2, 0, 0, 0, 0), 0);
	EXPECT_EQ(otp.OV_Opair(2, 0, 0, 0, 1), 1);
	EXPECT_EQ(otp.OV_Opair(2, 0, 1, 0, 0), 2);
	EXPECT_EQ(otp.OV_Opair(2, 0, 1, 0, 1), 3);
	EXPECT_EQ(otp.OV_Opair(2, 0, 2, 0, 0), 4);
	EXPECT_EQ(otp.OV_Opair(2, 0, 0, 1, 0), 5);
	EXPECT_EQ(otp.OV_Opair(2, 1, 0, 0, 0), 10);
	EXPECT_EQ(otp.OV_Opair(2, 1, 0, 1, 0), 15);
	EXPECT_EQ(otp.OV_Opair(2, 2, 0, 0, 0), 20);
	EXPECT_EQ(otp.OV_Opair(2, 2, 2, 0, 0), 24);

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


