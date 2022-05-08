#include "diago_elpa_utils.h"
#include "module_orbital/ORB_control.h"
#include "mpi.h"
#include "src_io/wf_local.h"
#include "src_lcao/local_orbital_wfc.h"
#include "src_pdiag/pdiag_double.h"
#include "src_pw/charge_broyden.h"
#include "src_pw/occupy.h"

#include "gtest/gtest.h"

#define PASSTHRESHOLD 1e-12
#define DETAILINFO    false
#define PRINT_HS      false
#define REPEATRUN     1

/************************************************
 *  unit test of LCAO diagonalization
 ***********************************************/

/**
 * Tested function:
 *  - ORB_control::setup_2d_division
 *  - Pdiag_Double::diago_double_begin (for Gamma_only)
 *  - Pdiag_Double::diago_complex_double_begin (for multiple K points)
 *
 */

// mock the uncorrected functions
Charge_Broyden::Charge_Broyden(){};
Charge_Broyden::~Charge_Broyden(){};
Charge_Pulay::Charge_Pulay(){};
Charge_Pulay::~Charge_Pulay(){};
Charge_Mixing::Charge_Mixing(){};
Charge_Mixing::~Charge_Mixing(){};
Charge::Charge(){};
Charge::~Charge(){};
namespace GlobalC {Charge_Broyden CHR;};
Local_Orbital_wfc::Local_Orbital_wfc(){};
Local_Orbital_wfc::~Local_Orbital_wfc(){};
void Local_Orbital_wfc::wfc_2d_to_grid(int out_wfc_lcao, const double *wfc_2d, double **wfc_grid){};
void Local_Orbital_wfc::wfc_2d_to_grid(int out_wfc_lcao,
                                       const std::complex<double> *wfc_2d,
                                       std::complex<double> **wfc_grid,
                                       int ik){};
Occupy::Occupy(){};
Occupy::~Occupy(){};
void Occupy::calculate_weights(){};
void WF_Local::write_lowf(const std::string &name, double **ctot){};
void WF_Local::write_lowf_complex(const std::string &name, std::complex<double> **ctot, const int &ik){};
void WF_Local::distri_lowf_complex(std::complex<double> **ctot, std::complex<double> **cc){};
void pdgseps(MPI_Comm comm_2D,
             int n,
             int nb,
             double *A,
             double *B,
             double *Z,
             double *eigen,
             LocalMatrix LM,
             char uplo,
             int &loc_size,
             int &loc_pos){};
void pzgseps(MPI_Comm comm_2D,
             int n,
             int nb,
             int &egnum,
             std::complex<double> *A,
             std::complex<double> *B,
             std::complex<double> *Z,
             double *eigen,
             LocalMatrix LM,
             char uplo,
             int &loc_size,
             int &loc_pos){};

// define class HSMatrix to store HS matrix
class HSMatrix
{
  public:
    std::vector<double> ev;
    std::vector<double> ev_lapack;

    std::vector<double> H_d;
    std::vector<double> S_d;
    std::vector<std::complex<double>> H_cd;
    std::vector<std::complex<double>> S_cd;

    std::vector<double> Hlocal_d;
    std::vector<double> Slocal_d;
    std::vector<double> Stmp_d;
    std::vector<std::complex<double>> Hlocal_cd;
    std::vector<std::complex<double>> Slocal_cd;
    std::vector<std::complex<double>> Stmp_cd;
};

class DiagoPrepare : public Pdiag_Double
{
  public:
    DiagoPrepare(int nlocal, int nbands, int nb2d, int sparsity, bool gamma_only, std::string ks_solver, std::string hfname, std::string sfname)
        : nlocal(nlocal), nbands(nbands), nb2d(nb2d), sparsity(sparsity), gamma_only(gamma_only), ks_solver(ks_solver), hfname(hfname), sfname(sfname)
    {
        if(strcmp("", hfname.c_str()) && strcmp("", sfname.c_str())) 
            readhs = true;
        MPI_Comm_size(MPI_COMM_WORLD, &dsize);
        MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    }

    int dsize, myrank;
    int nk = 1, ik = 0, nspin = 1, nb2d = 0, dcolor = 0, drank = 0;
    int nlocal, nbands, sparsity = 100;
    bool gamma_only, readhs = false;
    std::string calculation = "scf", ks_solver;
    HSMatrix hsm;
    Local_Orbital_wfc lowf;
    double hsolver_time = 0.0, lapack_time = 0.0;
    int na_row, na_col, icontext;
    std::string sfname="", hfname="";

    bool random_HS()
    {
        hsm.ev.resize(nlocal, 0.0);
        hsm.ev_lapack.resize(nlocal, 1.0);
        if (gamma_only)
        {
            hsm.H_d.resize(nlocal * nlocal);
            hsm.S_d.resize(nlocal * nlocal);
            LCAO_DIAGO_TEST::random_hs<double>(hsm.H_d.data(), hsm.S_d.data(), nlocal, sparsity);
        }
        else
        {
            hsm.H_cd.resize(nlocal * nlocal);
            hsm.S_cd.resize(nlocal * nlocal);
            LCAO_DIAGO_TEST::random_hs<std::complex<double>>(hsm.H_cd.data(), hsm.S_cd.data(), nlocal, sparsity);
        }
        return true;
    }

    bool read_HS()
    {
        bool readhfile = false;
        bool readsfile = false;
        if (this->myrank == 0)
        {
            int hdim, sdim;
            if (gamma_only)
            {
                readhfile = LCAO_DIAGO_TEST::read_hs<std::vector<double>>(hfname, hsm.H_d);
                readsfile = LCAO_DIAGO_TEST::read_hs<std::vector<double>>(sfname, hsm.S_d);
                hdim = sqrt(hsm.H_d.size());
                sdim = sqrt(hsm.S_d.size());
            }
            else
            {
                readhfile = LCAO_DIAGO_TEST::read_hs<std::vector<std::complex<double>>>(hfname,hsm.H_cd);
                readsfile = LCAO_DIAGO_TEST::read_hs<std::vector<std::complex<double>>>(sfname,hsm.S_cd);
                hdim = sqrt(hsm.H_cd.size());
                sdim = sqrt(hsm.S_cd.size());
            }
            if (hdim != sdim)
            {
                printf("Error: dimensions of H and S are not equal, %d, %d", hdim, sdim);
                readhfile = readsfile = false;
            }
            nlocal = hdim;
        }
        MPI_Bcast(&nlocal, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&readhfile, 1, MPI_C_BOOL, 0, MPI_COMM_WORLD);
        MPI_Bcast(&readsfile, 1, MPI_C_BOOL, 0, MPI_COMM_WORLD);
        hsm.ev.resize(nlocal, 0.0);
        hsm.ev_lapack.resize(nlocal, 1.0);
        nb2d = 0;
        GlobalV::NB2D = 0;
        nbands = nlocal/2;

        if (readhfile && readsfile)
            return true;
        else
            return false;
    }

    bool produce_HS()
    {
        if (readhs)
            return (this->read_HS());
        else
            return (this->random_HS());
    }

    void print_hs()
    {
        if (!PRINT_HS)
        {
            return;
        }
        if (myrank == 0)
        {
            std::ofstream fp("hmatrix.dat");
            if (gamma_only)
                LCAO_DIAGO_TEST::print_matrix(fp, hsm.H_d.data(), nlocal, nlocal, true);
            else
                LCAO_DIAGO_TEST::print_matrix(fp, hsm.H_cd.data(), nlocal, nlocal, true);
            fp.close();
            fp.open("smatrix.dat");
            if (gamma_only)
                LCAO_DIAGO_TEST::print_matrix(fp, hsm.S_d.data(), nlocal, nlocal, true);
            else
                LCAO_DIAGO_TEST::print_matrix(fp, hsm.S_cd.data(), nlocal, nlocal, true);
            fp.close();
        }
    }

    void pb2d(ORB_control &orb_con)
    {
        std::ofstream ofs_running, ofs_warning;
        orb_con.setup_2d_division(ofs_running, ofs_warning);
        na_row = orb_con.ParaV.nrow;
        na_col = orb_con.ParaV.ncol;
        icontext = orb_con.ParaV.blacs_ctxt;
        std::cout << "orb_con.ParaV.nb=" << orb_con.ParaV.nb << " " << GlobalV::NB2D<< std::endl;
        nb2d = orb_con.ParaV.nb;
    }

    void distribute_data()
    {
        if (gamma_only)
        {
            hsm.Hlocal_d.resize(lowf.ParaV->nrow * lowf.ParaV->ncol);
            hsm.Slocal_d.resize(lowf.ParaV->nrow * lowf.ParaV->ncol);
            hsm.Stmp_d.resize(lowf.ParaV->nrow * lowf.ParaV->ncol, 0.0);

            LCAO_DIAGO_TEST::distribute_data<double>(hsm.H_d.data(),
                                                     hsm.Hlocal_d.data(),
                                                     nlocal,
                                                     nb2d,
                                                     na_row,
                                                     na_col,
                                                     icontext);
            LCAO_DIAGO_TEST::distribute_data<double>(hsm.S_d.data(),
                                                     hsm.Slocal_d.data(),
                                                     nlocal,
                                                     nb2d,
                                                     na_row,
                                                     na_col,
                                                     icontext);
        }
        else
        {
            hsm.Hlocal_cd.resize(lowf.ParaV->nrow * lowf.ParaV->ncol);
            hsm.Slocal_cd.resize(lowf.ParaV->nrow * lowf.ParaV->ncol);
            hsm.Stmp_cd.resize(lowf.ParaV->nrow * lowf.ParaV->ncol, std::complex<double>{0.0, 0.0});
            LCAO_DIAGO_TEST::distribute_data<std::complex<double>>(hsm.H_cd.data(),
                                                                   hsm.Hlocal_cd.data(),
                                                                   nlocal,
                                                                   nb2d,
                                                                   na_row,
                                                                   na_col,
                                                                   icontext);
            LCAO_DIAGO_TEST::distribute_data<std::complex<double>>(hsm.S_cd.data(),
                                                                   hsm.Slocal_cd.data(),
                                                                   nlocal,
                                                                   nb2d,
                                                                   na_row,
                                                                   na_col,
                                                                   icontext);
        }
    }

    void set_env()
    {
        GlobalV::NLOCAL = nlocal;
        GlobalV::KS_SOLVER = ks_solver;
        GlobalV::NBANDS = nbands;
        GlobalV::DSIZE = dsize;
        lowf.wfc_gamma.resize(nspin);
        lowf.wfc_k.resize(nk);
        lowf.wfc_k_grid = new std::complex<double> **[nk];
        GlobalC::CHR.new_e_iteration = true;
    }

    void clean_env()
    {
        delete[] lowf.wfc_k_grid;
    }

    void diago()
    {
        ORB_control
            orb_con(gamma_only, nlocal, nbands, nspin, dsize, nb2d, dcolor, drank, myrank, calculation, ks_solver);
        this->pb2d(orb_con);
        lowf.ParaV = &(orb_con.ParaV);
        this->distribute_data();
        this->set_env();

        double starttime = 0.0, endtime = 0.0;
        starttime = MPI_Wtime();
        for(int i=0;i<REPEATRUN;i++){

            if (gamma_only)
            {
                std::vector<double> htmp = hsm.Hlocal_d; 
                std::vector<double> stmp = hsm.Slocal_d; 
                diago_double_begin(ik, lowf, htmp.data(), stmp.data(), hsm.Stmp_d.data(), hsm.ev.data());
            }
            else
            {
                std::vector<std::complex<double>> htmp = hsm.Hlocal_cd; 
                std::vector<std::complex<double>> stmp = hsm.Slocal_cd;
                diago_complex_begin(ik,lowf,htmp.data(), stmp.data(),hsm.Stmp_cd.data(),hsm.ev.data());
            }
        }
        endtime = MPI_Wtime();
        hsolver_time = (endtime - starttime)/REPEATRUN;
        this->clean_env();
    }

    void diago_lapack()
    {
        double starttime = 0.0, endtime = 0.0;
        starttime = MPI_Wtime();

        for(int i=0;i<REPEATRUN;i++){
            if (gamma_only)
                LCAO_DIAGO_TEST::lapack_diago(hsm.H_d.data(), hsm.S_d.data(), hsm.ev_lapack.data(), nlocal);
            else
                LCAO_DIAGO_TEST::lapack_diago(hsm.H_cd.data(), hsm.S_cd.data(), hsm.ev_lapack.data(), nlocal);
        }

        endtime = MPI_Wtime();
        lapack_time = (endtime - starttime)/REPEATRUN;
    }
};

class DiagoTest : public ::testing::TestWithParam<DiagoPrepare>
{
};

TEST_P(DiagoTest, LCAO)
{
    std::stringstream out_info;
    DiagoPrepare dp = GetParam();
    ASSERT_TRUE(dp.produce_HS());
    dp.print_hs();
    dp.diago();

    if (dp.myrank == 0)
    {
        dp.diago_lapack();
        double maxerror = 0.0;
        int iindex = 0;
        bool pass = true;
        for (int i = 0; i < dp.nbands; i++)
        {
            // EXPECT_NEAR(dp.hsm.ev_lapack[i], dp.hsm.ev[i], PASSTHRESHOLD);
            double error = abs(dp.hsm.ev_lapack[i] - dp.hsm.ev[i]);
            //std::cout << dp.hsm.ev_lapack[i] << " ";
            if (error > maxerror)
            {
                maxerror = error;
                iindex = i;
            }
            if (error > PASSTHRESHOLD)
                pass = false;
        }

        if (dp.readhs)
            out_info << "H/S matrix are read from files." << std::endl;
        else
            out_info << "H/S matrix are produced by random." << std::endl;
        out_info << "solver=" << dp.ks_solver << ", GammaOnly=" << dp.gamma_only;
        out_info << ", NLOCAL=" << dp.nlocal << ", nbands=" << dp.nbands << ", nb2d=" << dp.nb2d;
        if (!dp.readhs)
            out_info << ", Sparsity=" << dp.sparsity;
        out_info << ", solver time: " << dp.hsolver_time
                 << "s, LAPACK time(1 core):" << dp.lapack_time << "s" << std::endl;
        out_info << "Maximum difference between ks_hsolver and LAPACK is " << maxerror << " (" << iindex
                 << "-th eigenvalue), the pass threshold is " << PASSTHRESHOLD << std::endl;

        if (DETAILINFO)
        {
            std::cout << out_info.str();
            out_info.str("");
            out_info.clear();
        }
        EXPECT_TRUE(pass) << out_info.str();
    }
}

INSTANTIATE_TEST_SUITE_P(
    ElpaDoubleTest,
    DiagoTest,
    ::testing::Values( //int nlocal, int nbands, int nb2d, int sparsity, bool gamma_only, std::string ks_solver, bool readhs
        DiagoPrepare(0, 0, 1, 0, true, "genelpa", "H-GammaOnly-large.dat", "S-GammaOnly-large.dat")));

INSTANTIATE_TEST_SUITE_P(
    ElpaComplexDoubleTest,
    DiagoTest,
    ::testing::Values( //int nlocal, int nbands, int nb2d, int sparsity, bool gamma_only, std::string ks_solver, bool readhs
        DiagoPrepare(0, 0, 1, 0, false, "genelpa", "H-KPoints-large.dat", "S-KPoints-large.dat")));

INSTANTIATE_TEST_SUITE_P(
    ScalapackDoubleTest,
    DiagoTest,
    ::testing::Values( //int nlocal, int nbands, int nb2d, int sparsity, bool gamma_only, std::string ks_solver, bool readhs
        //DiagoPrepare(500, 500, 1, 7, true, "scalapack_gvx", "", ""),
        DiagoPrepare(0, 0, 1, 0, true, "scalapack_gvx", "H-GammaOnly.dat", "S-GammaOnly.dat"),
        DiagoPrepare(0, 0, 1, 0, true, "scalapack_gvx", "H-GammaOnly-large.dat", "S-GammaOnly-large.dat")       
        ));

INSTANTIATE_TEST_SUITE_P(
    ScalapackComplexDoubleTest,
    DiagoTest,
    ::testing::Values( //int nlocal, int nbands, int nb2d, int sparsity, bool gamma_only, std::string ks_solver, bool readhs
        //DiagoPrepare(500, 200, 1, 7, false, "scalapack_gvx", "", ""),
        DiagoPrepare(0, 0, 1, 0, false, "scalapack_gvx", "H-KPoints.dat", "S-KPoints.dat"),
        DiagoPrepare(0, 0, 1, 0, false, "scalapack_gvx", "H-KPoints-large.dat", "S-KPoints-large.dat")
        ));

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    testing::InitGoogleTest(&argc, argv);

    int mypnum, dsize;
    MPI_Comm_size(MPI_COMM_WORLD, &dsize);
    MPI_Comm_rank(MPI_COMM_WORLD, &mypnum);
    Parallel_Global::split_diag_world(dsize);
    ::testing::TestEventListeners &listeners = ::testing::UnitTest::GetInstance()->listeners();
    if (mypnum != 0)
    {
        delete listeners.Release(listeners.default_result_printer());
    }

    int result = RUN_ALL_TESTS();
    if (mypnum == 0 && result != 0)
    {
        std::cout << "ERROR:some tests are not passed" << std::endl;
        return result;
    }
    else
    {
        MPI_Finalize();
        return 0;
    }
}
