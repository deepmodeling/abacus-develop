#include "module_io/read_wfc_lcao.h"
#include "gtest/gtest.h"

TEST(ReadWfcLcaoTest, ReadAbacusLowfComplex) {
    
    // this test should only be executed on rank 0
#ifdef __MPI
    MPI_Barrier(MPI_COMM_WORLD);
    // printf("MPI environment detected, will use only rank 0\n");
    int iproc = 0, nprocs = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &iproc);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    if (iproc != 0) { GTEST_SKIP(); }
#endif
    int ik = 0;
    ModuleBase::Vector3<double> kvec_c;
    int nbands = 0;
    int nbasis = 0;
    std::vector<std::complex<double>> lowf;
    std::vector<double> ekb;
    std::vector<double> occ;
    double wk = -1.0;

    // first test
    std::string flowf = "./support/LOWF_K_1.txt";
    ModuleIO::read_abacus_lowf(flowf, ik, kvec_c, nbands, nbasis, lowf, ekb, occ, wk);
    EXPECT_EQ(1, ik);
    EXPECT_EQ(19, nbands);
    EXPECT_EQ(63, nbasis);
    EXPECT_EQ(-1.0, wk); // this is not read
    // kvec_c
    EXPECT_EQ(-0.1054038819688572770072454, kvec_c.x);
    EXPECT_EQ(-0.06085495962818950055339684, kvec_c.y);
    EXPECT_EQ(-0.04303095462162523365812206, kvec_c.z);
    // ekb
    EXPECT_EQ(ekb.size(), nbands);
    EXPECT_EQ(-6.0357194547014492158609755e-01, ekb[0]);
    EXPECT_EQ(-5.9803514165671245450539573e-01, ekb[1]);
    EXPECT_EQ(-5.9803514164885962500761707e-01, ekb[2]);
    EXPECT_EQ(1.5316129631211228279141778e+00, ekb[18]);
    // occ
    EXPECT_EQ(occ.size(), nbands);
    EXPECT_EQ(5.8309037900874851473309590e-03, occ[0]);
    EXPECT_EQ(5.8309037900874851473309590e-03, occ[1]);
    EXPECT_EQ(5.8309037900874851473309590e-03, occ[2]);
    EXPECT_EQ(0.0000000000000000000000000e+00, occ[18]);
    // lowf
    EXPECT_EQ(lowf.size(), nbands * nbasis);
    EXPECT_EQ(std::complex<double>(-6.7165115743031689188313749e-03, 2.2594638351419806671094292e-02), lowf[0]);
    EXPECT_EQ(std::complex<double>(1.4318012399152646452193949e-03, 1.4645148803525756542320835e-03), lowf[1]);
    EXPECT_EQ(std::complex<double>(2.3145203328431701583767222e-03, -1.1894969103781807828745798e-03), lowf[2]);
    EXPECT_EQ(std::complex<double>(1.8264875725444186787410628e-03, -2.1179988610690314576601168e-03), lowf[62]);

    // test reuse, expect to overwrite the previous values
    flowf = "./support/LOWF_K_2.txt";
    ModuleIO::read_abacus_lowf(flowf, ik, kvec_c, nbands, nbasis, lowf, ekb, occ, wk);
    EXPECT_EQ(2, ik);
    EXPECT_EQ(19, nbands);
    EXPECT_EQ(63, nbasis);
    EXPECT_EQ(-1.0, wk); // this is not read
    // kvec_c
    EXPECT_EQ(-0.07026925464590484671223436, kvec_c.x);
    EXPECT_EQ(-0.08113994617091932481933725, kvec_c.y);
    EXPECT_EQ(-0.05737460616216696895897087, kvec_c.z);
    // ekb
    EXPECT_EQ(ekb.size(), nbands);
    EXPECT_EQ(-6.0366727725523949832364679e-01, ekb[0]);
    EXPECT_EQ(-5.9786827612267734455286927e-01, ekb[1]);
    EXPECT_EQ(-5.9766242174005113074741757e-01, ekb[2]);
    EXPECT_EQ(1.7385703433295973674432844e+00, ekb[18]);
    // occ
    EXPECT_EQ(occ.size(), nbands);
    EXPECT_EQ(5.8309037900874851473309590e-03, occ[0]);
    EXPECT_EQ(5.8309037900874851473309590e-03, occ[1]);
    EXPECT_EQ(5.8309037900874851473309590e-03, occ[2]);
    EXPECT_EQ(0.0000000000000000000000000e+00, occ[18]);
    // lowf
    EXPECT_EQ(lowf.size(), nbands * nbasis);
    EXPECT_EQ(std::complex<double>(2.0993370513134518295927933e-02, 2.2061937126671427034096951e-03), lowf[0]);
    EXPECT_EQ(std::complex<double>(1.5245441697862502361537906e-03, -3.5413910543949378428862929e-04), lowf[1]);
    EXPECT_EQ(std::complex<double>(1.3119844351431805828944732e-03, -2.4487253841163000855907228e-03), lowf[2]);
    EXPECT_EQ(std::complex<double>(-1.1515848985910659071918438e-03, -1.7994003814514399237911579e-03), lowf[62]);

    // test reuse, the second time
    flowf = "./support/LOWF_K_3.txt";
    ModuleIO::read_abacus_lowf(flowf, ik, kvec_c, nbands, nbasis, lowf, ekb, occ, wk);
    EXPECT_EQ(3, ik);
    EXPECT_EQ(19, nbands);
    EXPECT_EQ(63, nbasis);
    EXPECT_EQ(-1.0, wk); // this is not read
    // kvec_c
    EXPECT_EQ(-0.03513462732295242335611718, kvec_c.x);
    EXPECT_EQ(-0.1014249327136491629630655, kvec_c.y);
    EXPECT_EQ(-0.0717182577027087320153953, kvec_c.z);
    // ekb
    EXPECT_EQ(-6.0466454411080960973379206e-01, ekb[0]);
    EXPECT_EQ(-5.9702547464581856573317964e-01, ekb[1]);
    EXPECT_EQ(-5.9687001897134039918313420e-01, ekb[2]);
    EXPECT_EQ(1.6979332331216445695076800e+00, ekb[18]);
    // occ
    EXPECT_EQ(5.8309037900874851473309590e-03, occ[0]);
    EXPECT_EQ(5.8309037900874851473309590e-03, occ[1]);
    EXPECT_EQ(5.8309037900874851473309590e-03, occ[2]);
    EXPECT_EQ(0.0000000000000000000000000e+00, occ[18]);
    // lowf
    EXPECT_EQ(std::complex<double>(-6.4441007210949141637001958e-03, 2.8610825297509245856986126e-03), lowf[0]);
    EXPECT_EQ(std::complex<double>(-5.8139841531629349313803345e-03, 4.0149570541960916125745484e-03), lowf[1]);
    EXPECT_EQ(std::complex<double>(-2.3215866618005558119630649e-03, 2.6254116611380343138115734e-03), lowf[2]);
    EXPECT_EQ(std::complex<double>(5.5896490231052415979806636e-04, 5.2186638990042212311176728e-04), lowf[62]);
}

TEST(ScatterLowfTest, ScatterLowfComplex)
{
#ifdef __MPI
    // this test should be run on all ranks
    int iproc = 0, nprocs = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &iproc);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    printf("run unittest ScatterLowfTest::ScatterLowfComplex with MPI env:\ntotal number of processes: %d\n", nprocs);
    // one get data on rank 0
    std::vector<std::complex<double>> lowf_glb;
    const int nbands = 10; // nrow_glb
    const int nbasis = 12; // ncol_glb
    if (iproc == 0) {
        lowf_glb.resize(nbands * nbasis);
        for (int i = 0; i < nbands * nbasis; i++) {
            const int j = i / nbands, k = i % nbands;
            lowf_glb[i] = std::complex<double>(k, j);
        }
    }
    /*
    (0,0) (0,1) (0,2) (0,3) (0,4) (0,5) (0,6) (0,7) (0,8) (0,9) (0,10) (0,11) 
    (1,0) (1,1) (1,2) (1,3) (1,4) (1,5) (1,6) (1,7) (1,8) (1,9) (1,10) (1,11) 
    (2,0) (2,1) (2,2) (2,3) (2,4) (2,5) (2,6) (2,7) (2,8) (2,9) (2,10) (2,11) 
    (3,0) (3,1) (3,2) (3,3) (3,4) (3,5) (3,6) (3,7) (3,8) (3,9) (3,10) (3,11)
    ...
    (9,0) (9,1) (9,2) (9,3) (9,4) (9,5) (9,6) (9,7) (9,8) (9,9) (9,10) (9,11)
    */
    // the other ranks get data
    MPI_Barrier(MPI_COMM_WORLD);
    std::vector<std::complex<double>> lowf_loc;
    Parallel_2D para2d_test; // an alternative to ParaV. But to use paraV, the desc would be desc_wfc instead of desc in Parallel_2D
    // initialize a para2d, as if it is paraV
    para2d_test.init(nbands, nbasis, 4, MPI_COMM_WORLD);
    ModuleIO::scatter_lowf(nbasis, nbands, lowf_glb, para2d_test.comm_2D, para2d_test.desc, para2d_test.blacs_ctxt, lowf_loc);
    // make lowf_loc row-major
    for(int i = 0; i < nprocs; i++)
    {
        if(iproc == i)
        {
            printf("rank %d: \n", iproc);
            for(int j = 0; j < lowf_loc.size(); j++)
            {
                const int j1 = j / para2d_test.nrow;
                const int j2 = j % para2d_test.nrow;
                const int j_ = j2 * para2d_test.nrow + j1;
                printf("(%.0f,%.0f)", lowf_loc[j_].real(), lowf_loc[j_].imag());
                if((j+1)%para2d_test.nrow == 0) { printf("\n"); }
            }
            printf("\n");
            usleep(10000);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    MPI_Barrier(MPI_COMM_WORLD);
#endif
}

TEST(ReadWfcLcaoTest, ReadAbacusLowfReal)
{
// not implemented yet
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    // print cwd
    char cwd[1024];
    if (getcwd(cwd, sizeof(cwd)) != nullptr) {
        printf("Current working dir: %s\n", cwd);
    } else {
        perror("getcwd() error");
        return 1;
    }
#ifdef __MPI
    MPI_Init(&argc, &argv);
#endif
    const int status = RUN_ALL_TESTS();
    printf("Unittest read_wfc_lcao_test exits with status %d\n", status);
#ifdef __MPI
    MPI_Finalize();
#endif
    return status;
}