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

TEST(ReadWfcLcaoTest, Pzgemr2dUseTest)
{
/*
(0,0) (0,1) (0,2) (0,3) (0,4) (0,5) (0,6) (0,7) (0,8) (0,9) (0,10) (0,11) 
(1,0) (1,1) (1,2) (1,3) (1,4) (1,5) (1,6) (1,7) (1,8) (1,9) (1,10) (1,11) 
(2,0) (2,1) (2,2) (2,3) (2,4) (2,5) (2,6) (2,7) (2,8) (2,9) (2,10) (2,11) 
(3,0) (3,1) (3,2) (3,3) (3,4) (3,5) (3,6) (3,7) (3,8) (3,9) (3,10) (3,11)
...
(9,0) (9,1) (9,2) (9,3) (9,4) (9,5) (9,6) (9,7) (9,8) (9,9) (9,10) (9,11)
*/
#ifdef __MPI
    // this test should be run on all ranks
    int iproc = 0, nprocs = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &iproc);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    if(nprocs != 4)
    {
        if(iproc == 0) { printf("Please run this unittest with nprocs = 4.\n"); } 
        GTEST_SKIP();
    }
    // one get data on rank 0
    std::vector<std::complex<double>> lowf_glb;
    const int nbands = 10; // nrow_glb
    const int nbasis = 12; // ncol_glb
    if (iproc == 0) {
        printf("Run unittest ScatterLowfTest::ScatterLowfComplex with MPI env:\n");
        printf("Total number of processes: %d\n\n", nprocs);
        printf("Row-major processor grid is used.\n\n");
        printf("First test the \"column-major\" matrix, which means for columns their memory\n");
        printf("are consecutive. The matrix in form of (i, j), where i runs over [0, 9] and \n");
        printf("j runs over [0, 11]: (0,0), (1,0), (2,0), ..., (9,0), (0,1), ...\n");
        lowf_glb.resize(nbands * nbasis);
        for (int i = 0; i < nbands * nbasis; i++) {
            const int j = i / nbands, k = i % nbands;
            lowf_glb[i] = std::complex<double>(k, j);
        }
    }
    // the other ranks get data
    MPI_Barrier(MPI_COMM_WORLD);
    std::vector<std::complex<double>> lowf_loc;
    Parallel_2D para2d_test; // an alternative to ParaV. But to use paraV, the desc would be desc_wfc instead of desc in Parallel_2D
    // initialize a para2d, as if it is paraV. 
    // BE CAREFUL! One must be careful about defining the dimension properly. Now the original 
    // matrix is made column-memory-consecutive, thus the vectorized data must be cut by "the 
    // number of rows", then by "the number of columns".
    // nbands is "the number of rows", nbasis is "the number of columns"
    para2d_test.init(nbands, nbasis, 4, MPI_COMM_WORLD);
    Parallel_2D para2d_glb;
    para2d_glb.init(nbands, nbasis, std::max(nbands, nbasis), MPI_COMM_WORLD);
    lowf_loc.resize(para2d_test.nrow * para2d_test.ncol); // the nrow and ncol are automatically
                                                          // calculated by Parallel_2D
    // wait, what is the meaning of nrow here? The "nrow" again is just the number to cut
    // the vectorized data into a matrix. The "nrow" is the number of rows in the matrix but
    // remember the matrix is column-memory-consecutive.
    const int one = 1;
    // the following function can do the scattering-gathering automatically.
    // a double counterpart is pdgemr2d_, int counterpart is pigemr2d_
    // Those in C style are Cpzgemr2d, Cdgemr2d, Cigemr2d...
    pzgemr2d_(&nbands, &nbasis, 
              lowf_glb.data(), &one, &one, para2d_glb.desc, 
              lowf_loc.data(), &one, &one, para2d_test.desc, 
              &(para2d_glb.blacs_ctxt));
    // what will happen if impose a row-major processor grid onto a column-major matrix?
    // you can get correct results, expect each block is column-major.
    // Have a look:
    MPI_Barrier(MPI_COMM_WORLD);
    std::vector<int> sizes_loc = {4*4*3, 4*4 + 4*2, 4*4*2, 4*4};
    std::vector<std::vector<double>> reals = {
        {0, 1, 2, 3, 8, 9},
        {0, 1, 2, 3, 8, 9},
        {4, 5, 6, 7},
        {4, 5, 6, 7}
    };
    std::vector<std::vector<double>> imags = {
        {0, 1, 2, 3, 8, 9, 10, 11},
        {4, 5, 6, 7},
        {0, 1, 2, 3, 8, 9, 10, 11},
        {4, 5, 6, 7}
    };
    for(int i = 0; i < nprocs; i++)
    {
        if(iproc == i)
        {
            EXPECT_EQ(lowf_loc.size(), sizes_loc[i]);
            printf(">>> rank %d: \n", iproc);
            printf("First print scattered matrix in the way that ELEMENTS WITH CONSECUTIVE\n");
            printf("MEMORIES ARE SHOWN (only shown) IN THE SAME LINE:\n");
            for(int j = 0; j < lowf_loc.size(); j++)
            {
                printf("(%2.0f,%2.0f)", lowf_loc[j].real(), lowf_loc[j].imag());
                if((j + 1)%para2d_test.nrow == 0) { printf("\n"); }
                const int k = j%para2d_test.nrow;
                EXPECT_EQ(lowf_loc[j].real(), reals[i][k]);
                const int l = j/para2d_test.nrow;
                EXPECT_EQ(lowf_loc[j].imag(), imags[i][l]);
            }
            printf("Or INDEX IT to show like \"row-major\":\n");
            // (i, j) -> (i', j') with i = j' and j = i'
            // x = i*ncol + j, x' = i'*ncol' + j' with ncol' = nrow and nrow' = ncol
            // i = x/ncol, j = x%ncol, x' = j*nrow + i = x%ncol*nrow + x/ncol
            for(int j = 0; j < lowf_loc.size(); j++)
            {
                const int x = j%para2d_test.ncol*para2d_test.nrow + j/para2d_test.ncol;
                printf("(%2.0f,%2.0f)", lowf_loc[x].real(), lowf_loc[x].imag());
                if((j + 1)%para2d_test.ncol == 0) { printf("\n"); }
            }
            printf("\n");
            usleep(10000);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    // test the other way around, the row-major matrix
    if (iproc == 0) {
        printf("Now test the \"row-major\" matrix, which means for rows their memory\n");
        printf("are consecutive. The matrix in form of (i, j), where i runs over [0, 9] and \n");
        printf("j runs over [0, 11]: (0,0), (0,1), (0,2), ..., (0,11), (1,0), ...\n");
        for (int i = 0; i < nbands * nbasis; i++) {
            const int irow = i / nbasis, icol = i % nbasis;
            lowf_glb[i] = std::complex<double>(irow, icol);
        }
    }
    // the other ranks get data
    MPI_Barrier(MPI_COMM_WORLD);
    // initialize a para2d, as if it is paraV.
    Parallel_2D para2d_test_prime;
    // note! this time the memory is first cut by "the number of columns", 
    // then by "the number of rows". Therefore the "nbasis" is put the first.
    // This is how ScaLAPCK defines a matrix: the first number defines the leading dimension.
    para2d_test_prime.init(nbasis, nbands, 4, MPI_COMM_WORLD);
    Parallel_2D para2d_glb_prime;
    para2d_glb_prime.init(nbasis, nbands, std::max(nbands, nbasis), MPI_COMM_WORLD);
    lowf_loc.resize(para2d_test_prime.nrow * para2d_test_prime.ncol);
    pzgemr2d_(&nbasis, &nbands, 
              lowf_glb.data(), &one, &one, para2d_glb_prime.desc, 
              lowf_loc.data(), &one, &one, para2d_test_prime.desc, 
              &(para2d_glb_prime.blacs_ctxt));
    MPI_Barrier(MPI_COMM_WORLD);
    sizes_loc = {4*4*3, 4*4*2, 4*4 + 4*2, 4*4};
    reals = {
        {0, 1, 2, 3, 8, 9},
        {4, 5, 6, 7},
        {0, 1, 2, 3, 8, 9},
        {4, 5, 6, 7}
    };
    imags = {
        {0, 1, 2, 3, 8, 9, 10, 11},
        {0, 1, 2, 3, 8, 9, 10, 11},
        {4, 5, 6, 7},
        {4, 5, 6, 7}
    };
    for(int i = 0; i < nprocs; i++)
    {
        if(iproc == i)
        {
            EXPECT_EQ(lowf_loc.size(), sizes_loc[i]);
            printf(">>> rank %d: \n", iproc);
            for(int j = 0; j < lowf_loc.size(); j++)
            {
                printf("(%2.0f,%2.0f)", lowf_loc[j].real(), lowf_loc[j].imag());
                if((j + 1)%para2d_test_prime.nrow == 0) { printf("\n"); }
                const int k = j/para2d_test_prime.nrow;
                EXPECT_EQ(lowf_loc[j].real(), reals[i][k]);
                const int l = j%para2d_test_prime.nrow;
                EXPECT_EQ(lowf_loc[j].imag(), imags[i][l]);
            }
            printf("\n");
            usleep(10000);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if(iproc == 0)
    {
        printf("BE CAREFUL!\n");
        printf("You note that the PROCESSOR GRID seems to be transposed. It is because\n");
        printf("in C/C++ it is always assumed memory in the same row is consecutive, while\n");
        printf("in FORTRAN or \"what ScaLAPACK supposes\" it is column-memory-consecutive.\n"); 
        usleep(10000);
    }
    MPI_Barrier(MPI_COMM_WORLD);
#else
    printf("Run unittest ScatterLowfTest::ScatterLowfComplex without MPI env:\n");
    printf("This test is not executed because MPI is not enabled.\n");
    GTEST_SKIP();
#endif
}

TEST(ReadWfcLcaoTest, ReadAbacusLowfReal)
{
    // this test should only be executed on rank 0
#ifdef __MPI
    MPI_Barrier(MPI_COMM_WORLD);
    // printf("MPI environment detected, will use only rank 0\n");
    int iproc = 0, nprocs = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &iproc);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    if (iproc != 0) { GTEST_SKIP(); }
#endif
    int ik = 1; // should be overwritten to 0
    ModuleBase::Vector3<double> kvec_c;
    int nbands = 0;
    int nbasis = 0;
    std::vector<double> lowf;
    std::vector<double> ekb;
    std::vector<double> occ;
    double wk = -1.0; // should be overwritten to 1.0

    // first test
    const std::string flowf = "./support/LOWF_GAMMA_S1.txt";
    ModuleIO::read_abacus_lowf(flowf, ik, kvec_c, nbands, nbasis, lowf, ekb, occ, wk);
    EXPECT_EQ(0, ik);
    EXPECT_EQ(24, nbands);
    EXPECT_EQ(31, nbasis);
    EXPECT_EQ(1.0, wk);
    // kvec_c, gamma point, 0, 0, 0
    EXPECT_EQ(0.0, kvec_c.x);
    EXPECT_EQ(0.0, kvec_c.y);
    EXPECT_EQ(0.0, kvec_c.z);
    // ekb
    EXPECT_EQ(ekb.size(), nbands);
    EXPECT_EQ(-1.22787155e+00, ekb[0]);
    EXPECT_EQ(-3.10595658e-01, ekb[1]);
    EXPECT_EQ(-3.00546690e-01, ekb[2]);
    EXPECT_EQ(3.17982970e+00, ekb[23]);
    // occ
    EXPECT_EQ(occ.size(), nbands);
    EXPECT_EQ(2.00000000e+00, occ[0]);
    EXPECT_EQ(2.00000000e+00, occ[1]);
    EXPECT_EQ(2.00000000e+00, occ[2]);
    EXPECT_EQ(0.00000000e+00, occ[23]);
    // lowf
    EXPECT_EQ(lowf.size(), nbands * nbasis);
    EXPECT_EQ(-1.51728369e-02, lowf[0]);
    EXPECT_EQ(-2.07808444e-03, lowf[1]);
    EXPECT_EQ(1.21298954e-17, lowf[2]);
    EXPECT_EQ(-5.44883791e-09, lowf[30]);
}

TEST(ReadWfcLcaoTest, Pdgemr2dUseTest)
{
    // you can find more information in unittest Pzgemr2dUseTest, present test
    // works identically to the previous one, but with real numbers.
#ifdef __MPI
    // this test should be run on all ranks
    int iproc = 0, nprocs = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &iproc);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    if(nprocs != 4)
    {
        if(iproc == 0) { printf("Please run this unittest with nprocs = 4.\n"); } 
        GTEST_SKIP();
    }
    std::vector<double> lowf_glb;
    const int nbands = 10;
    const int nbasis = 12;
    // still, the expected matrix is organized as row-memory-consecutive but here
    // just make the matrix column-memory-consecutive.
    // x = i*ncol + j, x' = j*nrow + i
    // i = x/ncol, j = x%ncol, x' = j*nrow + i = x%ncol*nrow + x/ncol
    if (iproc == 0) {
        lowf_glb.resize(nbands * nbasis);
        for (int i = 0; i < nbands * nbasis; i++) {
            lowf_glb[i] = i%nbasis*nbands + i/nbasis;
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    std::vector<double> lowf_loc;
    Parallel_2D para2d_test;
    para2d_test.init(nbasis, nbands, 4, MPI_COMM_WORLD);
    Parallel_2D para2d_glb;
    para2d_glb.init(nbasis, nbands, std::max(nbands, nbasis), MPI_COMM_WORLD);
    lowf_loc.resize(para2d_test.nrow * para2d_test.ncol);
    const int one = 1;
    pdgemr2d_(&nbasis, &nbands, 
              lowf_glb.data(), &one, &one, para2d_glb.desc, 
              lowf_loc.data(), &one, &one, para2d_test.desc, 
              &(para2d_glb.blacs_ctxt));
    MPI_Barrier(MPI_COMM_WORLD);
    std::vector<int> sizes_loc = {4*4*3, 4*4*2, 4*4 + 4*2, 4*4};
    for(int i = 0; i < nprocs; i++)
    {
        if(iproc == i) { EXPECT_EQ(lowf_loc.size(), sizes_loc[i]); }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    // test the other way around, the row-major matrix
    if (iproc == 0) {
        for (int i = 0; i < nbands * nbasis; i++) {
            lowf_glb[i] = i;
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    Parallel_2D para2d_test_prime;
    para2d_test_prime.init(nbands, nbasis, 4, MPI_COMM_WORLD);
    Parallel_2D para2d_glb_prime;
    para2d_glb_prime.init(nbands, nbasis, std::max(nbands, nbasis), MPI_COMM_WORLD);
    lowf_loc.resize(para2d_test_prime.nrow * para2d_test_prime.ncol);
    pdgemr2d_(&nbands, &nbasis, 
              lowf_glb.data(), &one, &one, para2d_glb_prime.desc, 
              lowf_loc.data(), &one, &one, para2d_test_prime.desc, 
              &(para2d_glb_prime.blacs_ctxt));
    MPI_Barrier(MPI_COMM_WORLD);
    sizes_loc = {4*4*3, 4*4 + 4*2, 4*4*2, 4*4};
    for(int i = 0; i < nprocs; i++)
    {
        if(iproc == i) { EXPECT_EQ(lowf_loc.size(), sizes_loc[i]); }
    }
    MPI_Barrier(MPI_COMM_WORLD);
#endif
}

TEST(ReadWfcLcaoTest, RestartFromFileParallel)
{
#ifdef __MPI
    MPI_Barrier(MPI_COMM_WORLD);
    // printf("MPI environment detected, will use only rank 0\n");
    int iproc = 0, nprocs = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &iproc);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    if(iproc == 0)
    {
        printf("Run unittest ReadWfcLcaoTest::RestartFromFileParallel with MPI env:\n");
        printf("Waiting for PR#4268 being merge, can this function be fully implemented.\n");
        usleep(10000);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    GTEST_SKIP();
#else
    printf("Run unittest ReadWfcLcaoTest::RestartFromFileParallel without MPI env:\n");
    printf("This test is not executed because MPI is not enabled.\n");
    GTEST_SKIP();
#endif
}

TEST(ReadWfcLcaoTest, RestartFromFileSerial)
{
#ifdef __MPI
    MPI_Barrier(MPI_COMM_WORLD);
    // printf("MPI environment detected, will use only rank 0\n");
    int iproc = 0, nprocs = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &iproc);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    if (iproc != 0) { GTEST_SKIP(); }
#endif
    const int nks = 4;
    int nbands = -1, nbasis = -1;
    std::vector<std::complex<double>> lowf;
    std::vector<double> ekb;
    std::vector<double> occ;
    std::vector<ModuleBase::Vector3<double>> kvec_c;
    std::vector<double> wk;
    const std::string out_dir = "./support";
    ModuleIO::restart_from_file(out_dir, nks, nbands, nbasis, lowf, ekb, occ, kvec_c, wk);
    EXPECT_EQ(19, nbands);
    EXPECT_EQ(63, nbasis);
    EXPECT_EQ(nks*nbands*nbasis, lowf.size());
    EXPECT_EQ(nks*nbands, ekb.size());
    EXPECT_EQ(nks*nbands, occ.size());
    EXPECT_EQ(nks, kvec_c.size());
    EXPECT_EQ(nks, wk.size());
    
    // test the first k-point
    int nbands_k0 = -1, nbasis_k0 = -1;
    int ik_k0;
    std::vector<std::complex<double>> lowf_k0;
    std::vector<double> ekb_k0;
    std::vector<double> occ_k0;
    ModuleBase::Vector3<double> kvec_c_k0;
    double wk_k0 = -1.0;
    ModuleIO::read_abacus_lowf("./support/LOWF_K_1.txt", ik_k0, kvec_c_k0, nbands_k0, nbasis_k0, lowf_k0, ekb_k0, occ_k0, wk_k0);

    EXPECT_EQ(1, ik_k0);
    EXPECT_EQ(19, nbands_k0);
    EXPECT_EQ(63, nbasis_k0);
    EXPECT_EQ(-1.0, wk_k0); // this is not read
    // kvec_c
    EXPECT_EQ(kvec_c_k0.x, kvec_c[0].x);
    EXPECT_EQ(kvec_c_k0.y, kvec_c[0].y);
    EXPECT_EQ(kvec_c_k0.z, kvec_c[0].z);
    // ekb
    for(int i = 0; i < nbands_k0; i++) { EXPECT_EQ(ekb_k0[i], ekb[i]); }
    // occ
    for(int i = 0; i < nbands_k0; i++) { EXPECT_EQ(occ_k0[i], occ[i]); }
    // lowf
    for(int i = 0; i < nbands_k0*nbasis_k0; i++) { EXPECT_EQ(lowf_k0[i], lowf[i]); }

    // test the second k-point
    int nbands_k1 = -1, nbasis_k1 = -1;
    int ik_k1;
    std::vector<std::complex<double>> lowf_k1;
    std::vector<double> ekb_k1;
    std::vector<double> occ_k1;
    ModuleBase::Vector3<double> kvec_c_k1;
    double wk_k1 = -1.0;
    ModuleIO::read_abacus_lowf("./support/LOWF_K_2.txt", ik_k1, kvec_c_k1, nbands_k1, nbasis_k1, lowf_k1, ekb_k1, occ_k1, wk_k1);

    EXPECT_EQ(2, ik_k1);
    EXPECT_EQ(19, nbands_k1);
    EXPECT_EQ(63, nbasis_k1);
    EXPECT_EQ(-1.0, wk_k1); // this is not read
    // kvec_c
    EXPECT_EQ(kvec_c_k1.x, kvec_c[1].x);
    EXPECT_EQ(kvec_c_k1.y, kvec_c[1].y);
    EXPECT_EQ(kvec_c_k1.z, kvec_c[1].z);
    // ekb
    for(int i = 0; i < nbands_k1; i++) { EXPECT_EQ(ekb_k1[i], ekb[i + nbands_k0]); }
    // occ
    for(int i = 0; i < nbands_k1; i++) { EXPECT_EQ(occ_k1[i], occ[i + nbands_k0]); }
    // lowf
    for(int i = 0; i < nbands_k1*nbasis_k1; i++) { EXPECT_EQ(lowf_k1[i], lowf[i + nbands_k0*nbasis_k0]); }
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