#include "module_io/read_wfc_lcao.h"

#include "gtest/gtest.h"

TEST(ReadWfcLcaoTest, ReadAbacusLowfComplex)
{

    // this test should only be executed on rank 0
#ifdef __MPI
    MPI_Barrier(MPI_COMM_WORLD);
    // printf("MPI environment detected, will use only rank 0\n");
    int iproc = 0;
    int nprocs = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &iproc);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    if (iproc != 0)
    {
        GTEST_SKIP();
    }
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
    std::string flowf = "./support/WFC_NAO_K1.txt";
    ModuleIO::read_abacus_lowf(flowf, ik, kvec_c, nbands, nbasis, lowf, ekb, occ, wk);
    EXPECT_EQ(1, ik);
    EXPECT_EQ(3, nbands);
    EXPECT_EQ(63, nbasis);
    EXPECT_EQ(1.0, wk); // this is overwritten with 1.0
    // kvec_c
    EXPECT_NEAR(-0.10540388, kvec_c.x, 1e-7);
    EXPECT_NEAR(-0.060854959, kvec_c.y, 1e-7);
    EXPECT_NEAR(-0.043030954, kvec_c.z, 1e-7);
    // ekb
    EXPECT_EQ(ekb.size(), nbands);
    EXPECT_NEAR(-6.03571945e-01, ekb[0], 1e-7);
    EXPECT_NEAR(-5.98035141e-01, ekb[1], 1e-7);
    EXPECT_NEAR(-5.98035141e-01, ekb[2], 1e-7);
    // occ
    EXPECT_EQ(occ.size(), nbands);
    EXPECT_NEAR(5.83090379e-03, occ[0], 1e-7);
    EXPECT_NEAR(5.83090379e-03, occ[1], 1e-7);
    EXPECT_NEAR(5.83090379e-03, occ[2], 1e-7);
    // lowf
    EXPECT_EQ(lowf.size(), nbands * nbasis);
    EXPECT_NEAR(lowf[0].real(), -6.71651157e-03, 1e-7);
    EXPECT_NEAR(lowf[0].imag(), 2.25946383e-02, 1e-7);
    EXPECT_NEAR(lowf[1].real(), 1.43180123e-03, 1e-7);
    EXPECT_NEAR(lowf[1].imag(), 1.46451488e-03, 1e-7);
    EXPECT_NEAR(lowf[2].real(), 2.31452033e-03, 1e-7);
    EXPECT_NEAR(lowf[2].imag(), -1.18949691e-03, 1e-7);
    EXPECT_NEAR(lowf[62].real(), 1.82648757e-03, 1e-7);
    EXPECT_NEAR(lowf[62].imag(), -2.11799886e-03, 1e-7);
    // test reuse, expect to overwrite the previous values
    flowf = "./support/WFC_NAO_K2.txt";
    ModuleIO::read_abacus_lowf(flowf, ik, kvec_c, nbands, nbasis, lowf, ekb, occ, wk);
    EXPECT_EQ(2, ik);
    EXPECT_EQ(3, nbands);
    EXPECT_EQ(63, nbasis);
    EXPECT_EQ(1.0, wk); // this is overwritten with 1.0
    // kvec_c
    EXPECT_NEAR(-0.070269254, kvec_c.x, 1e-7);
    EXPECT_NEAR(-0.081139946, kvec_c.y, 1e-7);
    EXPECT_NEAR(-0.057374606, kvec_c.z, 1e-7);
    // ekb
    EXPECT_EQ(ekb.size(), nbands);
    EXPECT_NEAR(-6.03667277e-01, ekb[0], 1e-7);
    EXPECT_NEAR(-5.97868276e-01, ekb[1], 1e-7);
    EXPECT_NEAR(-5.97662421e-01, ekb[2], 1e-7);
    // occ
    EXPECT_EQ(occ.size(), nbands);
    EXPECT_NEAR(5.83090379e-03, occ[0], 1e-7);
    EXPECT_NEAR(5.83090379e-03, occ[1], 1e-7);
    EXPECT_NEAR(5.83090379e-03, occ[2], 1e-7);
    // lowf
    EXPECT_EQ(lowf.size(), nbands * nbasis);
    EXPECT_NEAR(2.09933705e-02, lowf[0].real(), 1e-7);
    EXPECT_NEAR(2.20619371e-03, lowf[0].imag(), 1e-7);
    EXPECT_NEAR(1.52454416e-03, lowf[1].real(), 1e-7);
    EXPECT_NEAR(-3.54139105e-04, lowf[1].imag(), 1e-7);
    EXPECT_NEAR(1.31198443e-03, lowf[2].real(), 1e-7);
    EXPECT_NEAR(-2.44872538e-03, lowf[2].imag(), 1e-7);
    EXPECT_NEAR(lowf[62].real(), -1.15158489e-03, 1e-7);
    EXPECT_NEAR(lowf[62].imag(), -1.79940038e-03, 1e-7);
    // test reuse, the second time
    flowf = "./support/WFC_NAO_K3.txt";
    ModuleIO::read_abacus_lowf(flowf, ik, kvec_c, nbands, nbasis, lowf, ekb, occ, wk);
    EXPECT_EQ(3, ik);
    EXPECT_EQ(3, nbands);
    EXPECT_EQ(63, nbasis);
    EXPECT_EQ(1.0, wk); // this is overwritten with 1.0
    // kvec_c
    EXPECT_NEAR(-0.035134627, kvec_c.x, 1e-7);
    EXPECT_NEAR(-0.10142493, kvec_c.y, 1e-7);
    EXPECT_NEAR(-0.07171825, kvec_c.z, 1e-7);
    // ekb
    EXPECT_NEAR(-6.04664544e-01, ekb[0], 1e-7);
    EXPECT_NEAR(-5.97025474e-01, ekb[1], 1e-7);
    EXPECT_NEAR(-5.96870018e-01, ekb[2], 1e-7);
    // occ
    EXPECT_NEAR(5.83090379e-03, occ[0], 1e-7);
    EXPECT_NEAR(5.83090379e-03, occ[1], 1e-7);
    EXPECT_NEAR(5.83090379e-03, occ[2], 1e-7);
    // lowf
    EXPECT_NEAR(-6.44410072e-03, lowf[0].real(), 1e-7);
    EXPECT_NEAR(2.86108252e-03, lowf[0].imag(), 1e-7);
    EXPECT_NEAR(-5.81398415e-03, lowf[1].real(), 1e-7);
    EXPECT_NEAR(4.01495705e-03, lowf[1].imag(), 1e-7);
    EXPECT_NEAR(-2.32158666e-03, lowf[2].real(), 1e-7);
    EXPECT_NEAR(2.62541166e-03, lowf[2].imag(), 1e-7);
    EXPECT_NEAR(lowf[62].real(), 5.58964902e-04, 1e-7);
    EXPECT_NEAR(lowf[62].imag(), 5.21866389e-04, 1e-7);
}

TEST(ReadWfcLcaoTest, Cpzgemr2dUseTest)
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
    int iproc = 0;
    int nprocs = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &iproc);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    if (nprocs != 4)
    {
        if (iproc == 0)
        {
            printf("Please run this unittest with nprocs = 4.\n");
        }
        GTEST_SKIP();
    }
    // one get data on rank 0
    std::vector<std::complex<double>> lowf_glb;
    const int nbands = 10; // nrow_glb
    const int nbasis = 12; // ncol_glb
    if (iproc == 0)
    {
        printf("Run unittest ScatterLowfTest::ScatterLowfComplex with MPI env:\n");
        printf("Total number of processes: %d\n\n", nprocs);
        printf("Row-major processor grid is used.\n\n");
        printf("First test the \"column-major\" matrix, which means for columns their memory\n");
        printf("are consecutive. The matrix in form of (i, j), where i runs over [0, 9] and \n");
        printf("j runs over [0, 11]: (0,0), (1,0), (2,0), ..., (9,0), (0,1), ...\n");
        lowf_glb.resize(nbands * nbasis);
        for (int i = 0; i < nbands * nbasis; i++)
        {
            const int j = i / nbands, k = i % nbands;
            lowf_glb[i] = std::complex<double>(k, j);
        }
    }
    // the other ranks get data
    MPI_Barrier(MPI_COMM_WORLD);
    std::vector<std::complex<double>> lowf_loc;
    Parallel_2D para2d_test; // an alternative to ParaV. But to use paraV, the desc would be desc_wfc instead of desc in
                             // Parallel_2D
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

    // the following function can do the scattering-gathering automatically.
    // a double counterpart is pdgemr2d_, int counterpart is pigemr2d_
    // Those in C style are Cpzgemr2d, Cdgemr2d, Cigemr2d...
    Cpxgemr2d(nbands, nbasis, lowf_glb.data(), 1, 1, const_cast<int*>(para2d_glb.desc), lowf_loc.data(), 1, 1,
              const_cast<int*>(para2d_test.desc), para2d_glb.blacs_ctxt);
    // what will happen if impose a row-major processor grid onto a column-major matrix?
    // you can get correct results, expect each block is column-major.
    // Have a look:
    MPI_Barrier(MPI_COMM_WORLD);
    std::vector<int> sizes_loc = {4 * 4 * 3, 4 * 4 + 4 * 2, 4 * 4 * 2, 4 * 4};
    std::vector<std::vector<double>> reals = {{0, 1, 2, 3, 8, 9}, {0, 1, 2, 3, 8, 9}, {4, 5, 6, 7}, {4, 5, 6, 7}};
    std::vector<std::vector<double>> imags
        = {{0, 1, 2, 3, 8, 9, 10, 11}, {4, 5, 6, 7}, {0, 1, 2, 3, 8, 9, 10, 11}, {4, 5, 6, 7}};
    for (int i = 0; i < nprocs; i++)
    {
        if (iproc == i)
        {
            EXPECT_EQ(lowf_loc.size(), sizes_loc[i]);
            printf(">>> rank %d: \n", iproc);
            printf("First print scattered matrix in the way that ELEMENTS WITH CONSECUTIVE\n");
            printf("MEMORIES ARE SHOWN (only shown) IN THE SAME LINE:\n");
            for (int j = 0; j < lowf_loc.size(); j++)
            {
                printf("(%2.0f,%2.0f)", lowf_loc[j].real(), lowf_loc[j].imag());
                if ((j + 1) % para2d_test.nrow == 0)
                {
                    printf("\n");
                }
                const int k = j % para2d_test.nrow;
                EXPECT_NEAR(lowf_loc[j].real(), reals[i][k], 1e-7);
                const int l = j / para2d_test.nrow;
                EXPECT_NEAR(lowf_loc[j].imag(), imags[i][l], 1e-7);
            }
            printf("Or INDEX IT to show like \"row-major\":\n");
            // (i, j) -> (i', j') with i = j' and j = i'
            // x = i*ncol + j, x' = i'*ncol' + j' with ncol' = nrow and nrow' = ncol
            // i = x/ncol, j = x%ncol, x' = j*nrow + i = x%ncol*nrow + x/ncol
            for (int j = 0; j < lowf_loc.size(); j++)
            {
                const int x = j % para2d_test.ncol * para2d_test.nrow + j / para2d_test.ncol;
                printf("(%2.0f,%2.0f)", lowf_loc[x].real(), lowf_loc[x].imag());
                if ((j + 1) % para2d_test.ncol == 0)
                {
                    printf("\n");
                }
            }
            printf("\n");
            usleep(10000);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    // test the other way around, the row-major matrix
    if (iproc == 0)
    {
        printf("Now test the \"row-major\" matrix, which means for rows their memory\n");
        printf("are consecutive. The matrix in form of (i, j), where i runs over [0, 9] and \n");
        printf("j runs over [0, 11]: (0,0), (0,1), (0,2), ..., (0,11), (1,0), ...\n");
        for (int i = 0; i < nbands * nbasis; i++)
        {
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
    Cpxgemr2d(nbasis, nbands, lowf_glb.data(), 1, 1, const_cast<int*>(para2d_glb_prime.desc), lowf_loc.data(), 1, 1,
              const_cast<int*>(para2d_test_prime.desc), para2d_glb_prime.blacs_ctxt);
    MPI_Barrier(MPI_COMM_WORLD);
    sizes_loc = {4 * 4 * 3, 4 * 4 * 2, 4 * 4 + 4 * 2, 4 * 4};
    reals = {{0, 1, 2, 3, 8, 9}, {4, 5, 6, 7}, {0, 1, 2, 3, 8, 9}, {4, 5, 6, 7}};
    imags = {{0, 1, 2, 3, 8, 9, 10, 11}, {0, 1, 2, 3, 8, 9, 10, 11}, {4, 5, 6, 7}, {4, 5, 6, 7}};
    for (int i = 0; i < nprocs; i++)
    {
        if (iproc == i)
        {
            EXPECT_EQ(lowf_loc.size(), sizes_loc[i]);
            printf(">>> rank %d: \n", iproc);
            for (int j = 0; j < lowf_loc.size(); j++)
            {
                printf("(%2.0f,%2.0f)", lowf_loc[j].real(), lowf_loc[j].imag());
                if ((j + 1) % para2d_test_prime.nrow == 0)
                {
                    printf("\n");
                }
                const int k = j / para2d_test_prime.nrow;
                EXPECT_NEAR(lowf_loc[j].real(), reals[i][k], 1e-7);
                const int l = j % para2d_test_prime.nrow;
                EXPECT_NEAR(lowf_loc[j].imag(), imags[i][l], 1e-7);
            }
            printf("\n");
            usleep(10000);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if (iproc == 0)
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
    int iproc = 0;
    int nprocs = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &iproc);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    if (iproc != 0)
    {
        GTEST_SKIP();
    }
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
    const std::string flowf = "./support/WFC_NAO_GAMMA1.txt";
    ModuleIO::read_abacus_lowf(flowf, ik, kvec_c, nbands, nbasis, lowf, ekb, occ, wk);
    EXPECT_EQ(0, ik);
    EXPECT_EQ(3, nbands);
    EXPECT_EQ(31, nbasis);
    EXPECT_EQ(1.0, wk);
    // kvec_c, gamma point, 0, 0, 0
    EXPECT_NEAR(0.0, kvec_c.x, 1e-7);
    EXPECT_NEAR(0.0, kvec_c.y, 1e-7);
    EXPECT_NEAR(0.0, kvec_c.z, 1e-7);
    // ekb
    EXPECT_EQ(ekb.size(), nbands);
    EXPECT_NEAR(-1.22787155e+00, ekb[0], 1e-7);
    EXPECT_NEAR(-3.10595658e-01, ekb[1], 1e-7);
    EXPECT_NEAR(-3.00546690e-01, ekb[2], 1e-7);
    // occ
    EXPECT_EQ(occ.size(), nbands);
    EXPECT_NEAR(2.00000000e+00, occ[0], 1e-7);
    EXPECT_NEAR(2.00000000e+00, occ[1], 1e-7);
    EXPECT_NEAR(2.00000000e+00, occ[2], 1e-7);
    // lowf
    EXPECT_EQ(lowf.size(), nbands * nbasis);
    EXPECT_NEAR(-1.51728369e-02, lowf[0], 1e-7);
    EXPECT_NEAR(-2.07808444e-03, lowf[1], 1e-7);
    EXPECT_NEAR(1.21298954e-17, lowf[2], 1e-7);
    EXPECT_NEAR(-5.44883791e-09, lowf[30], 1e-7);
}

TEST(ReadWfcLcaoTest, Cpdgemr2dUseTest)
{
    // you can find more information in unittest Pzgemr2dUseTest, present test
    // works identically to the previous one, but with real numbers.
#ifdef __MPI
    // this test should be run on all ranks
    int iproc = 0;
    int nprocs = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &iproc);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    if (nprocs != 4)
    {
        if (iproc == 0)
        {
            printf("Please run this unittest with nprocs = 4.\n");
        }
        GTEST_SKIP();
    }
    std::vector<double> lowf_glb;
    const int nbands = 10;
    const int nbasis = 12;
    // still, the expected matrix is organized as row-memory-consecutive but here
    // just make the matrix column-memory-consecutive.
    // x = i*ncol + j, x' = j*nrow + i
    // i = x/ncol, j = x%ncol, x' = j*nrow + i = x%ncol*nrow + x/ncol
    if (iproc == 0)
    {
        lowf_glb.resize(nbands * nbasis);
        for (int i = 0; i < nbands * nbasis; i++)
        {
            lowf_glb[i] = i % nbasis * nbands + i / nbasis;
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    std::vector<double> lowf_loc;
    Parallel_2D para2d_test;
    para2d_test.init(nbasis, nbands, 4, MPI_COMM_WORLD);
    Parallel_2D para2d_glb;
    para2d_glb.init(nbasis, nbands, std::max(nbands, nbasis), MPI_COMM_WORLD);
    lowf_loc.resize(para2d_test.nrow * para2d_test.ncol);
    Cpxgemr2d(nbasis, nbands, lowf_glb.data(), 1, 1, const_cast<int*>(para2d_glb.desc), lowf_loc.data(), 1, 1,
              const_cast<int*>(para2d_test.desc), para2d_glb.blacs_ctxt);
    MPI_Barrier(MPI_COMM_WORLD);
    std::vector<int> sizes_loc = {4 * 4 * 3, 4 * 4 * 2, 4 * 4 + 4 * 2, 4 * 4};
    for (int i = 0; i < nprocs; i++)
    {
        if (iproc == i)
        {
            EXPECT_EQ(lowf_loc.size(), sizes_loc[i]);
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    // test the other way around, the row-major matrix
    if (iproc == 0)
    {
        for (int i = 0; i < nbands * nbasis; i++)
        {
            lowf_glb[i] = i;
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    Parallel_2D para2d_test_prime;
    para2d_test_prime.init(nbands, nbasis, 4, MPI_COMM_WORLD);
    Parallel_2D para2d_glb_prime;
    para2d_glb_prime.init(nbands, nbasis, std::max(nbands, nbasis), MPI_COMM_WORLD);
    lowf_loc.resize(para2d_test_prime.nrow * para2d_test_prime.ncol);
    Cpxgemr2d(nbands, nbasis, lowf_glb.data(), 1, 1, const_cast<int*>(para2d_glb_prime.desc), lowf_loc.data(), 1, 1,
              const_cast<int*>(para2d_test_prime.desc), para2d_glb_prime.blacs_ctxt);
    MPI_Barrier(MPI_COMM_WORLD);
    sizes_loc = {4 * 4 * 3, 4 * 4 + 4 * 2, 4 * 4 * 2, 4 * 4};
    for (int i = 0; i < nprocs; i++)
    {
        if (iproc == i)
        {
            EXPECT_EQ(lowf_loc.size(), sizes_loc[i]);
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
#endif
}

TEST(ReadWfcLcaoTest, RestartFromFileParallel)
{
#ifdef __MPI
    MPI_Barrier(MPI_COMM_WORLD);
    // printf("MPI environment detected, will use only rank 0\n");
    int iproc = 0;
    int nprocs = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &iproc);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    if (nprocs != 4)
    {
        if (iproc == 0)
        {
            printf("Please run this unittest with nprocs = 4.\n");
        }
        GTEST_SKIP();
    }
    Parallel_2D para2d_test;
    // this time the nbands and nbasis are given as priori, from reading of INPUT and of orb files.
    para2d_test.init(63, 3, 2, MPI_COMM_WORLD);
    // therefore it is a nrows x ncols matrix, where nr = 3, nc = 63
    // scatter it with block size 2 x 2, on 4-processors grid. From testcase ReadWfcLcaoTest::Cpzgemr2dTest
    // it is known the processor grid will be transposed to column-major, thus
    // the global matrix will be split to 2 x 2 blocks, one 2 x 1 block, 1 x 2 blocks and 1 x 1 block
    // number of elements:
    // rank0: ceil(63/2)*2 = 64
    // rank1: ceil(63/2)*1 = 32
    // rank2: floor(63/2)*2 = 62
    // rank3: floor(63/2)*1 = 31
    // total size check
    if (iproc == 0)
    {
        EXPECT_EQ(64, para2d_test.nrow * para2d_test.ncol);
    }
    if (iproc == 1)
    {
        EXPECT_EQ(32, para2d_test.nrow * para2d_test.ncol);
    }
    if (iproc == 2)
    {
        EXPECT_EQ(62, para2d_test.nrow * para2d_test.ncol);
    }
    if (iproc == 3)
    {
        EXPECT_EQ(31, para2d_test.nrow * para2d_test.ncol);
    }
    // nrows check
    if (iproc == 0)
    {
        EXPECT_EQ(2, para2d_test.ncol);
    }
    if (iproc == 1)
    {
        EXPECT_EQ(1, para2d_test.ncol);
    }
    if (iproc == 2)
    {
        EXPECT_EQ(2, para2d_test.ncol);
    }
    if (iproc == 3)
    {
        EXPECT_EQ(1, para2d_test.ncol);
    }

    const std::string out_dir = "./support";
    const int nks = 4;
    int nbands = -1;
    int nbasis = -1;
    std::vector<std::complex<double>> lowf;
    std::vector<double> ekb;
    std::vector<double> occ;
    std::vector<ModuleBase::Vector3<double>> kvec_c;
    std::vector<double> wk;
    ModuleIO::restart_from_file(out_dir, para2d_test, nks, nbands, nbasis, lowf, ekb, occ, kvec_c, wk);
    EXPECT_EQ(3, nbands);
    EXPECT_EQ(63, nbasis);
    // ekb, occ, kvec_c, wk will have size irrelevant to scatter
    EXPECT_EQ(nks * nbands, ekb.size());
    EXPECT_EQ(nks * nbands, occ.size());
    EXPECT_EQ(nks, kvec_c.size());
    EXPECT_EQ(nks, wk.size());
    // value test
    const std::vector<double> ekb_ref
        = {-6.03571945e-01, -5.98035141e-01, -5.98035141e-01, -6.03667277e-01, -5.97868276e-01, -5.97662421e-01,
           -6.04664544e-01, -5.97025474e-01, -5.96870018e-01, -6.05615293e-01, -5.96302906e-01, -5.96302906e-01};
    const std::vector<double> occ_ref
        = {5.83090379e-03, 5.83090379e-03, 5.83090379e-03, 5.83090379e-03, 5.83090379e-03, 5.83090379e-03,
           5.83090379e-03, 5.83090379e-03, 5.83090379e-03, 5.83090379e-03, 5.83090379e-03, 5.83090379e-03};
    const std::vector<ModuleBase::Vector3<double>> kvec_c_ref
        = {ModuleBase::Vector3<double>(-0.10540388, -0.060854959, -0.043030954),
           ModuleBase::Vector3<double>(-0.070269254, -0.081139946, -0.057374606),
           ModuleBase::Vector3<double>(-0.035134627, -0.10142493, -0.07171825),
           ModuleBase::Vector3<double>(0.00000000, -0.12170991, -0.086061909)};
    const std::vector<double> wk_ref = {1.0, 1.0, 1.0, 1.0};
    for (int i = 0; i < nprocs; i++)
    {
        if (iproc == i)
        {
            // ekb
            for (int j = 0; j < nks * nbands; j++)
            {
                EXPECT_NEAR(ekb_ref[j], ekb[j], 1e-7);
            }
            // occ
            for (int j = 0; j < nks * nbands; j++)
            {
                EXPECT_NEAR(occ_ref[j], occ[j], 1e-7);
            }
            // kvec_c
            for (int j = 0; j < nks; j++)
            {
                EXPECT_NEAR(kvec_c_ref[j].x, kvec_c[j].x, 1e-7);
                EXPECT_NEAR(kvec_c_ref[j].y, kvec_c[j].y, 1e-7);
                EXPECT_NEAR(kvec_c_ref[j].z, kvec_c[j].z, 1e-7);
            }
            // wk
            for (int j = 0; j < nks; j++)
            {
                EXPECT_NEAR(wk_ref[j], wk[j], 1e-7);
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    // lowf will be the result of scattering
    // rank0: 64, rank1: 32, rank2: 62, rank3: 31, each should multiply the number of k-points
    if (iproc == 0)
    {
        EXPECT_EQ(nks * 64, lowf.size());
    }
    if (iproc == 1)
    {
        EXPECT_EQ(nks * 32, lowf.size());
    }
    if (iproc == 2)
    {
        EXPECT_EQ(nks * 62, lowf.size());
    }
    if (iproc == 3)
    {
        EXPECT_EQ(nks * 31, lowf.size());
    }
    // value test on lowf
    // rank0
    if (iproc == 0)
    {
        EXPECT_NEAR(lowf[0].real(), -6.71651157e-03, 1e-7);
        EXPECT_NEAR(lowf[0].imag(), 2.25946383e-02, 1e-7);
        EXPECT_NEAR(lowf[1].real(), 1.43180123e-03, 1e-7);
        EXPECT_NEAR(lowf[1].imag(), 1.46451488e-03, 1e-7);
        EXPECT_NEAR(lowf[2].real(), -9.45760546e-03, 1e-7);
        EXPECT_NEAR(lowf[2].imag(), 1.29911511e-02, 1e-7);
        EXPECT_NEAR(lowf[3].real(), -5.46035106e-03, 1e-7);
        EXPECT_NEAR(lowf[3].imag(), 7.50044462e-03, 1e-7);
        EXPECT_NEAR(lowf[lowf.size() - 1].real(), -1.39799597e-03, 1e-7);
        EXPECT_NEAR(lowf[lowf.size() - 1].imag(), -1.68192980e-03, 1e-7);
    }
    else if (iproc == 1)
    {
        EXPECT_NEAR(lowf[0].real(), 9.86470874e-13, 1e-7);
        EXPECT_NEAR(lowf[0].imag(), -5.95387122e-12, 1e-7);
        EXPECT_NEAR(lowf[1].real(), 4.82573453e-13, 1e-7);
        EXPECT_NEAR(lowf[1].imag(), 5.54264959e-12, 1e-7);
        EXPECT_NEAR(lowf[2].real(), 5.20920946e-04, 1e-7);
        EXPECT_NEAR(lowf[2].imag(), -2.33076310e-03, 1e-7);
        EXPECT_NEAR(lowf[3].real(), -8.35155442e-04, 1e-7);
        EXPECT_NEAR(lowf[3].imag(), 3.75083842e-03, 1e-7);
        EXPECT_NEAR(lowf[lowf.size() - 1].real(), -6.18118243e-05, 1e-7);
        EXPECT_NEAR(lowf[lowf.size() - 1].imag(), -7.43658388e-05, 1e-7);
    }
    else if (iproc == 2)
    {
        EXPECT_NEAR(lowf[0].real(), 2.31452033e-03, 1e-7);
        EXPECT_NEAR(lowf[0].imag(), -1.18949691e-03, 1e-7);
        EXPECT_NEAR(lowf[1].real(), 3.86105126e-03, 1e-7);
        EXPECT_NEAR(lowf[1].imag(), -5.30361525e-03, 1e-7);
        EXPECT_NEAR(lowf[2].real(), 1.07440727e-03, 1e-7);
        EXPECT_NEAR(lowf[2].imag(), -1.52230629e-03, 1e-7);
        EXPECT_NEAR(lowf[3].real(), -2.63174959e-03, 1e-7);
        EXPECT_NEAR(lowf[3].imag(), 3.72887365e-03, 1e-7);
        EXPECT_NEAR(lowf[lowf.size() - 1].real(), 1.19038759e-04, 1e-7);
        EXPECT_NEAR(lowf[lowf.size() - 1].imag(), 1.17824924e-04, 1e-7);
    }
    else if (iproc == 3)
    {
        EXPECT_NEAR(lowf[0].real(), 3.66087151e-13, 1e-7);
        EXPECT_NEAR(lowf[0].imag(), 1.96386245e-13, 1e-7);
        EXPECT_NEAR(lowf[1].real(), 9.49023673e-05, 1e-7);
        EXPECT_NEAR(lowf[1].imag(), -4.04693771e-04, 1e-7);
        EXPECT_NEAR(lowf[2].real(), 1.33229060e-04, 1e-7);
        EXPECT_NEAR(lowf[2].imag(), 9.69176971e-04, 1e-7);
        EXPECT_NEAR(lowf[3].real(), 8.23664081e-04, 1e-7);
        EXPECT_NEAR(lowf[3].imag(), 5.56014508e-03, 1e-7);
        EXPECT_NEAR(lowf[lowf.size() - 1].real(), -2.69229582e-03, 1e-7);
        EXPECT_NEAR(lowf[lowf.size() - 1].imag(), -2.66484241e-03, 1e-7);
    }
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
    int iproc = 0;
    int nprocs = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &iproc);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    if (iproc != 0)
    {
        GTEST_SKIP();
    }
#endif
    const int nks = 4;
    int nbands = -1;
    int nbasis = -1;
    std::vector<std::complex<double>> lowf;
    std::vector<double> ekb;
    std::vector<double> occ;
    std::vector<ModuleBase::Vector3<double>> kvec_c;
    std::vector<double> wk;
    const std::string out_dir = "./support";
    ModuleIO::restart_from_file(out_dir, nks, nbands, nbasis, lowf, ekb, occ, kvec_c, wk);
    EXPECT_EQ(3, nbands);
    EXPECT_EQ(63, nbasis);
    EXPECT_EQ(nks * nbands * nbasis, lowf.size());
    EXPECT_EQ(nks * nbands, ekb.size());
    EXPECT_EQ(nks * nbands, occ.size());
    EXPECT_EQ(nks, kvec_c.size());
    EXPECT_EQ(nks, wk.size());

    // test the first k-point
    int nbands_k0 = -1;
    int nbasis_k0 = -1;
    int ik_k0;
    std::vector<std::complex<double>> lowf_k0;
    std::vector<double> ekb_k0;
    std::vector<double> occ_k0;
    ModuleBase::Vector3<double> kvec_c_k0;
    double wk_k0 = -1.0;
    ModuleIO::read_abacus_lowf("./support/WFC_NAO_K1.txt", ik_k0, kvec_c_k0, nbands_k0, nbasis_k0, lowf_k0, ekb_k0,
                               occ_k0, wk_k0);

    EXPECT_EQ(1, ik_k0);
    EXPECT_EQ(3, nbands_k0);
    EXPECT_EQ(63, nbasis_k0);
    EXPECT_EQ(1.0, wk_k0); // this is overwritten with 1.0
    // kvec_c
    EXPECT_NEAR(kvec_c_k0.x, kvec_c[0].x, 1e-7);
    EXPECT_NEAR(kvec_c_k0.y, kvec_c[0].y, 1e-7);
    EXPECT_NEAR(kvec_c_k0.z, kvec_c[0].z, 1e-7);
    // ekb
    for (int i = 0; i < nbands_k0; i++)
    {
        EXPECT_NEAR(ekb_k0[i], ekb[i], 1e-7);
    }
    // occ
    for (int i = 0; i < nbands_k0; i++)
    {
        EXPECT_NEAR(occ_k0[i], occ[i], 1e-7);
    }
    // lowf
    for (int i = 0; i < nbands_k0 * nbasis_k0; i++)
    {
        EXPECT_NEAR(lowf_k0[i].real(), lowf[i].real(), 1e-7);
        EXPECT_NEAR(lowf_k0[i].imag(), lowf[i].imag(), 1e-7);
    }

    // test the second k-point
    int nbands_k1 = -1;
    int nbasis_k1 = -1;
    int ik_k1;
    std::vector<std::complex<double>> lowf_k1;
    std::vector<double> ekb_k1;
    std::vector<double> occ_k1;
    ModuleBase::Vector3<double> kvec_c_k1;
    double wk_k1 = -1.0;
    ModuleIO::read_abacus_lowf("./support/WFC_NAO_K2.txt", ik_k1, kvec_c_k1, nbands_k1, nbasis_k1, lowf_k1, ekb_k1,
                               occ_k1, wk_k1);

    EXPECT_EQ(2, ik_k1);
    EXPECT_EQ(3, nbands_k1);
    EXPECT_EQ(63, nbasis_k1);
    EXPECT_EQ(1.0, wk_k1); // this is overwritten with 1.0
    // kvec_c
    EXPECT_NEAR(kvec_c_k1.x, kvec_c[1].x, 1e-7);
    EXPECT_NEAR(kvec_c_k1.y, kvec_c[1].y, 1e-7);
    EXPECT_NEAR(kvec_c_k1.z, kvec_c[1].z, 1e-7);
    // ekb
    for (int i = 0; i < nbands_k1; i++)
    {
        EXPECT_NEAR(ekb_k1[i], ekb[i + nbands_k0], 1e-7);
    }
    // occ
    for (int i = 0; i < nbands_k1; i++)
    {
        EXPECT_NEAR(occ_k1[i], occ[i + nbands_k0], 1e-7);
    }
    // lowf
    for (int i = 0; i < nbands_k1 * nbasis_k1; i++)
    {
        EXPECT_NEAR(lowf_k1[i].real(), lowf[i + nbands_k0 * nbasis_k0].real(), 1e-7);
        EXPECT_NEAR(lowf_k1[i].imag(), lowf[i + nbands_k0 * nbasis_k0].imag(), 1e-7);
    }
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    // print cwd
    char cwd[1024];
    if (getcwd(cwd, sizeof(cwd)) != nullptr)
    {
        printf("Current working dir: %s\n", cwd);
    }
    else
    {
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