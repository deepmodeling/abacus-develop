#include "../kernels/math_kernel_op.h"
#include "../para_gemm.h"

#include <gtest/gtest.h>
#include <iostream>
#include <mpi.h>
#include <vector>

void random_data(std::vector<double>& A_global,
                 std::vector<double>& B_global,
                 std::vector<double>& Cref_global,
                 std::vector<double>& C_global,
                 double& alpha,
                 double& beta)
{
    for (auto& val: A_global)
    {
        val = std::rand() / (RAND_MAX + 1.0);
    }
    for (auto& val: B_global)
    {
        val = std::rand() / (RAND_MAX + 1.0);
    }
    for (auto& val: Cref_global)
    {
        val = std::rand() / (RAND_MAX + 1.0);
    }
    C_global = Cref_global;

    alpha = std::rand() / (RAND_MAX + 1.0);
    beta = std::rand() / (RAND_MAX + 1.0);
}
void random_data(std::vector<std::complex<double>>& A_global,
                 std::vector<std::complex<double>>& B_global,
                 std::vector<std::complex<double>>& Cref_global,
                 std::vector<std::complex<double>>& C_global,
                 std::complex<double>& alpha,
                 std::complex<double>& beta)
{
    for (auto& val: A_global)
    {
        val = std::complex<double>(std::rand() / (RAND_MAX + 1.0), std::rand() / (RAND_MAX + 1.0));
    }
    for (auto& val: B_global)
    {
        val = std::complex<double>(std::rand() / (RAND_MAX + 1.0), std::rand() / (RAND_MAX + 1.0));
    }
    for (auto& val: Cref_global)
    {
        val = std::complex<double>(std::rand() / (RAND_MAX + 1.0), std::rand() / (RAND_MAX + 1.0));
    }
    C_global = Cref_global;

    alpha = std::complex<double>(std::rand() / (RAND_MAX + 1.0), std::rand() / (RAND_MAX + 1.0));
    beta = std::complex<double>(std::rand() / (RAND_MAX + 1.0), std::rand() / (RAND_MAX + 1.0));
}
double get_double(std::complex<double>& val)
{
    return val.real() + val.imag();
}
double get_double(double& val)
{
    return val;
}

void scatterv_data(const double* sendbuf,
                   const int* sendcounts,
                   const int* displs,
                   double* recvbuf,
                   const int recvcount,
                   MPI_Comm comm)
{
    MPI_Scatterv(sendbuf, sendcounts, displs, MPI_DOUBLE, recvbuf, recvcount, MPI_DOUBLE, 0, comm);
}
void scatterv_data(const std::complex<double>* sendbuf,
                   const int* sendcounts,
                   const int* displs,
                   std::complex<double>* recvbuf,
                   const int recvcount,
                   MPI_Comm comm)
{
    MPI_Scatterv(sendbuf, sendcounts, displs, MPI_DOUBLE_COMPLEX, recvbuf, recvcount, MPI_DOUBLE_COMPLEX, 0, comm);
}
template <typename T>
class PgemmTest : public ::testing::Test
{
  protected:
    void SetUp() override
    {
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    }
    void TearDown() override
    {
        MPI_Comm_free(&col_world);
        MPI_Comm_free(&row_world);
    }

  public:
    void decide_ngroup(const int& willing_ncolgroup, const int& willing_nrowgroup)
    {
        ncolgroup = willing_ncolgroup;
        nrowgroup = willing_nrowgroup;
        if (nproc % (ncolgroup * nrowgroup) != 0)
        {
            ncolgroup = nproc;
            nrowgroup = 1;
        }
        else
        {
            nrowgroup = nproc / ncolgroup;
        }

        MPI_Comm_split(MPI_COMM_WORLD, rank % nrowgroup, rank / nrowgroup, &col_world);
        MPI_Comm_split(MPI_COMM_WORLD, rank / nrowgroup, rank % nrowgroup, &row_world);
        MPI_Comm_rank(col_world, &rank_col);
        MPI_Comm_rank(row_world, &rank_row);
        MPI_Comm_size(col_world, &nproc_col);
        MPI_Comm_size(row_world, &nproc_row);
    }
    void randomize_initialization()
    {
        random_data(A_global, B_global, Cref_global, C_global, alpha, beta);
    }

    void prepare(const int& ncolA_global,
                 const int& ncolB_global,
                 const int& nrow_global,
                 const int& LDA_global,
                 const int& LDB_global,
                 const int& LDC_global)
    {
        A_global = std::vector<T>(LDA_global * ncolA_global, 0.0);
        B_global = std::vector<T>(LDB_global * ncolB_global, 0.0);
        C_global = std::vector<T>(LDC_global * ncolB_global, 0.0);
        Cref_global = std::vector<T>(LDC_global * ncolB_global, 0.0);
        if (rank == 0)
        {

            this->randomize_initialization();
            const base_device::DEVICE_CPU* ctx = {};
            char transC = 'C';
            char transN = 'N';
            ModuleBase::gemm_op<T, base_device::DEVICE_CPU>()(ctx,
                                                              transC,
                                                              transN,
                                                              ncolA_global,
                                                              ncolB_global,
                                                              nrow_global,
                                                              &alpha,
                                                              A_global.data(),
                                                              LDA_global,
                                                              B_global.data(),
                                                              LDB_global,
                                                              &beta,
                                                              Cref_global.data(),
                                                              LDC_global);
        }

        if (std::is_same<T, double>::value)
        {
            MPI_Bcast(A_global.data(), A_global.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
            MPI_Bcast(B_global.data(), B_global.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
            MPI_Bcast(C_global.data(), C_global.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
            MPI_Bcast(Cref_global.data(), Cref_global.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
            MPI_Bcast(&alpha, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            MPI_Bcast(&beta, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        }
        else if (std::is_same<T, std::complex<double>>::value)
        {
            MPI_Bcast(A_global.data(), A_global.size(), MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
            MPI_Bcast(B_global.data(), B_global.size(), MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
            MPI_Bcast(C_global.data(), C_global.size(), MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
            MPI_Bcast(Cref_global.data(), Cref_global.size(), MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
            MPI_Bcast(&alpha, 1, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
            MPI_Bcast(&beta, 1, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
        }

        // Broadcast A_global and B_global to all ranks
        getncol_and_row(ncolA_global, ncolB_global, nrow_global);
        LDA = nrow + 1;
        LDB = nrow + 2;

        A_local = std::vector<T>(LDA * ncolA, 0.0);
        B_local = std::vector<T>(LDB * ncolB, 0.0);

        scatter_matrix(ncolA_global, ncolB_global, nrow_global, LDA_global, LDB_global);
    }

    void getncol_and_row(const int& ncolA_global, const int& ncolB_global, const int& nrow_global)
    {
        ncolA = ncolA_global / ncolgroup;
        if (ncolA_global % ncolgroup > rank_col)
        {
            ncolA += 1;
        }
        ncolB = ncolB_global / ncolgroup;
        if (ncolB_global % ncolgroup > rank_col)
        {
            ncolB += 1;
        }

        nrow = nrow_global / nrowgroup;
        if (nrow_global % nrowgroup > rank_row)
        {
            nrow += 1;
        }

        ncolA_ip.resize(nproc_col);
        ncolB_ip.resize(nproc_col);
        nrow_ip.resize(nproc_row);
        MPI_Allgather(&ncolA, 1, MPI_INT, ncolA_ip.data(), 1, MPI_INT, col_world);
        MPI_Allgather(&ncolB, 1, MPI_INT, ncolB_ip.data(), 1, MPI_INT, col_world);
        if (row_world != MPI_COMM_NULL)
        {
            MPI_Allgather(&nrow, 1, MPI_INT, nrow_ip.data(), 1, MPI_INT, row_world);
        }
    }

    void scatter_matrix(const int& ncolA_global,
                        const int& ncolB_global,
                        const int& nrow_global,
                        const int& LDA_global,
                        const int& LDB_global)
    {
        std::vector<T> A_semiglobal(ncolA * LDA_global, 0.0);
        std::vector<T> B_semiglobal(ncolB * LDB_global, 0.0);

        // Scatter A_global and B_global to A_semiglobal and B_semiglobal
        std::vector<int> sendcounts(nproc_col, 0);
        std::vector<int> displs(nproc_col, 0);
        for (int i = 0; i < nproc_col; i++)
        {
            sendcounts[i] = ncolA_ip[i] * LDA_global;
        }
        displs[0] = 0;
        for (int i = 1; i < nproc_col; i++)
        {
            displs[i] = displs[i - 1] + sendcounts[i - 1];
        }
        scatterv_data(A_global.data(),
                      sendcounts.data(),
                      displs.data(),
                      A_semiglobal.data(),
                      ncolA * LDA_global,
                      col_world);

        for (int i = 0; i < nproc_col; i++)
        {
            sendcounts[i] = ncolB_ip[i] * LDB_global;
        }
        displs[0] = 0;
        for (int i = 1; i < nproc_col; i++)
        {
            displs[i] = displs[i - 1] + sendcounts[i - 1];
        }
        scatterv_data(B_global.data(),
                      sendcounts.data(),
                      displs.data(),
                      B_semiglobal.data(),
                      ncolB * LDB_global,
                      col_world);

        // Scatter A_semiglobal and B_semiglobal to A_local and B_local
        sendcounts.resize(nproc_row, 0);
        displs.resize(nproc_row, 0);
        for (int i = 0; i < nproc_row; i++)
        {
            sendcounts[i] = nrow_ip[i];
        }
        displs[0] = 0;
        for (int i = 1; i < nproc_row; i++)
        {
            displs[i] = displs[i - 1] + sendcounts[i - 1];
        }
        for (int i = 0; i < ncolA; i++)
        {
            scatterv_data(A_semiglobal.data() + i * LDA_global,
                          sendcounts.data(),
                          displs.data(),
                          A_local.data() + i * LDA,
                          nrow,
                          row_world);
        }

        for (int i = 0; i < ncolB; i++)
        {
            scatterv_data(B_semiglobal.data() + i * LDB_global,
                          sendcounts.data(),
                          displs.data(),
                          B_local.data() + i * LDB,
                          nrow,
                          row_world);
        }
    }

    void compare_result(const int& nrowC_global, const int& ncolC_global, const int& LDC_global)
    {
        for (int i = 0; i < ncolC_global; i++)
        {
            for (int j = 0; j < nrowC_global; j++)
            {
                EXPECT_NEAR(get_double(Cref_global[i * LDC_global + j]),
                            get_double(C_global[i * LDC_global + j]),
                            1e-10);
            }
        }
    }

    int rank = 0, nproc = 0;
    T alpha = 0, beta = 0;
    std::vector<T> A_global, B_global, Cref_global, C_global;
    std::vector<T> A_local, B_local;
    int ncolA = 0, ncolB = 0, nrow = 0, LDA = 0, LDB = 0;
    int ncolgroup = 1, nrowgroup = 1;
    int rank_col = 0, rank_row = 0;
    int nproc_col = 0, nproc_row = 0;
    ModuleBase::PGemmCN<T> pgemm;
    MPI_Comm col_world;
    MPI_Comm row_world;
    std::vector<int> ncolA_ip, ncolB_ip, nrow_ip;
};

typedef ::testing::Types<double, std::complex<double>> MyTypes;

TYPED_TEST_SUITE(PgemmTest, MyTypes);

TYPED_TEST(PgemmTest, even_case)
{
    const int ncolA_global = 16, ncolB_global = 8, nrow_global = 12;
    const int LDA_global = 17, LDB_global = 18, LDC_global = 19;

    this->decide_ngroup(2, 2);
    this->prepare(ncolA_global, ncolB_global, nrow_global, LDA_global, LDB_global, LDC_global);

    this->pgemm.set_dimension(this->col_world,
                              this->row_world,
                              this->ncolA,
                              this->LDA,
                              this->ncolB,
                              this->LDB,
                              this->nrow,
                              LDC_global);
    this->pgemm.multiply(this->alpha, this->A_local.data(), this->B_local.data(), this->beta, this->C_global.data());

    this->compare_result(ncolA_global, ncolB_global, LDC_global);
}

TYPED_TEST(PgemmTest, odd_case)
{
    const int ncolA_global = 17, ncolB_global = 7, nrow_global = 13;
    const int LDA_global = 17, LDB_global = 18, LDC_global = 19;

    this->decide_ngroup(2, 2);
    this->prepare(ncolA_global, ncolB_global, nrow_global, LDA_global, LDB_global, LDC_global);

    this->pgemm.set_dimension(this->col_world,
                              this->row_world,
                              this->ncolA,
                              this->LDA,
                              this->ncolB,
                              this->LDB,
                              this->nrow,
                              LDC_global);
    this->pgemm.multiply(this->alpha, this->A_local.data(), this->B_local.data(), this->beta, this->C_global.data());

    this->compare_result(ncolA_global, ncolB_global, LDC_global);
}

TYPED_TEST(PgemmTest, odd_case_not_gather)
{
    const int ncolA_global = 17, ncolB_global = 7, nrow_global = 13;
    const int LDA_global = 17, LDB_global = 18, LDC_global = 19;

    this->decide_ngroup(2, 2);
    this->prepare(ncolA_global, ncolB_global, nrow_global, LDA_global, LDB_global, LDC_global);
    std::vector<int> colB_loc(this->nproc_col);
    MPI_Allgather(&this->ncolB, 1, MPI_INT, colB_loc.data(), 1, MPI_INT, this->col_world);
    std::vector<int> displs(this->nproc_col);
    displs[0] = 0;
    for (int i = 1; i < this->nproc_col; i++)
    {
        displs[i] = (displs[i - 1] + colB_loc[i - 1]) * LDC_global;
    }
    int start = displs[this->rank_col];

    this->pgemm.set_dimension(this->col_world,
                              this->row_world,
                              this->ncolA,
                              this->LDA,
                              this->ncolB,
                              this->LDB,
                              this->nrow,
                              LDC_global,
                              false);
    this->pgemm.multiply(this->alpha, this->A_local.data(), this->B_local.data(), this->beta, this->C_global.data()+ start);

    

    for (int i = 0; i < this->ncolB; i++)
    {
        for (int j = 0; j < ncolA_global; j++)
        {
            EXPECT_NEAR(get_double(this->Cref_global[i * LDC_global + start + j]),
                        get_double(this->C_global[i * LDC_global + start + j]),
                        1e-10);
        }
    }
}

TYPED_TEST(PgemmTest, row_parallel)
{
    const int ncolA_global = 17, ncolB_global = 7, nrow_global = 13;
    const int LDA_global = 17, LDB_global = 18, LDC_global = 19;

    this->decide_ngroup(1, 4);
    this->prepare(ncolA_global, ncolB_global, nrow_global, LDA_global, LDB_global, LDC_global);

    this->pgemm.set_dimension(this->col_world,
                              this->row_world,
                              this->ncolA,
                              this->LDA,
                              this->ncolB,
                              this->LDB,
                              this->nrow,
                              LDC_global);
    this->pgemm.multiply(this->alpha, this->A_local.data(), this->B_local.data(), this->beta, this->C_global.data());

    this->compare_result(ncolA_global, ncolB_global, LDC_global);
}

TYPED_TEST(PgemmTest, col_parallel)
{
    const int ncolA_global = 17, ncolB_global = 7, nrow_global = 13;
    const int LDA_global = 17, LDB_global = 18, LDC_global = 19;

    this->decide_ngroup(4, 1);
    this->prepare(ncolA_global, ncolB_global, nrow_global, LDA_global, LDB_global, LDC_global);

    this->pgemm.set_dimension(this->col_world,
                              this->row_world,
                              this->ncolA,
                              this->LDA,
                              this->ncolB,
                              this->LDB,
                              this->nrow,
                              LDC_global);
    this->pgemm.multiply(this->alpha, this->A_local.data(), this->B_local.data(), this->beta, this->C_global.data());

    this->compare_result(ncolA_global, ncolB_global, LDC_global);
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    MPI_Init(&argc, &argv);

    int RANK, NPROC;
    MPI_Comm_rank(MPI_COMM_WORLD, &RANK);
    MPI_Comm_size(MPI_COMM_WORLD, &NPROC);

    int result = RUN_ALL_TESTS();

    MPI_Finalize();
    return result;
}