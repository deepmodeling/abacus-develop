#include "../para_linear_transform.h"

#include <gtest/gtest.h>
#ifdef __MPI
#include <mpi.h>
#endif

void random_data(std::vector<double>& A_global,
                 std::vector<double>& B_global,
                 std::vector<double>& U_global,
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
    for (auto& val: U_global)
    {
        val = std::rand() / (RAND_MAX + 1.0);
    }
    alpha = std::rand() / (RAND_MAX + 1.0);
    beta = std::rand() / (RAND_MAX + 1.0);
}
void random_data(std::vector<std::complex<double>>& A_global,
                 std::vector<std::complex<double>>& B_global,
                 std::vector<std::complex<double>>& U_global,
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
    for (auto& val: U_global)
    {
        val = std::complex<double>(std::rand() / (RAND_MAX + 1.0), std::rand() / (RAND_MAX + 1.0));
    }
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

template <typename T>
class ParaLinearTransformTest : public ::testing::Test
{
  protected:
    void SetUp() override
    {
    }

    void TearDown() override
    {
    }
    void prepare(const int nrow, const int ncolA_glo, const int ncolB_glo, const int LDA)
    {
        int rank = 0;
        int nproc = 1;
        int colA_start = 0;
        int colB_start = 0;
        this->ncolA_glo = ncolA_glo;
        this->ncolB_glo = ncolB_glo;
        this->ncolA_loc = ncolA_glo;
        this->ncolB_loc = ncolB_glo;
#ifdef __MPI
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &nproc);
        this->ncolA_loc = ncolA_glo / nproc;
        this->ncolB_loc = ncolB_glo / nproc;
        if (rank < ncolA_glo % nproc)
        {
            ncolA_loc++;
            ncolB_loc++;
        }
        std::vector<int> ncolA_ip(nproc);
        std::vector<int> ncolB_ip(nproc);
        MPI_Allgather(&ncolA_loc, 1, MPI_INT, ncolA_ip.data(), 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Allgather(&ncolB_loc, 1, MPI_INT, ncolB_ip.data(), 1, MPI_INT, MPI_COMM_WORLD);
        for (int i = 0; i < rank; ++i)
        {
            colA_start += ncolA_ip[i];
            colB_start += ncolB_ip[i];
        }
#endif
        A_global.resize(LDA * ncolA_glo);
        B_global.resize(LDA * ncolB_glo);
        B_global_ref.resize(LDA * ncolB_glo);
        U_global.resize(ncolA_glo * ncolB_glo);
        if (rank == 0)
        {
            random_data(A_global, B_global, U_global, alpha, beta);
            B_global_ref = B_global;
            const base_device::DEVICE_CPU* ctx = {};
            ModuleBase::gemm_op<T, base_device::DEVICE_CPU>()(ctx,
                                                              'N',
                                                              'N',
                                                              nrow,
                                                              ncolB_glo,
                                                              ncolA_glo,
                                                              &alpha,
                                                              A_global.data(),
                                                              LDA,
                                                              U_global.data(),
                                                              ncolA_glo,
                                                              &beta,
                                                              B_global_ref.data(),
                                                              LDA);
        }
        if (std::is_same<T, double>::value)
        {
#ifdef __MPI
            MPI_Bcast(A_global.data(), A_global.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
            MPI_Bcast(B_global.data(), B_global.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
            MPI_Bcast(U_global.data(), U_global.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
            MPI_Bcast(B_global_ref.data(), B_global_ref.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
            MPI_Bcast(&alpha, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            MPI_Bcast(&beta, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif
        }
        else if (std::is_same<T, std::complex<double>>::value)
        {
#ifdef __MPI
            MPI_Bcast(A_global.data(), A_global.size(), MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
            MPI_Bcast(B_global.data(), B_global.size(), MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
            MPI_Bcast(U_global.data(), U_global.size(), MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
            MPI_Bcast(B_global_ref.data(), B_global_ref.size(), MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
            MPI_Bcast(&alpha, 1, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
            MPI_Bcast(&beta, 1, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
#endif
        }

        A.resize(LDA * ncolA_loc);
        B.resize(LDA * ncolB_loc);
        B_ref.resize(LDA * ncolB_loc);
        for (int i = 0; i < LDA * ncolA_loc; ++i)
        {
            A[i] = A_global[colA_start * LDA + i];
        }
        for (int i = 0; i < LDA * ncolB_loc; ++i)
        {
            B[i] = B_global[colB_start * LDA + i];
            B_ref[i] = B_global_ref[colB_start * LDA + i];
        }
    }
    std::vector<T> A, B;
    std::vector<T> B_ref;
    std::vector<T> A_global;
    std::vector<T> B_global;
    std::vector<T> U_global;
    std::vector<T> B_global_ref;
    int ncolA_glo = 1;
    int ncolB_glo = 1;
    int ncolA_loc = 1;
    int ncolB_loc = 1;
    T alpha;
    T beta;
    hsolver::PLinearTransform<T, base_device::DEVICE_CPU> lt;
};

typedef ::testing::Types<double, std::complex<double>> MyTypes;
TYPED_TEST_SUITE(ParaLinearTransformTest, MyTypes);

TYPED_TEST(ParaLinearTransformTest, globalU)
{
    const int nrowA = 7;
    const int ncolA_glo = 13;
    const int ncolB_glo = 11;
    const int LDA = 9;

    this->prepare(nrowA, ncolA_glo, ncolB_glo, LDA);
    int rank_col = 0, nproc_col = 1;
#ifdef __MPI
    MPI_Comm col_world = MPI_COMM_WORLD;
    MPI_Comm_rank(col_world, &rank_col);
    MPI_Comm_size(col_world, &nproc_col);
#endif

    this->lt.set_dimension(nrowA,
                           this->ncolA_loc,
                           this->ncolB_loc,
                           LDA,
#ifdef __MPI
                           col_world,
#endif
                           false);
    this->lt.act(this->alpha, this->A.data(), this->U_global.data(), this->beta, this->B.data());

    for (int i = 0; i < this->ncolB_loc; ++i)
    {
        for (int j = 0; j < nrowA; ++j)
        {
            EXPECT_NEAR(get_double(this->B[j + i * LDA]), get_double(this->B_ref[j + i * LDA]), 1e-10);
        }
    }
}
#ifdef __MPI
TYPED_TEST(ParaLinearTransformTest, localU)
{
    const int nrowA = 7;
    const int ncolA_glo = 13;
    const int ncolB_glo = 11;
    const int LDA = 9;

    this->prepare(nrowA, ncolA_glo, ncolB_glo, LDA);
    int rank_col = 0, nproc_col = 1;

    MPI_Comm col_world = MPI_COMM_WORLD;
    MPI_Comm_rank(col_world, &rank_col);
    MPI_Comm_size(col_world, &nproc_col);
    std::vector<int> ncolB_ip(nproc_col);
    std::vector<int> start_colB(nproc_col);
    MPI_Allgather(&this->ncolB_loc, 1, MPI_INT, ncolB_ip.data(), 1, MPI_INT, col_world);
    start_colB[0] = 0;
    for (int i = 1; i < nproc_col; ++i)
    {
        start_colB[i] = start_colB[i - 1] + ncolB_ip[i - 1];
    }
    int start = start_colB[rank_col];

    this->lt.set_dimension(nrowA, this->ncolA_loc, this->ncolB_loc, LDA, col_world, true);

    this->lt.act(this->alpha, this->A.data(), this->U_global.data() + start * ncolA_glo, this->beta, this->B.data());

    for (int i = 0; i < this->ncolB_loc; ++i)
    {
        for (int j = 0; j < nrowA; ++j)
        {
            EXPECT_NEAR(get_double(this->B[j + i * LDA]), get_double(this->B_ref[j + i * LDA]), 1e-10);
        }
    }
}
#endif

int main(int argc, char** argv)
{
#ifdef __MPI
    MPI_Init(&argc, &argv);
#endif
    ::testing::InitGoogleTest(&argc, argv);
    int result = RUN_ALL_TESTS();
#ifdef __MPI
    MPI_Finalize();
#endif
    return result;
}