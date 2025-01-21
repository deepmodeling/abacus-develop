#include "../para_linear_transform.h"

#include <gtest/gtest.h>
#ifdef __MPI
#include <mpi.h>
#endif

void random_data(std::vector<double>& A_global, std::vector<double>& U_global, double& alpha, double& beta)
{
    for (auto& val: A_global)
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
                 std::vector<std::complex<double>>& U_global,
                 std::complex<double>& alpha,
                 std::complex<double>& beta)
{
    for (auto& val: A_global)
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
    void prepare(const int nrow, const int ncol_glo, const int LDA)
    {
        int rank = 0;
        int nproc = 1;
        int colA_start = 0;
        this->ncol_glo = ncol_glo;
        this->ncol_loc = ncol_glo;
#ifdef __MPI
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &nproc);
        this->ncol_loc = ncol_glo / nproc;
        if (rank < ncol_glo % nproc)
        {
            ncol_loc++;
        }
        std::vector<int> ncolA_ip(nproc);
        MPI_Allgather(&ncol_loc, 1, MPI_INT, ncolA_ip.data(), 1, MPI_INT, MPI_COMM_WORLD);
        for (int i = 0; i < rank; ++i)
        {
            colA_start += ncolA_ip[i];
        }
#endif
        A_global.resize(LDA * ncol_glo);
        A_global_ref.resize(LDA * ncol_glo);
        U_global.resize(ncol_glo * ncol_glo);
        if (rank == 0)
        {
            random_data(A_global, U_global, alpha, beta);
            A_global_ref = A_global;
            std::vector<T> A_global_tmp = A_global;
            const base_device::DEVICE_CPU* ctx = {};
            ModuleBase::gemm_op<T, base_device::DEVICE_CPU>()(ctx,
                                                              'N',
                                                              'N',
                                                              nrow,
                                                              ncol_glo,
                                                              ncol_glo,
                                                              &alpha,
                                                              A_global_tmp.data(),
                                                              LDA,
                                                              U_global.data(),
                                                              ncol_glo,
                                                              &beta,
                                                              A_global_ref.data(),
                                                              LDA);
        }
        if (std::is_same<T, double>::value)
        {
#ifdef __MPI
            MPI_Bcast(A_global.data(), A_global.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
            MPI_Bcast(U_global.data(), U_global.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
            MPI_Bcast(A_global_ref.data(), A_global_ref.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
            MPI_Bcast(&alpha, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            MPI_Bcast(&beta, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif
        }
        else if (std::is_same<T, std::complex<double>>::value)
        {
#ifdef __MPI
            MPI_Bcast(A_global.data(), A_global.size(), MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
            MPI_Bcast(U_global.data(), U_global.size(), MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
            MPI_Bcast(A_global_ref.data(), A_global_ref.size(), MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
            MPI_Bcast(&alpha, 1, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
            MPI_Bcast(&beta, 1, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
#endif
        }

        A.resize(LDA * ncol_loc);
        A_ref.resize(LDA * ncol_loc);
        for (int i = 0; i < LDA * ncol_loc; ++i)
        {
            A[i] = A_global[colA_start * LDA + i];
            A_ref[i] = A_global_ref[colA_start * LDA + i];
        }
    }
    std::vector<T> A;
    std::vector<T> A_ref;
    std::vector<T> A_global;
    std::vector<T> U_global;
    std::vector<T> A_global_ref;
    int ncol_glo = 1;
    int ncol_loc = 1;
    T alpha;
    T beta;
};

typedef ::testing::Types<double, std::complex<double>> MyTypes;
TYPED_TEST_SUITE(ParaLinearTransformTest, MyTypes);

TYPED_TEST(ParaLinearTransformTest, cpucase)
{
    const int nrow = 7;
    const int ncol_glo = 13;
    const int LDA = 9;

    this->prepare(nrow, ncol_glo, LDA);
    int rank_col = 0, nproc_col = 1;
#ifdef __MPI
    MPI_Comm col_world = MPI_COMM_WORLD;
    MPI_Comm_rank(col_world, &rank_col);
    MPI_Comm_size(col_world, &nproc_col);
#endif

    hsolver::para_linear_transform_op<TypeParam, base_device::DEVICE_CPU>()(this->A.data(),
                                                                            this->alpha,
                                                                            this->beta,
                                                                            this->U_global.data(),
                                                                            nrow,
                                                                            LDA,
                                                                            this->ncol_loc,
                                                                            ncol_glo,
#ifdef __MPI
                                                                            col_world,
#endif
                                                                            rank_col,
                                                                            nproc_col);

    for (int i = 0; i < this->ncol_loc; ++i)
    {
        for (int j = 0; j < nrow; ++j)
        {
            EXPECT_NEAR(get_double(this->A[j + i * LDA]), get_double(this->A_ref[j + i * LDA]), 1e-10);
        }
    }
}

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