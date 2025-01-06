#include <gtest/gtest.h>
#ifdef __MPI
#include "mpi.h"
#endif
#include "../AX.h"

#include "module_lr/utils/lr_util.h"

#define rand01 (static_cast<double>(rand()) / static_cast<double>(RAND_MAX) - 0.5 )
struct matsize
{
    int nks = 1;
    int naos = 2;
    int nocc = 1;
    int nvirt = 1;
    int nb = 1;
    matsize(int nks, int naos, int nocc, int nvirt, int nb = 1)
        :nks(nks), naos(naos), nocc(nocc), nvirt(nvirt), nb(nb) {
        assert(nocc + nvirt <= naos);
    };
};

class AXTest : public testing::Test
{
public:
    std::vector<matsize> sizes{
        // {2, 3, 2, 1},
        {2, 13, 7, 4},
        // {2, 14, 8, 5}
    };
    int nstate = 2;
    std::ofstream ofs_running;
    int my_rank = 0;
#ifdef __MPI
    void SetUp() override
    {
        MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
        this->ofs_running.open("log" + std::to_string(my_rank) + ".txt");
        ofs_running << "my_rank = " << my_rank << std::endl;
    }
    void TearDown() override
    {
        ofs_running.close();
    }
#endif

    void set_ones(double* data, int size) { for (int i = 0;i < size;++i) { data[i] = 1.0; } };
    void set_int(double* data, int size) { for (int i = 0;i < size;++i) { data[i] = static_cast<double>(i + 1); } };
    void set_int(std::complex<double>* data, int size) { for (int i = 0;i < size;++i) { data[i] = std::complex<double>(i + 1, -i - 1); } };
    void set_rand(double* data, int size) { for (int i = 0;i < size;++i) { data[i] = rand01 * 10; } };
    void set_rand(std::complex<double>* data, int size) { for (int i = 0;i < size;++i) { data[i] = std::complex<double>(rand01 * 10, rand01 * 10); } };
    void check_eq(double* data1, double* data2, int size) { for (int i = 0;i < size;++i) { EXPECT_NEAR(data1[i], data2[i], 1e-10); } };
    void check_eq(std::complex<double>* data1, std::complex<double>* data2, int size)
    {
        for (int i = 0;i < size;++i)
        {
            EXPECT_NEAR(data1[i].real(), data2[i].real(), 1e-10);
            EXPECT_NEAR(data1[i].imag(), data2[i].imag(), 1e-10);
        }
    };
};

TEST_F(AXTest, DoubleSerial)
{
    for (auto s : this->sizes)
    {
        psi::Psi<double> AX_for(s.nks, nstate, s.nocc * s.nvirt, nullptr, false);
        psi::Psi<double> AX_blas(s.nks, nstate, s.nocc * s.nvirt, nullptr, false);
        int size_c = s.nks * (s.nocc + s.nvirt) * s.naos;
        int size_v = s.naos * s.naos;
        for (int istate = 0;istate < nstate;++istate)
        {
            std::vector<int> temp(s.nks, s.naos);
            psi::Psi<double> c(s.nks, s.nocc + s.nvirt, s.naos, temp.data(), true);
            std::vector<container::Tensor> V(s.nks, container::Tensor(DAT::DT_DOUBLE, DEV::CpuDevice, { s.naos, s.naos }));
            set_rand(&c(0, 0, 0), size_c);
            for (auto& v : V) { set_rand(v.data<double>(), size_v); }
            LR::cal_AX_forloop_serial(V, c, s.nocc, s.nvirt, &AX_for(istate, 0, 0));
            LR::cal_AX_blas(V, c, s.nocc, s.nvirt, &AX_blas(istate, 0, 0), false);
        }
        check_eq(&AX_for(0, 0, 0), &AX_blas(0, 0, 0), nstate * s.nks * s.nocc * s.nvirt);
    }
}

TEST_F(AXTest, ComplexSerial)
{
    for (auto s : this->sizes)
    {
        psi::Psi<std::complex<double>> AX_for(s.nks, nstate, s.nocc * s.nvirt, nullptr, false);
        psi::Psi<std::complex<double>> AX_blas(s.nks, nstate, s.nocc * s.nvirt, nullptr, false);
        int size_c = s.nks * (s.nocc + s.nvirt) * s.naos;
        int size_v = s.naos * s.naos;
        for (int istate = 0;istate < nstate;++istate)
        {
            std::vector<int> temp(s.nks, s.naos);
            psi::Psi<std::complex<double>> c(s.nks, s.nocc + s.nvirt, s.naos, temp.data(), true);
            std::vector<container::Tensor> V(s.nks, container::Tensor(DAT::DT_COMPLEX_DOUBLE, DEV::CpuDevice, { s.naos, s.naos }));
            set_rand(&c(0, 0, 0), size_c);
            for (auto& v : V) { set_rand(v.data<std::complex<double>>(), size_v); }
            LR::cal_AX_forloop_serial(V, c, s.nocc, s.nvirt, &AX_for(istate, 0, 0));
            LR::cal_AX_blas(V, c, s.nocc, s.nvirt, &AX_blas(istate, 0, 0), false);
        }
        check_eq(&AX_for(0, 0, 0), &AX_blas(0, 0, 0), nstate * s.nks * s.nocc * s.nvirt);
    }
}
#ifdef __MPI
TEST_F(AXTest, DoubleParallel)
{
    for (auto s : this->sizes)
    {
        // c: nao*nbands in para2d, nbands*nao in psi  (row-para and constructed: nao)
        // X: nvirt*nocc in para2d, nocc*nvirt in psi (row-para and constructed: nvirt)
        Parallel_2D pV;
        LR_Util::setup_2d_division(pV, s.nb, s.naos, s.naos);
        std::vector<container::Tensor> V(s.nks, container::Tensor(DAT::DT_DOUBLE, DEV::CpuDevice, { pV.get_col_size(), pV.get_row_size() }));
        Parallel_2D pc;
        LR_Util::setup_2d_division(pc, s.nb, s.naos, s.nocc + s.nvirt, pV.blacs_ctxt);
        
        std::vector<int> ngk_temp(s.nks, pc.get_row_size());
        psi::Psi<double> c(s.nks, pc.get_col_size(), pc.get_row_size(), ngk_temp.data(), true);
        Parallel_2D px;
        LR_Util::setup_2d_division(px, s.nb, s.nvirt, s.nocc, pV.blacs_ctxt);

        EXPECT_EQ(pV.dim0, pc.dim0);
        EXPECT_EQ(pV.dim1, pc.dim1);
        EXPECT_GE(s.nvirt, px.dim0);
        EXPECT_GE(s.nocc, px.dim1);
        EXPECT_GE(s.naos, pc.dim0);
        psi::Psi<double> AX_pblas_loc(s.nks, nstate, px.get_local_size(), nullptr, false);
        psi::Psi<double> AX_gather(s.nks, nstate, s.nocc * s.nvirt, nullptr, false);
        for (int istate = 0;istate < nstate;++istate)
        {
            for (int isk = 0;isk < s.nks;++isk)
            {
                set_rand(V.at(isk).data<double>(), pV.get_local_size());
                set_rand(&c(isk, 0, 0), pc.get_local_size());
            }
            LR::cal_AX_pblas(V, pV, c, pc, s.naos, s.nocc, s.nvirt, px, &AX_pblas_loc(istate, 0, 0), false);
            // gather AX and output
            for (int isk = 0;isk < s.nks;++isk)
            {
                LR_Util::gather_2d_to_full(px, &AX_pblas_loc(istate, isk, 0), &AX_gather(istate, isk, 0), false/*pblas: row first*/, s.nvirt, s.nocc);
            }
            // compare to global AX
            std::vector<container::Tensor> V_full(s.nks, container::Tensor(DAT::DT_DOUBLE, DEV::CpuDevice, { s.naos, s.naos }));
            
            std::vector<int> ngk_temp_1(s.nks, s.naos);
            psi::Psi<double> c_full(s.nks, s.nocc + s.nvirt, s.naos, ngk_temp_1.data(), true);
            for (int isk = 0;isk < s.nks;++isk)
            {
                LR_Util::gather_2d_to_full(pV, V.at(isk).data<double>(), V_full.at(isk).data<double>(), false, s.naos, s.naos);
                LR_Util::gather_2d_to_full(pc, &c(isk, 0, 0), &c_full(isk, 0, 0), false, s.naos, s.nocc + s.nvirt);
            }
            if (my_rank == 0)
            {
                psi::Psi<double>  AX_full_istate(s.nks, 1, s.nocc * s.nvirt, nullptr, false);
                LR::cal_AX_blas(V_full, c_full, s.nocc, s.nvirt, &AX_full_istate(0, 0, 0), false);
                check_eq(&AX_full_istate(0, 0, 0), &AX_gather(istate, 0, 0), s.nks * s.nocc * s.nvirt);
            }
        }
    }
}
TEST_F(AXTest, ComplexParallel)
{
    for (auto s : this->sizes)
    {
        // c: nao*nbands in para2d, nbands*nao in psi  (row-para and constructed: nao)
        // X: nvirt*nocc in para2d, nocc*nvirt in psi (row-para and constructed: nvirt)
        Parallel_2D pV;
        LR_Util::setup_2d_division(pV, s.nb, s.naos, s.naos);
        std::vector<container::Tensor> V(s.nks, container::Tensor(DAT::DT_COMPLEX_DOUBLE, DEV::CpuDevice, { pV.get_col_size(), pV.get_row_size() }));
        Parallel_2D pc;
        LR_Util::setup_2d_division(pc, s.nb, s.naos, s.nocc + s.nvirt, pV.blacs_ctxt);

        std::vector<int> ngk_temp_1(s.nks, pc.get_row_size());
        psi::Psi<std::complex<double>> c(s.nks, pc.get_col_size(), pc.get_row_size(), ngk_temp_1.data(), true);
        Parallel_2D px;
        LR_Util::setup_2d_division(px, s.nb, s.nvirt, s.nocc, pV.blacs_ctxt);

        psi::Psi<std::complex<double>> AX_pblas_loc(s.nks, nstate, px.get_local_size(), nullptr, false);
        psi::Psi<std::complex<double>> AX_gather(s.nks, nstate, s.nocc * s.nvirt, nullptr, false);
        for (int istate = 0;istate < nstate;++istate)
        {
            for (int isk = 0;isk < s.nks;++isk)
            {
                set_rand(V.at(isk).data<std::complex<double>>(), pV.get_local_size());
                set_rand(&c(isk, 0, 0), pc.get_local_size());
            }
            LR::cal_AX_pblas(V, pV, c, pc, s.naos, s.nocc, s.nvirt, px, &AX_pblas_loc(istate, 0, 0), false);

            // gather AX and output
            for (int isk = 0;isk < s.nks;++isk)
            {
                LR_Util::gather_2d_to_full(px, &AX_pblas_loc(istate, isk, 0), &AX_gather(istate, isk, 0), false/*pblas: row first*/, s.nvirt, s.nocc);
            }
            // compare to global AX
            std::vector<container::Tensor> V_full(s.nks, container::Tensor(DAT::DT_COMPLEX_DOUBLE, DEV::CpuDevice, { s.naos, s.naos }));


            std::vector<int> ngk_temp_2(s.nks, s.naos);
            psi::Psi<std::complex<double>> c_full(s.nks, s.nocc + s.nvirt, s.naos, ngk_temp_2.data(), true);
            for (int isk = 0;isk < s.nks;++isk)
            {
                LR_Util::gather_2d_to_full(pV, V.at(isk).data<std::complex<double>>(), V_full.at(isk).data<std::complex<double>>(), false, s.naos, s.naos);
                LR_Util::gather_2d_to_full(pc, &c(isk, 0, 0), &c_full(isk, 0, 0), false, s.naos, s.nocc + s.nvirt);
            }
            if (my_rank == 0)
            {
                psi::Psi<std::complex<double>>  AX_full_istate(s.nks, 1, s.nocc * s.nvirt, nullptr, false);
                LR::cal_AX_blas(V_full, c_full, s.nocc, s.nvirt, &AX_full_istate(0, 0, 0), false);
                check_eq(&AX_full_istate(0, 0, 0), &AX_gather(istate, 0, 0), s.nks * s.nocc * s.nvirt);
            }
        }
    }
}
#endif


int main(int argc, char** argv)
{
    srand(time(nullptr));  // for random number generator
#ifdef __MPI
    MPI_Init(&argc, &argv);
#endif
    testing::InitGoogleTest(&argc, argv);
    int result = RUN_ALL_TESTS();
#ifdef __MPI
    MPI_Finalize();
#endif
    return result;
}
