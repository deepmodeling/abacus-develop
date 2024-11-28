#include <gtest/gtest.h>
#include <fstream>
#include "../dm_trans.h"
#ifdef __MPI
#include "mpi.h"
#include "module_lr/utils/lr_util.h"
#endif
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

class DMTransTest : public testing::Test
{
public:
    std::vector<matsize> sizes{
        {2, 14, 9, 4},
        {2, 20, 10, 7}
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

    void set_ones(double* data, int size) { for (int i = 0;i < size;++i) data[i] = 1.0; };
    void set_int(double* data, int size) { for (int i = 0;i < size;++i) data[i] = static_cast<double>(i + 1); };
    void set_int(std::complex<double>* data, int size) { for (int i = 0;i < size;++i) data[i] = std::complex<double>(i + 1, -i - 1); };
    void set_rand(double* data, int size) { for (int i = 0;i < size;++i) data[i] = double(rand()) / double(RAND_MAX) * 10.0 - 5.0; };
    void set_rand(std::complex<double>* data, int size) { for (int i = 0;i < size;++i) data[i] = std::complex<double>(rand(), rand()) / double(RAND_MAX) * 10.0 - 5.0; };
    void check_eq(double* data1, double* data2, int size) { for (int i = 0;i < size;++i) EXPECT_NEAR(data1[i], data2[i], 1e-10); };
    void check_eq(std::complex<double>* data1, std::complex<double>* data2, int size)
    {
        for (int i = 0;i < size;++i)
        {
            EXPECT_NEAR(data1[i].real(), data2[i].real(), 1e-10);
            EXPECT_NEAR(data1[i].imag(), data2[i].imag(), 1e-10);
        }
    };
};

TEST_F(DMTransTest, DoubleSerial)
{
    for (auto s : this->sizes)
    {
        psi::Psi<double> X(s.nks, nstate, s.nocc * s.nvirt, nullptr, false);
        set_rand(X.get_pointer(), nstate * s.nks * s.nocc * s.nvirt);
        for (int istate = 0;istate < nstate;++istate)
        {
            int size_c = s.nks * (s.nocc + s.nvirt) * s.naos;
            psi::Psi<double> c(s.nks, s.nocc + s.nvirt, s.naos);
            set_rand(c.get_pointer(), size_c);
            X.fix_b(istate);
            const std::vector<container::Tensor>& dm_for = LR::cal_dm_trans_forloop_serial(X.get_pointer(), c, s.nocc, s.nvirt);
            const std::vector<container::Tensor>& dm_blas = LR::cal_dm_trans_blas(X.get_pointer(), c, s.nocc, s.nvirt);
            for (int isk = 0;isk < s.nks;++isk) check_eq(dm_for[isk].data<double>(), dm_blas[isk].data<double>(), s.naos * s.naos);
        }

    }
}
TEST_F(DMTransTest, ComplexSerial)
{
    for (auto s : this->sizes)
    {
        psi::Psi<std::complex<double>> X(s.nks, nstate, s.nocc * s.nvirt, nullptr, false);
        set_rand(X.get_pointer(), nstate * s.nks * s.nocc * s.nvirt);
        for (int istate = 0;istate < nstate;++istate)
        {
            int size_c = s.nks * (s.nocc + s.nvirt) * s.naos;
            psi::Psi<std::complex<double>> c(s.nks, s.nocc + s.nvirt, s.naos);
            set_rand(c.get_pointer(), size_c);
            X.fix_b(istate);
            const std::vector<container::Tensor>& dm_for = LR::cal_dm_trans_forloop_serial(X.get_pointer(), c, s.nocc, s.nvirt);
            const std::vector<container::Tensor>& dm_blas = LR::cal_dm_trans_blas(X.get_pointer(), c, s.nocc, s.nvirt);
            for (int isk = 0;isk < s.nks;++isk) check_eq(dm_for[isk].data<std::complex<double>>(), dm_blas[isk].data<std::complex<double>>(), s.naos * s.naos);
        }

    }
}

#ifdef __MPI
TEST_F(DMTransTest, DoubleParallel)
{
    for (auto s : this->sizes)
    {
        // c: nao*nbands in para2d, nbands*nao in psi  (row-para and constructed: nao)
        // X: nvirt*nocc in para2d, nocc*nvirt in psi (row-para and constructed: nvirt)
        Parallel_2D px;
        LR_Util::setup_2d_division(px, s.nb, s.nvirt, s.nocc);
        psi::Psi<double> X(s.nks, nstate, px.get_local_size(), nullptr, false);
        Parallel_2D pc;
        LR_Util::setup_2d_division(pc, s.nb, s.naos, s.nocc + s.nvirt, px.blacs_ctxt);
        psi::Psi<double> c(s.nks, pc.get_col_size(), pc.get_row_size());
        Parallel_2D pmat;
        LR_Util::setup_2d_division(pmat, s.nb, s.naos, s.naos, px.blacs_ctxt);

        EXPECT_EQ(px.dim0, pc.dim0);
        EXPECT_EQ(px.dim1, pc.dim1);
        EXPECT_GE(s.nvirt, px.dim0);
        EXPECT_GE(s.nocc, px.dim1);
        EXPECT_GE(s.naos, pc.dim0);

        set_rand(X.get_pointer(), nstate * s.nks * px.get_local_size());        //set X and X_full
        psi::Psi<double> X_full(s.nks, nstate, s.nocc * s.nvirt, nullptr, false);        // allocate X_full
        for (int istate = 0;istate < nstate;++istate)
        {
            X.fix_b(istate);
            X_full.fix_b(istate);
            for (int isk = 0;isk < s.nks;++isk)
            {
                X.fix_k(isk);
                X_full.fix_k(isk);
                LR_Util::gather_2d_to_full(px, X.get_pointer(), X_full.get_pointer(), false, s.nvirt, s.nocc);
            }
        }
        for (int istate = 0;istate < nstate;++istate)
        {
            c.fix_k(0);
            set_rand(c.get_pointer(), s.nks * pc.get_local_size()); // set c 

            X.fix_b(istate);
            X_full.fix_b(istate);

            std::vector<container::Tensor> dm_pblas_loc = LR::cal_dm_trans_pblas(X.get_pointer(), px, c, pc, s.naos, s.nocc, s.nvirt, pmat);

            // gather dm and output
            std::vector<container::Tensor> dm_gather(s.nks, container::Tensor(DAT::DT_DOUBLE, DEV::CpuDevice, { s.naos, s.naos }));
            for (int isk = 0;isk < s.nks;++isk)
                LR_Util::gather_2d_to_full(pmat, dm_pblas_loc[isk].data<double>(), dm_gather[isk].data<double>(), false, s.naos, s.naos);

            // compare to global matrix
            psi::Psi<double> c_full(s.nks, s.nocc + s.nvirt, s.naos);
            for (int isk = 0;isk < s.nks;++isk)
            {
                c.fix_k(isk);
                c_full.fix_k(isk);
                LR_Util::gather_2d_to_full(pc, c.get_pointer(), c_full.get_pointer(), false, s.naos, s.nocc + s.nvirt);
            }
            if (my_rank == 0)
            {
                const std::vector<container::Tensor>& dm_full = LR::cal_dm_trans_blas(X_full.get_pointer(), c_full, s.nocc, s.nvirt);
                for (int isk = 0;isk < s.nks;++isk) check_eq(dm_full[isk].data<double>(), dm_gather[isk].data<double>(), s.naos * s.naos);
            }
        }
    }
}
TEST_F(DMTransTest, ComplexParallel)
{
    for (auto s : this->sizes)
    {
        // c: nao*nbands in para2d, nbands*nao in psi  (row-para and constructed: nao)
        // X: nvirt*nocc in para2d, nocc*nvirt in psi (row-para and constructed: nvirt)
        Parallel_2D px;
        LR_Util::setup_2d_division(px, s.nb, s.nvirt, s.nocc);
        psi::Psi<std::complex<double>> X(s.nks, nstate, px.get_local_size(), nullptr, false);
        Parallel_2D pc;
        LR_Util::setup_2d_division(pc, s.nb, s.naos, s.nocc + s.nvirt, px.blacs_ctxt);
        psi::Psi<std::complex<double>> c(s.nks, pc.get_col_size(), pc.get_row_size());
        Parallel_2D pmat;
        LR_Util::setup_2d_division(pmat, s.nb, s.naos, s.naos, px.blacs_ctxt);

        set_rand(X.get_pointer(), nstate * s.nks * px.get_local_size());        //set X and X_full
        psi::Psi<std::complex<double>> X_full(s.nks, nstate, s.nocc * s.nvirt, nullptr, false);        // allocate X_full
        for (int istate = 0;istate < nstate;++istate)
        {
            X.fix_b(istate);
            X_full.fix_b(istate);
            for (int isk = 0;isk < s.nks;++isk)
            {
                X.fix_k(isk);
                X_full.fix_k(isk);
                LR_Util::gather_2d_to_full(px, X.get_pointer(), X_full.get_pointer(), false, s.nvirt, s.nocc);
            }
        }
        for (int istate = 0;istate < nstate;++istate)
        {
            c.fix_k(0);
            set_rand(c.get_pointer(), s.nks * pc.get_local_size()); // set c 

            X.fix_b(istate);
            X_full.fix_b(istate);

            std::vector<container::Tensor> dm_pblas_loc = LR::cal_dm_trans_pblas(X.get_pointer(), px, c, pc, s.naos, s.nocc, s.nvirt, pmat);

            // gather dm and output
            std::vector<container::Tensor> dm_gather(s.nks, container::Tensor(DAT::DT_COMPLEX_DOUBLE, DEV::CpuDevice, { s.naos, s.naos }));
            for (int isk = 0;isk < s.nks;++isk)
                LR_Util::gather_2d_to_full(pmat, dm_pblas_loc[isk].data<std::complex<double>>(), dm_gather[isk].data<std::complex<double>>(), false, s.naos, s.naos);

            // compare to global matrix
            psi::Psi<std::complex<double>> c_full(s.nks, s.nocc + s.nvirt, s.naos);
            for (int isk = 0;isk < s.nks;++isk)
            {
                c.fix_k(isk);
                c_full.fix_k(isk);
                LR_Util::gather_2d_to_full(pc, c.get_pointer(), c_full.get_pointer(), false, s.naos, s.nocc + s.nvirt);
            }
            if (my_rank == 0)
            {
                std::vector<container::Tensor> dm_full = LR::cal_dm_trans_blas(X_full.get_pointer(), c_full, s.nocc, s.nvirt);
                for (int isk = 0;isk < s.nks;++isk) check_eq(dm_full[isk].data<std::complex<double>>(), dm_gather[isk].data<std::complex<double>>(), s.naos * s.naos);
            }
        }
    }
}
#endif


int main(int argc, char** argv)
{
    srand(time(NULL));  // for random number generator
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
