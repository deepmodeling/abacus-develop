#include <gtest/gtest.h>

#include "../lr_util.h"
#include "../lr_util_print.h"

inline void check_double_eq(double* data1, double* data2, int size)
{
    for (int i = 0;i < size;++i)
        EXPECT_NEAR(data1[i], data2[i], 1e-10);
};
inline void check_double_eq(std::complex<double>* data1, std::complex<double>* data2, int size)
{
    for (int i = 0;i < size;++i)
    {
        EXPECT_NEAR(data1[i].real(), data2[i].real(), 1e-10);
        EXPECT_NEAR(data1[i].imag(), data2[i].imag(), 1e-10);
    }
};
inline void check_norm_eq(double* data1, double* data2, int size)
{
    for (int i = 0;i < size;++i)
    {
        EXPECT_NEAR(std::norm(data1[i]), std::norm(data2[i]), 1e-10);
    }
};
inline void check_norm_eq(std::complex<double>* data1, std::complex<double>* data2, int size)
{
    for (int i = 0;i < size;++i)
    {
        EXPECT_NEAR(std::norm(data1[i]), std::norm(data2[i]), 1e-10);
    }
}

TEST(LR_Util, PsiWrapper)
{
    int nk = 2;
    int nbands = 5;
    int nbasis = 6;

    psi::Psi<float> k1(1, nbands, nk * nbasis, nk * nbasis, true);
    for (int i = 0;i < nbands * nk * nbasis;++i)k1.get_pointer()[i] = i;

    k1.fix_b(2);
    psi::Psi<float> bf = LR_Util::k1_to_bfirst_wrapper(k1, nk, nbasis);
    EXPECT_EQ(k1.get_current_k(), 0);
    EXPECT_EQ(k1.get_current_b(), 2); // invariance after wrapper
    EXPECT_EQ(bf.get_current_k(), 0);
    EXPECT_EQ(bf.get_current_b(), 0);

    bf.fix_kb(1, 3);
    psi::Psi<float> kb = LR_Util::bfirst_to_k1_wrapper(bf);
    EXPECT_EQ(bf.get_current_k(), 1);
    EXPECT_EQ(bf.get_current_b(), 3);
    EXPECT_EQ(kb.get_current_k(), 0);
    EXPECT_EQ(kb.get_current_b(), 0);


    EXPECT_EQ(bf.get_k_first(), false);
    EXPECT_EQ(bf.get_nk(), nk);
    EXPECT_EQ(bf.get_nbands(), nbands);
    EXPECT_EQ(bf.get_nbasis(), nbasis);

    EXPECT_EQ(kb.get_k_first(), true);
    EXPECT_EQ(kb.get_nk(), 1);
    EXPECT_EQ(kb.get_nbands(), nbands);
    EXPECT_EQ(kb.get_nbasis(), nk * nbasis);

    k1.fix_b(0);
    bf.fix_kb(0, 0);
    EXPECT_EQ(bf.get_pointer(), k1.get_pointer());
    EXPECT_EQ(bf.get_pointer(), kb.get_pointer());
    for (int ik = 0; ik < nk; ik++)
    {
        for (int ib = 0; ib < nbands; ib++)
        {
            bf.fix_kb(ik, ib);
            kb.fix_b(ib);
            k1.fix_b(ib);
            for (int ibasis = 0; ibasis < nbasis; ibasis++)
            {
                int ikb = ik * nbasis + ibasis;
                EXPECT_EQ(kb(ikb), bf(ibasis));
                EXPECT_EQ(k1(ikb), kb(ikb));
            }
        }
    }
}
#ifdef __MPI
void set_rand(double* data, int size) { for (int i = 0;i < size;++i) data[i] = double(rand()) / double(RAND_MAX) * 10.0 - 5.0; };
TEST(LR_Util, MatSymDouble)
{
    int n = 7;
    std::vector<double> din(n * n);
    set_rand(din.data(), n * n);
    std::vector<double> dref(n * n, 0.0);
    LR_Util::matsym(din.data(), n, dref.data());

    Parallel_2D pmat;
    LR_Util::setup_2d_division(pmat, 1, n, n);
    std::vector<double> din_local(pmat.get_local_size(), 0.0);
    for (int i = 0;i < pmat.get_row_size();++i)
        for (int j = 0;j < pmat.get_col_size();++j)
            din_local[j * pmat.get_row_size() + i] = din[pmat.local2global_col(j) * n + pmat.local2global_row(i)];

    std::vector<double> dout_local(pmat.get_local_size(), 0.0);
    LR_Util::matsym(din_local.data(), n, pmat, dout_local.data());
    for (int i = 0;i < pmat.get_row_size();++i)
        for (int j = 0;j < pmat.get_col_size();++j)
            EXPECT_DOUBLE_EQ(dout_local[j * pmat.get_row_size() + i], dref[pmat.local2global_col(j) * n + pmat.local2global_row(i)]);

    //in-place version
    LR_Util::matsym(din.data(), n);
    for (int i = 0;i < n * n;++i)
        EXPECT_DOUBLE_EQ(din[i], dref[i]);

    LR_Util::matsym(din_local.data(), n, pmat);
    for (int i = 0;i < pmat.get_local_size();++i)
        EXPECT_DOUBLE_EQ(din_local[i], dout_local[i]);
}

void set_rand(std::complex<double>* data, int size) { for (int i = 0;i < size;++i) data[i] = std::complex<double>(rand(), rand()) / double(RAND_MAX) * 10.0 - 5.0; };
TEST(LR_Util, MatSymComplex)
{
    int n = 5;
    std::vector<std::complex<double>> din(n * n);
    set_rand(din.data(), n * n);
    std::vector<std::complex<double>> dref(n * n, std::complex<double>(0.0, 0.0));
    LR_Util::matsym(din.data(), n, dref.data());

    Parallel_2D pmat;
    LR_Util::setup_2d_division(pmat, 1, n, n);
    std::vector<std::complex<double>> din_local(pmat.get_local_size(), std::complex<double>(0.0, 0.0));
    for (int i = 0;i < pmat.get_row_size();++i)
        for (int j = 0;j < pmat.get_col_size();++j)
            din_local[j * pmat.get_row_size() + i] = din[pmat.local2global_col(j) * n + pmat.local2global_row(i)];

    std::vector<std::complex<double>> dout_local(pmat.get_local_size(), std::complex<double>(0.0, 0.0));
    LR_Util::matsym(din_local.data(), n, pmat, dout_local.data());
    for (int i = 0;i < pmat.get_row_size();++i)
        for (int j = 0;j < pmat.get_col_size();++j)
        {
            EXPECT_DOUBLE_EQ(dout_local[j * pmat.get_row_size() + i].real(), dref[pmat.local2global_col(j) * n + pmat.local2global_row(i)].real());
            EXPECT_DOUBLE_EQ(dout_local[j * pmat.get_row_size() + i].imag(), dref[pmat.local2global_col(j) * n + pmat.local2global_row(i)].imag());
        }

    //in-place version
    LR_Util::matsym(din.data(), n);
    for (int i = 0;i < n * n;++i)
    {
        EXPECT_DOUBLE_EQ(din[i].real(), dref[i].real());
        EXPECT_DOUBLE_EQ(din[i].imag(), dref[i].imag());
    }

    LR_Util::matsym(din_local.data(), n, pmat);
    for (int i = 0;i < pmat.get_local_size();++i)
    {
        EXPECT_DOUBLE_EQ(din_local[i].real(), dout_local[i].real());
        EXPECT_DOUBLE_EQ(din_local[i].imag(), dout_local[i].imag());
    }
}

TEST(LR_Util, RWValue)
{
    const std::string file = "RWValue.txt";
    std::ofstream ofs(file);
    std::vector<int> vec(2 * 3 * 4 * 5, 0);
    for (int i = 0;i < vec.size();++i) vec[i] = i;
    LR_Util::write_value(ofs, vec.data(), 2, 3, 4, 5);
    ofs.close();

    std::vector<int> vec1(2 * 3 * 4 * 5, 0);
    std::ifstream ifs1(file);
    EXPECT_EQ(LR_Util::read_value(ifs1, vec1.data(), 2, 3, 4, 5), 120);
    ifs1.close();
    for (int i = 0;i < vec1.size();++i) { EXPECT_EQ(vec1[i], vec[i]); };
    std::vector<int> vec2(2 * 3 * 4 * 5, 0);
    std::ifstream ifs2(file);
    EXPECT_EQ(LR_Util::read_value(ifs2, vec2.data(), 2 * 3, 4 * 5), 120);
    ifs2.close();
    for (int i = 0;i < vec2.size();++i) { EXPECT_EQ(vec2[i], vec[i]); };
}

TEST(LR_Util, DiagScaLapackDouble)
{
    // setup the matrix
    const int dim = 14;
    std::vector<double> mat(dim * dim);
    set_rand(mat.data(), dim * dim);
    LR_Util::matsym(mat.data(), dim);
    Parallel_2D pmat;
    LR_Util::setup_2d_division(pmat, 1, dim, dim);
    std::vector<double> mat_local(pmat.get_local_size(), 0.0);
    LR_Util::set_local_from_global(pmat, mat.data(), mat_local.data());

    // serial
    std::vector<double> eig(dim);
    LR_Util::diag_lapack(dim, mat.data(), eig.data());

    // parallel
    std::vector<double> eig_para(dim);
    std::vector<double> eigvec_para(pmat.get_local_size());
    LR_Util::diag_scalapack(dim, mat_local.data(), eig_para.data(), eigvec_para.data(), pmat.desc);

    // compare
    check_double_eq(eig_para.data(), eig.data(), dim);
    std::vector<double> eigvec_serial_local(pmat.get_local_size());
    LR_Util::set_local_from_global(pmat, mat.data(), eigvec_serial_local.data());
    check_norm_eq(eigvec_para.data(), eigvec_serial_local.data(), pmat.get_local_size());
}


TEST(LR_Util, DiagScaLapackGeneralComplex)
{
    // setup the matrix
    const int dim = 15;
    std::vector<std::complex<double>> mat(dim * dim);
    set_rand(mat.data(), dim * dim);
    LR_Util::matsym(mat.data(), dim);
    Parallel_2D pmat;
    LR_Util::setup_2d_division(pmat, 1, dim, dim);
    std::vector<std::complex<double>> hmat_local(pmat.get_local_size());
    LR_Util::set_local_from_global(pmat, mat.data(), hmat_local.data());
    std::vector<std::complex<double>> smat_local(pmat.get_local_size(), 0.0);
    for (int lj = 0;lj < pmat.get_col_size();++lj)
        for (int li = 0;li < pmat.get_row_size();++li)
            if (pmat.local2global_row(li) == pmat.local2global_col(lj))  // diagonal elements
                smat_local[li * pmat.get_row_size() + lj] = std::complex<double>(1.0, 0.0);

    // serial
    std::vector<double> eig(dim);
    LR_Util::diag_lapack(dim, mat.data(), eig.data());

    // parallel
    std::vector<double> eig_para(dim);
    std::vector<std::complex<double>> eigvec_para(pmat.get_local_size());
    LR_Util::diag_scalapack(dim, hmat_local.data(), smat_local.data(), eig_para.data(), eigvec_para.data(), pmat.desc);

    // compare
    check_double_eq(eig_para.data(), eig.data(), dim);
    std::vector<std::complex<double>> eigvec_serial_local(pmat.get_local_size());
    LR_Util::set_local_from_global(pmat, mat.data(), eigvec_serial_local.data());
    check_norm_eq(eigvec_para.data(), eigvec_serial_local.data(), pmat.get_local_size());
}

TEST(LR_Util, DiagScaLapackComplex)
{
    // setup the matrix
    const int dim = 15;
    std::vector<std::complex<double>> mat(dim * dim);
    set_rand(mat.data(), dim * dim);
    LR_Util::matsym(mat.data(), dim);
    Parallel_2D pmat;
    LR_Util::setup_2d_division(pmat, 1, dim, dim);
    std::vector<std::complex<double>> mat_local(pmat.get_local_size());
    LR_Util::set_local_from_global(pmat, mat.data(), mat_local.data());

    // serial
    std::vector<double> eig(dim);
    LR_Util::diag_lapack(dim, mat.data(), eig.data());

    // parallel
    std::vector<double> eig_para(dim);
    std::vector<std::complex<double>> eigvec_para(pmat.get_local_size());
    LR_Util::diag_scalapack(dim, mat_local.data(), eig_para.data(), eigvec_para.data(), pmat.desc);

    // compare
    check_double_eq(eig_para.data(), eig.data(), dim);
    std::vector<std::complex<double>> eigvec_serial_local(pmat.get_local_size());
    LR_Util::set_local_from_global(pmat, mat.data(), eigvec_serial_local.data());
    check_norm_eq(eigvec_para.data(), eigvec_serial_local.data(), pmat.get_local_size());
}

int main(int argc, char** argv)
{
    srand(time(NULL));  // for random number generator
    MPI_Init(&argc, &argv);
    testing::InitGoogleTest(&argc, argv);
    int result = RUN_ALL_TESTS();
    MPI_Finalize();
    return result;
}
#endif