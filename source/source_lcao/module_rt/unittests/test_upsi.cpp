#include "source_base/global_variable.h"
#include "source_base/matrix2d_view.h"
#include "source_base/parallel_blkcyc.h"
#include "source_lcao/module_rt/upsi.h"

#include <gtest/gtest.h>

#include <complex>
#include <vector>

extern int ictxt;

namespace
{
constexpr double threshold = 1.0e-8;
} // namespace

TEST(UpsiTest, testUpsi1)
{
    using Complex = std::complex<double>;
    const int nband = 2;
    const int nlocal = 2;
    Parallel_Blkcyclic matrix_layout;
    Parallel_Blkcyclic wfc_layout;
    ASSERT_EQ(matrix_layout.set(nlocal, nlocal, 1, ictxt), 0);
    ASSERT_EQ(wfc_layout.set(nlocal, nband, 1, ictxt), 0);

    std::vector<Complex> U_operator(matrix_layout.size(), Complex(1.0, 0.0));
    std::vector<Complex> psi_k_laststep(wfc_layout.size(), Complex(1.0, 0.0));
    std::vector<Complex> psi_k(wfc_layout.size(), Complex(0.0, 0.0));
    const Matrix2dView<const Complex> U_operator_view(
        U_operator.data(),
        matrix_layout.nrloc(),
        matrix_layout.ncloc(),
        matrix_layout.desc());
    const Matrix2dView<const Complex> psi_last_view(
        psi_k_laststep.data(),
        wfc_layout.nrloc(),
        wfc_layout.ncloc(),
        wfc_layout.desc());
    const Matrix2dView<Complex> psi_view(
        psi_k.data(),
        wfc_layout.nrloc(),
        wfc_layout.ncloc(),
        wfc_layout.desc());

    const bool print_matrix = false;
    module_rt::upsi(U_operator_view,
                    psi_last_view,
                    psi_view,
                    GlobalV::ofs_running,
                    print_matrix);

    for (int i = 0; i < psi_view.nrloc(); ++i)
    {
        for (int j = 0; j < psi_view.ncloc(); ++j)
        {
            EXPECT_NEAR(psi_view.at_local(i, j).real(), nlocal, threshold);
            EXPECT_NEAR(psi_view.at_local(i, j).imag(), 0.0, threshold);
        }
    }
}

TEST(UpsiTest, testUpsi2)
{
    using Complex = std::complex<double>;
    const int nband = 3;
    const int nlocal = 4;
    Parallel_Blkcyclic matrix_layout;
    Parallel_Blkcyclic wfc_layout;
    ASSERT_EQ(matrix_layout.set(nlocal, nlocal, 1, ictxt), 0);
    ASSERT_EQ(wfc_layout.set(nlocal, nband, 1, ictxt), 0);

    std::vector<Complex> U_operator(matrix_layout.size(), Complex(1.0, 0.0));
    std::vector<Complex> psi_k_laststep(wfc_layout.size(), Complex(2.0, 0.0));
    std::vector<Complex> psi_k(wfc_layout.size(), Complex(0.0, 0.0));
    const Matrix2dView<const Complex> U_operator_view(
        U_operator.data(),
        matrix_layout.nrloc(),
        matrix_layout.ncloc(),
        matrix_layout.desc());
    const Matrix2dView<const Complex> psi_last_view(
        psi_k_laststep.data(),
        wfc_layout.nrloc(),
        wfc_layout.ncloc(),
        wfc_layout.desc());
    const Matrix2dView<Complex> psi_view(
        psi_k.data(),
        wfc_layout.nrloc(),
        wfc_layout.ncloc(),
        wfc_layout.desc());

    const bool print_matrix = false;
    module_rt::upsi(U_operator_view,
                    psi_last_view,
                    psi_view,
                    GlobalV::ofs_running,
                    print_matrix);

    for (int i = 0; i < psi_view.nrloc(); ++i)
    {
        for (int j = 0; j < psi_view.ncloc(); ++j)
        {
            EXPECT_NEAR(psi_view.at_local(i, j).real(),
                        2.0 * nlocal,
                        threshold);
            EXPECT_NEAR(psi_view.at_local(i, j).imag(), 0.0, threshold);
        }
    }
}

TEST(UpsiTest, testUpsi3)
{
    using Complex = std::complex<double>;
    const int nband = 3;
    const int nlocal = 4;
    Parallel_Blkcyclic matrix_layout;
    Parallel_Blkcyclic wfc_layout;
    ASSERT_EQ(matrix_layout.set(nlocal, nlocal, 1, ictxt), 0);
    ASSERT_EQ(wfc_layout.set(nlocal, nband, 1, ictxt), 0);

    std::vector<Complex> U_operator(matrix_layout.size());
    std::vector<Complex> psi_k_laststep(wfc_layout.size());
    std::vector<Complex> psi_k(wfc_layout.size());
    Matrix2dView<Complex> U_operator_view(
        U_operator.data(),
        matrix_layout.nrloc(),
        matrix_layout.ncloc(),
        matrix_layout.desc());
    Matrix2dView<Complex> psi_last_view(
        psi_k_laststep.data(),
        wfc_layout.nrloc(),
        wfc_layout.ncloc(),
        wfc_layout.desc());
    Matrix2dView<Complex> psi_view(
        psi_k.data(),
        wfc_layout.nrloc(),
        wfc_layout.ncloc(),
        wfc_layout.desc());

    for (int i = 0; i < U_operator_view.nrloc(); ++i)
    {
        for (int j = 0; j < U_operator_view.ncloc(); ++j)
        {
            U_operator_view.at_local(i, j)
                = i == j ? Complex(1.0, 0.0) : Complex(0.0, 0.0);
        }
    }
    for (int i = 0; i < psi_last_view.nrloc(); ++i)
    {
        for (int j = 0; j < psi_last_view.ncloc(); ++j)
        {
            psi_last_view.at_local(i, j) = Complex(i * nband + j, 0.0);
            psi_view.at_local(i, j) = Complex(0.0, 0.0);
        }
    }

    const bool print_matrix = false;
    module_rt::upsi(U_operator_view,
                    psi_last_view,
                    psi_view,
                    GlobalV::ofs_running,
                    print_matrix);

    for (int i = 0; i < psi_view.nrloc(); ++i)
    {
        for (int j = 0; j < psi_view.ncloc(); ++j)
        {
            EXPECT_NEAR(psi_view.at_local(i, j).real(),
                        i * nband + j,
                        threshold);
            EXPECT_NEAR(psi_view.at_local(i, j).imag(), 0.0, threshold);
        }
    }
}
