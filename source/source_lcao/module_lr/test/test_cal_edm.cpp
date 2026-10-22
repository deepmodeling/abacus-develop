#include <gtest/gtest.h>
#include "../cal_edm.h"
#include <complex>
#include <vector>

namespace
{
template<typename T> T value(double real, double imag);
template <>
double value<double>(double real, double)
{
    return real;
}
template<> std::complex<double> value<std::complex<double>>(double real, double imag)
{
    return std::complex<double>(real, imag);
}
double conjugate(double x)
{
    return x;
}
std::complex<double> conjugate(std::complex<double> x)
{
    return std::conj(x);
}

template<typename T> T coefficient(int row, int band)
{
    const double real = 0.3 * (row + 1) - 0.2 * band;
    const double imag = 0.1 * (row - band);
    return value<T>(real, imag);
}
template<typename T> T amplitude(int virt, int occ)
{
    const double real = 0.4 * (virt + 1) + 0.15 * occ;
    const double imag = 0.2 * virt - 0.3 * (occ + 1);
    return value<T>(real, imag);
}

template<typename T> void check_contractions()
{
    const int naos = 5;
    const int nocc = 2;
    const int nvirt = 3;
    const int nbands = nocc + nvirt;
    const int block = 1;
    Parallel_2D pc;
    Parallel_2D px;
    Parallel_2D pao;
    Parallel_2D pmat;
    LR_Util::setup_2d_division(pc, block, naos, nbands);
#ifdef __MPI
    LR_Util::setup_2d_division(px, block, nvirt, nocc, pc.blacs_ctxt);
    LR_Util::setup_2d_division(pao, block, naos, nocc, pc.blacs_ctxt);
    LR_Util::setup_2d_division(pmat, block, naos, naos, pc.blacs_ctxt);
#else
    LR_Util::setup_2d_division(px, block, nvirt, nocc);
    LR_Util::setup_2d_division(pao, block, naos, nocc);
    LR_Util::setup_2d_division(pmat, block, naos, naos);
#endif
    std::vector<T> c(pc.get_local_size());
    std::vector<T> x(px.get_local_size());
    for (int col = 0; col < pc.get_col_size(); ++col)
    {
        const int band = pc.local2global_col(col);
        for (int row = 0; row < pc.get_row_size(); ++row)
        {
            const int ao = pc.local2global_row(row);
            c[col * pc.get_row_size() + row] = coefficient<T>(ao, band);
        }
    }
    for (int col = 0; col < px.get_col_size(); ++col)
    {
        const int occ = px.local2global_col(col);
        for (int row = 0; row < px.get_row_size(); ++row)
        {
            const int virt = px.local2global_row(row);
            x[col * px.get_row_size() + row] = amplitude<T>(virt, occ);
        }
    }
    std::vector<T> ao_occ(pao.get_local_size(), T(100.0));
    LR::cal_X_ao_occ(x.data(), px, c.data(), pc, ao_occ.data(), pao);
    for (int col = 0; col < pao.get_col_size(); ++col)
    {
        const int occ = pao.local2global_col(col);
        for (int row = 0; row < pao.get_row_size(); ++row)
        {
            const int ao = pao.local2global_row(row);
            T expected = T(0.0);
            for (int a = 0; a < nvirt; ++a)
            {
                expected += coefficient<T>(ao, nocc + a) * amplitude<T>(a, occ);
            }
            const T actual = ao_occ[col * pao.get_row_size() + row];
            EXPECT_NEAR(std::abs(actual - expected), 0.0, 1e-12);
        }
    }
    // Different left/right factors detect both operand order and conjugation mistakes.
    std::vector<T> right(pao.get_local_size());
    for (int col = 0; col < pao.get_col_size(); ++col)
    {
        const int band = pao.local2global_col(col);
        for (int row = 0; row < pao.get_row_size(); ++row)
        {
            const int ao = pao.local2global_row(row);
            right[col * pao.get_row_size() + row] = coefficient<T>(ao, band);
        }
    }
    std::vector<T> dm(pmat.get_local_size(), T(100.0));
    LR::matdot(ao_occ.data(), right.data(), pao, dm.data(), pmat);
    for (int col = 0; col < pmat.get_col_size(); ++col)
    {
        const int nu = pmat.local2global_col(col);
        for (int row = 0; row < pmat.get_row_size(); ++row)
        {
            const int mu = pmat.local2global_row(row);
            T expected = T(0.0);
            for (int i = 0; i < nocc; ++i)
            {
                T left = T(0.0);
                for (int a = 0; a < nvirt; ++a)
                {
                    left += coefficient<T>(mu, nocc + a) * amplitude<T>(a, i);
                }
                expected += left * conjugate(coefficient<T>(nu, i));
            }
            const T actual = dm[col * pmat.get_row_size() + row];
            EXPECT_NEAR(std::abs(actual - expected), 0.0, 1e-12);
        }
    }
}
}

TEST(CalEdm, RealColumnMajorContractions)
{
    check_contractions<double>();
}
TEST(CalEdm, ComplexAdjointContractions)
{
    check_contractions<std::complex<double>>();
}
