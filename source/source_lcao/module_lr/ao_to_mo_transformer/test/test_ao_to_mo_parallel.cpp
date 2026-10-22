#include <gtest/gtest.h>
#include "../ao_to_mo.h"
#include "../cvcx.h"
#include "../../dm_trans/dm_trans.h"
#include "../../dm_trans/dm_diff.h"
#include "../../utils/lr_util.h"
#include <algorithm>
#include <complex>

#ifdef __MPI
namespace
{
template<typename T> T coefficient(const int ao, const int band)
{
    return T(0.001 * (ao + band + 1));
}
template<> std::complex<double> coefficient(const int ao, const int band)
{
    return std::complex<double>(0.001 * (ao + band + 1), 0.0007 * (ao - band));
}
double conjugate(const double value)
{
    return value;
}
std::complex<double> conjugate(const std::complex<double> value)
{
    return std::conj(value);
}
template <typename T> void check_empty_blocks(const int naos, const int nocc, const int nvirt, const int block)
{
    // The virtual window is smaller than the KS block; its first band is
    // also inside a block rather than at a block boundary.
    Parallel_2D pao;
    Parallel_2D pc;
    Parallel_2D px;
    LR_Util::setup_2d_division(pao, block, naos, naos);
    LR_Util::setup_2d_division(pc, block, naos, nocc + nvirt, pao.blacs_ctxt);
    LR_Util::setup_2d_division(px, block, nvirt, nocc, pao.blacs_ctxt);
    const int csize = std::max<int>(1, pc.get_local_size());
    const int xsize = std::max<int>(1, px.get_local_size());
    std::vector<T> coefficients(csize, T(0));
    std::vector<T> result(xsize, T(0));
    for (int col = 0; col < pc.get_col_size(); ++col)
    {
        const int band = pc.local2global_col(col);
        for (int row = 0; row < pc.get_row_size(); ++row)
        {
            const int ao = pc.local2global_row(row);
            coefficients[col * pc.get_row_size() + row] = coefficient<T>(ao, band);
        }
    }
    psi::Psi<T> c(coefficients.data(), 1, pc.get_col_size(), pc.get_row_size(), true);
    auto identity = LR_Util::newTensor<T>({pao.get_col_size(), pao.get_row_size()});
    identity.zero();
    for (int col = 0; col < pao.get_col_size(); ++col)
    {
        for (int row = 0; row < pao.get_row_size(); ++row)
        {
            if (pao.local2global_row(row) == pao.local2global_col(col))
            {
                identity.template data<T>()[col * pao.get_row_size() + row] = T(1);
            }
        }
    }
    std::vector<T> amplitudes(xsize, T(0));
    std::vector<T> full_amplitudes(nvirt * nocc);
    psi::Psi<T> full_c(1, nocc + nvirt, naos, naos, true);
    for (int band = 0; band < nocc + nvirt; ++band)
    {
        for (int ao = 0; ao < naos; ++ao)
        {
            full_c(0, band, ao) = coefficient<T>(ao, band);
        }
    }
    for (int occ = 0; occ < nocc; ++occ)
    {
        for (int virt = 0; virt < nvirt; ++virt)
        {
            full_amplitudes[occ * nvirt + virt] = T(0.002 * (occ + virt + 1));
        }
    }
    for (int col = 0; col < px.get_col_size(); ++col)
    {
        for (int row = 0; row < px.get_row_size(); ++row)
        {
            const int occ = px.local2global_col(col);
            const int virt = px.local2global_row(row);
            amplitudes[col * px.get_row_size() + row] = full_amplitudes[occ * nvirt + virt];
        }
    }
    auto full_identity = LR_Util::newTensor<T>({naos, naos});
    full_identity.zero();
    for (int ao = 0; ao < naos; ++ao)
    {
        full_identity.template data<T>()[ao * naos + ao] = T(1);
    }
    const std::vector<container::Tensor> full_potential = {full_identity};
    const auto spin_view = LR_Util::get_psi_spin(c, 0, 1);
    EXPECT_EQ(spin_view.get_nbands(), c.get_nbands());
    EXPECT_EQ(spin_view.get_nbasis(), c.get_nbasis());
    const std::vector<container::Tensor> potential = {identity};
    const T factor(1);
    T* const local_result = px.get_local_size() == 0 ? nullptr : result.data();
    const T* const local_amplitudes = px.get_local_size() == 0 ? nullptr : amplitudes.data();
    LR::ao_to_mo_pblas(potential, pao, c, pc, naos, nocc, nvirt, px,
        local_result, false, LR_Util::MO_TYPE::VO, factor);
    for (int col = 0; col < px.get_col_size(); ++col)
    {
        const int occ = px.local2global_col(col);
        for (int row = 0; row < px.get_row_size(); ++row)
        {
            const int virt = nocc + px.local2global_row(row);
            T expected(0);
            for (int ao = 0; ao < naos; ++ao)
            {
                expected += conjugate(coefficient<T>(ao, virt)) * coefficient<T>(ao, occ);
            }
            EXPECT_NEAR(std::abs(result[col * px.get_row_size() + row] - T(expected)), 0, 1e-9);
        }
    }
    std::vector<T> full_result(nvirt * nocc);
    for (const bool occupied : {true, false})
    {
        if (occupied)
        {
            LR::CVCX_occ_pblas(potential, pao, c, pc, local_amplitudes, px,
                naos, nocc, nvirt, local_result, false, factor);
            LR::CVCX_occ_blas(full_potential, full_c, full_amplitudes.data(),
                naos, nocc, nvirt, full_result.data(), false, factor);
        }
        else
        {
            LR::CVCX_virt_pblas(potential, pao, c, pc, local_amplitudes, px,
                naos, nocc, nvirt, local_result, false, factor);
            LR::CVCX_virt_blas(full_potential, full_c, full_amplitudes.data(),
                naos, nocc, nvirt, full_result.data(), false, factor);
        }
        for (int col = 0; col < px.get_col_size(); ++col)
        {
            for (int row = 0; row < px.get_row_size(); ++row)
            {
                const int global = px.local2global_col(col) * nvirt + px.local2global_row(row);
                const int local = col * px.get_row_size() + row;
                EXPECT_NEAR(std::abs(result[local] - full_result[global]), 0, 1e-9);
            }
        }
    }
    for (const bool transition : {true, false})
    {
        std::vector<container::Tensor> distributed;
        std::vector<container::Tensor> serial;
        if (transition)
        {
            distributed = LR::cal_dm_trans_pblas(local_amplitudes, px, c, pc,
                naos, nocc, nvirt, pao, factor, LR_Util::MO_TYPE::VO);
            serial = LR::cal_dm_trans_blas(full_amplitudes.data(), full_c,
                nocc, nvirt, factor, LR_Util::MO_TYPE::VO);
        }
        else
        {
            distributed = LR::cal_dm_diff_pblas(local_amplitudes, px, c, pc,
                naos, nocc, nvirt, pao, false, 1);
            serial = LR::cal_dm_diff_blas(full_amplitudes.data(), full_c,
                naos, nocc, nvirt, false, 1);
        }
        for (int col = 0; col < pao.get_col_size(); ++col)
        {
            for (int row = 0; row < pao.get_row_size(); ++row)
            {
                const int global = pao.local2global_col(col) * naos + pao.local2global_row(row);
                const int local = col * pao.get_row_size() + row;
                const T difference = distributed[0].template data<T>()[local] - serial[0].template data<T>()[global];
                EXPECT_NEAR(std::abs(difference), 0, 1e-9);
            }
        }
    }
}
}
TEST(AOToMO, EmptyLocalBlocksReal)
{
    check_empty_blocks<double>(16, 10, 3, 32);
    check_empty_blocks<double>(64, 20, 3, 32);
    check_empty_blocks<double>(512, 250, 10, 32);
}
TEST(AOToMO, EmptyLocalBlocksComplex)
{
    check_empty_blocks<std::complex<double>>(16, 10, 3, 32);
    check_empty_blocks<std::complex<double>>(64, 20, 3, 32);
    check_empty_blocks<std::complex<double>>(512, 250, 10, 32);
}
#endif
