#ifndef ABACUS_LR_STATE_OVLP_TEST_H
#define ABACUS_LR_STATE_OVLP_TEST_H

#include "source_base/parallel_2d.h"
#include "source_base/parallel_reduce.h"
#include "source_lcao/module_lr/state_ovlp.h"
#include <gtest/gtest.h>
#include <algorithm>
#include <cmath>
#include <complex>
#include <vector>

// Independent contraction oracles and distribution helpers for the overlap tests.
namespace
{
template <typename T>
std::vector<T> overlap_block(const std::vector<T>& full, int bands, int start, int count)
{
    std::vector<T> block(count * count);
    for (int i = 0; i < count; ++i)
        for (int j = 0; j < count; ++j)
            block[i * count + j] = full[(start + i) * bands + start + j];
    return block;
}

void setup(Parallel_2D& target, int rows, int cols, const Parallel_2D& grid)
{
#ifdef __MPI
    const int block_size = grid.get_block_size();
    target.set(rows, cols, block_size, grid.blacs_ctxt);
#else
    target.set_serial(rows, cols);
#endif
}

template <typename T>
std::vector<T> scatter(const std::vector<T>& full, const Parallel_2D& layout)
{
    std::vector<T> local(layout.get_local_size(), T(0));
    const int columns = layout.get_global_col_size();
    for (int col = 0; col < layout.get_col_size(); ++col)
    {
        const int global_col = layout.local2global_col(col);
        for (int row = 0; row < layout.get_row_size(); ++row)
        {
            const int global_row = layout.local2global_row(row);
            local[col * layout.get_row_size() + row] = full[global_row * columns + global_col];
        }
    }
    return local;
}

// Independent test oracles: direct formula contractions, with no BLAS/PBLAS calls.
template <typename T>
T conjugate(T value)
{
    return value;
}
template <>
std::complex<double> conjugate(std::complex<double> value)
{
    return std::conj(value);
}

template <typename T>
std::vector<T> reference_mo(const std::vector<T>& old, const std::vector<double>& ao,
    const std::vector<T>& current, int naos, int bands)
{
    std::vector<T> result(bands * bands, T(0));
    for (int i = 0; i < bands; ++i)
    {
        for (int j = 0; j < bands; ++j)
        {
            for (int mu = 0; mu < naos; ++mu)
            {
                for (int nu = 0; nu < naos; ++nu)
                {
                    result[i * bands + j] += conjugate(old[i * naos + mu])
                        * ao[mu * naos + nu] * current[j * naos + nu];
                }
            }
        }
    }
    return result;
}

template <typename T>
std::vector<T> reference_projection(const std::vector<T>& previous,
    const std::vector<T>& overlap, int no, int nv)
{
    const int bands = no + nv;
    std::vector<T> result(no * nv, T(0));
    for (int j = 0; j < no; ++j)
    {
        for (int b = 0; b < nv; ++b)
        {
            for (int i = 0; i < no; ++i)
            {
                for (int a = 0; a < nv; ++a)
                {
                    result[j * nv + b] += overlap[i * bands + j] * previous[i * nv + a]
                        * conjugate(overlap[(no + a) * bands + no + b]);
                }
            }
        }
    }
    return result;
}

template <typename T>
std::vector<T> gather(const std::vector<T>& local, const Parallel_2D& layout)
{
    const int columns = layout.get_global_col_size();
    const int count = layout.get_global_row_size() * columns;
    std::vector<T> result(count, T(0));
    EXPECT_EQ(local.size(), layout.get_local_size());
    for (int col = 0; col < layout.get_col_size(); ++col)
    {
        const int global_col = layout.local2global_col(col);
        for (int row = 0; row < layout.get_row_size(); ++row)
        {
            const int global_row = layout.local2global_row(row);
            result[global_row * columns + global_col] = local[col * layout.get_row_size() + row];
        }
    }
    Parallel_Reduce::reduce_all(result.data(), count);
    return result;
}

template <typename T>
std::vector<T> transform_for_test(const std::vector<T>& old, const std::vector<double>& ao,
    const std::vector<T>& current, int naos, int bands, int start, int count, int block_size)
{
    Parallel_2D pc;
#ifdef __MPI
    pc.init(naos, bands, block_size, MPI_COMM_WORLD);
#else
    pc.set_serial(naos, bands);
#endif
    Parallel_2D ps;
    Parallel_2D po;
    setup(ps, naos, naos, pc);
    setup(po, count, count, pc);
    std::vector<T> old_rows(old.size());
    std::vector<T> current_rows(current.size());
    for (int i = 0; i < naos; ++i)
    {
        for (int band = 0; band < bands; ++band)
        {
            old_rows[i * bands + band] = old[band * naos + i];
            current_rows[i * bands + band] = current[band * naos + i];
        }
    }
    const std::vector<T> typed_ao(ao.begin(), ao.end());
    const auto local_old = scatter(old_rows, pc);
    const auto local_current = scatter(current_rows, pc);
    const auto local_ao = scatter(typed_ao, ps);
    const auto actual = LR::mo_overlap_dist(local_old, local_ao, local_current, pc, ps, po, start);
    return gather(actual, po);
}

template <typename T>
std::vector<T> project(const std::vector<T>& previous, const std::vector<T>& overlap, int no, int nv)
{
    const int bands = no + nv;
    Parallel_2D po;
#ifdef __MPI
    po.init(bands, bands, 1, MPI_COMM_WORLD);
#else
    po.set_serial(bands, bands);
#endif
    Parallel_2D px;
    setup(px, nv, no, po);
    std::vector<T> previous_rows(previous.size());
    for (int i = 0; i < no; ++i)
    {
        for (int a = 0; a < nv; ++a)
        {
            previous_rows[a * no + i] = previous[i * nv + a];
        }
    }
    const auto local_previous = scatter(previous_rows, px);
    Parallel_2D pi;
    Parallel_2D pv;
    setup(pi, no, no, po);
    setup(pv, nv, nv, po);
    const auto occupied = overlap_block(overlap, bands, 0, no);
    const auto virtuals = overlap_block(overlap, bands, no, nv);
    const auto local_occupied = scatter(occupied, pi);
    const auto local_virtuals = scatter(virtuals, pv);
    const auto actual = LR::project_reference_dist(local_previous, local_occupied, local_virtuals, pi, pv, px);
    const auto rows = gather(actual, px);
    std::vector<T> result(previous.size());
    for (int i = 0; i < no; ++i)
    {
        for (int a = 0; a < nv; ++a)
        {
            result[i * nv + a] = rows[a * no + i];
        }
    }
    return result;
}

template <typename T>
void check_mo(int naos, int bands, T phase, int start, int count, int block_size)
{
    std::vector<double> ao(naos * naos);
    std::vector<T> old(naos * bands);
    std::vector<T> current(naos * bands);
    for (int i = 0; i < naos; ++i)
    {
        for (int j = 0; j < naos; ++j)
        {
            ao[i * naos + j] = (2 * i - j + 1) / 7.0;
        }
        for (int band = 0; band < bands; ++band)
        {
            old[band * naos + i] = phase * ((i + band + 1) / 11.0);
            current[band * naos + i] = T((i - 2 * band + 2) / 13.0);
        }
    }
    const auto actual = transform_for_test(old, ao, current, naos, bands, start, count, block_size);
    const auto full = reference_mo(old, ao, current, naos, bands);
    const auto expected = overlap_block(full, bands, start, count);
    ASSERT_EQ(actual.size(), expected.size());
    for (std::size_t i = 0; i < actual.size(); ++i)
    {
        EXPECT_NEAR(std::abs(actual[i] - expected[i]), 0.0, 1e-13);
    }
}
template <typename T>
void check_projection(int no, int nv, T phase, int block_size)
{
    const int bands = no + nv;
    Parallel_2D po;
#ifdef __MPI
    po.init(bands, bands, block_size, MPI_COMM_WORLD);
#else
    po.set_serial(bands, bands);
#endif
    Parallel_2D px;
    setup(px, nv, no, po);
    std::vector<T> overlap(bands * bands);
    std::vector<T> previous(no * nv);
    for (int i = 0; i < bands; ++i)
    {
        for (int j = 0; j < bands; ++j)
        {
            overlap[i * bands + j] = T((i - 2 * j + 1) / 17.0) + phase * ((i + j + 1) / 19.0);
        }
    }
    for (int i = 0; i < no * nv; ++i)
    {
        previous[i] = phase * ((i + 1) / 23.0);
    }
    // Convert occupied-major X into row-major A=X^T before scattering into px.
    std::vector<T> previous_rows(previous.size());
    for (int i = 0; i < no; ++i)
    {
        for (int a = 0; a < nv; ++a)
        {
            previous_rows[a * no + i] = previous[i * nv + a];
        }
    }
    const auto local_previous = scatter(previous_rows, px);
    Parallel_2D pi;
    Parallel_2D pv;
    setup(pi, no, no, po);
    setup(pv, nv, nv, po);
    const auto occupied = overlap_block(overlap, bands, 0, no);
    const auto virtuals = overlap_block(overlap, bands, no, nv);
    const auto local_occupied = scatter(occupied, pi);
    const auto local_virtuals = scatter(virtuals, pv);
    const auto actual = LR::project_reference_dist(local_previous, local_occupied, local_virtuals, pi, pv, px);
    const auto expected = reference_projection(previous, overlap, no, nv);
    std::vector<T> expected_rows(expected.size());
    for (int i = 0; i < no; ++i)
    {
        for (int a = 0; a < nv; ++a)
        {
            expected_rows[a * no + i] = expected[i * nv + a];
        }
    }
    const auto local_expected = scatter(expected_rows, px);
    ASSERT_EQ(actual.size(), local_expected.size());
    for (std::size_t i = 0; i < actual.size(); ++i)
    {
        // Keep the original unit-block tolerance; larger new oracle contractions
        // allow roundoff proportional to the independently computed result.
        const double scaled_tolerance = 1e-12 * std::max(1.0, std::abs(local_expected[i]));
        const double tolerance = block_size == 1 ? 1e-13 : scaled_tolerance;
        EXPECT_NEAR(std::abs(actual[i] - local_expected[i]), 0.0, tolerance);
    }
}

}

#endif
