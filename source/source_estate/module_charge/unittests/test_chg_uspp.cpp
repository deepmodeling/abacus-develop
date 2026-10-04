#include "gtest/gtest.h"

#include "source_estate/module_charge/chg_uspp.h"

#include <complex>
#include <vector>

/************************************************
 *  unit test of module_charge/chg_uspp.cpp
 ***********************************************/

/**
 * - Tested Functions:
 *   - split_dgrid: split dense reciprocal data into smooth and
 *     high-frequency parts on the USPP double grid
 *     - normal split with nspin=1 and nspin=2
 *     - boundary: npw_smooth == 0 (all high-frequency)
 *     - boundary: npw_dense == npw_smooth (no high-frequency)
 *     - multi-spin channel isolation
 *     - abort on invalid inputs (null pointer, bad nspin/npw, size mismatch)
 *   - merge_dgrid: merge smooth and high-frequency parts back into dense
 *     - round-trip with split_dgrid reproduces the original data
 *     - abort on invalid inputs
 */

class ChgUsppTest : public ::testing::Test
{
  protected:
    // build a dense buffer of shape [nspin * npw_dense] with distinct
    // per-element values so split/merge correctness is easy to verify
    static std::vector<std::complex<double>> make_dense(int nspin, int npw_dense)
    {
        std::vector<std::complex<double>> buf(nspin * npw_dense);
        for (int is = 0; is < nspin; ++is)
        {
            for (int ig = 0; ig < npw_dense; ++ig)
            {
                const double v = static_cast<double>(is * 1000 + ig);
                buf[is * npw_dense + ig] = std::complex<double>(v, v + 0.5);
            }
        }
        return buf;
    }
};

TEST_F(ChgUsppTest, SplitDgridNormalNspin1)
{
    const int nspin = 1;
    const int npw_smooth = 3;
    const int npw_dense = 5;
    auto data_d = make_dense(nspin, npw_dense);

    std::vector<std::complex<double>> data_s(nspin * npw_smooth);
    std::vector<std::complex<double>> data_hf(nspin * (npw_dense - npw_smooth));

    module_charge::split_dgrid(data_d.data(), data_s, data_hf, nspin, npw_smooth, npw_dense);

    // smooth part == first npw_smooth entries
    for (int ig = 0; ig < npw_smooth; ++ig)
    {
        EXPECT_EQ(data_s[ig], data_d[ig]);
    }
    // high-frequency part == remaining entries
    for (int ig = 0; ig < npw_dense - npw_smooth; ++ig)
    {
        EXPECT_EQ(data_hf[ig], data_d[npw_smooth + ig]);
    }
}

TEST_F(ChgUsppTest, SplitDgridNormalNspin2)
{
    const int nspin = 2;
    const int npw_smooth = 2;
    const int npw_dense = 4;
    auto data_d = make_dense(nspin, npw_dense);

    std::vector<std::complex<double>> data_s(nspin * npw_smooth);
    std::vector<std::complex<double>> data_hf(nspin * (npw_dense - npw_smooth));

    module_charge::split_dgrid(data_d.data(), data_s, data_hf, nspin, npw_smooth, npw_dense);

    // each spin channel is split independently
    for (int is = 0; is < nspin; ++is)
    {
        for (int ig = 0; ig < npw_smooth; ++ig)
        {
            EXPECT_EQ(data_s[is * npw_smooth + ig], data_d[is * npw_dense + ig]);
        }
        for (int ig = 0; ig < npw_dense - npw_smooth; ++ig)
        {
            EXPECT_EQ(data_hf[is * (npw_dense - npw_smooth) + ig],
                      data_d[is * npw_dense + npw_smooth + ig]);
        }
    }
}

TEST_F(ChgUsppTest, SplitDgridSmoothIsZero)
{
    // npw_smooth == 0: the whole dense buffer is high-frequency
    const int nspin = 1;
    const int npw_smooth = 0;
    const int npw_dense = 3;
    auto data_d = make_dense(nspin, npw_dense);

    std::vector<std::complex<double>> data_s(0);
    std::vector<std::complex<double>> data_hf(nspin * npw_dense);

    module_charge::split_dgrid(data_d.data(), data_s, data_hf, nspin, npw_smooth, npw_dense);

    EXPECT_TRUE(data_s.empty());
    for (int ig = 0; ig < npw_dense; ++ig)
    {
        EXPECT_EQ(data_hf[ig], data_d[ig]);
    }
}

TEST_F(ChgUsppTest, SplitDgridDenseEqualsSmooth)
{
    // npw_dense == npw_smooth: no high-frequency tail, data_hf is empty
    const int nspin = 2;
    const int npw_smooth = 3;
    const int npw_dense = 3;
    auto data_d = make_dense(nspin, npw_dense);

    std::vector<std::complex<double>> data_s(nspin * npw_smooth);
    std::vector<std::complex<double>> data_hf(0);

    module_charge::split_dgrid(data_d.data(), data_s, data_hf, nspin, npw_smooth, npw_dense);

    EXPECT_TRUE(data_hf.empty());
    for (int i = 0; i < nspin * npw_smooth; ++i)
    {
        EXPECT_EQ(data_s[i], data_d[i]);
    }
}

TEST_F(ChgUsppTest, MergeDgridRoundTripNspin1)
{
    const int nspin = 1;
    const int npw_smooth = 3;
    const int npw_dense = 5;
    auto data_d = make_dense(nspin, npw_dense);

    std::vector<std::complex<double>> data_s(nspin * npw_smooth);
    std::vector<std::complex<double>> data_hf(nspin * (npw_dense - npw_smooth));
    module_charge::split_dgrid(data_d.data(), data_s, data_hf, nspin, npw_smooth, npw_dense);

    std::vector<std::complex<double>> merged(nspin * npw_dense);
    module_charge::merge_dgrid(merged.data(), data_s, data_hf, nspin, npw_smooth, npw_dense);

    for (int i = 0; i < nspin * npw_dense; ++i)
    {
        EXPECT_EQ(merged[i], data_d[i]);
    }
}

TEST_F(ChgUsppTest, MergeDgridRoundTripNspin2)
{
    const int nspin = 2;
    const int npw_smooth = 2;
    const int npw_dense = 5;
    auto data_d = make_dense(nspin, npw_dense);

    std::vector<std::complex<double>> data_s(nspin * npw_smooth);
    std::vector<std::complex<double>> data_hf(nspin * (npw_dense - npw_smooth));
    module_charge::split_dgrid(data_d.data(), data_s, data_hf, nspin, npw_smooth, npw_dense);

    std::vector<std::complex<double>> merged(nspin * npw_dense);
    module_charge::merge_dgrid(merged.data(), data_s, data_hf, nspin, npw_smooth, npw_dense);

    for (int i = 0; i < nspin * npw_dense; ++i)
    {
        EXPECT_EQ(merged[i], data_d[i]);
    }
}

TEST_F(ChgUsppTest, SplitDgridNullDataAborts)
{
    const int nspin = 1;
    const int npw_smooth = 2;
    const int npw_dense = 4;
    std::vector<std::complex<double>> data_s(nspin * npw_smooth);
    std::vector<std::complex<double>> data_hf(nspin * (npw_dense - npw_smooth));
    EXPECT_DEATH(module_charge::split_dgrid(nullptr, data_s, data_hf, nspin, npw_smooth, npw_dense),
                 "");
}

TEST_F(ChgUsppTest, SplitDgridBadNspinAborts)
{
    const int nspin = 0;
    const int npw_smooth = 2;
    const int npw_dense = 4;
    auto data_d = make_dense(1, npw_dense);
    std::vector<std::complex<double>> data_s(npw_smooth);
    std::vector<std::complex<double>> data_hf(npw_dense - npw_smooth);
    EXPECT_DEATH(module_charge::split_dgrid(data_d.data(), data_s, data_hf, nspin, npw_smooth, npw_dense),
                 "");
}

TEST_F(ChgUsppTest, SplitDgridBadNpwAborts)
{
    const int nspin = 1;
    // npw_dense < npw_smooth is invalid
    const int npw_smooth = 5;
    const int npw_dense = 3;
    auto data_d = make_dense(nspin, npw_dense);
    std::vector<std::complex<double>> data_s(nspin * npw_smooth);
    std::vector<std::complex<double>> data_hf(nspin * (npw_smooth - npw_dense));
    EXPECT_DEATH(module_charge::split_dgrid(data_d.data(), data_s, data_hf, nspin, npw_smooth, npw_dense),
                 "");
}

TEST_F(ChgUsppTest, SplitDgridSizeMismatchAborts)
{
    const int nspin = 1;
    const int npw_smooth = 2;
    const int npw_dense = 4;
    auto data_d = make_dense(nspin, npw_dense);
    // data_s is too small
    std::vector<std::complex<double>> data_s(1);
    std::vector<std::complex<double>> data_hf(npw_dense - npw_smooth);
    EXPECT_DEATH(module_charge::split_dgrid(data_d.data(), data_s, data_hf, nspin, npw_smooth, npw_dense),
                 "");
}

TEST_F(ChgUsppTest, MergeDgridNullDataAborts)
{
    const int nspin = 1;
    const int npw_smooth = 2;
    const int npw_dense = 4;
    std::vector<std::complex<double>> data_s(nspin * npw_smooth);
    std::vector<std::complex<double>> data_hf(nspin * (npw_dense - npw_smooth));
    EXPECT_DEATH(module_charge::merge_dgrid(nullptr, data_s, data_hf, nspin, npw_smooth, npw_dense),
                 "");
}
