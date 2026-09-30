#include "source_io/module_wf/exx_source_io.h"
#include "source_base/matrix.h"
#include "source_cell/klist.h"

#include "gtest/gtest.h"

#include <sstream>
#include <stdexcept>
#include <string>

namespace
{
std::string checkpoint(int nspin)
{
    const int nks = nspin == 2 ? 4 : 2;
    K_Vectors points;
    points.kvec_d.resize(nks);
    points.wk.resize(nks);
    points.isk.resize(nks);
    points.nmp[0] = 2;
    points.nmp[1] = 1;
    points.nmp[2] = 1;
    ModuleBase::matrix weights(nks, 3);
    for (int iq = 0; iq < nks; ++iq)
    {
        points.kvec_d[iq].x = 0.5 * (iq % 2);
        points.wk[iq] = nspin == 1 ? 1.0 : 0.5;
        points.isk[iq] = iq / 2;
        if (nspin == 1)
        {
            points.isk[iq] = 0;
        }
        weights(iq, 0) = points.wk[iq];
        weights(iq, 1) = points.wk[iq] * 0.125;
    }
    std::ostringstream out;
    ModuleIO::write_exx_source(out, points, weights, nspin, 15, 15, 15, "HSE:10:40");
    return out.str();
}
}

TEST(ExxSourceIO, PreservesFractionalOccupationsAndSpinOrder)
{
    for (int nspin = 1; nspin <= 2; ++nspin)
    {
        const std::string saved = checkpoint(nspin);
        std::istringstream in(saved);
        K_Vectors points;
        ModuleBase::matrix weights;
        ModuleIO::read_exx_source(in, points, weights, nspin, 15, 15, 15, "HSE:10:40");
        const int nks = nspin == 1 ? 2 : 4;
        EXPECT_EQ(points.get_nks(), nks);
        EXPECT_EQ(points.get_nkstot_nospin(), 2);
        EXPECT_EQ(points.get_spin_mult(), nspin);
        EXPECT_EQ(weights.nc, 3);
        EXPECT_DOUBLE_EQ(weights(1, 1), points.wk[1] * 0.125);
        EXPECT_DOUBLE_EQ(weights(1, 2), 0.0);
        EXPECT_DOUBLE_EQ(points.kvec_d[1].x, 0.5);
        EXPECT_EQ(points.isk.back(), nspin - 1);
    }
}

TEST(ExxSourceIO, RejectsIncompatibleConfigurationAndGrid)
{
    const std::string saved = checkpoint(1);
    K_Vectors points;
    ModuleBase::matrix weights;
    std::istringstream changed_functional(saved);
    EXPECT_THROW(ModuleIO::read_exx_source(changed_functional, points, weights, 1, 15, 15, 15, "PBE0:10:40"), std::runtime_error);
    std::istringstream changed_grid(saved);
    EXPECT_THROW(ModuleIO::read_exx_source(changed_grid, points, weights, 1, 16, 15, 15, "HSE:10:40"), std::runtime_error);
    std::istringstream changed_spin(saved);
    EXPECT_THROW(ModuleIO::read_exx_source(changed_spin, points, weights, 2, 15, 15, 15, "HSE:10:40"), std::runtime_error);
}

TEST(ExxSourceIO, RejectsTruncationInvalidOccupationsAndTrailingData)
{
    const std::string saved = checkpoint(1);
    K_Vectors points;
    ModuleBase::matrix weights;
    const std::string truncated_text = saved.substr(0, saved.size() - 6);
    std::istringstream truncated(truncated_text);
    EXPECT_THROW(ModuleIO::read_exx_source(truncated, points, weights, 1, 15, 15, 15, "HSE:10:40"), std::runtime_error);
    std::istringstream trailing(saved + "unexpected");
    EXPECT_THROW(ModuleIO::read_exx_source(trailing, points, weights, 1, 15, 15, 15, "HSE:10:40"), std::runtime_error);
    std::string invalid = saved;
    const auto occupation = invalid.find("0.125");
    ASSERT_NE(occupation, std::string::npos);
    invalid.replace(occupation, 5, "-0.25");
    std::istringstream negative(invalid);
    EXPECT_THROW(ModuleIO::read_exx_source(negative, points, weights, 1, 15, 15, 15, "HSE:10:40"), std::runtime_error);
}
