#include "source_lcao/module_ri/sternheimer_fd_nonlocal_projector.h"

#include <complex>
#include <gtest/gtest.h>
#include <stdexcept>
#include <vector>

namespace
{

using Projector = ModuleRI::SternheimerFDNonlocalProjector;
using Complex = Projector::Complex;
using Vector = Projector::Vector;

} // namespace

TEST(SternheimerFDNonlocalProjector, SingleProjectorActsAsLowRankPotential)
{
    const double volume_element = 0.5;
    const Vector beta = {Complex(1.0, 0.0), Complex(1.0, 0.0)};
    Projector::ProjectorBlock block;
    block.projectors = {beta};
    block.d_matrix = {{Complex(3.0, 0.0)}};

    Projector projector(2, volume_element, {block});

    Vector output;
    projector.apply(beta, output);

    ASSERT_EQ(output.size(), beta.size());
    EXPECT_NEAR(output[0].real(), 3.0, 1.0e-12);
    EXPECT_NEAR(output[0].imag(), 0.0, 1.0e-12);
    EXPECT_NEAR(output[1].real(), 3.0, 1.0e-12);
    EXPECT_NEAR(output[1].imag(), 0.0, 1.0e-12);

    const Vector orthogonal = {Complex(1.0, 0.0), Complex(-1.0, 0.0)};
    projector.apply(orthogonal, output);
    EXPECT_NEAR(std::abs(output[0]), 0.0, 1.0e-12);
    EXPECT_NEAR(std::abs(output[1]), 0.0, 1.0e-12);
}

TEST(SternheimerFDNonlocalProjector, CoupledProjectorsUseDMatrix)
{
    const double volume_element = 1.0;
    Projector::ProjectorBlock block;
    block.projectors = {{Complex(1.0, 0.0), Complex(0.0, 0.0)}, {Complex(0.0, 0.0), Complex(1.0, 0.0)}};
    block.d_matrix = {{Complex(2.0, 0.0), Complex(0.5, -0.25)}, {Complex(0.5, 0.25), Complex(4.0, 0.0)}};

    Projector projector(2, volume_element, {block});
    const Vector psi = {Complex(1.0, 2.0), Complex(-3.0, 0.5)};

    Vector output;
    projector.apply(psi, output);

    ASSERT_EQ(output.size(), psi.size());
    EXPECT_NEAR(output[0].real(), 2.0 * 1.0 + (0.5 * -3.0 + 0.25 * 0.5), 1.0e-12);
    EXPECT_NEAR(output[0].imag(), 2.0 * 2.0 + (0.5 * 0.5 - 0.25 * -3.0), 1.0e-12);
    EXPECT_NEAR(output[1].real(), (0.5 * 1.0 - 0.25 * 2.0) + 4.0 * -3.0, 1.0e-12);
    EXPECT_NEAR(output[1].imag(), (0.5 * 2.0 + 0.25 * 1.0) + 4.0 * 0.5, 1.0e-12);
}

TEST(SternheimerFDNonlocalProjector, AddToAccumulatesOnExistingVector)
{
    const double volume_element = 1.0;
    Projector::ProjectorBlock block;
    block.projectors = {{Complex(1.0, 0.0), Complex(0.0, 0.0)}};
    block.d_matrix = {{Complex(2.0, 0.0)}};
    Projector projector(2, volume_element, {block});

    const Vector psi = {Complex(3.0, 0.0), Complex(4.0, 0.0)};
    Vector output = {Complex(10.0, 0.0), Complex(20.0, 0.0)};
    projector.add_to(psi, output);

    EXPECT_NEAR(output[0].real(), 16.0, 1.0e-12);
    EXPECT_NEAR(output[0].imag(), 0.0, 1.0e-12);
    EXPECT_NEAR(output[1].real(), 20.0, 1.0e-12);
    EXPECT_NEAR(output[1].imag(), 0.0, 1.0e-12);
}

TEST(SternheimerFDNonlocalProjector, SparseBatchMatchesRepeatedScalarApplications)
{
    const double volume_element = 0.25;
    Projector::ProjectorBlock block;
    block.projectors = {{Complex(0.0, 0.0),
                         Complex(1.0, -0.5),
                         Complex(0.0, 0.0),
                         Complex(0.0, 0.0),
                         Complex(-0.25, 0.75),
                         Complex(0.0, 0.0)},
                        {Complex(0.5, 0.25),
                         Complex(0.0, 0.0),
                         Complex(0.0, 0.0),
                         Complex(-1.0, 0.0),
                         Complex(0.0, 0.0),
                         Complex(0.0, 0.0)}};
    block.d_matrix = {{Complex(2.0, 0.0), Complex(0.25, -0.5)}, {Complex(0.25, 0.5), Complex(3.0, 0.0)}};
    Projector projector(6, volume_element, {block});

    EXPECT_EQ(projector.projector_value_count(), 12U);
    EXPECT_EQ(projector.nonzero_projector_value_count(), 4U);

    const Projector::Matrix input = {{Complex(1.0, 0.0),
                                      Complex(2.0, -1.0),
                                      Complex(-3.0, 0.5),
                                      Complex(0.0, 2.0),
                                      Complex(1.5, -0.25),
                                      Complex(-0.5, 0.0)},
                                     {Complex(-1.0, 1.0),
                                      Complex(0.0, 0.5),
                                      Complex(2.0, 0.0),
                                      Complex(1.0, -1.0),
                                      Complex(-2.0, 0.75),
                                      Complex(0.25, 0.5)},
                                     {Complex(0.5, -0.5),
                                      Complex(-1.0, 0.0),
                                      Complex(0.0, 0.0),
                                      Complex(3.0, 0.25),
                                      Complex(0.75, 1.0),
                                      Complex(-2.0, -1.0)}};

    Projector::Matrix batch_output;
    projector.apply_batch(input, batch_output);
    ASSERT_EQ(batch_output.size(), input.size());
    for (std::size_t column = 0; column != input.size(); ++column)
    {
        Vector scalar_output;
        projector.apply(input[column], scalar_output);
        ASSERT_EQ(batch_output[column].size(), scalar_output.size());
        for (std::size_t ir = 0; ir != scalar_output.size(); ++ir)
        {
            EXPECT_NEAR(std::abs(batch_output[column][ir] - scalar_output[ir]), 0.0, 1.0e-13);
        }
    }

    Projector::Matrix accumulated(input.size(), Vector(6, Complex(0.75, -0.25)));
    projector.add_to_batch(input, accumulated);
    for (std::size_t column = 0; column != input.size(); ++column)
    {
        Vector scalar_accumulated(6, Complex(0.75, -0.25));
        projector.add_to(input[column], scalar_accumulated);
        for (std::size_t ir = 0; ir != scalar_accumulated.size(); ++ir)
        {
            EXPECT_NEAR(std::abs(accumulated[column][ir] - scalar_accumulated[ir]), 0.0, 1.0e-13);
        }
    }
}

TEST(SternheimerFDNonlocalProjector, RejectsInvalidBlocks)
{
    Projector::ProjectorBlock block;
    block.projectors = {{Complex(1.0, 0.0), Complex(0.0, 0.0)}};
    block.d_matrix = {{Complex(1.0, 0.0)}};

    EXPECT_THROW(Projector(2, 0.0, {block}), std::invalid_argument);

    block.projectors[0].resize(1);
    EXPECT_THROW(Projector(2, 1.0, {block}), std::invalid_argument);

    block.projectors = {{Complex(1.0, 0.0), Complex(0.0, 0.0)}};
    block.d_matrix = {};
    EXPECT_THROW(Projector(2, 1.0, {block}), std::invalid_argument);
}
