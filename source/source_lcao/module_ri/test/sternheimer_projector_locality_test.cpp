#include "source_lcao/module_ri/sternheimer_rpa.h"
#include <gtest/gtest.h>
#include <cmath>

TEST(SternheimerProjectorLocality, LargeComplexBorrowedViewsMatchExplicitProjection)
{
    using Projector = ModuleRI::SternheimerSubspaceProjector;
    using Vector = Projector::Vector;
    using Complex = Projector::Complex;
    constexpr int grid = 257;
    constexpr double weight = 0.125;
    const double pi = std::acos(-1.0);
    for (const int count : {20, 196, 216})
    {
        std::vector<Vector> basis(count, Vector(grid));
        for (int b = 0; b < count-1; ++b)
            for (int g = 0; g < grid; ++g)
            {
                const double phase = 2.0*pi*(b+1)*g/grid;
                basis[b][g] = Complex(std::cos(phase), std::sin(phase));
            }
        std::vector<const Vector*> views;
        for (const auto& vector : basis) views.push_back(&vector);
        const Projector projector(views, weight, true);
        const auto retained = projector.storage_bytes();
        for (const int width : {1, 2, 3})
        {
            std::vector<Vector> input(width, Vector(grid));
            for (int c = 0; c < width; ++c)
                for (int g = 0; g < grid; ++g)
                    input[c][g] = Complex(std::sin(0.07*(g+c)), std::cos(0.11*(g+c)));
            auto expected = input;
            for (int c = 0; c < width; ++c)
            {
                std::vector<Complex> coefficient(count, Complex(0.0, 0.0));
                for (int b = 0; b < count; ++b)
                {
                    Complex norm(0.0, 0.0);
                    for (int g = 0; g < grid; ++g)
                    {
                        norm += std::conj(basis[b][g])*basis[b][g];
                        coefficient[b] += std::conj(basis[b][g])*input[c][g];
                    }
                    coefficient[b] = std::abs(norm) == 0.0 ? Complex(0.0, 0.0) : coefficient[b]/norm;
                }
                for (int g = 0; g < grid; ++g)
                {
                    Complex correction(0.0, 0.0);
                    for (int b = 0; b < count; ++b) correction += coefficient[b]*basis[b][g];
                    expected[c][g] -= correction;
                }
            }
            projector.project_batch(input);
            EXPECT_EQ(projector.storage_bytes(), retained);
            for (int c = 0; c < width; ++c)
                for (int g = 0; g < grid; ++g)
                    EXPECT_NEAR(std::abs(input[c][g]-expected[c][g]), 0.0, 2.0e-13);
        }
    }
}
