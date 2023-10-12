#include "module_basis/module_representation/representation.h"
#include <gtest/gtest.h>

class RepresentationTest : public ::testing::Test
{
    protected:
        void SetUp() override {
            Representation<std::complex<double>, psi::DEVICE_CPU> rep;
        }
        void TearDown() override {}
};

TEST_F(RepresentationTest, Constructor)
{

}

TEST_F(RepresentationTest, Configure)
{

}

TEST_F(RepresentationTest, SetRepIn)
{

}

TEST_F(RepresentationTest, AddRepOut)
{

}

TEST_F(RepresentationTest, AddTransformPair)
{

}

TEST_F(RepresentationTest, CleanRepresentations)
{

}

TEST_F(RepresentationTest, Transform)
{
    
}