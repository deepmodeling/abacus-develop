#include <gtest/gtest.h>

#include <RI/global/Tensor.h>
#include <array>
#include <vector>
#include <fstream>
#include "../write_ri_cv.hpp"

using TC = std::array<int, 3>;
using TAC = std::pair<int, TC>;
template <typename T>
using TLRI = std::map<int, std::map<TAC, RI::Tensor<T>>>;

TEST(LRI_CV_Tools, ReadCs)
{
    const TLRI<std::complex<double>>& Cs = LRI_CV_Tools::read_Cs_ao<std::complex<double>>("./support/Cs", 1e-10);
    LRI_CV_Tools::write_Cs_ao(Cs, "./support/Cs_out");
    EXPECT_EQ(Cs.size(), 2);
    EXPECT_EQ(Cs.at(0).size(), 2);
    EXPECT_EQ(Cs.at(1).size(), 2);
    EXPECT_EQ(Cs.at(1).at({ 1, {0, 0, 0} }).shape.size(), 3);
    EXPECT_EQ(Cs.at(1).at({ 1, {0, 0, 0} }).shape[0], 15);
    EXPECT_EQ(Cs.at(1).at({ 1, {0, 0, 0} }).shape[1], 5);
    EXPECT_EQ(Cs.at(1).at({ 1, {0, 0, 0} }).shape[2], 5);
    EXPECT_DOUBLE_EQ(Cs.at(1).at({ 1, {0, 0, 0} })(0, 0, 0).real(), -0.00444197);


    const TLRI<double>& Vs = LRI_CV_Tools::read_Vs_abf<double>("./support/Vs", 1e-10);
    LRI_CV_Tools::write_Vs_abf(Vs, "./support/Vs_out");
    EXPECT_EQ(Vs.size(), 2);
    EXPECT_EQ(Vs.at(0).size(), 6);
    EXPECT_EQ(Vs.at(1).at({ 1, {1, 0, 0} }).shape.size(), 2);

}