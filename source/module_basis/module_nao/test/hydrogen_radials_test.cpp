#include <gtest/gtest.h>
#include <iostream>
#include "module_basis/module_nao/hydrogen_radials.h"
#include "module_base/math_integral.h"

#ifdef __MPI
#include <mpi.h>
#endif

class HydrogenRadialsTest : public ::testing::Test
{
    protected:
        virtual void SetUp()
        {
            // set up the test case
            itype_ = 1;
            charge_ = 1.0;
            nmax_ = 3;
            rcut_ = 20.0;
            dr_ = 0.01;
            rank_ = 0;
            ptr_log_ = NULL;
        }

        virtual void TearDown()
        {
            // tear down the test case
        }

        int itype_;
        double charge_;
        int nmax_;
        double rcut_;
        double dr_;
        int rank_;
        std::ofstream* ptr_log_;
};

TEST_F(HydrogenRadialsTest, UnzipStrategy)
{
    HydrogenRadials hr;
    std::vector<std::pair<int, int>> nl_pairs;
    // minimal, 1s
    nl_pairs = hr.unzip_strategy(1, "minimal");
    EXPECT_EQ(nl_pairs.size(), 1);
    EXPECT_EQ(nl_pairs[0].first, 1);
    EXPECT_EQ(nl_pairs[0].second, 0);
    // minimal, 1s, 2p, 3d, 4f
    nl_pairs = hr.unzip_strategy(4, "minimal");
    EXPECT_EQ(nl_pairs.size(), 4);
    EXPECT_EQ(nl_pairs[0].first, 1);
    EXPECT_EQ(nl_pairs[0].second, 0);
    EXPECT_EQ(nl_pairs[1].first, 2);
    EXPECT_EQ(nl_pairs[1].second, 1);
    EXPECT_EQ(nl_pairs[2].first, 3);
    EXPECT_EQ(nl_pairs[2].second, 2);
    EXPECT_EQ(nl_pairs[3].first, 4);
    EXPECT_EQ(nl_pairs[3].second, 3);
    // full, 1s
    nl_pairs = hr.unzip_strategy(1, "full");
    EXPECT_EQ(nl_pairs.size(), 1);
    EXPECT_EQ(nl_pairs[0].first, 1);
    EXPECT_EQ(nl_pairs[0].second, 0);
    // full, 1s, 2s, 2p, 3s, 3p, 3d, 4s, 4p, 4d, 4f
    nl_pairs = hr.unzip_strategy(4, "full");
    EXPECT_EQ(nl_pairs.size(), 10);
    EXPECT_EQ(nl_pairs[0].first, 1);
    EXPECT_EQ(nl_pairs[0].second, 0);
    EXPECT_EQ(nl_pairs[1].first, 2);
    EXPECT_EQ(nl_pairs[1].second, 0);
    EXPECT_EQ(nl_pairs[2].first, 2);
    EXPECT_EQ(nl_pairs[2].second, 1);
    EXPECT_EQ(nl_pairs[3].first, 3);
    EXPECT_EQ(nl_pairs[3].second, 0);
    EXPECT_EQ(nl_pairs[4].first, 3);
    EXPECT_EQ(nl_pairs[4].second, 1);
    EXPECT_EQ(nl_pairs[5].first, 3);
    EXPECT_EQ(nl_pairs[5].second, 2);
    EXPECT_EQ(nl_pairs[6].first, 4);
    EXPECT_EQ(nl_pairs[6].second, 0);
    EXPECT_EQ(nl_pairs[7].first, 4);
    EXPECT_EQ(nl_pairs[7].second, 1);
    EXPECT_EQ(nl_pairs[8].first, 4);
    EXPECT_EQ(nl_pairs[8].second, 2);
    EXPECT_EQ(nl_pairs[9].first, 4);
    EXPECT_EQ(nl_pairs[9].second, 3);
    nl_pairs = hr.unzip_strategy(1, "energy"); // H, 1s1 -> 1s -> [1s]
    EXPECT_EQ(nl_pairs.size(), 1);
    EXPECT_EQ(nl_pairs[0].first, 1);
    EXPECT_EQ(nl_pairs[0].second, 0);
    nl_pairs = hr.unzip_strategy(6, "energy"); // C, 1s1 2s2 2p2 -> 1s2s2p -> [1s][2s][2p3s]
    EXPECT_EQ(nl_pairs.size(), 4);
    EXPECT_EQ(nl_pairs[0].first, 1); // 1s
    EXPECT_EQ(nl_pairs[0].second, 0);
    EXPECT_EQ(nl_pairs[1].first, 2); // 2s
    EXPECT_EQ(nl_pairs[1].second, 0);
    EXPECT_EQ(nl_pairs[2].first, 2); // 2p
    EXPECT_EQ(nl_pairs[2].second, 1);
    EXPECT_EQ(nl_pairs[3].first, 3); // 3s
    EXPECT_EQ(nl_pairs[3].second, 0);
    nl_pairs = hr.unzip_strategy(29, "energy"); // Cu, 1s1 2s2 2p6 3s2 3p6 3d10 4s1 -> 1s2s2p3s3p4s3d -> [1s][2s][2p3s][3p4s][3d4p5s]
    EXPECT_EQ(nl_pairs.size(), 9);
    EXPECT_EQ(nl_pairs[0].first, 1); // 1s
    EXPECT_EQ(nl_pairs[0].second, 0);
    EXPECT_EQ(nl_pairs[1].first, 2); // 2s
    EXPECT_EQ(nl_pairs[1].second, 0);
    EXPECT_EQ(nl_pairs[2].first, 2); // 2p
    EXPECT_EQ(nl_pairs[2].second, 1);
    EXPECT_EQ(nl_pairs[3].first, 3); // 3s
    EXPECT_EQ(nl_pairs[3].second, 0);
    EXPECT_EQ(nl_pairs[4].first, 3); // 3p
    EXPECT_EQ(nl_pairs[4].second, 1);
    EXPECT_EQ(nl_pairs[5].first, 4); // 4s
    EXPECT_EQ(nl_pairs[5].second, 0);
    EXPECT_EQ(nl_pairs[6].first, 3); // 3d
    EXPECT_EQ(nl_pairs[6].second, 2);
    EXPECT_EQ(nl_pairs[7].first, 4); // 4p
    EXPECT_EQ(nl_pairs[7].second, 1);
    EXPECT_EQ(nl_pairs[8].first, 5); // 5s
    EXPECT_EQ(nl_pairs[8].second, 0);
}

TEST_F(HydrogenRadialsTest, RadialNorm)
{
    HydrogenRadials hr;
    std::vector<double> r;
    std::vector<double> f;
    double dr = 0.01;
    double rmax = 10.0;
    for (double r_ = 0.0; r_ <= rmax; r_ += dr)
    {
        r.push_back(r_);
        f.push_back(r_);
    }
    // radial norm computes the integral of r^2*f^2
    std::vector<double> r2f2;
    for (int i = 0; i < r.size(); i++)
    {
        r2f2.push_back(r[i]*r[i]*f[i]*f[i]);
    }
    double norm = hr.radial_norm(r, f);
    EXPECT_EQ(norm, sqrt(ModuleBase::Integral::simpson(r.size(), r2f2.data(), dr)));
}

TEST_F(HydrogenRadialsTest, MappingNLLZeta)
{
    HydrogenRadials hr;
    std::map<std::pair<int, int>, std::pair<int, int>> nl_lzeta_map;
    int l;
    int lzeta;
    // minimal, 1s, map (1, 0) to (0, 0)
    nl_lzeta_map = hr.mapping_nl_lzeta(1, "minimal");
    EXPECT_EQ(nl_lzeta_map.size(), 1);
    l = nl_lzeta_map[std::make_pair(1, 0)].first;
    lzeta = nl_lzeta_map[std::make_pair(1, 0)].second;
    EXPECT_EQ(l, 0);
    EXPECT_EQ(lzeta, 0);
    // minimal, 1s, 2p, 3d, 4f, map (1, 0) to (0, 0), (2, 1) to (1, 0), (3, 2) to (2, 0), (4, 3) to (3, 0)
    nl_lzeta_map = hr.mapping_nl_lzeta(4, "minimal");
    EXPECT_EQ(nl_lzeta_map.size(), 4);
    l = nl_lzeta_map[std::make_pair(1, 0)].first;
    lzeta = nl_lzeta_map[std::make_pair(1, 0)].second;
    EXPECT_EQ(l, 0);
    EXPECT_EQ(lzeta, 0);
    l = nl_lzeta_map[std::make_pair(2, 1)].first;
    lzeta = nl_lzeta_map[std::make_pair(2, 1)].second;
    EXPECT_EQ(l, 1);
    EXPECT_EQ(lzeta, 0);
    l = nl_lzeta_map[std::make_pair(3, 2)].first;
    lzeta = nl_lzeta_map[std::make_pair(3, 2)].second;
    EXPECT_EQ(l, 2);
    EXPECT_EQ(lzeta, 0);
    l = nl_lzeta_map[std::make_pair(4, 3)].first;
    lzeta = nl_lzeta_map[std::make_pair(4, 3)].second;
    EXPECT_EQ(l, 3);
    EXPECT_EQ(lzeta, 0);
    // full, 1s, map (1, 0) to (0, 0)
    nl_lzeta_map = hr.mapping_nl_lzeta(1, "full");
    EXPECT_EQ(nl_lzeta_map.size(), 1);
    l = nl_lzeta_map[std::make_pair(1, 0)].first;
    lzeta = nl_lzeta_map[std::make_pair(1, 0)].second;
    EXPECT_EQ(l, 0);
    EXPECT_EQ(lzeta, 0);
    // full, 1s, 2s, 2p, 3s, 3p, 3d, 4s, 4p, 4d, 4f, 
    // map (1, 0), (2, 0), (3, 0), (4, 0) to (0, 0), (0, 1), (0, 2), (0, 3)
    //     (2, 1), (3, 1), (4, 1) to (1, 0), (1, 1), (1, 2)
    //     (3, 2), (4, 2) to (2, 0), (2, 1)
    //     (4, 3) to (3, 0)
    nl_lzeta_map = hr.mapping_nl_lzeta(4, "full");
    EXPECT_EQ(nl_lzeta_map.size(), 10);
    l = nl_lzeta_map[std::make_pair(1, 0)].first;
    lzeta = nl_lzeta_map[std::make_pair(1, 0)].second;
    EXPECT_EQ(l, 0);
    EXPECT_EQ(lzeta, 0);
    l = nl_lzeta_map[std::make_pair(2, 0)].first;
    lzeta = nl_lzeta_map[std::make_pair(2, 0)].second;
    EXPECT_EQ(l, 0);
    EXPECT_EQ(lzeta, 1);
    l = nl_lzeta_map[std::make_pair(3, 0)].first;
    lzeta = nl_lzeta_map[std::make_pair(3, 0)].second;
    EXPECT_EQ(l, 0);
    EXPECT_EQ(lzeta, 2);
    l = nl_lzeta_map[std::make_pair(4, 0)].first;
    lzeta = nl_lzeta_map[std::make_pair(4, 0)].second;
    EXPECT_EQ(l, 0);
    EXPECT_EQ(lzeta, 3);
    l = nl_lzeta_map[std::make_pair(2, 1)].first;
    lzeta = nl_lzeta_map[std::make_pair(2, 1)].second;
    EXPECT_EQ(l, 1);
    EXPECT_EQ(lzeta, 0);
    l = nl_lzeta_map[std::make_pair(3, 1)].first;
    lzeta = nl_lzeta_map[std::make_pair(3, 1)].second;
    EXPECT_EQ(l, 1);
    EXPECT_EQ(lzeta, 1);
    l = nl_lzeta_map[std::make_pair(4, 1)].first;
    lzeta = nl_lzeta_map[std::make_pair(4, 1)].second;
    EXPECT_EQ(l, 1);
    EXPECT_EQ(lzeta, 2);
    l = nl_lzeta_map[std::make_pair(3, 2)].first;
    lzeta = nl_lzeta_map[std::make_pair(3, 2)].second;
    EXPECT_EQ(l, 2);
    EXPECT_EQ(lzeta, 0);
    l = nl_lzeta_map[std::make_pair(4, 2)].first;
    lzeta = nl_lzeta_map[std::make_pair(4, 2)].second;
    EXPECT_EQ(l, 2);
    EXPECT_EQ(lzeta, 1);
    l = nl_lzeta_map[std::make_pair(4, 3)].first;
    lzeta = nl_lzeta_map[std::make_pair(4, 3)].second;
    EXPECT_EQ(l, 3);
    EXPECT_EQ(lzeta, 0);
}

TEST_F(HydrogenRadialsTest, GenerateHydrogenRadialToconv)
{
    HydrogenRadials hr;
    std::vector<double> r;
    std::vector<double> Rnl;

    double rmax_chg1_n1l0 = hr.generate_hydrogen_radial_toconv(
        1.0,
        1,
        0,
        1e-7,
        0,
        r,
        Rnl
    );
    std::vector<double> r2Rnl2;
    for (int i = 0; i < r.size(); i++)
    {
        r2Rnl2.push_back(r[i]*r[i]*Rnl[i]*Rnl[i]);
    }
    double norm = ModuleBase::Integral::simpson(r.size(), r2Rnl2.data(), 0.01);
    EXPECT_NEAR(norm, 1.0, 1e-6);

    double rmax_chg4_n2l1 = hr.generate_hydrogen_radial_toconv(
        4.0,
        2,
        1,
        1e-7,
        0,
        r,
        Rnl
    );
    r2Rnl2.clear();
    for (int i = 0; i < r.size(); i++)
    {
        r2Rnl2.push_back(r[i]*r[i]*Rnl[i]*Rnl[i]);
    }
    norm = ModuleBase::Integral::simpson(r.size(), r2Rnl2.data(), 0.01);
    EXPECT_NEAR(norm, 1.0, 1e-6);

    EXPECT_NE(rmax_chg1_n1l0, rmax_chg4_n2l1);
}

TEST_F(HydrogenRadialsTest, Build)
{
    HydrogenRadials hr;
    // build 1s 2p 3d
    hr.build(
        itype_,
        charge_,
        nmax_,
        rcut_,
        dr_,
        1e-6,
        rank_,
        "H",
        "minimal",
        ptr_log_
    );
    // nmax = 1, minimal, yields 1s orbital
    EXPECT_EQ(hr.lmax(), 2);
    EXPECT_EQ(hr.nzeta(0), 1);
    EXPECT_EQ(hr.nzeta_max(), 1);
    EXPECT_EQ(hr.nchi(), 3);
    // build 1s 2s 2p 3s 3p 3d 4s 4p 4d 4f
    hr.build(
        itype_,
        charge_,
        4,
        rcut_,
        dr_,
        1e-6,
        rank_,
        "H",
        "full",
        ptr_log_
    );
    // nmax = 4, full, yields 1s 2s 2p 3s 3p 3d 4s 4p 4d 4f orbitals
    EXPECT_EQ(hr.lmax(), 3);
    EXPECT_EQ(hr.nzeta(0), 4);
    EXPECT_EQ(hr.nzeta(1), 3);
    EXPECT_EQ(hr.nzeta(2), 2);
    EXPECT_EQ(hr.nzeta(3), 1);
    EXPECT_EQ(hr.nzeta_max(), 4);
    EXPECT_EQ(hr.nchi(), 10);
}

int main(int argc, char** argv)
{
    testing::InitGoogleTest(&argc, argv);
    int result = RUN_ALL_TESTS();
    return result;
}