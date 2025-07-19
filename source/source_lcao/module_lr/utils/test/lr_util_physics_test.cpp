#include <gtest/gtest.h>

#include "../lr_util.h"
struct Atom_pseudo_Test
{
    size_t zv;
};
struct Atom_Test
{
    size_t na;
    Atom_pseudo_Test ncpp;
};
struct UnitCell_Test
{
    size_t ntype;
    Atom_Test* atoms;
};

TEST(LR_Util, cal_nelec)
{
    UnitCell_Test ucell;
    ucell.ntype = 2;
    ucell.atoms = new Atom_Test[2];
    ucell.atoms[0].ncpp.zv = 1;
    ucell.atoms[0].na = 1;
    ucell.atoms[1].ncpp.zv = 2;
    ucell.atoms[1].na = 2;
    size_t nelec = LR_Util::cal_nelec(ucell);
    EXPECT_EQ(nelec, 5);
}

TEST(LR_Util, cal_nocc)
{
    size_t nelec = 5;
    size_t nocc = LR_Util::cal_nocc(nelec);
    EXPECT_EQ(nocc, 3);
}

TEST(LR_Util, set_ix_map_diagonal)
{

    int nc, nv;
    std::pair<ModuleBase::matrix, std::vector<std::pair<int, int>>> res;

    auto check_result = [&nc, &nv, &res](std::vector<int>& cv2x_list) -> void
        {
            ModuleBase::matrix iciv2ix = std::get<0>(res);
            std::vector<std::pair<int, int>> ix2iciv = std::get<1>(res);
            EXPECT_EQ(iciv2ix.nr, nc);
            EXPECT_EQ(iciv2ix.nc, nv);
            EXPECT_EQ(ix2iciv.size(), nc * nv);
            for (int ic = 0;ic < nc;++ic)
                for (int iv = 0;iv < nv;++iv)
                {
                    EXPECT_EQ(iciv2ix(ic, iv), cv2x_list[ic * nv + iv]);
                    EXPECT_EQ(std::get<0>(ix2iciv[cv2x_list[ic * nv + iv]]), ic);
                    EXPECT_EQ(std::get<1>(ix2iciv[cv2x_list[ic * nv + iv]]), iv);
                }
        };

    // case1: 1*3, 3*1
    // ic\iv   0 1 2
    //  0       0 1 2
    nc = 1, nv = 3;
    std::vector<int> cv2x_list_1_0{ 0, 1, 2 };
    res = LR_Util::set_ix_map_diagonal(0, nc, nv);
    check_result(cv2x_list_1_0);
    res = LR_Util::set_ix_map_diagonal(1, nc, nv);
    check_result(cv2x_list_1_0);
    nc = 3;nv = 1;
    std::vector<int> cv2x_list_1_1{ 2, 1, 0 };
    res = LR_Util::set_ix_map_diagonal(0, nc, nv);
    check_result(cv2x_list_1_1);
    res = LR_Util::set_ix_map_diagonal(1, nc, nv);
    check_result(cv2x_list_1_1);


    // case 2: 2*4
    nc = 2;
    nv = 4;
    // mode 0: 
    // ic\iv  0 1 2 3
    // 1        0 2 4 6
    // 0       1 3 5 7
    res = LR_Util::set_ix_map_diagonal(0, nc, nv);
    std::vector<int> cv2x_list_2_0{ 1, 3, 5, 7, 0, 2, 4, 6 };
    check_result(cv2x_list_2_0);
    // mode 1: 
    // ic\iv  0 1 2 3
    // 1        0 1 3 5
    // 0       2 4 6 7
    res = LR_Util::set_ix_map_diagonal(1, nc, nv);
    std::vector<int> cv2x_list_2_1{ 2, 4, 6, 7, 0, 1, 3, 5 };
    check_result(cv2x_list_2_1);


    // case 3: 5*3
    nc = 5;
    nv = 3;
    // mode 0: 
    // ic\iv  0 1 2
    // 4        0 2 5
    // 3        1 4 8
    // 2        3 7 11
    // 1        6 10 13
    // 0        9 12 14
    res = LR_Util::set_ix_map_diagonal(0, nc, nv);
    std::vector<int> cv2x_list_3_0{ 9, 12, 14, 6, 10, 13, 3, 7, 11, 1, 4, 8, 0, 2, 5 };
    check_result(cv2x_list_3_0);
    // mode 1: 
    // ic\iv  0 1 2
    // 4        0 1 3
    // 3        2 4 6
    // 2        5 7 9
    // 1        8 10 12
    // 0        11 13 14
    res = LR_Util::set_ix_map_diagonal(1, nc, nv);
    std::vector<int> cv2x_list_3_1{ 11, 13, 14, 8, 10, 12, 5, 7, 9, 2, 4, 6, 0, 1, 3 };
    check_result(cv2x_list_3_1);

    //case 4: 4*4
    nc = 4;
    nv = 4;
    // mode 0:
    // ic\iv  0 1 2 3
    // 3        0 2 5 9
    // 2        1 4 8 12
    // 1        3 7 11 14
    // 0        6 10 13 15
    res = LR_Util::set_ix_map_diagonal(0, nc, nv);
    std::vector<int> cv2x_list_4_0{ 6, 10, 13, 15, 3, 7, 11, 14, 1, 4, 8, 12, 0, 2, 5, 9 };
    check_result(cv2x_list_4_0);
    // mode 1:
    // ic\iv  0 1 2 3
    // 3        0 1 3 6
    // 2        2 4 7 10
    // 1        5 8 11 13
    // 0        9 12 14 15
    res = LR_Util::set_ix_map_diagonal(1, nc, nv);
    std::vector<int> cv2x_list_4_1{ 9, 12, 14, 15, 5, 8, 11, 13, 2, 4, 7, 10, 0, 1, 3, 6 };
    check_result(cv2x_list_4_1);
}