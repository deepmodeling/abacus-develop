#include "symmetry_test.h"
// cases
/**
    switch(ibrav)
    {
        case 1: return "01. Cubic P (simple)";
        case 2: return "02. Cubic I (body-centered)";
        case 3: return "03. Cubic F (face-centered)";
        case 4: return "04. Hexagonal cell";
        case 5: return "05. Tetrogonal P (simple)";
        case 6: return "06. Tetrogonal I (body-centered)";
        case 7: return "07. Rhombohedral (Trigonal) cell";
        case 8: return "08. Orthorhombic P(simple)";
        case 9: return "09. Orthorhombic I (body-centered)";
        case 10: return "10. Orthorhombic F (face-centered)";
        case 11: return "11. Orthorhombic C (base-centered)";
        case 12: return "12. Monoclinic P (simple)";
        case 13: return "13. Monoclinic A (base-center)";
        case 14: return "14. Triclinic cell";
        case 15: return "wrong !! ";
    }
*/

std::vector<stru_> stru_lib{
    stru_{1,
          "O_h",
          "m-3m",
          "Pm-3m",
          std::vector<double>{1., 0., 0., 0., 1., 0., 0., 0., 1.},
          std::vector<atomtype_>{atomtype_{"C",
                                           std::vector<std::vector<double>>{
                                               {0., 0., 0.},
                                           }}}, "C"},
    stru_{2,
          "O_h",
          "m-3m",
          "Im-3m",
          std::vector<double>{-0.5, 0.5, 0.5, 0.5, -0.5, 0.5, 0.5, 0.5, -0.5},
          std::vector<atomtype_>{atomtype_{"C",
                                           std::vector<std::vector<double>>{
                                               {0., 0., 0.},
                                           }}}, "C"},
    stru_{3,
          "O_h",
          "m-3m",
          "Fm-3m",
          std::vector<double>{0., 0.5, 0.5, 0.5, 0., 0.5, 0.5, 0.5, 0.},
          std::vector<atomtype_>{atomtype_{"C",
                                           std::vector<std::vector<double>>{
                                               {0., 0., 0.},
                                           }}}, "C"},
    stru_{4,
          "D_6h",
          "6/mmm",
          "P6/mmm",
          std::vector<double>{1., 0., 0., -0.5, 0.8660254, 0., 0., 0., 2.},
          std::vector<atomtype_>{atomtype_{"C",
                                           std::vector<std::vector<double>>{
                                               {0., 0., 0.},
                                           }}}, "C"},
    stru_{5,
          "D_4h",
          "4/mmm",
          "P4/mmm",
          std::vector<double>{1., 0., 0., 0., 1., 0., 0., 0., 2.},
          std::vector<atomtype_>{atomtype_{"C",
                                           std::vector<std::vector<double>>{
                                               {0., 0., 0.},
                                           }}}, "C"},
    stru_{6,
          "D_4h",
          "4/mmm",
          "I4/mmm",
          std::vector<double>{-0.35355339, 0.35355339, 1., 0.35355339, -0.35355339, 1., 0.35355339, 0.35355339, -1.},
          std::vector<atomtype_>{atomtype_{"C",
                                           std::vector<std::vector<double>>{
                                               {0., 0., 0.},
                                           }}}, "C"},
    stru_{7,
          "D_3d",
          "-3m",
          "R-3m",
          std::vector<double>{0.57357644,
                              0.33115451,
                              0.74923078,
                              -0.57357644,
                              0.33115451,
                              0.74923078,
                              0.,
                              -0.66230902,
                              0.74923078},
          std::vector<atomtype_>{atomtype_{"C",
                                           std::vector<std::vector<double>>{
                                               {-0., 0., 0.},
                                           }}}, "C"},
    stru_{8,
          "D_2h",
          "mmm",
          "Pmmm",
          std::vector<double>{1., 0., 0., 0., 2., 0., 0., 0., 3.},
          std::vector<atomtype_>{atomtype_{"C",
                                           std::vector<std::vector<double>>{
                                               {0., 0., 0.},
                                           }}}, "C"},
    stru_{9,
          "D_2h",
          "mmm",
          "Immm",
          std::vector<double>{-0.25, 0.75, 1., 0.25, -0.75, 1., 0.25, 0.75, -1.},
          std::vector<atomtype_>{atomtype_{"C",
                                           std::vector<std::vector<double>>{
                                               {0., 0., 0.},
                                           }}}, "C"},
    stru_{10,
          "D_2h",
          "mmm",
          "Fmmm",
          std::vector<double>{0., 1., 1.5, 0.5, 0., 1.5, 0.5, 1., 0.},
          std::vector<atomtype_>{atomtype_{"C",
                                           std::vector<std::vector<double>>{
                                               {0., 0., 0.},
                                           }}}, "C"},
    stru_{11,
          "D_2h",
          "mmm",
          "Cmmm",
          std::vector<double>{0.5, -1.5, 0., 0.5, 1.5, 0., 0., 0., 2.},
          std::vector<atomtype_>{atomtype_{"C",
                                           std::vector<std::vector<double>>{
                                               {0., 0., 0.},
                                           }}}, "C"},
    stru_{12,
          "C_2h",
          "2/m",
          "P2/m",
          std::vector<double>{1., 0., 0., 0., 2., 0., -0.02606043, 0., 2.81907786},
          std::vector<atomtype_>{atomtype_{"C",
                                           std::vector<std::vector<double>>{
                                               {0., 0., 0.},
                                           }}}, "C"},
   stru_{13,
         "C_2h",
         "2/m",
         "C2/m",
         std::vector<double>{0.5, -1., 0., 0.5, 1., 0., -0.40192379, 0., 1.5},
         std::vector<atomtype_>{atomtype_{"C",
                                          std::vector<std::vector<double>>{
                                              {0., 0., 0.},
                                          }}}, "C"},
   stru_{14,
         "S_2",
         "-1",
         "P-1",
         std::vector<double>{1., 0., 0., -0.28989928, 1.53691386, 0., -0.31595971, -0.66789914, 1.75670135},
         std::vector<atomtype_>{atomtype_{"C",
                                          std::vector<std::vector<double>>{
                                              {0., 0., 0.},
                                          }}}, "C"},

};

// test cases for space group and primitive cell analysis
// ibrav here means the number of primitive cells 
std::vector<stru_> supercell_lib{
    // bcc, 2 primitive cells
    stru_{2,
        "O_h",
        "",
        "O_h",
        std::vector<double>{1., 0., 0., 0., 1., 0., 0., 0., 1.},
        std::vector<atomtype_>{atomtype_{"C",
        std::vector<std::vector<double>>{
            {0., 0., 0.}, {0.5, 0.5, 0.5}}}}, "D",
        std::vector<int>{0, 1},
        {},
        {},
        {{0, 1}, {0, 2}, {1, 0}, {1, 2}, {2, 0}, {2, 1}},
        {{{0, 0}, {1, 1}, {2, 2}}}},
    // bct, 2
    stru_{2,
        "D_4h",
        "",
        "D_4h",
        std::vector<double>{1.2, 0., 0., 0., 1., 0., 0., 0., 1.},
        std::vector<atomtype_>{atomtype_{"C",
        std::vector<std::vector<double>>{
            {0., 0., 0.}, {0.5, 0.5, 0.5},}}} , "D",
        std::vector<int>{0, 1},
        {},
        {},
        {{0, 1}, {0, 2}, {1, 0}, {1, 2}, {2, 0}, {2, 1}},
        {{{1, 1}, {2, 2}}}},
    //bct, 2
    stru_{2,
        "D_2h",
        "",
        "D_2h",
        std::vector<double>{1.2, 0., 0., 0., 1.1, 0., 0., 0., 1.},
        std::vector<atomtype_>{atomtype_{"C",
        std::vector<std::vector<double>>{
            {0., 0., 0.}, {0.5, 0.5, 0.5}}}} , "D",
        std::vector<int>{0, 1},
        {},
        {},
        {{0, 1}, {0, 2}, {1, 0}, {1, 2}, {2, 0}, {2, 1}},
        {}},
    //fcc, 4
    stru_{4,
        "O_h",
        "",
        "O_h",
        std::vector<double>{1., 0., 0., 0., 1., 0., 0., 0., 1.},
        std::vector<atomtype_>{atomtype_{"C",
        std::vector<std::vector<double>>{
            {0., 0., 0.}, {0.5, 0.5, 0.},{0.5, 0., 0.5},{0., 0.5, 0.5}}}} , "D",
        std::vector<int>{0, 1, 2, 3},
        {},
        {},
        {{0, 1}, {0, 2}, {1, 0}, {1, 2}, {2, 0}, {2, 1}},
        {{{0, 0}, {1, 1}, {2, 2}}}},
    //fco, 4
    stru_{4,
        "D_2h",
        "",
        "D_2h",
        std::vector<double>{1.2, 0., 0., 0., 1.1, 0., 0., 0., 1.},
        std::vector<atomtype_>{atomtype_{"C",
        std::vector<std::vector<double>>{
            {0., 0., 0.}, {0.5, 0.5, 0.},{0.5, 0., 0.5},{0., 0.5, 0.5}}}} , "D",
        std::vector<int>{0, 1, 2, 3},
        {},
        {},
        {{0, 1}, {0, 2}, {1, 0}, {1, 2}, {2, 0}, {2,1}},
        {}},
    //3 in x
    stru_{3,
        "C_1h",
        "",
        "D_2h",
        std::vector<double>{3., -3., -3., -1., 1., -1., -1., -1., 1.},
        std::vector<atomtype_>{atomtype_{"C",
        std::vector<std::vector<double>>{
            {0., 0.1, 0.5}, {1. / 3., 0.1, 0.5},{2. / 3., 0.1, 0.5}}}} , "D",
        std::vector<int>{1},
        std::map<int,int>{{0, 2}},
        {},
        {{0, 1}, {0, 2}, {1, 0}, {2, 0}},
        {{{1, 2}, {2, 1}}, {{1, 1}, {2, 2}}}},
    //3 in y
    stru_{3,
        "C_1h",
        "",
        "D_2h",
        std::vector<double>{1., -1., -1., -3., 3., -3., -1., -1., 1.},
        std::vector<atomtype_>{atomtype_{"C",
        std::vector<std::vector<double>>{
            {0.4, 0., 0.1}, {0.4, 1. / 3., 0.1},{0.4, 2. / 3., 0.1}}}} , "D",
        std::vector<int>{2},
        std::map<int,int>{{0, 1}},
        {},
        {{0, 1}, {1, 0}, {1, 2}, {2, 1}},
        {{{0, 2}, {2, 0}}, {{0, 0}, {2, 2}}}},
    //3 in z
    stru_{3,
        "C_1h",
        "",
        "D_2h",
        std::vector<double>{1., -1., -1., -1., 1., -1., -3., -3., 3.},
        std::vector<atomtype_>{atomtype_{"C",
        std::vector<std::vector<double>>{
            {0.3, 0.1, 0.}, {0.3, 0.1, 1. / 3.},{0.3, 0.1, 2. / 3.}}}} , "D",
        std::vector<int>{1},
        std::map<int,int>{{0, 2}},
        {},
        {{0, 2}, {1, 2}, {2, 0}, {2, 1}},
        {{{0, 0}, {1, 1}}, {{0, 1}, {1, 0}}} },
    //6 in xy
    stru_{6,
        "C_1",
        "",
        "S_2",
        std::vector<double>{2., -2., -2., -3., 3., -3., -1., -1., 1.},
        std::vector<atomtype_>{atomtype_{"C",
        std::vector<std::vector<double>>{
            {0., 0., 0.1}, {0., 1. / 3., 0.1},{0., 2. / 3., 0.1},{0.5, 0., 0.1}, {0.5, 1. / 3., 0.1},{0.5, 2. / 3., 0.1}}}} , "D",
        std::vector<int>{0, 3},
        std::map<int,int>{{1, 2}, {4, 5}},
        {} },
    //6 in yz
    stru_{6,
        "C_1",
        "",
        "S_2",
        std::vector<double>{1., -1., -1., -2., 2., -2., -3., -3., 3.},
        std::vector<atomtype_>{atomtype_{"C",
        std::vector<std::vector<double>>{
            {0.1, 0., 0.}, {0.1, 0., 1. / 3.},{0.1, 0., 2. / 3.},{0.1, 0.5, 0.}, {0.1, 0.5, 1. / 3.},{0.1, 0.5, 2. / 3.}}}} , "D",
        std::vector<int>{2, 5},
        std::map<int,int>{{0, 1}, {3, 4}},
        {} },
    //6 in zx
    stru_{6,
        "C_1",
        "",
        "S_2",
        std::vector<double>{3., -3., -3., -1., 1., -1., -2., -2., 2.},
        std::vector<atomtype_>{atomtype_{"C",
        std::vector<std::vector<double>>{
            {0., 0.1, 0.}, {1. / 3., 0.1, 0.},{2. / 3., 0.1, 0.},  {0., 0.1, 0.5}, {1. / 3., 0.1, 0.5},{2. / 3., 0.1, 0.5}}}} , "D",
        std::vector<int>{0, 3},
        std::map<int,int>{{1, 2}, {4, 5}},
        {} },
    //hex: 3 in a1 - 231
    stru_{3,
        "C_1h",
        "",
        "C_2v",
        std::vector<double>{0., 1.59516, 2.76289, 20., 0., 0., 0., 9.57096, 0.},
        std::vector<atomtype_>{atomtype_{"Mo",
        std::vector<std::vector<double>>{
            {2./3., 0.1859875, 0.22222222}, {2./3., 0.1859875, 0.55555555},{2./3., 0.1859875, 0.88888888},
    }},
        atomtype_{"S",
        std::vector<std::vector<double>>{
            {1./3., 0.2642317, 0.11111111}, {1./3., 0.1077433, 0.11111111},{1./3., 0.2642317, 0.44444444},
            {1./3., 0.1077433, 0.44444444}, {1./3., 0.2642317, 0.77777777},{1./3., 0.1077433, 0.77777777},
    }}}
    , "D",
        {},
        {},
        std::vector<std::vector<int>>{{3, 4, 0, 1, 0}, {5, 6, 0, 1, 0}, {7, 8, 0, 1, 0}},
        {{0, 1}, {0, 2}, {1, 0}, {2, 0}},
        {} },
};