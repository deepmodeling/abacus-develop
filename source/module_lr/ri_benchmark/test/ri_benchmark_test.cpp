#include <gtest/gtest.h>
#include "../ri_benchmark.h"
pseudo::pseudo() {}
pseudo::~pseudo() {}
Atom_pseudo::Atom_pseudo() {}
Atom_pseudo::~Atom_pseudo() {}
InfoNonlocal::InfoNonlocal() {}
InfoNonlocal::~InfoNonlocal() {}
Magnetism::Magnetism() {}
Magnetism::~Magnetism() {}
Atom::Atom() { this->nw = 2; }
Atom::~Atom() {}
UnitCell::UnitCell() {
    atoms = new Atom[1];
    iat2it = new int[1]; iat2it[0] = 0;
    iat2iwt.resize(1, 0);
}
UnitCell::~UnitCell() {
    delete[] atoms;
}
// inline const int* UnitCell::get_iat2iwt(int iat) { return iat2iwt; }

TEST(RI_Benchmark, SlicePsi)
{
    const int nk = 1, nbands = 2, nbasis = 3;
    psi::Psi<double> psi(nk, nbands, nbasis, nbasis, true);
    for (int i = 0; i < nk * nbands * nbasis; i++) {
        psi.get_pointer()[i] = i;
}
    std::vector<double> psi_slice = RI_Benchmark::slice_psi(psi, 1, 1, 1, 2);
    EXPECT_DOUBLE_EQ(psi_slice[0], 4);
    EXPECT_DOUBLE_EQ(psi_slice[1], 5);
}
TEST(RI_Benchmark, CalCV)
{
    const int nabf1 = 2, nabf2 = 3, nw = 1;
    RI_Benchmark::TLRI<double> Cs_a{ {0, {{{0, {0, 0, 0}}, RI::Tensor<double>({nabf1, nw, nw})}}} };
    for (int i = 0;i < nabf1;++i) { Cs_a[0][{0, { 0, 0, 0 }}](i, 0, 0) = static_cast<double>(i); }
    RI_Benchmark::TLRI<double> Vs{ {0, {{{0, {0, 0, 0}}, RI::Tensor<double>({nabf1, nabf2})}}} };
    for (int i = 0;i < nabf1;++i) { for (int j = 0;j < nabf2;++j) { Vs[0][{0, { 0, 0, 0 }}](i, j) = static_cast<double>(i + j); } }
    RI_Benchmark::TLRI<double> CV = RI_Benchmark::cal_CV(Cs_a, Vs);
    EXPECT_EQ(CV.size(), 1);
    EXPECT_EQ(CV[0].size(), 1);
    for (int i = 0;i < nabf2;++i) { EXPECT_DOUBLE_EQ((CV[0][{0, { 0, 0, 0 }}](i, 0, 0)), i + 1); }
}
TEST(RI_Benchmark, CalCsMO)
{
    const int nabf = 2, nocc = 1, nvirt = 1, nao = 2;
    RI_Benchmark::TLRI<double> Cs_ao{ {0, {{{0, {0, 0, 0}}, RI::Tensor<double>({nabf, nao, nao })}}} };
    for (int i = 0;i < nabf * nao * nao;++i) { Cs_ao[0][{0, { 0, 0, 0 }}].ptr()[i] = static_cast<double>(i); }

    const UnitCell ucell;
    psi::Psi<double> psi_ks(1, 2, 2, 2, true);
    for (int i = 0;i < 4;++i) { psi_ks.get_pointer()[i] = static_cast<double>(i); }
    RI_Benchmark::TLRI<double> Cs_a_mo = RI_Benchmark::cal_Cs_mo(ucell, Cs_ao, psi_ks, nocc, nvirt, false);
    std::vector<double> Cs_a_mo_ref = { 11,31 };
    for (int i = 0;i < 2;++i) { EXPECT_DOUBLE_EQ(Cs_a_mo.at(0).at({ 0, { 0, 0, 0 } })(i, 0, 0), Cs_a_mo_ref[i]); }
    RI_Benchmark::TLRI<double> Cs_b_mo = RI_Benchmark::cal_Cs_mo(ucell, Cs_ao, psi_ks, nocc, nvirt, true);
    std::vector<double> Cs_b_mo_ref = { 13,33 };
    for (int i = 0;i < 2;++i) { EXPECT_DOUBLE_EQ(Cs_b_mo.at(0).at({ 0, { 0, 0, 0 } })(i, 0, 0), Cs_b_mo_ref[i]); }
}

int main(int argc, char** argv)
{
    srand(time(nullptr));  // for random number generator
    MPI_Init(&argc, &argv);
    testing::InitGoogleTest(&argc, argv);
    int result = RUN_ALL_TESTS();
    MPI_Finalize();
    return result;
}