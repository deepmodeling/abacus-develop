
#include "../module_mulliken/output_mulliken.h"

#include "source_cell/cell_index.h"
#include "source_io/module_mulliken/output_dmk.h"
#include "source_io/module_mulliken/output_sk.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <algorithm>
#include <cmath>
#include <complex>

template <typename TK>
class OutputMullikenTest : public testing::Test
{
  protected:
    std::vector<std::string> atomLabels = {"Si"};
    std::vector<int> atomCounts = {1};
    std::vector<std::vector<int>> lnchiCounts = {{2, 2, 1}};
    Parallel_Orbitals paraV;
    int nrow;
    int ncol;
};

using MyTypes = ::testing::Types<double, std::complex<double>>;
TYPED_TEST_SUITE(OutputMullikenTest, MyTypes);

TYPED_TEST(OutputMullikenTest, OrbInfo)
{
    CellIndex cell_index = CellIndex(this->atomLabels, this->atomCounts, this->lnchiCounts, 1);
    cell_index.write_orb_info("./");
    std::ifstream ifs("./Orbital");
    std::string str((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    EXPECT_THAT(str, testing::HasSubstr("#io    spec    l    m    z  sym"));
    EXPECT_THAT(str, testing::HasSubstr("0      Si    2    3    1       dx^2-y^2"));
    remove("./Orbital");
}

#ifdef __MPI
TYPED_TEST(OutputMullikenTest, nspin1)
{
    this->nrow = 13;
    this->ncol = 13;
    this->paraV.init(this->nrow, this->ncol, 1, MPI_COMM_WORLD, 0);
    auto cell_index = CellIndex(this->atomLabels, this->atomCounts, this->lnchiCounts, 1);
    auto out_s_k = ModuleIO::Output_Sk<TypeParam>(nullptr, &this->paraV, 1, 1);
    auto out_dm_k = ModuleIO::Output_DMK<TypeParam>(nullptr, &this->paraV, 1, 1);
    auto mulp = ModuleIO::Output_Mulliken<TypeParam>(&(out_s_k), &(out_dm_k), &(this->paraV), &(cell_index), {0}, 1);
    mulp.write(0, "./");
    std::vector<double> tot_chg = mulp.get_tot_chg();
    EXPECT_NEAR(tot_chg[0], 4.0, 1e-5);
    std::ifstream ifs("./mulliken.txt");
    std::string str((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    EXPECT_THAT(str, testing::HasSubstr(" Total charge 4"));
    EXPECT_THAT(str, testing::HasSubstr("total charge    on atom 1     4.0000"));
    remove("./mulliken.txt");
}

TYPED_TEST(OutputMullikenTest, nspin2)
{
    this->nrow = 13;
    this->ncol = 13;
    this->paraV.init(this->nrow, this->ncol, 1, MPI_COMM_WORLD, 0);
    auto cell_index = CellIndex(this->atomLabels, this->atomCounts, this->lnchiCounts, 2);
    auto out_s_k = ModuleIO::Output_Sk<TypeParam>(nullptr, &this->paraV, 2, 1);
    auto out_dm_k = ModuleIO::Output_DMK<TypeParam>(nullptr, &this->paraV, 2, 1);
    auto mulp = ModuleIO::Output_Mulliken<TypeParam>(&(out_s_k), &(out_dm_k), &(this->paraV), &(cell_index), {0, 1}, 2);
    mulp.write(0, "./");
    std::vector<double> tot_chg = mulp.get_tot_chg();
    EXPECT_NEAR(tot_chg[0], 3.0, 1e-5);
    EXPECT_NEAR(tot_chg[1], 1.0, 1e-5);
    std::ifstream ifs("./mulliken.txt");
    std::string str((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    EXPECT_THAT(str, testing::HasSubstr(" Total charge 4"));
    EXPECT_THAT(str, testing::HasSubstr(" Total charge of spin1 3"));
    EXPECT_THAT(str, testing::HasSubstr(" Total charge of spin2 1"));
    EXPECT_THAT(str, testing::HasSubstr("total charge    on atom 1     4.0000"));
    EXPECT_THAT(str, testing::HasSubstr("total magnetism on atom 1     2.0000"));
    remove("./mulliken.txt");
}

TYPED_TEST(OutputMullikenTest, nspin4)
{
    this->nrow = 26;
    this->ncol = 26;
    this->paraV.init(this->nrow, this->ncol, 1, MPI_COMM_WORLD, 0);
    auto cell_index = CellIndex(this->atomLabels, this->atomCounts, this->lnchiCounts, 4);
    auto out_s_k = ModuleIO::Output_Sk<std::complex<double>>(nullptr, &this->paraV, 4, 1);
    auto out_dm_k = ModuleIO::Output_DMK<std::complex<double>>(nullptr, &this->paraV, 4, 1);
    auto mulp
        = ModuleIO::Output_Mulliken<std::complex<double>>(&(out_s_k), &(out_dm_k), &(this->paraV), &(cell_index), {0}, 4);
    mulp.write(0, "./");
    std::vector<double> tot_chg = mulp.get_tot_chg();
    EXPECT_NEAR(tot_chg[0], 4.0, 1e-5);
    EXPECT_NEAR(tot_chg[1], 0.0, 1e-5);
    EXPECT_NEAR(tot_chg[2], 0.0, 1e-5);
    EXPECT_NEAR(tot_chg[3], 2.0, 1e-5);
    std::ifstream ifs("./mulliken.txt");
    std::string str((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    EXPECT_THAT(str, testing::HasSubstr(" Total charge 4"));
    EXPECT_THAT(str, testing::HasSubstr("total charge    on atom 1     4.0000"));
    EXPECT_THAT(str, testing::HasSubstr("total magnetism on atom 1     0.0000     0.0000     2.0000"));
    remove("./mulliken.txt");
}

namespace
{
using cd = std::complex<double>;

/// spinor of the +1 eigenstate of mhat.sigma for a unit vector mhat
void spinor_from_direction(const double mhat[3], cd c[2])
{
    double cos_theta = mhat[2];
    cos_theta = std::max(-1.0, cos_theta);
    cos_theta = std::min(1.0, cos_theta);
    const double theta = std::acos(cos_theta);
    const double phi = std::atan2(mhat[1], mhat[0]);
    const double half_theta = 0.5 * theta;
    const double cos_half = std::cos(half_theta);
    const double sin_half = std::sin(half_theta);
    const double cos_phi = std::cos(phi);
    const double sin_phi = std::sin(phi);
    c[0] = cd(cos_half, 0.0);
    c[1] = sin_half * cd(cos_phi, sin_phi);
}

/// <c|sigma_i|c> with the standard Pauli matrices (sigma_y = [[0,-i],[i,0]])
double pauli_expectation(const cd c[2], const int component)
{
    const cd sigma[3][2][2] = {
        {{cd(0.0, 0.0), cd(1.0, 0.0)}, {cd(1.0, 0.0), cd(0.0, 0.0)}},
        {{cd(0.0, 0.0), cd(0.0, -1.0)}, {cd(0.0, 1.0), cd(0.0, 0.0)}},
        {{cd(1.0, 0.0), cd(0.0, 0.0)}, {cd(0.0, 0.0), cd(-1.0, 0.0)}},
    };
    cd value(0.0, 0.0);
    for (int a = 0; a < 2; ++a)
    {
        for (int b = 0; b < 2; ++b)
        {
            value += std::conj(c[a]) * sigma[component][a][b] * c[b];
        }
    }
    return value.real();
}
} // namespace

// Regression test for the nspin=4 Mulliken M_y sign.
// The stored DM is conj-first (DM_{ab} = conj(c_a) c_b, see dm_from_psi.h), so
// the physical Pauli components of mud = DM.S (S = 1 here) are
//   Mx = Re(mud_{up,down} + mud_{down,up})
//   My = Im(mud_{up,down} - mud_{down,up})
//   Mz = Re(mud_{up,up} - mud_{down,down}).
// The pre-fix code accumulated -Im(mud_{up,down}) + Im(mud_{down,up}) = -My,
// so this test fails before the fix and passes after it.
TEST(OutputMullikenSpinTest, Nspin4MyFollowsPauliConvention)
{
    const double dirs[4][3] = {
        {0.0, 1.0, 0.0},    // pure +y: transverse channel, maximal |My|
        {0.0, -1.0, 0.0},   // pure -y
        {0.6, 0.8, 0.0},    // in-plane moment with My > 0
        {0.36, 0.48, -0.8}, // general three-dimensional direction
    };

    for (const auto& mhat : dirs)
    {
        cd c[2];
        spinor_from_direction(mhat, c);

        // DM in the conj-first storage convention of dm_from_psi
        cd dm[4];
        dm[0] = std::conj(c[0]) * c[0];
        dm[1] = std::conj(c[0]) * c[1];
        dm[2] = std::conj(c[1]) * c[0];
        dm[3] = std::conj(c[1]) * c[1];

        // one spatial orbital, one k-point, overlap S = 1: mud = DM * S^T = DM.
        // mud is a ScaLAPACK column-major buffer; collect_MW_nspin4 accesses
        // mud(ic, ir) as the (ir, ic) element of the logical matrix.
        const int nw = 2;
        Parallel_Orbitals pv;
        pv.set_serial(nw, nw);
        ModuleBase::ComplexMatrix mud(pv.ncol, pv.nrow, true);
        for (int ir = 0; ir < nw; ++ir)
        {
            for (int ic = 0; ic < nw; ++ic)
            {
                const int idx = ir * 2 + ic;
                mud(ic, ir) = dm[idx];
            }
        }

        ModuleBase::matrix MecMulP(4, 1, true);
        ModuleIO::collect_MW_nspin4(MecMulP, mud, nw, &pv);

        const double expected_charge = 1.0;
        const double expected_mx = pauli_expectation(c, 0);
        const double expected_my = pauli_expectation(c, 1);
        const double expected_mz = pauli_expectation(c, 2);
        EXPECT_NEAR(MecMulP(0, 0), expected_charge, 1e-12);
        EXPECT_NEAR(MecMulP(1, 0), expected_mx, 1e-12)
            << "Mx wrong for direction (" << mhat[0] << "," << mhat[1] << "," << mhat[2] << ")";
        EXPECT_NEAR(MecMulP(2, 0), expected_my, 1e-12)
            << "My sign/value wrong (conj-first DM) for direction (" << mhat[0] << "," << mhat[1] << ","
            << mhat[2] << ")";
        EXPECT_NEAR(MecMulP(3, 0), expected_mz, 1e-12)
            << "Mz wrong for direction (" << mhat[0] << "," << mhat[1] << "," << mhat[2] << ")";
    }
}

#include "mpi.h"
int main(int argc, char** argv)
{

    MPI_Init(&argc, &argv);
    testing::InitGoogleTest(&argc, argv);

    int nprocs=0;
    int myrank=0;

    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    int result = RUN_ALL_TESTS();

    MPI_Finalize();

    return result;
}
#endif
