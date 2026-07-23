#include "mpi.h"
#include "../symm_rotation.h"
#include  "gtest/gtest.h"
#include <cmath>
#define DOUBLETHRESHOLD 1e-8

/*

tested functions:
- wigner_d
- wigner_D
- ovlp_Ylm_Slm
- get_euler_angle
- cal_rotmat_Slm
- get_return_lattice

untested functions:
- cal_Ms (depending on UnitCell,K_Vectors, Parallel_2D)
- restore_dm  (depending on Ms)
- rot_matrix_ao (depending on Ms)

*/

// mock the useless functions
pseudo::pseudo() {}
pseudo::~pseudo() {}
Atom::Atom() {}
Atom::~Atom() {}
Atom_pseudo::Atom_pseudo() {}
Atom_pseudo::~Atom_pseudo() {}
UnitCell::UnitCell() {}
UnitCell::~UnitCell() {}

Magnetism::Magnetism() {}
Magnetism::~Magnetism() {}
SepPot::SepPot(){}
SepPot::~SepPot(){}
Sep_Cell::Sep_Cell() noexcept {}
Sep_Cell::~Sep_Cell() noexcept {}

class SymmetryRotationTest : public testing::Test
{
protected:
    void SetUp() override
    {
        //init pv
        int myrank, dsize;
        MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
        MPI_Comm_size(MPI_COMM_WORLD, &dsize);
        pv.init(matsize, matsize, 1, MPI_COMM_WORLD);
    }
    ModuleBase::Matrix3 C41 = ModuleBase::Matrix3(0, 1, 0, -1, 0, 0, 0, 0, 1);
    std::vector<std::complex<double>> wigerD_p_C41_ref = { ModuleBase::IMAG_UNIT, 0, 0, 0, 1, 0, 0, 0, -ModuleBase::IMAG_UNIT };
    std::vector<std::complex<double>> cmm_p_C41_ref = { -ModuleBase::IMAG_UNIT / sqrt(2), 0, -1 / sqrt(2), 0, 1, 0, -ModuleBase::IMAG_UNIT / sqrt(2), 0, 1 / sqrt(2) };
    std::vector<std::complex<double>> c_dagger_D_c_C41_ref = { 1, 0, 0, 0, 0, -1, 0, 1, 0 };
    ModuleSymmetry::Symmetry_rotation symrot;
    Parallel_2D pv;
    const int matsize = 5;  //2s1p
};

// inline void outmat(RI::Tensor<std::complex<double>>& mat, int size, std::string name)
// {
//     std::cout << name << std::endl;
//     for (int i = 0;i < size;++i)
//     {
//         for (int j = 0;j < size;++j)std::cout << mat(i, j) << " ";
//         std::cout << std::endl;
//     }
// }
TEST_F(SymmetryRotationTest, Wignerd)
{
    EXPECT_NEAR(symrot.wigner_d(0, 1, 0, 0), 1, DOUBLETHRESHOLD);
    EXPECT_NEAR(symrot.wigner_d(0, 1, 1, 1), 1, DOUBLETHRESHOLD);
    EXPECT_NEAR(symrot.wigner_d(0, 1, -1, -1), 1, DOUBLETHRESHOLD);
}
TEST_F(SymmetryRotationTest, WignerD)
{
    RI::Tensor<std::complex<double>> wignerD_p_C41({ 3, 3 });
    int l = 1;
    for (int m1 = -l;m1 <= l;++m1)
        for (int m2 = -l;m2 <= l;++m2)
        {
            int i = m1 + l, j = m2 + l;
            wignerD_p_C41(i, j) = symrot.wigner_D(ModuleBase::Vector3<double>(0, 0, ModuleBase::PI / 2), 1, m1, m2, false);
            EXPECT_NEAR(wignerD_p_C41(i, j).real(), wigerD_p_C41_ref[i * 3 + j].real(), DOUBLETHRESHOLD);
            EXPECT_NEAR(wignerD_p_C41(i, j).imag(), wigerD_p_C41_ref[i * 3 + j].imag(), DOUBLETHRESHOLD);
            // alpha and gamma are the same when beta = 0
            wignerD_p_C41(i, j) = symrot.wigner_D(ModuleBase::Vector3<double>(ModuleBase::PI / 2, 0, 0), 1, m1, m2, false);
            EXPECT_NEAR(wignerD_p_C41(i, j).real(), wigerD_p_C41_ref[i * 3 + j].real(), DOUBLETHRESHOLD);
            EXPECT_NEAR(wignerD_p_C41(i, j).imag(), wigerD_p_C41_ref[i * 3 + j].imag(), DOUBLETHRESHOLD);
        }
    // outmat(wignerD_p_C41, 3, "wignerD_p_C41_cal");

}

TEST_F(SymmetryRotationTest, EulerAngle)
{
    ModuleBase::Vector3<double> euler_angle = symrot.get_euler_angle(C41);
    EXPECT_NEAR(euler_angle.x + euler_angle.z, ModuleBase::PI / 2, DOUBLETHRESHOLD);
    EXPECT_NEAR(euler_angle.y, 0, DOUBLETHRESHOLD);
}

TEST_F(SymmetryRotationTest, OvlpYS)
{
    RI::Tensor<std::complex<double>> c_mm_p({ 3, 3 });
    int l = 1;
    for (int m1 = -l;m1 <= l;++m1)
        for (int m2 = -l;m2 <= l;++m2)
        {
            int i = m1 + l, j = m2 + l;
            c_mm_p(i, j) = symrot.ovlp_Ylm_Slm(l, m1, m2);
            EXPECT_NEAR(c_mm_p(i, j).real(), cmm_p_C41_ref[i * 3 + j].real(), DOUBLETHRESHOLD);
            EXPECT_NEAR(c_mm_p(i, j).imag(), cmm_p_C41_ref[i * 3 + j].imag(), DOUBLETHRESHOLD);
        }
}

TEST_F(SymmetryRotationTest, RotMat)
{
    symrot.cal_rotmat_Slm(&C41, 1, -1);
    ModuleBase::ComplexMatrix& rotmat = symrot.get_rotmat_Slm()[0][1];
    int l = 1;
    for (int m1 = -l;m1 <= l;++m1)
        for (int m2 = -l;m2 <= l;++m2)
        {
            int i = m1 + l, j = m2 + l;
            EXPECT_NEAR(rotmat(i, j).real(), c_dagger_D_c_C41_ref[i * 3 + j].real(), DOUBLETHRESHOLD);
            EXPECT_NEAR(rotmat(i, j).imag(), c_dagger_D_c_C41_ref[i * 3 + j].imag(), DOUBLETHRESHOLD);
        }
}

TEST_F(SymmetryRotationTest, RotMatHighLIdentityFinite)
{
    ModuleBase::Matrix3 identity(1, 0, 0, 0, 1, 0, 0, 0, 1);
    symrot.cal_rotmat_Slm(&identity, 8);
    RI::Tensor<std::complex<double>>& rotmat = symrot.get_rotmat_Slm()[0][8];
    const int dim = 2 * 8 + 1;
    for (int i = 0; i < dim; ++i)
    {
        for (int j = 0; j < dim; ++j)
        {
            EXPECT_TRUE(std::isfinite(rotmat(i, j).real()));
            EXPECT_TRUE(std::isfinite(rotmat(i, j).imag()));
            const double target = (i == j) ? 1.0 : 0.0;
            EXPECT_NEAR(rotmat(i, j).real(), target, 1e-8);
            EXPECT_NEAR(rotmat(i, j).imag(), 0.0, 1e-8);
        }
    }
}

TEST_F(SymmetryRotationTest, GetReturnLattice)
{
    ModuleBase::Vector3<double> posd_a1(1. / 3., 1. / 3., 0.2);
    ModuleBase::Vector3<double> posd_a2(1. / 3., 1. / 3., -0.2);
    ModuleBase::Vector3<double> gtransd(0, 0, 0);
    ModuleBase::Matrix3 gmatd(-1, 1, 0, -1, 0, 0, 0, 0, -1);
    ModuleBase::Vector3<double> return_lattice = symrot.get_return_lattice(ModuleSymmetry::Symmetry(), gmatd, gtransd, posd_a1, posd_a2);
    EXPECT_NEAR(return_lattice.x, -1, DOUBLETHRESHOLD);
    EXPECT_NEAR(return_lattice.y, 0, DOUBLETHRESHOLD);
    EXPECT_NEAR(return_lattice.z, -1, DOUBLETHRESHOLD);
}

TEST_F(SymmetryRotationTest, SetBlockToMat2d)
{
    std::vector<std::complex<double>> obj_mat(pv.get_local_size());
    for (int j = 0;j < pv.get_col_size();++j)
        for (int i = 0;i < pv.get_row_size();++i)
            obj_mat[j * pv.get_row_size() + i] = std::complex<double>(static_cast<double>(pv.local2global_row(i)), static_cast<double>(pv.local2global_col(j)));
    ModuleBase::ComplexMatrix block(2, 2);
    block(0, 0) = 0; block(0, 1) = -1; block(1, 0) = -2; block(1, 1) = -3;
    symrot.set_block_to_mat2d(2, 3, block, obj_mat, pv);
    for (int i = 2;i < 4;++i)
        for (int j = 3;j < 5;++j)
        {
            int local_index = pv.global2local_col(j) * pv.get_row_size() + pv.global2local_row(i);
            if (pv.in_this_processor(i, j))
            {
                EXPECT_NEAR(obj_mat[local_index].real(), block(j - 3, i - 2).real(), DOUBLETHRESHOLD);
                EXPECT_NEAR(obj_mat[local_index].imag(), block(j - 3, i - 2).imag(), DOUBLETHRESHOLD);
            }
        }
}

TEST_F(SymmetryRotationTest, RotateAOCoefficientsUsesTransposeConvention)
{
    const std::vector<std::complex<double>> matrix = {
        {0.0, 0.0}, {0.0, 1.0},
        {1.0, 0.0}, {0.0, 0.0}};
    const std::vector<std::complex<double>> coefficients = {{1.0, 2.0}, {3.0, -1.0}};

    const auto rotated = ModuleSymmetry::rotate_ao_coefficients_dense(matrix, coefficients);

    ASSERT_EQ(rotated.size(), 2);
    EXPECT_NEAR(rotated[0].real(), 1.0, DOUBLETHRESHOLD);
    EXPECT_NEAR(rotated[0].imag(), 3.0, DOUBLETHRESHOLD);
    EXPECT_NEAR(rotated[1].real(), 1.0, DOUBLETHRESHOLD);
    EXPECT_NEAR(rotated[1].imag(), 2.0, DOUBLETHRESHOLD);
}

TEST_F(SymmetryRotationTest, RotateAOCoefficientsAppliesTimeReversalBeforeRotation)
{
    const std::vector<std::complex<double>> matrix = {
        {0.0, 0.0}, {0.0, 1.0},
        {1.0, 0.0}, {0.0, 0.0}};
    const std::vector<std::complex<double>> coefficients = {{1.0, 2.0}, {3.0, -1.0}};

    const auto rotated = ModuleSymmetry::rotate_ao_coefficients_dense(matrix, coefficients, true);

    ASSERT_EQ(rotated.size(), 2);
    EXPECT_NEAR(rotated[0].real(), -1.0, DOUBLETHRESHOLD);
    EXPECT_NEAR(rotated[0].imag(), 3.0, DOUBLETHRESHOLD);
    EXPECT_NEAR(rotated[1].real(), 1.0, DOUBLETHRESHOLD);
    EXPECT_NEAR(rotated[1].imag(), -2.0, DOUBLETHRESHOLD);
}

TEST_F(SymmetryRotationTest, RotatedCoefficientsMatchDensityMatrixConvention)
{
    const double inv_sqrt_two = 1.0 / std::sqrt(2.0);
    const std::vector<std::complex<double>> matrix = {
        {inv_sqrt_two, 0.0}, {0.0, inv_sqrt_two},
        {0.0, inv_sqrt_two}, {inv_sqrt_two, 0.0}};
    const std::vector<std::complex<double>> coefficients = {{1.0, 2.0}, {3.0, -1.0}};

    for (const bool time_reversal : {false, true})
    {
        const auto rotated =
            ModuleSymmetry::rotate_ao_coefficients_dense(matrix, coefficients, time_reversal);
        for (int i = 0; i < 2; ++i)
        {
            for (int j = 0; j < 2; ++j)
            {
                std::complex<double> expected = 0.0;
                for (int a = 0; a < 2; ++a)
                {
                    for (int b = 0; b < 2; ++b)
                    {
                        const auto density = time_reversal
                                                 ? std::conj(coefficients[a]) * coefficients[b]
                                                 : coefficients[a] * std::conj(coefficients[b]);
                        expected += matrix[i * 2 + a] * density
                                    * std::conj(matrix[j * 2 + b]);
                    }
                }
                EXPECT_NEAR((rotated[i] * std::conj(rotated[j])).real(),
                            expected.real(),
                            DOUBLETHRESHOLD);
                EXPECT_NEAR((rotated[i] * std::conj(rotated[j])).imag(),
                            expected.imag(),
                            DOUBLETHRESHOLD);
            }
        }
    }
}

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    testing::InitGoogleTest(&argc, argv);
    int result = RUN_ALL_TESTS();
    MPI_Finalize();
    return result;
}
