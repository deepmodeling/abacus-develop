#include "mpi.h"
#include "../symmetry_rotation.h"
#include  "gtest/gtest.h"
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
InfoNonlocal::InfoNonlocal() {}
InfoNonlocal::~InfoNonlocal() {}
Magnetism::Magnetism() {}
Magnetism::~Magnetism() {}

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
    symrot.cal_rotmat_Slm(&C41, 1);
    RI::Tensor<std::complex<double>>& rotmat = symrot.get_rotmat_Slm()[0][1];
    int l = 1;
    for (int m1 = -l;m1 <= l;++m1)
        for (int m2 = -l;m2 <= l;++m2)
        {
            int i = m1 + l, j = m2 + l;
            EXPECT_NEAR(rotmat(i, j).real(), c_dagger_D_c_C41_ref[i * 3 + j].real(), DOUBLETHRESHOLD);
            EXPECT_NEAR(rotmat(i, j).imag(), c_dagger_D_c_C41_ref[i * 3 + j].imag(), DOUBLETHRESHOLD);
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
    RI::Tensor<std::complex<double>> block({ 2, 2 });
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

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    testing::InitGoogleTest(&argc, argv);
    int result = RUN_ALL_TESTS();
    MPI_Finalize();
    return result;
}
