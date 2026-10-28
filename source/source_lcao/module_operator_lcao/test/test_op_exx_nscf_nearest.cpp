#include "../op_exx_lcao.h"
#include "source_io/module_restart/restart_exx_csr.h"

#include "gtest/gtest.h"
#include <array>
#include <cstdlib>
#include <fstream>
#include <string>
#include <vector>

namespace
{

std::string make_test_dir()
{
    const std::string dir = "./op_exx_nscf_nearest_tmp";
    std::system(("rm -rf " + dir).c_str());
    std::system(("mkdir -p " + dir).c_str());
    return dir;
}

// The csrFileReader on develop expects the full sparse-file header (title,
// spin lines, dimension and R counts, lattice block) before the CSR blocks,
// in the shape produced by the SR.csr fixture of source_base tests.
void write_csr_header(std::ofstream& ofs)
{
    ofs << " --- Ionic Step 0 ---\n";
    ofs << " # minimal HexxR fixture for the cell-nearest operator test\n";
    ofs << " 1 # number of spin directions\n";
    ofs << " 1 # spin index\n";
    ofs << " 2 # number of localized basis\n";
    ofs << " 2 # number of Bravais lattice vector R\n";
    ofs << "\n";
    ofs << " user_defined_lattice\n";
    ofs << " 1\n";
    ofs << " 1 0 0\n";
    ofs << " 0 1 0\n";
    ofs << " 0 0 1\n";
    ofs << " X\n";
    ofs << " 2\n";
    ofs << " Direct\n";
    ofs << " 0 0 0\n";
    ofs << " 0.9 0 0\n";
}

void write_csr_r_block(std::ofstream& ofs, const int rx, const int ry, const int rz, const double value)
{
    ofs << " " << rx << " " << ry << " " << rz << " 1\n";
    ofs << " # CSR values\n " << value << "\n";
    ofs << " # CSR column indices\n 1\n";
    ofs << " # CSR row pointers\n 0 1 1\n";
    ofs << "\n";
}

void write_minimal_hexx_csr(const std::string& dir)
{
    std::system(("mkdir -p " + dir).c_str());
    std::ofstream ofs(dir + "/HexxR0_0.csr");
    write_csr_header(ofs);
    write_csr_r_block(ofs, 0, 0, 0, 1.0);
    write_csr_r_block(ofs, 1, 0, 0, 2.0);
}

void write_centered_hexx_csr(const std::string& dir)
{
    std::system(("mkdir -p " + dir).c_str());
    std::ofstream ofs(dir + "/HexxR0_0.csr");
    write_csr_header(ofs);
    write_csr_r_block(ofs, -1, 0, 0, 1.0);
    write_csr_r_block(ofs, 0, 0, 0, 2.0);
}

} // namespace

class OperatorExxNscfCellNearestTest : public ::testing::Test
{
  protected:
    void SetUp() override
    {
#ifdef __MPI
        MPI_Comm_size(MPI_COMM_WORLD, &world_size_);
        MPI_Comm_rank(MPI_COMM_WORLD, &world_rank_);
#endif

        dir_ = make_test_dir();
        write_minimal_hexx_csr(dir_);

        ucell_.ntype = 1;
        ucell_.nat = 2;
        ucell_.atoms = new Atom[ucell_.ntype];
        ucell_.iat2it.resize(ucell_.nat);
        ucell_.iat2ia.resize(ucell_.nat);
        ucell_.itia2iat.create(ucell_.ntype, ucell_.nat);
        ucell_.atoms[0].na = 2;
        ucell_.atoms[0].nw = 1;
        ucell_.atoms[0].tau.resize(ucell_.nat);
        ucell_.atoms[0].tau[0] = ModuleBase::Vector3<double>(0.0, 0.0, 0.0);
        ucell_.atoms[0].tau[1] = ModuleBase::Vector3<double>(0.9, 0.0, 0.0);
        ucell_.a1 = ModuleBase::Vector3<double>(1.0, 0.0, 0.0);
        ucell_.a2 = ModuleBase::Vector3<double>(0.0, 1.0, 0.0);
        ucell_.a3 = ModuleBase::Vector3<double>(0.0, 0.0, 1.0);
        for (int iat = 0; iat < ucell_.nat; ++iat)
        {
            ucell_.iat2it[iat] = 0;
            ucell_.iat2ia[iat] = iat;
            ucell_.itia2iat(0, iat) = iat;
        }
        ucell_.set_iat2iwt(1);
        ucell_.iwt2iat = {0, 1};
        ucell_.iwt2iw = {0, 0};

        paraV_ = new Parallel_Orbitals();
#ifdef __MPI
        paraV_->init(2, 2, 1, MPI_COMM_WORLD);
        paraV_->set_atomic_trace(ucell_.get_iat2iwt(), ucell_.nat, 2);
#endif

        hR_ = new hamilt::HContainer<double>(paraV_);

    }

    void TearDown() override
    {
        // tmp_mocks have partial UnitCell ownership; process cleanup avoids double free.
    }

    int world_size_ = 1;
    int world_rank_ = 0;
    std::string dir_;
    UnitCell ucell_;
    Parallel_Orbitals* paraV_ = nullptr;
    hamilt::HContainer<double>* hR_ = nullptr;
};

TEST_F(OperatorExxNscfCellNearestTest, ReallocatesHexxRWithNearestImageWhenPeriodIsInferredFromFile)
{
    std::vector<std::map<int, std::map<hamilt::TAC, RI::Tensor<double>>>> Hexxd;
    ModuleIO::read_Hexxs_csr(dir_ + "/HexxR0", ucell_, 1, 2, Hexxd);

    const auto inferred = hamilt::infer_complete_Rs_period_from_Hexxs(Hexxd);
    const std::array<int, 3> expected_period = {2, 1, 1};
    ASSERT_TRUE(inferred.first);
    EXPECT_EQ(expected_period, inferred.second);

    const std::array<int, 3> raw_R = {1, 0, 0};
    const auto cell_nearest = hamilt::init_cell_nearest(ucell_, inferred.second);
    const auto nearest = cell_nearest.get_cell_nearest_discrete(0, 1, raw_R);
    ASSERT_TRUE(nearest != raw_R);
    hamilt::reallocate_hcontainer(Hexxd, hR_, &cell_nearest);

    EXPECT_NE(nullptr, hR_->find_matrix(0, 1, nearest[0], nearest[1], nearest[2]));
    EXPECT_EQ(nullptr, hR_->find_matrix(0, 1, 1, 0, 0));
}

TEST_F(OperatorExxNscfCellNearestTest, InfersPeriodFromCenteredCompleteRsSet)
{
    const std::string centered_dir = dir_ + "/centered";
    write_centered_hexx_csr(centered_dir);

    std::vector<std::map<int, std::map<hamilt::TAC, RI::Tensor<double>>>> Hexxd;
    ModuleIO::read_Hexxs_csr(centered_dir + "/HexxR0", ucell_, 1, 2, Hexxd);

    const auto inferred = hamilt::infer_complete_Rs_period_from_Hexxs(Hexxd);
    const std::array<int, 3> expected_period = {2, 1, 1};
    ASSERT_TRUE(inferred.first);
    EXPECT_EQ(expected_period, inferred.second);
}

TEST_F(OperatorExxNscfCellNearestTest, WriteHexxsCsrRoundTripsThroughReader)
{
    // read_Hexxs_csr and write_Hexxs_csr must agree on the file layout; the
    // reader parses the sparse-file format with the lattice block, so the
    // writer has to emit that layout, not the bare save_sparse header.
    std::vector<std::map<int, std::map<hamilt::TAC, RI::Tensor<double>>>> Hexxd;
    ModuleIO::read_Hexxs_csr(dir_ + "/HexxR0", ucell_, 1, 2, Hexxd);

    // fields only touched by the writer's lattice block
    ucell_.latName = "user_defined_lattice";
    ucell_.lat0 = 1.0;
    ucell_.latvec = ModuleBase::Matrix3(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0);
    ucell_.atoms[0].label = "X";
    ucell_.atoms[0].stapos_wf = 0;
    ucell_.atoms[0].taud.resize(ucell_.nat);
    ucell_.atoms[0].taud[0] = ModuleBase::Vector3<double>(0.0, 0.0, 0.0);
    ucell_.atoms[0].taud[1] = ModuleBase::Vector3<double>(0.9, 0.0, 0.0);

    const std::string out_dir = dir_ + "/roundtrip";
    std::system(("mkdir -p " + out_dir).c_str());
    ModuleIO::write_Hexxs_csr(out_dir + "/HexxR0", ucell_, Hexxd);

    std::vector<std::map<int, std::map<hamilt::TAC, RI::Tensor<double>>>> Hexxd2;
    ModuleIO::read_Hexxs_csr(out_dir + "/HexxR0", ucell_, 1, 2, Hexxd2);

    ASSERT_EQ(1u, Hexxd2.size());
    ASSERT_EQ(1u, Hexxd2[0].count(0));
    EXPECT_DOUBLE_EQ(1.0, Hexxd2[0].at(0).at({1, {0, 0, 0}})(0, 0));
    EXPECT_DOUBLE_EQ(2.0, Hexxd2[0].at(0).at({1, {1, 0, 0}})(0, 0));
}

TEST(OperatorExxNscfCellNearestPolicyTest, UsesInferredPeriodIndependentlyOfPathOffset)
{
    EXPECT_TRUE(hamilt::can_remap_wigner_seitz_for_nscf(hamilt::Add_Hexx_Type::R, false, false));
    EXPECT_TRUE(hamilt::can_remap_wigner_seitz_for_nscf(hamilt::Add_Hexx_Type::R, true, true));
    EXPECT_FALSE(hamilt::can_remap_wigner_seitz_for_nscf(hamilt::Add_Hexx_Type::R, true, false));
    EXPECT_FALSE(hamilt::can_remap_wigner_seitz_for_nscf(hamilt::Add_Hexx_Type::k, false, true));
}

TEST_F(OperatorExxNscfCellNearestTest, FindsAllEqualDistanceImagesInZnOPrimitiveCell)
{
    ucell_.a1 = ModuleBase::Vector3<double>(2.292, 2.292, 0.0);
    ucell_.a2 = ModuleBase::Vector3<double>(2.292, 0.0, 2.292);
    ucell_.a3 = ModuleBase::Vector3<double>(0.0, 2.292, 2.292);

    const auto images = hamilt::find_nearest_bvk_cells(ucell_, 0, 0, {4, 0, 0}, {8, 8, 8});

    ASSERT_EQ(2, images.size());
    EXPECT_NE(images.end(), std::find(images.begin(), images.end(), std::array<int, 3>{-4, 0, 0}));
    EXPECT_NE(images.end(), std::find(images.begin(), images.end(), std::array<int, 3>{4, 0, 0}));
}

TEST_F(OperatorExxNscfCellNearestTest, SplitsHexxEquallyAcrossWignerSeitzBoundaryImages)
{
    ucell_.a1 = ModuleBase::Vector3<double>(2.292, 2.292, 0.0);
    ucell_.a2 = ModuleBase::Vector3<double>(2.292, 0.0, 2.292);
    ucell_.a3 = ModuleBase::Vector3<double>(0.0, 2.292, 2.292);

    std::vector<std::map<int, std::map<hamilt::TAC, RI::Tensor<double>>>> Hexxd(1);
    RI::Tensor<double> block({1, 1});
    block(0, 0) = 2.0;
    Hexxd[0][0][{0, {4, 0, 0}}] = block;

    const auto stats = hamilt::remap_Hexxs_wigner_seitz(ucell_, {8, 8, 8}, Hexxd);

    EXPECT_EQ(1, stats.remapped_blocks);
    EXPECT_EQ(1, stats.split_blocks);
    EXPECT_EQ(2, stats.max_images);
    ASSERT_EQ(2, Hexxd[0][0].size());
    EXPECT_DOUBLE_EQ(1.0, Hexxd[0][0].at({0, {-4, 0, 0}})(0, 0));
    EXPECT_DOUBLE_EQ(1.0, Hexxd[0][0].at({0, {4, 0, 0}})(0, 0));
}

int main(int argc, char** argv)
{
#ifdef __MPI
    MPI_Init(&argc, &argv);
#endif
    ::testing::InitGoogleTest(&argc, argv);
    const int result = RUN_ALL_TESTS();
#ifdef __MPI
    MPI_Finalize();
#endif
    return result;
}
