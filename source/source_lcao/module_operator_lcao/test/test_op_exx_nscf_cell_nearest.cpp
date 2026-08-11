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
    const std::string dir = "./op_exx_nscf_cell_nearest_tmp";
    std::system(("rm -rf " + dir).c_str());
    std::system(("mkdir -p " + dir).c_str());
    return dir;
}

void write_minimal_hexx_csr(const std::string& dir)
{
    std::ofstream ofs(dir + "/HexxR0_0.csr");
    ofs << "STEP: 0\n";
    ofs << "Matrix Dimension of Hexxs_0(R): 2\n";
    ofs << "Matrix number of Hexxs_0(R): 2\n";

    ofs << "0 0 0 1\n";
    ofs << "1.0\n";
    ofs << "1\n";
    ofs << "0 1 1\n";

    ofs << "1 0 0 1\n";
    ofs << "2.0\n";
    ofs << "1\n";
    ofs << "0 1 1\n";
}

void write_centered_hexx_csr(const std::string& dir)
{
    std::system(("mkdir -p " + dir).c_str());
    std::ofstream ofs(dir + "/HexxR0_0.csr");
    ofs << "STEP: 0\n";
    ofs << "Matrix Dimension of Hexxs_0(R): 2\n";
    ofs << "Matrix number of Hexxs_0(R): 2\n";

    ofs << "-1 0 0 1\n";
    ofs << "1.0\n";
    ofs << "1\n";
    ofs << "0 1 1\n";

    ofs << "0 0 0 1\n";
    ofs << "2.0\n";
    ofs << "1\n";
    ofs << "0 1 1\n";
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
        ucell_.iat2it = new int[ucell_.nat];
        ucell_.iat2ia = new int[ucell_.nat];
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
        ucell_.iwt2iat = new int[2]{0, 1};
        ucell_.iwt2iw = new int[2]{0, 0};

        paraV_ = new Parallel_Orbitals();
#ifdef __MPI
        paraV_->init(2, 2, 1, MPI_COMM_WORLD);
        paraV_->set_atomic_trace(ucell_.get_iat2iwt(), ucell_.nat, 2);
#endif

        hR_ = new hamilt::HContainer<double>(paraV_);

        kv_.set_nks(1);
        kv_.set_nkstot(1);
        kv_.set_nkstot_full(1);
        kv_.kvec_d.resize(1, ModuleBase::Vector3<double>(0.0, 0.0, 0.0));
        kv_.nmp[0] = 0;
        kv_.nmp[1] = 0;
        kv_.nmp[2] = 0;
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
    K_Vectors kv_;
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
