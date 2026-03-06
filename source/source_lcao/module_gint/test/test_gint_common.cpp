#include "gtest/gtest.h"

#include "../gint_common.h"
#include "../gint_info.h"

#include <vector>

class GintCommonTest : public ::testing::Test
{
  protected:
    void SetUp() override
    {
        ucell_.ntype = 1;
        ucell_.nat = 2;
        ucell_.atoms = new Atom[ucell_.ntype];
        ucell_.iat2it = new int[ucell_.nat];
        ucell_.iat2ia = new int[ucell_.nat];

        ucell_.atoms[0].na = 2;
        ucell_.atoms[0].nw = 2;

        ucell_.iat2it[0] = 0;
        ucell_.iat2it[1] = 0;
        ucell_.iat2ia[0] = 0;
        ucell_.iat2ia[1] = 1;
        ucell_.set_iat2iwt(1);

        // 2 atom-pairs, each with a single R = (0, 0, 0)
        ijr_info_ = {
            2,
            0, 0, 1, 0, 0, 0,
            0, 1, 1, 0, 0, 0
        };
        gint_info_ = ModuleGint::GintInfo::make_test_instance_ptr(ucell_, ijr_info_);
    }

    void TearDown() override
    {
    }

    template<typename T>
    static void fill_values(hamilt::HContainer<T>& hr)
    {
        T* values = hr.get_wrapper();
        ASSERT_NE(values, nullptr);
        for (size_t i = 0; i < hr.get_nnr(); ++i)
        {
            values[i] = static_cast<T>(i + 1) / static_cast<T>(8);
        }
    }

    UnitCell ucell_;
    std::vector<int> ijr_info_;
    ModuleGint::GintInfo* gint_info_ = nullptr;
};

TEST_F(GintCommonTest, GetHrFloatBuildsExpectedShape)
{
    auto hr = gint_info_->get_hr<float>();

    EXPECT_EQ(hr.size_atom_pairs(), 2);
    EXPECT_EQ(hr.get_nnr(), 8);
    EXPECT_FALSE(hr.is_gamma_only());
    ASSERT_NE(hr.find_pair(0, 0), nullptr);
    ASSERT_NE(hr.find_pair(0, 1), nullptr);
}

TEST_F(GintCommonTest, CastHcontainerValuesPreservesLayoutAndValues)
{
    auto src = gint_info_->get_hr<double>();
    auto dst = gint_info_->get_hr<float>();
    fill_values(src);

    ModuleGint::cast_hcontainer_values(src, dst);

    EXPECT_EQ(src.get_ijr_info(), dst.get_ijr_info());
    ASSERT_NE(dst.get_wrapper(), nullptr);
    for (size_t i = 0; i < src.get_nnr(); ++i)
    {
        EXPECT_FLOAT_EQ(dst.get_wrapper()[i], static_cast<float>(src.get_wrapper()[i]));
    }
}

TEST_F(GintCommonTest, MakeCastHcontainerBuildsNewTypedCopy)
{
    auto src = gint_info_->get_hr<double>();
    fill_values(src);

    auto dst = ModuleGint::make_cast_hcontainer<float>(src);

    EXPECT_EQ(src.get_ijr_info(), dst.get_ijr_info());
    EXPECT_EQ(src.get_nnr(), dst.get_nnr());
    for (size_t i = 0; i < src.get_nnr(); ++i)
    {
        EXPECT_FLOAT_EQ(dst.get_wrapper()[i], static_cast<float>(src.get_wrapper()[i]));
    }
}

TEST_F(GintCommonTest, TransferDm2dToGintSupportsDoubleToFloat)
{
    auto dm = gint_info_->get_hr<double>();
    fill_values(dm);

    std::vector<hamilt::HContainer<double>*> dm_vec{&dm};
    std::vector<hamilt::HContainer<float>> dm_gint;
    dm_gint.push_back(gint_info_->get_hr<float>());

    ModuleGint::transfer_dm_2d_to_gint(*gint_info_, dm_vec, dm_gint);

    ASSERT_EQ(dm_gint.size(), 1);
    EXPECT_EQ(dm.get_ijr_info(), dm_gint[0].get_ijr_info());
    for (size_t i = 0; i < dm.get_nnr(); ++i)
    {
        EXPECT_FLOAT_EQ(dm_gint[0].get_wrapper()[i], static_cast<float>(dm.get_wrapper()[i]));
    }
}
