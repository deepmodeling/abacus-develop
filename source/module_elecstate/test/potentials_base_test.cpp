#include "gtest/gtest.h"
#include "module_elecstate/potentials/pot_base.h"

ModuleBase::matrix::~matrix(){}

/***************************************************************
 *  unit test of functions in pot_base.h
 ****************************************************************/

/**
 * - Tested functions:
 *   - elecstate::PotBase::cal_v_eff()
 *   - elecstate::PotBase::cal_fixed_v()
 */

namespace elecstate
{
class MockPot : public PotBase
{
  public:
    bool get_fixed_mode() const
    {
        return fixed_mode;
    }

    bool get_dynamic_mode() const
    {
        return dynamic_mode;
    }
};

} // namespace elecstate

class PotBaseTest : public ::testing::Test
{
  protected:
    void SetUp() override
    {
        pot_base_ = new elecstate::MockPot;
    }

    void TearDown() override
    {
        delete pot_base_;
    }

    elecstate::MockPot* pot_base_;
};

TEST_F(PotBaseTest, CalVeff)
{
    ModuleBase::matrix v_eff;
    pot_base_->fixed_mode = true;
    if (pot_base_->get_fixed_mode())
    {
        EXPECT_NO_THROW(pot_base_->cal_v_eff(nullptr, nullptr, v_eff));
    }
}

TEST_F(PotBaseTest, CalFixedV)
{
    pot_base_->dynamic_mode = true;
    if (pot_base_->get_dynamic_mode())
    {
        EXPECT_NO_THROW(pot_base_->cal_fixed_v(nullptr));
    }
}