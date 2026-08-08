#include "source_estate/module_pot/potential_new.h"

#include "source_hamilt/module_xc/xc_functional.h"
#include "source_io/module_parameter/parameter.h"

#include "gtest/gtest.h"
#include <memory>
#include <string>
#include <vector>

class TestParameters
{
  public:
    static void reset()
    {
        PARAM.input.nspin = 1;
        PARAM.input.basis_type = "pw";
        PARAM.input.device = "cpu";
        PARAM.input.precision = "double";
        PARAM.sys.has_float_data = false;
        PARAM.sys.has_double_data = true;
    }

    static void use_cpu_single()
    {
        PARAM.input.device = "cpu";
        PARAM.input.precision = "single";
        PARAM.sys.has_float_data = true;
        PARAM.sys.has_double_data = false;
    }
};

Structure_Factor::Structure_Factor()
{
}
Structure_Factor::~Structure_Factor()
{
}
UnitCell::UnitCell()
{
}
UnitCell::~UnitCell()
{
}
Magnetism::Magnetism()
{
}
Magnetism::~Magnetism()
{
}
SepPot::SepPot()
{
}
SepPot::~SepPot()
{
}
Sep_Cell::Sep_Cell() noexcept
{
}
Sep_Cell::~Sep_Cell() noexcept
{
}
Charge::Charge()
{
}
Charge::~Charge()
{
}
surchem::surchem()
{
}
surchem::~surchem()
{
}

int XC_Functional::func_type = 1;
bool XC_Functional::ked_flag = false;

void XC_Functional::set_xc_type(const std::string xc_func_in)
{
    if (xc_func_in == "meta")
    {
        func_type = 3;
        ked_flag = true;
    }
    else
    {
        func_type = 1;
        ked_flag = false;
    }
}

namespace elecstate
{

class MockPotComponent : public PotBase
{
  public:
    MockPotComponent(const std::string& type, const int grid_size) : type_(type), grid_size_(grid_size)
    {
        this->fixed_mode = (type == "fixed");
        this->dynamic_mode = (type == "dynamic" || type == "ramp");
        ++created;
    }

    ~MockPotComponent() override
    {
        ++destroyed;
    }

    void cal_fixed_v(double* vl_pseudo) override
    {
        ++fixed_calls;
        for (int ir = 0; ir < grid_size_; ++ir)
        {
            vl_pseudo[ir] += 2.0;
        }
    }

    void cal_v_eff(const Charge* const chg, const UnitCell* const ucell, ModuleBase::matrix& v_eff) override
    {
        ++dynamic_calls;
        for (int is = 0; is < v_eff.nr; ++is)
        {
            for (int ir = 0; ir < v_eff.nc; ++ir)
            {
                v_eff(is, ir) += (type_ == "ramp" ? is + ir : dynamic_calls);
            }
        }
    }

    static void reset()
    {
        created = 0;
        destroyed = 0;
        fixed_calls = 0;
        dynamic_calls = 0;
    }

    static int created;
    static int destroyed;
    static int fixed_calls;
    static int dynamic_calls;

  private:
    std::string type_;
    int grid_size_;
};

int MockPotComponent::created = 0;
int MockPotComponent::destroyed = 0;
int MockPotComponent::fixed_calls = 0;
int MockPotComponent::dynamic_calls = 0;

PotBase* Potential::get_pot_type(const std::string& pot_type)
{
    const int grid_size = this->get_rho_basis() == nullptr ? 0 : this->get_rho_basis()->nrxx;
    return new MockPotComponent(pot_type, grid_size);
}

} // namespace elecstate

class PotentialNewTest : public ::testing::Test
{
  protected:
    void SetUp() override
    {
        TestParameters::reset();
        XC_Functional::set_xc_type("lda");
        elecstate::MockPotComponent::reset();

        smooth_basis.reset(new ModulePW::PW_Basis);
        dense_basis.reset(new ModulePW::PW_Basis_Sup);
        ucell.reset(new UnitCell);
        vloc.reset(new ModuleBase::matrix);
        structure_factors.reset(new Structure_Factor);
        solvent.reset(new surchem);
        etxc.reset(new double(0.0));
        vtxc.reset(new double(0.0));
    }

    void TearDown() override
    {
        potential.reset();
        XC_Functional::set_xc_type("lda");
        TestParameters::reset();
        elecstate::MockPotComponent::reset();
    }

    void create_potential(const ModulePW::PW_Basis* dense, const ModulePW::PW_Basis* smooth)
    {
        potential.reset(new elecstate::Potential(dense,
                                                 smooth,
                                                 ucell.get(),
                                                 vloc.get(),
                                                 structure_factors.get(),
                                                 solvent.get(),
                                                 etxc.get(),
                                                 vtxc.get()));
    }

    void setup_smooth_basis()
    {
        smooth_basis->initgrids(4, ModuleBase::Matrix3(1, 0, 0, 0, 1, 0, 0, 0, 1), 4);
        smooth_basis->initparameters(false, 4);
        smooth_basis->setuptransform();
        smooth_basis->collect_local_pw();
    }

    void setup_dense_basis()
    {
        dense_basis->initgrids(4, ModuleBase::Matrix3(1, 0, 0, 0, 1, 0, 0, 0, 1), 6);
        dense_basis->initparameters(false, 6);
        dense_basis->setuptransform(smooth_basis.get());
        dense_basis->collect_local_pw();
    }

    std::unique_ptr<ModulePW::PW_Basis> smooth_basis;
    std::unique_ptr<ModulePW::PW_Basis_Sup> dense_basis;
    std::unique_ptr<UnitCell> ucell;
    std::unique_ptr<ModuleBase::matrix> vloc;
    std::unique_ptr<Structure_Factor> structure_factors;
    std::unique_ptr<surchem> solvent;
    std::unique_ptr<double> etxc;
    std::unique_ptr<double> vtxc;
    std::unique_ptr<elecstate::Potential> potential;
};

TEST_F(PotentialNewTest, ConstructorAndPublicGetters)
{
    smooth_basis->nrxx = 100;
    create_potential(smooth_basis.get(), smooth_basis.get());

    EXPECT_TRUE(potential->fixed_mode);
    EXPECT_TRUE(potential->dynamic_mode);
    EXPECT_EQ(potential->get_eff_v().nr, 1);
    EXPECT_EQ(potential->get_eff_v().nc, 100);
    EXPECT_NE(potential->get_fixed_v(), nullptr);
    EXPECT_EQ(potential->get_rho_basis(), smooth_basis.get());
    EXPECT_EQ(potential->get_ucell(), ucell.get());
    EXPECT_EQ(potential->get_vloc(), vloc.get());
    EXPECT_EQ(potential->get_veff_smooth_data<double>(), potential->get_veff_smooth().c);
    EXPECT_EQ(potential->get_veff_smooth_data<float>(), nullptr);

    const elecstate::Potential& const_potential = *potential;
    EXPECT_EQ(const_potential.get_eff_v().nr, 1);
    EXPECT_EQ(const_potential.get_eff_v(0), potential->get_eff_v(0));
    EXPECT_EQ(const_potential.get_fixed_v(), potential->get_fixed_v());
    EXPECT_EQ(const_potential.get_veff_smooth().c, potential->get_veff_smooth().c);
}

TEST_F(PotentialNewTest, EmptyPotentialReturnsNullData)
{
    potential.reset(new elecstate::Potential);
    const elecstate::Potential& const_potential = *potential;

    EXPECT_EQ(potential->get_eff_v(0), nullptr);
    EXPECT_EQ(const_potential.get_eff_v(0), nullptr);
    EXPECT_EQ(potential->get_eff_vofk(0), nullptr);
    EXPECT_EQ(const_potential.get_eff_vofk(0), nullptr);
    EXPECT_EQ(potential->get_veff_smooth_data<float>(), nullptr);
    EXPECT_EQ(potential->get_veff_smooth_data<double>(), nullptr);
    EXPECT_EQ(potential->get_vofk_smooth_data<float>(), nullptr);
    EXPECT_EQ(potential->get_vofk_smooth_data<double>(), nullptr);
}

TEST_F(PotentialNewTest, ConstructorCPUSingle)
{
    TestParameters::use_cpu_single();
    smooth_basis->nrxx = 100;
    create_potential(smooth_basis.get(), smooth_basis.get());

    EXPECT_NE(potential->get_veff_smooth_data<float>(), nullptr);
    EXPECT_EQ(potential->get_veff_smooth_data<double>(), nullptr);
}

TEST_F(PotentialNewTest, MetaPotentialPublicGetters)
{
    XC_Functional::set_xc_type("meta");
    smooth_basis->nrxx = 100;
    create_potential(smooth_basis.get(), smooth_basis.get());

    EXPECT_EQ(potential->get_eff_vofk().nr, 1);
    EXPECT_EQ(potential->get_eff_vofk().nc, 100);
    EXPECT_EQ(potential->get_vofk_smooth().nr, 1);
    EXPECT_EQ(potential->get_vofk_smooth().nc, 100);
    EXPECT_EQ(potential->get_vofk_smooth_data<double>(), potential->get_vofk_smooth().c);

    potential->get_eff_vofk()(0, 0) = 7.0;
    EXPECT_DOUBLE_EQ(potential->get_eff_vofk(0)[0], 7.0);

    const elecstate::Potential& const_potential = *potential;
    EXPECT_EQ(const_potential.get_eff_vofk(0), potential->get_eff_vofk(0));
    EXPECT_EQ(const_potential.get_vofk_smooth().c, potential->get_vofk_smooth().c);
}

TEST_F(PotentialNewTest, PotRegisterReplacesAndDestroysComponents)
{
    potential.reset(new elecstate::Potential);
    potential->pot_register({"fixed", "dynamic"});
    EXPECT_EQ(elecstate::MockPotComponent::created, 2);
    EXPECT_EQ(elecstate::MockPotComponent::destroyed, 0);

    potential->pot_register({"dynamic"});
    EXPECT_EQ(elecstate::MockPotComponent::created, 3);
    EXPECT_EQ(elecstate::MockPotComponent::destroyed, 2);

    potential.reset();
    EXPECT_EQ(elecstate::MockPotComponent::destroyed, 3);
}

TEST_F(PotentialNewTest, PublicUpdateFlowsScheduleComponents)
{
    smooth_basis->nrxx = 8;
    create_potential(smooth_basis.get(), smooth_basis.get());
    potential->pot_register({"fixed", "dynamic"});
    Charge charge;

    potential->update_from_charge(&charge, ucell.get());
    EXPECT_EQ(elecstate::MockPotComponent::fixed_calls, 1);
    EXPECT_EQ(elecstate::MockPotComponent::dynamic_calls, 1);
    for (int ir = 0; ir < smooth_basis->nrxx; ++ir)
    {
        EXPECT_DOUBLE_EQ(potential->get_fixed_v()[ir], 2.0);
        EXPECT_DOUBLE_EQ(potential->get_eff_v()(0, ir), 3.0);
        EXPECT_DOUBLE_EQ(potential->get_veff_smooth()(0, ir), 3.0);
    }

    potential->update_from_charge(&charge, ucell.get());
    EXPECT_EQ(elecstate::MockPotComponent::fixed_calls, 1);
    EXPECT_EQ(elecstate::MockPotComponent::dynamic_calls, 2);
    EXPECT_DOUBLE_EQ(potential->get_eff_v()(0, 0), 4.0);

    potential->init_pot(&charge);
    EXPECT_EQ(elecstate::MockPotComponent::fixed_calls, 2);
    EXPECT_EQ(elecstate::MockPotComponent::dynamic_calls, 3);
    EXPECT_DOUBLE_EQ(potential->get_eff_v()(0, 0), 5.0);

    ModuleBase::matrix vnew;
    potential->get_vnew(&charge, vnew);
    EXPECT_EQ(elecstate::MockPotComponent::fixed_calls, 2);
    EXPECT_EQ(elecstate::MockPotComponent::dynamic_calls, 4);
    ASSERT_EQ(vnew.nr, 1);
    ASSERT_EQ(vnew.nc, smooth_basis->nrxx);
    for (int ir = 0; ir < vnew.nc; ++ir)
    {
        EXPECT_DOUBLE_EQ(vnew(0, ir), 1.0);
    }
}

TEST_F(PotentialNewTest, DifferentBasisObjectsTriggerInterpolation)
{
    XC_Functional::set_xc_type("meta");
    setup_smooth_basis();
    setup_dense_basis();
    create_potential(dense_basis.get(), smooth_basis.get());
    potential->pot_register({"ramp"});

    for (int ir = 0; ir < potential->get_eff_vofk().nc; ++ir)
    {
        potential->get_eff_vofk()(0, ir) = 2.0 * ir;
    }

    Charge charge;
    potential->update_from_charge(&charge, ucell.get());

    const std::vector<double> expected_veff
        = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26};
    const std::vector<double> expected_vofk
        = {0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 50, 52};
    ASSERT_EQ(potential->get_veff_smooth().nc, expected_veff.size());
    ASSERT_EQ(potential->get_vofk_smooth().nc, expected_vofk.size());
    for (int ir = 0; ir < potential->get_veff_smooth().nc; ++ir)
    {
        EXPECT_DOUBLE_EQ(potential->get_veff_smooth()(0, ir), expected_veff[ir]);
        EXPECT_DOUBLE_EQ(potential->get_vofk_smooth()(0, ir), expected_vofk[ir]);
    }
}

TEST_F(PotentialNewTest, SameBasisObjectCopiesPotentialDirectly)
{
    XC_Functional::set_xc_type("meta");
    setup_smooth_basis();
    create_potential(smooth_basis.get(), smooth_basis.get());
    potential->pot_register({"ramp"});

    for (int ir = 0; ir < potential->get_eff_vofk().nc; ++ir)
    {
        potential->get_eff_vofk()(0, ir) = 2.0 * ir;
    }

    Charge charge;
    potential->update_from_charge(&charge, ucell.get());

    for (int ir = 0; ir < potential->get_veff_smooth().nc; ++ir)
    {
        EXPECT_DOUBLE_EQ(potential->get_veff_smooth()(0, ir), ir);
        EXPECT_DOUBLE_EQ(potential->get_vofk_smooth()(0, ir), 2.0 * ir);
    }
}

TEST_F(PotentialNewTest, DifferentBasisObjectsRequireMatchingGammaOnly)
{
    setup_smooth_basis();
    setup_dense_basis();
    dense_basis->gamma_only = !smooth_basis->gamma_only;
    create_potential(dense_basis.get(), smooth_basis.get());
    potential->pot_register({"ramp"});
    Charge charge;

    EXPECT_EXIT(potential->update_from_charge(&charge, ucell.get()), ::testing::ExitedWithCode(1), "");
}
