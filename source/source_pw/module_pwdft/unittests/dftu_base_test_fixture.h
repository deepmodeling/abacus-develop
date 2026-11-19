#ifndef DFTU_BASE_TEST_FIXTURE_H
#define DFTU_BASE_TEST_FIXTURE_H

#include "source_pw/module_pwdft/dftu_types.h"
#include "source_cell/unitcell.h"
#include "gtest/gtest.h"

#include <string>
#include <vector>

class Plus_U_Base;

/// Give each file-based test a private directory and remove its fixtures.
class DFTUTestDirectory
{
  public:
    DFTUTestDirectory();
    ~DFTUTestDirectory();
    const std::string& path() const { return path_; }
    std::string file(const std::string& name);

  private:
    std::string path_;
    std::vector<std::string> files_;
};

class DFTUBaseTest : public testing::Test
{
  protected:
    UnitCell ucell;

    void SetUp() override;
    void TearDown() override;
    void init_dftu(Plus_U_Base& dftu, const bool yukawa_potential);
    void init_dftu_for_spin(Plus_U_Base& dftu,
                            const bool yukawa_potential,
                            const int nspin,
                            const int npol,
                            const DFTU_BASE::OccmatSocLayout& soc_layout);
};

#endif // DFTU_BASE_TEST_FIXTURE_H
