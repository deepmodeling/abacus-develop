#include <gtest/gtest.h>
#include "for_test.h"
#include "module_relax/relax_old/bfgs.h"
#include "module_cell/unitcell.h"
#include "module_base/matrix.h"
#include "module_relax/relax_old/ions_move_basic.h"
#include "module_relax/relax_old/matrix_methods.h"

class BFGSTest : public ::testing::Test {
protected:
    BFGS bfgs;
    UnitCell ucell;
    std::vector<std::vector<double>> force;

    void SetUp() override {
        int size = 10; 
        bfgs.allocate(size);

        ucell.ntype = 2;
        ucell.lat0 = 1.0;
        ucell.nat = 10;
        ucell.atoms = new Atom[ucell.ntype];
        for (int i = 0; i < ucell.ntype; i++) {
            ucell.atoms[i].na = 5; 
            ucell.atoms[i].tau = std::vector<ModuleBase::Vector3<double>>(5);
            ucell.atoms[i].taud = std::vector<ModuleBase::Vector3<double>>(5);
            ucell.atoms[i].mbl = std::vector<ModuleBase::Vector3<int>>(5, {1, 1, 1});
        }

        force = std::vector<std::vector<double>>(size, std::vector<double>(3, 0.0));
        for (int i = 0; i < force.size(); ++i) {
            for (int j = 0; j < 3; ++j) {
                force[i][j] =  -0.1 * (i + 1);
            }
        }
    }
};

TEST_F(BFGSTest, PrepareStep) {
    bfgs.PrepareStep(force, bfgs.pos, bfgs.H, bfgs.pos0, bfgs.force0, bfgs.steplength, bfgs.dpos, ucell);
    EXPECT_EQ(bfgs.steplength.size(), 10);
    for (int i = 0; i < 10; ++i) {
        EXPECT_GT(bfgs.steplength[i], 0);
    }
}


TEST_F(BFGSTest, AllocateTest) {
    BFGS bfgs;
    int size = 5;
    bfgs.allocate(size);


    EXPECT_EQ(bfgs.steplength.size(), size);
    EXPECT_EQ(bfgs.force0.size(), 3*size);
    EXPECT_EQ(bfgs.H.size(), 3*size);
    for (const auto& row : bfgs.H) {
        EXPECT_EQ(row.size(), 3*size);
    }
}

TEST_F(BFGSTest, FullStepTest) 
{ 
    BFGS bfgs; 
    UnitCell ucell; 
    ModuleBase::matrix force(3, 3); 
    int size = 3; 
    bfgs.allocate(size);  
    force(0, 0)=-0.5; 
    force(1, 1)=-0.3; 
    force(2, 2)=0.1; 
    EXPECT_EQ(bfgs.force.size(), size); 
    EXPECT_EQ(bfgs.pos.size(), size); 
}