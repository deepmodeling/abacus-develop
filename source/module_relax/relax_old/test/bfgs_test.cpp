#include <gtest/gtest.h>
#include "for_test.h"
#include "module_relax/relax_old/bfgs.h"
#include "module_cell/unitcell.h"
#include "module_base/matrix.h"
#include "module_relax/relax_old/ions_move_basic.h"


TEST(BFGSTest, AllocateTest) {
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


/*TEST(BFGSTest, RelaxStepTest) {
    BFGS bfgs;
    UnitCell ucell;
    ModuleBase::matrix force(3, 3,true);  
    int size = 3;

    bfgs.allocate(size);

    force(0, 0)=0.1;
    force(1, 1)=0.2;
    force(2, 2)=0.3;

    ASSERT_NO_THROW(bfgs.relax_step(force, ucell));  


    EXPECT_EQ(bfgs.pos.size(), size);
}

TEST(BFGSTest, PrepareStepIntegrationTest) {
    BFGS bfgs;
    int size = 3;
    bfgs.allocate(size);

    std::vector<std::vector<double>> force = {{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}, {7.0, 8.0, 9.0}};
    std::vector<std::vector<double>> pos = {{0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}, {2.0, 2.0, 2.0}};
    std::vector<std::vector<double>> H = {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};
    std::vector<double> pos0(size, 0.0);
    std::vector<double> force0(size, 0.0);
    std::vector<double> steplength(size, 0.0);
    std::vector<std::vector<double>> dpos(size, std::vector<double>(size, 0.0));

    bfgs.PrepareStep(force, pos, H, pos0, force0, steplength, dpos);

    for (double step : steplength) {
        EXPECT_GT(step, 0.0);
    }
}*/


TEST(BFGSTest, FullStepTest) 
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