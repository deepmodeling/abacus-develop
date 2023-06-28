#include "symmetry_test_cases.h"
#include "mpi.h"

/************************************************
 *  unit test of class Symmetry
 * 4. function: `force_symmetry`
 *
***********************************************/
// mock the useless functions
void output::printM3(std::ofstream& ofs, const std::string& description, const ModuleBase::Matrix3& m) {}
pseudo_nc::pseudo_nc() {}
pseudo_nc::~pseudo_nc() {}
Atom::Atom() {}
Atom::~Atom() {}
Atom_pseudo::Atom_pseudo() {}
Atom_pseudo::~Atom_pseudo() {}
UnitCell::UnitCell() {}
UnitCell::~UnitCell() {}
Magnetism::Magnetism() {}
Magnetism::~Magnetism() {}

TEST_F(SymmetryTest, ForceSymmetry)
{
    auto check_force = [](stru_& conf, ModuleBase::matrix& force)
    {
        // 1. check zeros  
        for (auto iat : conf.force_zero_iat)
            for (int j = 0; j < 3; ++j)
                EXPECT_NEAR(force(iat, j), 0.0, DOUBLETHRESHOLD);
        // 2. check opposites
        for (auto oppo_pair : conf.force_oppo_iat)
            for (int j = 0; j < 3; ++j)
                EXPECT_NEAR(force(oppo_pair.first, j), -force(oppo_pair.second, j), DOUBLETHRESHOLD);
        for (auto oppo_xyz : conf.force_oppo_iat_xyz)
            for (int j = 0;j < 3;++j)
                if (oppo_xyz[j + 2] == 1)
                    EXPECT_NEAR(force(oppo_xyz[0], j), -force(oppo_xyz[1], j), DOUBLETHRESHOLD);
                else
                    EXPECT_NEAR(force(oppo_xyz[0], j), force(oppo_xyz[1], j), DOUBLETHRESHOLD);
    };

    for (int stru = 0; stru < supercell_lib.size(); ++stru)
    {
        ModuleSymmetry::Symmetry symm;
        construct_ucell(supercell_lib[stru]);
        symm.analy_sys(ucell, ofs_running);

        ModuleBase::matrix force(ucell.nat, 3, true);
        //generate random number for force and restrict to [-100,100)
        srand(time(NULL));
        for (int i = 0;i < ucell.nat;++i)
            for (int j = 0;j < 3;++j)
                force(i, j) = double(rand()) / double(RAND_MAX) * 200 - 100;

        // allocate pos
        std::vector<double> pos(ucell.nat * 3, 0.0);
        int iat = 0;
        for (int it = 0; it < ucell.ntype; it++)
        {
            for (int ia = 0; ia < ucell.atoms[it].na; ia++)
            {
                pos[3 * iat] = ucell.atoms[it].taud[ia].x;
                pos[3 * iat + 1] = ucell.atoms[it].taud[ia].y;
                pos[3 * iat + 2] = ucell.atoms[it].taud[ia].z;
                for (int k = 0; k < 3; ++k)
                {
                    symm.check_translation(pos[iat * 3 + k], -floor(pos[iat * 3 + k]));
                    symm.check_boundary(pos[iat * 3 + k]);
                }
                iat++;
            }
        }
        symm.force_symmetry(force, pos.data(), ucell);
        check_force(supercell_lib[stru], force);
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
