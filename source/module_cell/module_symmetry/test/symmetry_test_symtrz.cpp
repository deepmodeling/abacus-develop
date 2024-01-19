#include "symmetry_test_cases.h"
#include "mpi.h"

/************************************************
 *  unit test of class Symmetry
 * 4. function: `symmetrize_vec3_nat`
 * 5. function `symmetrize_mat3`
 *
***********************************************/
// mock the useless functions
void output::printM3(std::ofstream& ofs, const std::string& description, const ModuleBase::Matrix3& m) {}
pseudo::pseudo() {}
pseudo::~pseudo() {}
Atom::Atom() {}
Atom::~Atom() {}
Atom_pseudo::Atom_pseudo() {}
Atom_pseudo::~Atom_pseudo() {}
UnitCell::UnitCell() {}
UnitCell::~UnitCell() {}
Magnetism::Magnetism() {}
Magnetism::~Magnetism() {}

inline std::vector<double> allocate_pos(ModuleSymmetry::Symmetry& symm, UnitCell& ucell)
{
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
            ++iat;
        }
    }
    return pos;
}

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
        symm.analy_sys(ucell.lat, ucell.st, ucell.atoms, ofs_running);

        ModuleBase::matrix force(ucell.nat, 3, true);
        //generate random number for force and restrict to [-100,100)
        for (int i = 0;i < ucell.nat;++i)
            for (int j = 0;j < 3;++j)
                force(i, j) = double(rand()) / double(RAND_MAX) * 200 - 100;

        std::vector<double> pos = allocate_pos(symm, ucell);
        symm.symmetrize_vec3_nat(force.c);
        check_force(supercell_lib[stru], force);
    }
}

TEST_F(SymmetryTest, StressSymmetry)
{
    auto check_stress = [](stru_& conf, ModuleBase::matrix& stress)
        {
            // 1. check zeros
            for (auto elm : conf.stress_zero)
                EXPECT_NEAR(stress(elm.first, elm.second), 0.0, DOUBLETHRESHOLD);
            // 2. check equals
            for (auto eq_set : conf.stress_eq)
                for (int i = 1;i < eq_set.size();++i)
                    EXPECT_NEAR(stress(eq_set[i].first, eq_set[i].second), stress(eq_set[0].first, eq_set[0].second), DOUBLETHRESHOLD);
        };

    for (int stru = 0; stru < supercell_lib.size(); ++stru)
    {
        ModuleSymmetry::Symmetry symm;
        construct_ucell(supercell_lib[stru]);
        symm.analy_sys(ucell.lat, ucell.st, ucell.atoms, ofs_running);

        ModuleBase::matrix stress(3, 3, true);
        //generate random number for stress and restrict to [-1e5,1e5)
        for (int i = 0;i < 3;++i)
            for (int j = 0;j < 3;++j)
                stress(i, j) = double(rand()) / double(RAND_MAX) * 2e5 - 1e5;

        symm.symmetrize_mat3(stress, ucell.lat);
        check_stress(supercell_lib[stru], stress);

    }
}

int main(int argc, char** argv)
{
    srand(time(NULL));  // for random number generator
    MPI_Init(&argc, &argv);
    testing::InitGoogleTest(&argc, argv);
    int result = RUN_ALL_TESTS();
    MPI_Finalize();
    return result;
}
