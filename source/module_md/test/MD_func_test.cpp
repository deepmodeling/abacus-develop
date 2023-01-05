#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include "setcell.h"
#include "module_md/MD_func.h"
#include "module_esolver/esolver_lj.h"

#define doublethreshold 1e-12

class MD_func_test : public testing::Test
{
protected:
    UnitCell ucell;
    double *allmass;                     // atom mass 
    ModuleBase::Vector3<double> *pos;    // atom position
    ModuleBase::Vector3<double> *vel;    // atom velocity
    ModuleBase::Vector3<int> *ionmbl;    // atom is frozen or not
    ModuleBase::Vector3<double> *force;  // atom force
    ModuleBase::matrix virial;           // virial for this lattice
    ModuleBase::matrix stress;           // stress for this lattice
    double potential;                    // potential energy
    int natom;                           // atom number
    double temperature;                  // temperature
    int frozen_freedom;                  // frozen_freedom

    void SetUp()
    {
        Setcell::setupcell(ucell);
        Setcell::parameters();
        natom = ucell.nat;
        allmass = new double [natom];
        pos = new ModuleBase::Vector3<double> [natom];
        ionmbl = new ModuleBase::Vector3<int> [natom];
        vel = new ModuleBase::Vector3<double> [natom];
        force = new ModuleBase::Vector3<double> [natom];
        stress.create(3,3);
        virial.create(3,3);
    }

    void TearDown()
    {
        delete[] allmass;
        delete[] pos;
        delete[] vel;
        delete[] ionmbl;
        delete[] force;
    }
};

TEST_F(MD_func_test, gaussrand)
{
    EXPECT_DOUBLE_EQ(MD_func::gaussrand(), 1.1122716058967226);
    EXPECT_DOUBLE_EQ(MD_func::gaussrand(), -0.34532367182326629);
    EXPECT_DOUBLE_EQ(MD_func::gaussrand(), 0.60805637857480721);
}

TEST_F(MD_func_test, initpos)
{
    MD_func::InitPos(ucell, pos);
    
    EXPECT_DOUBLE_EQ(pos[0].x, 0.0);
    EXPECT_DOUBLE_EQ(pos[0].y, 0.0);
    EXPECT_DOUBLE_EQ(pos[0].z, 0.0);
    EXPECT_DOUBLE_EQ(pos[1].x, 5.2);
    EXPECT_DOUBLE_EQ(pos[1].y, 5.2);
    EXPECT_DOUBLE_EQ(pos[1].z, 0.0);
    EXPECT_DOUBLE_EQ(pos[2].x, 5.1);
    EXPECT_DOUBLE_EQ(pos[2].y, 0.0);
    EXPECT_DOUBLE_EQ(pos[2].z, 5.0);
    EXPECT_DOUBLE_EQ(pos[3].x, 0.0);
    EXPECT_DOUBLE_EQ(pos[3].y, 5.3);
    EXPECT_DOUBLE_EQ(pos[3].z, 5.0);
}

TEST_F(MD_func_test, randomvel)
{
    ucell.init_vel = 0;
    temperature = 300 / ModuleBase::Hartree_to_K;
    MD_func::InitVel(ucell, temperature, allmass, frozen_freedom, ionmbl, vel);
    
    EXPECT_NEAR(vel[0].x, 9.9105892783200826e-06, doublethreshold);
    EXPECT_NEAR(vel[0].y, -3.343699576563167e-05, doublethreshold);
    EXPECT_NEAR(vel[0].z, 9.385130426808701e-05, doublethreshold);
    EXPECT_NEAR(vel[1].x, -0.00017919300771203808, doublethreshold);
    EXPECT_NEAR(vel[1].y, 5.7074002254799079e-05, doublethreshold);
    EXPECT_NEAR(vel[1].z, -3.1088136026582953e-05, doublethreshold);
    EXPECT_NEAR(vel[2].x, 0.000141316492668737, doublethreshold);
    EXPECT_NEAR(vel[2].y, -0.00015841124290501442, doublethreshold);
    EXPECT_NEAR(vel[2].z, 1.900921882689748e-05, doublethreshold);
    EXPECT_NEAR(vel[3].x, 2.7965925764981002e-05, doublethreshold);
    EXPECT_NEAR(vel[3].y, 0.00013477423641584702, doublethreshold);
    EXPECT_NEAR(vel[3].z, -8.177238706840153e-05, doublethreshold);
}

TEST_F(MD_func_test, getmassmbl)
{
    ucell.init_vel = 0;
    temperature = 300 / ModuleBase::Hartree_to_K;
    MD_func::InitVel(ucell, temperature, allmass, frozen_freedom, ionmbl, vel);
    
    for(int i=0; i<natom; ++i)
    {
        EXPECT_DOUBLE_EQ(allmass[i], 39.948 / ModuleBase::AU_to_MASS);
        EXPECT_TRUE(ionmbl[i].x == 1);
        EXPECT_TRUE(ionmbl[i].y == 1);
        EXPECT_TRUE(ionmbl[i].z == 1);
    }
    
    EXPECT_TRUE(frozen_freedom == 3);
}

TEST_F(MD_func_test, readvel)
{
    MD_func::ReadVel(ucell, vel);
    
    EXPECT_DOUBLE_EQ(vel[0].x, -0.0001320807363640);
    EXPECT_DOUBLE_EQ(vel[0].y, 7.13429429835e-05);
    EXPECT_DOUBLE_EQ(vel[0].z, -1.40179977966e-05);
    EXPECT_DOUBLE_EQ(vel[1].x, 0.000153039878532);
    EXPECT_DOUBLE_EQ(vel[1].y, -0.000146533266608);
    EXPECT_DOUBLE_EQ(vel[1].z, 9.64491480698e-05);
    EXPECT_DOUBLE_EQ(vel[2].x, -0.000133789480226);
    EXPECT_DOUBLE_EQ(vel[2].y, -3.0451038112e-06);
    EXPECT_DOUBLE_EQ(vel[2].z, -5.40998380137e-05);
    EXPECT_DOUBLE_EQ(vel[3].x, 0.000112830338059);
    EXPECT_DOUBLE_EQ(vel[3].y, 7.82354274358e-05);
    EXPECT_DOUBLE_EQ(vel[3].z, -2.83313122596e-05);
}

TEST_F(MD_func_test, compute_stress)
{
    MD_func::InitVel(ucell, temperature, allmass, frozen_freedom, ionmbl, vel);
    MD_func::compute_stress(ucell, vel, allmass, virial, stress);
    EXPECT_DOUBLE_EQ(stress(0,0), 5.2064533063673623e-06);
    EXPECT_DOUBLE_EQ(stress(0,1), -1.6467487572445481e-06);
    EXPECT_DOUBLE_EQ(stress(0,2), 1.5039983732220751e-06);
    EXPECT_DOUBLE_EQ(stress(1,0), -1.6467487572445481e-06);
    EXPECT_DOUBLE_EQ(stress(1,1), 2.3806464376131247e-06);
    EXPECT_DOUBLE_EQ(stress(1,2), -1.251414906590483e-06);
    EXPECT_DOUBLE_EQ(stress(2,0), 1.5039983732220751e-06);
    EXPECT_DOUBLE_EQ(stress(2,1), -1.251414906590483e-06);
    EXPECT_DOUBLE_EQ(stress(2,2), 9.6330189688582584e-07);
}

TEST_F(MD_func_test, MDdump)
{
    MD_func::InitPos(ucell, pos);
    ModuleESolver::ESolver *p_esolver = new ModuleESolver::ESolver_LJ();
    p_esolver->Init(INPUT, ucell);
    MD_func::force_virial(p_esolver, 0, ucell, potential, force, virial);

    MD_func::MDdump(0, ucell, virial, force);
    std::ifstream ifs("MD_dump");
    std::string output_str;
    getline(ifs,output_str);
    EXPECT_THAT(output_str,testing::HasSubstr("MDSTEP:  0"));
    getline(ifs,output_str);
    EXPECT_THAT(output_str,testing::HasSubstr("LATTICE_CONSTANT: 1.000000000000"));
    getline(ifs,output_str);
    EXPECT_THAT(output_str,testing::HasSubstr("LATTICE_VECTORS"));
    getline(ifs,output_str);
    EXPECT_THAT(output_str,testing::HasSubstr("  10.000000000000  0.000000000000  0.000000000000"));
    getline(ifs,output_str);
    EXPECT_THAT(output_str,testing::HasSubstr("  0.000000000000  10.000000000000  0.000000000000"));
    getline(ifs,output_str);
    EXPECT_THAT(output_str,testing::HasSubstr("  0.000000000000  0.000000000000  10.000000000000"));
    getline(ifs,output_str);
    EXPECT_THAT(output_str,testing::HasSubstr("VIRIAL (KBAR)"));
    getline(ifs,output_str);
    EXPECT_THAT(output_str,testing::HasSubstr("  0.236428176219  0.050626986749  -0.000000000000"));
    getline(ifs,output_str);
    EXPECT_THAT(output_str,testing::HasSubstr("  0.050626986749  0.312766562817  0.000000000000"));
    getline(ifs,output_str);
    EXPECT_THAT(output_str,testing::HasSubstr("  -0.000000000000  0.000000000000  0.189105034410"));
    getline(ifs,output_str);
    EXPECT_THAT(output_str,testing::HasSubstr("INDEX    LABEL    POSITIONS    FORCE (eV/Angstrom)"));
    getline(ifs,output_str);
    EXPECT_THAT(output_str,testing::HasSubstr("  0  Ar  0.000000000000  0.000000000000  0.000000000000  0.025617329374  0.042288127307  -0.000000000000"));
    getline(ifs,output_str);
    EXPECT_THAT(output_str,testing::HasSubstr("  1  Ar  5.200000000000  5.200000000000  0.000000000000  -0.033300246168  -0.034414246659  -0.000000000000"));
    getline(ifs,output_str);
    EXPECT_THAT(output_str,testing::HasSubstr("  2  Ar  5.100000000000  0.000000000000  5.000000000000  -0.018209198893  0.041606045042  -0.000000000000"));
    getline(ifs,output_str);
    EXPECT_THAT(output_str,testing::HasSubstr("  3  Ar  0.000000000000  5.300000000000  5.000000000000  0.025892115687  -0.049479925690  0.000000000000"));
    ifs.close();

    // append
    MD_func::MDdump(1, ucell, virial, force);
    std::ifstream ifs2("MD_dump");
    getline(ifs2,output_str);
    EXPECT_THAT(output_str,testing::HasSubstr("MDSTEP:  0"));
    getline(ifs2,output_str);
    EXPECT_THAT(output_str,testing::HasSubstr("LATTICE_CONSTANT: 1.000000000000"));
    getline(ifs2,output_str);
    EXPECT_THAT(output_str,testing::HasSubstr("LATTICE_VECTORS"));
    getline(ifs2,output_str);
    EXPECT_THAT(output_str,testing::HasSubstr("  10.000000000000  0.000000000000  0.000000000000"));
    getline(ifs2,output_str);
    EXPECT_THAT(output_str,testing::HasSubstr("  0.000000000000  10.000000000000  0.000000000000"));
    getline(ifs2,output_str);
    EXPECT_THAT(output_str,testing::HasSubstr("  0.000000000000  0.000000000000  10.000000000000"));
    getline(ifs2,output_str);
    EXPECT_THAT(output_str,testing::HasSubstr("VIRIAL (KBAR)"));
    getline(ifs2,output_str);
    EXPECT_THAT(output_str,testing::HasSubstr("  0.236428176219  0.050626986749  -0.000000000000"));
    getline(ifs2,output_str);
    EXPECT_THAT(output_str,testing::HasSubstr("  0.050626986749  0.312766562817  0.000000000000"));
    getline(ifs2,output_str);
    EXPECT_THAT(output_str,testing::HasSubstr("  -0.000000000000  0.000000000000  0.189105034410"));
    getline(ifs2,output_str);
    EXPECT_THAT(output_str,testing::HasSubstr("INDEX    LABEL    POSITIONS    FORCE (eV/Angstrom)"));
    getline(ifs2,output_str);
    EXPECT_THAT(output_str,testing::HasSubstr("  0  Ar  0.000000000000  0.000000000000  0.000000000000  0.025617329374  0.042288127307  -0.000000000000"));
    getline(ifs2,output_str);
    EXPECT_THAT(output_str,testing::HasSubstr("  1  Ar  5.200000000000  5.200000000000  0.000000000000  -0.033300246168  -0.034414246659  -0.000000000000"));
    getline(ifs2,output_str);
    EXPECT_THAT(output_str,testing::HasSubstr("  2  Ar  5.100000000000  0.000000000000  5.000000000000  -0.018209198893  0.041606045042  -0.000000000000"));
    getline(ifs2,output_str);
    EXPECT_THAT(output_str,testing::HasSubstr("  3  Ar  0.000000000000  5.300000000000  5.000000000000  0.025892115687  -0.049479925690  0.000000000000"));
    getline(ifs2,output_str);
    getline(ifs2,output_str);
    getline(ifs2,output_str);
    EXPECT_THAT(output_str,testing::HasSubstr("MDSTEP:  1"));
    getline(ifs2,output_str);
    EXPECT_THAT(output_str,testing::HasSubstr("LATTICE_CONSTANT: 1.000000000000"));
    getline(ifs2,output_str);
    EXPECT_THAT(output_str,testing::HasSubstr("LATTICE_VECTORS"));
    getline(ifs2,output_str);
    EXPECT_THAT(output_str,testing::HasSubstr("  10.000000000000  0.000000000000  0.000000000000"));
    getline(ifs2,output_str);
    EXPECT_THAT(output_str,testing::HasSubstr("  0.000000000000  10.000000000000  0.000000000000"));
    getline(ifs2,output_str);
    EXPECT_THAT(output_str,testing::HasSubstr("  0.000000000000  0.000000000000  10.000000000000"));
    getline(ifs2,output_str);
    EXPECT_THAT(output_str,testing::HasSubstr("VIRIAL (KBAR)"));
    getline(ifs2,output_str);
    EXPECT_THAT(output_str,testing::HasSubstr("  0.236428176219  0.050626986749  -0.000000000000"));
    getline(ifs2,output_str);
    EXPECT_THAT(output_str,testing::HasSubstr("  0.050626986749  0.312766562817  0.000000000000"));
    getline(ifs2,output_str);
    EXPECT_THAT(output_str,testing::HasSubstr("  -0.000000000000  0.000000000000  0.189105034410"));
    getline(ifs2,output_str);
    EXPECT_THAT(output_str,testing::HasSubstr("INDEX    LABEL    POSITIONS    FORCE (eV/Angstrom)"));
    getline(ifs2,output_str);
    EXPECT_THAT(output_str,testing::HasSubstr("  0  Ar  0.000000000000  0.000000000000  0.000000000000  0.025617329374  0.042288127307  -0.000000000000"));
    getline(ifs2,output_str);
    EXPECT_THAT(output_str,testing::HasSubstr("  1  Ar  5.200000000000  5.200000000000  0.000000000000  -0.033300246168  -0.034414246659  -0.000000000000"));
    getline(ifs2,output_str);
    EXPECT_THAT(output_str,testing::HasSubstr("  2  Ar  5.100000000000  0.000000000000  5.000000000000  -0.018209198893  0.041606045042  -0.000000000000"));
    getline(ifs2,output_str);
    EXPECT_THAT(output_str,testing::HasSubstr("  3  Ar  0.000000000000  5.300000000000  5.000000000000  0.025892115687  -0.049479925690  0.000000000000"));
    ifs2.close();

    remove("MD_dump");
}

TEST_F(MD_func_test, outStress)
{
    MD_func::InitPos(ucell, pos);
    MD_func::ReadVel(ucell, vel);
    ModuleESolver::ESolver *p_esolver = new ModuleESolver::ESolver_LJ();
    p_esolver->Init(INPUT, ucell);
    MD_func::force_virial(p_esolver, 0, ucell, potential, force, virial);
    MD_func::compute_stress(ucell, vel, allmass, virial, stress);
    GlobalV::ofs_running.open("running.log");
    MD_func::outStress(virial, stress);

    std::ifstream ifs("running.log");
    std::string output_str;
    getline(ifs,output_str);
    getline(ifs,output_str);
    EXPECT_THAT(output_str,testing::HasSubstr("output Pressure for check!"));
    getline(ifs,output_str);
    EXPECT_THAT(output_str,testing::HasSubstr("Virtual Pressure is 0.2461 Kbar "));
    getline(ifs,output_str);
    EXPECT_THAT(output_str,testing::HasSubstr("Virial Term is 0.2461 Kbar "));
    getline(ifs,output_str);
    EXPECT_THAT(output_str,testing::HasSubstr("Kenetic Term is 2.66338e-10 Kbar "));
    getline(ifs,output_str);
    getline(ifs,output_str);
    getline(ifs,output_str);
    EXPECT_THAT(output_str,testing::HasSubstr(" ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><"));
    getline(ifs,output_str);
    getline(ifs,output_str);
    EXPECT_THAT(output_str,testing::HasSubstr(" MD STRESS (KBAR)"));
    getline(ifs,output_str);
    getline(ifs,output_str);
    EXPECT_THAT(output_str,testing::HasSubstr(" ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><"));
    getline(ifs,output_str);
    getline(ifs,output_str);
    EXPECT_THAT(output_str,testing::HasSubstr("     0.23642818    0.050626987 -1.7378254e-10"));
    getline(ifs,output_str);
    EXPECT_THAT(output_str,testing::HasSubstr("    0.050626987     0.31276656 -1.0047995e-10"));
    getline(ifs,output_str);
    EXPECT_THAT(output_str,testing::HasSubstr(" -1.7378254e-10 -1.0047995e-10     0.18910503"));

    ifs.close();
    remove("running.log");
}