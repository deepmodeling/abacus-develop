#include "module_cell/unitcell.h"
#include "module_cell/setup_nonlocal.h"
#include "module_io/input.h"
#include "module_base/mathzone.h"
#include "module_base/vector3.h"
#include"gtest/gtest.h"
#include"gmock/gmock.h"
#include "mpi.h"
#define private public
#include "module_hamilt_general/module_vdw/vdwd2_parameters.h"
#include "module_hamilt_general/module_vdw/vdwd3_parameters.h"
#include "module_hamilt_general/module_vdw/vdwd2.h"
#include "module_hamilt_general/module_vdw/vdwd3.h"
#include "module_hamilt_general/module_vdw/vdw.h"


/************************************************
*  unit test of class VDW and related functions
***********************************************/

/**
* - Tested functions:
*   - vdw::make_vdw():
*       Based on the value of INPUT.vdw_method, construct 
*       Vdwd2 or Vdwd3 class, and do the initialization.
*   - vdw::get_energy()/vdw::get_force()/vdw::get_stress():
*       Calculate the VDW (d2, d3_0 and d3_bj types) enerygy, force, stress.    
*   - Vdwd2Parameters::initial_parameters()
*   - Vdwd3Parameters::initial_parameters()
*/

pseudo_nc::pseudo_nc()
{
}
pseudo_nc::~pseudo_nc()
{
}
Atom::Atom()
{
}
Atom::~Atom()
{
}
Atom_pseudo::Atom_pseudo()
{
}
Atom_pseudo::~Atom_pseudo()
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
InfoNonlocal::InfoNonlocal(){}
InfoNonlocal::~InfoNonlocal(){}

struct atomtype_
{
    std::string atomname;
    std::vector<std::vector<double>> coordinate;
};

struct stru_
{
    std::vector<double> cell;
    std::vector<atomtype_> all_type;
};

void construct_ucell(stru_ &stru, UnitCell &ucell)
{
    std::vector<atomtype_> coord = stru.all_type;
    ucell.a1 = ModuleBase::Vector3<double>(stru.cell[0], stru.cell[1], stru.cell[2]);
    ucell.a2 = ModuleBase::Vector3<double>(stru.cell[3], stru.cell[4], stru.cell[5]);
    ucell.a3 = ModuleBase::Vector3<double>(stru.cell[6], stru.cell[7], stru.cell[8]);
    ucell.latvec.e11 = stru.cell[0]; ucell.latvec.e12 = stru.cell[1]; ucell.latvec.e13 = stru.cell[2];
    ucell.latvec.e21 = stru.cell[3]; ucell.latvec.e22 = stru.cell[4]; ucell.latvec.e23 = stru.cell[5];
    ucell.latvec.e31 = stru.cell[6]; ucell.latvec.e32 = stru.cell[7]; ucell.latvec.e33 = stru.cell[8];
    ucell.lat0 = 10.2;
    ucell.ntype = stru.all_type.size();
    ucell.atoms = new Atom[ucell.ntype];
    ucell.nat = 0;
    int nmax = 0;
    for (int i = 0; i < coord.size(); i++)
    {
        ucell.atoms[i].label = coord[i].atomname;
        ucell.atoms[i].ncpp.psd = coord[i].atomname;
        ucell.atoms[i].na = coord[i].coordinate.size();
        ucell.atoms[i].tau = new ModuleBase::Vector3<double>[ucell.atoms[i].na];
        ucell.atoms[i].taud = new ModuleBase::Vector3<double>[ucell.atoms[i].na];
        for (int j = 0; j < ucell.atoms[i].na; j++)
        {
            std::vector<double> this_atom = coord[i].coordinate[j];
            ucell.atoms[i].tau[j] = ModuleBase::Vector3<double>(this_atom[0], this_atom[1], this_atom[2]);
            ModuleBase::Mathzone::Cartesian_to_Direct(ucell.atoms[i].tau[j].x,
                                                      ucell.atoms[i].tau[j].y,
                                                      ucell.atoms[i].tau[j].z,
                                                      ucell.a1.x,
                                                      ucell.a1.y,
                                                      ucell.a1.z,
                                                      ucell.a2.x,
                                                      ucell.a2.y,
                                                      ucell.a2.z,
                                                      ucell.a3.x,
                                                      ucell.a3.y,
                                                      ucell.a3.z,
                                                      ucell.atoms[i].taud[j].x,
                                                      ucell.atoms[i].taud[j].y,
                                                      ucell.atoms[i].taud[j].z);
        }
        ucell.nat += ucell.atoms[i].na;
        if (ucell.atoms[i].na > nmax)
        {
            nmax = ucell.atoms[i].na;
        }
    }

    ucell.omega = std::abs(ucell.latvec.Det()) * ucell.lat0 * ucell.lat0 * ucell.lat0;
    ucell.itia2iat.create(ucell.ntype, nmax);
	int iat=0;
	for(int it = 0;it < ucell.ntype;it++)
	{
		for(int ia=0; ia<ucell.atoms[it].na; ia++)
		{
			ucell.itia2iat(it, ia) = iat;
            ++iat;
		}	
	}
}

void ClearUcell(UnitCell &ucell)
{
    for (int i = 0; i < ucell.ntype; i++)
    {
        delete[] ucell.atoms[i].tau;
        delete[] ucell.atoms[i].taud;
    }
    delete[] ucell.atoms;
}

class vdwd2Test: public testing::Test
{
    protected:
    UnitCell ucell;
    Input input;

    void SetUp(){
        stru_ structure{std::vector<double>{0.5,0.5,0.0,0.5,0.0,0.5,0.0,0.5,0.5},
                        std::vector<atomtype_>{atomtype_{"Si",
                                                         std::vector<std::vector<double>>{
                                                             {0., 0., 0.},
                                                             {0.3, 0.25, 0.25}
                                                         }}}};
        construct_ucell(structure,ucell);

        input.vdw_method = "d2";
        input.vdw_s6 = "0.75";
        input.vdw_d  = 20;
        input.vdw_C6_file = "default";
        input.vdw_C6_unit = "Jnm6/mol";
        input.vdw_R0_file = "default";
        input.vdw_R0_unit = "A";
        input.vdw_cutoff_type = "radius";
        input.vdw_cutoff_radius = "56.6918";
        input.vdw_radius_unit = "Bohr";
    }

    void TearDown(){
        ClearUcell(ucell);
    }
};

TEST_F(vdwd2Test, D2Default)
{
    std::unique_ptr<vdw::Vdwd2> vdwd2_test = std::make_unique<vdw::Vdwd2>(ucell);

    vdwd2_test->parameter().initial_parameters(input);
    vdwd2_test->parameter().initset(ucell);
    EXPECT_EQ(vdwd2_test->parameter().scaling(), 0.75);
    EXPECT_EQ(vdwd2_test->parameter().damping(), input.vdw_d);
    EXPECT_EQ(vdwd2_test->parameter().model(), input.vdw_cutoff_type);
    EXPECT_EQ(vdwd2_test->parameter().radius_, 56.6918);
    double Si_C6 = 9.23*1e6 / (ModuleBase::ELECTRONVOLT_SI * ModuleBase::NA) / pow(ModuleBase::BOHR_TO_A, 6)/ ModuleBase::Ry_to_eV;
    EXPECT_NEAR(vdwd2_test->parameter().C6_["Si"], Si_C6,1e-13);
    EXPECT_EQ(vdwd2_test->parameter().R0_["Si"], 1.716/ModuleBase::BOHR_TO_A);
    EXPECT_EQ(vdwd2_test->parameter().period().x, 2 * ceil(56.6918 / ucell.lat0 / sqrt(ucell.a1.norm2())) + 1);
    EXPECT_EQ(vdwd2_test->parameter().period().y, 2 * ceil(56.6918 / ucell.lat0 / sqrt(ucell.a2.norm2())) + 1);
    EXPECT_EQ(vdwd2_test->parameter().period().z, 2 * ceil(56.6918 / ucell.lat0 / sqrt(ucell.a3.norm2())) + 1);
}

TEST_F(vdwd2Test, WrongVdwType)
{
    input.vdw_method = "d2d3";
    EXPECT_EQ(vdw::make_vdw(ucell, input),nullptr);
}

TEST_F(vdwd2Test, OneAtomWarning)
{
    UnitCell ucell1;
    stru_ structure1{std::vector<double>{0.5, 0.5, 0.0, 0.5, 0.0, 0.5, 0.0, 0.5, 0.5},
                     std::vector<atomtype_>{atomtype_{"Si", std::vector<std::vector<double>>{{0., 0., 0.}}}}};
    construct_ucell(structure1,ucell1);
    GlobalV::ofs_warning.open("warning.log");
    std::ifstream ifs;
    std::string output;
 
    std::unique_ptr<vdw::Vdw> vdw_test = vdw::make_vdw(ucell1, input);

    GlobalV::ofs_warning.close();
	ifs.open("warning.log");
	getline(ifs,output);
	EXPECT_THAT(output,testing::HasSubstr("warning"));
    EXPECT_EQ(vdw_test,nullptr);

    ifs.close();
    ClearUcell(ucell1);
}

TEST_F(vdwd2Test, D2ReadFile)
{
    input.vdw_C6_file = "c6.txt";
    input.vdw_R0_file = "r0.txt";
    std::unique_ptr<vdw::Vdwd2> vdwd2_test = std::make_unique<vdw::Vdwd2>(ucell);
    
    vdwd2_test->parameter().initial_parameters(input);
    double Si_C6 = 9.13*1e6 / (ModuleBase::ELECTRONVOLT_SI * ModuleBase::NA) / pow(ModuleBase::BOHR_TO_A, 6)/ ModuleBase::Ry_to_eV;
    EXPECT_NEAR(vdwd2_test->parameter().C6_["Si"], Si_C6,1e-13);
    EXPECT_EQ(vdwd2_test->parameter().R0_["Si"], 1.626/ModuleBase::BOHR_TO_A);
}

TEST_F(vdwd2Test, D2ReadFileError)
{
    input.vdw_C6_file = "c6_wrong.txt";
    input.vdw_R0_file = "r0_wrong.txt";
    std::unique_ptr<vdw::Vdwd2> vdwd2_test = std::make_unique<vdw::Vdwd2>(ucell);
    
    testing::internal::CaptureStdout();
    EXPECT_EXIT(vdwd2_test->parameter().C6_input(input.vdw_C6_file, input.vdw_C6_unit), ::testing::ExitedWithCode(0), "");
    EXPECT_EXIT(vdwd2_test->parameter().R0_input(input.vdw_R0_file, input.vdw_R0_unit), ::testing::ExitedWithCode(0), "");
    std::string output = testing::internal::GetCapturedStdout();
}

TEST_F(vdwd2Test, D2c6UniteVA6)
{
    input.vdw_C6_unit = "eVA6";
    std::unique_ptr<vdw::Vdwd2> vdwd2_test = std::make_unique<vdw::Vdwd2>(ucell);

    vdwd2_test->parameter().initial_parameters(input);
    double Si_C6 = 9.23 / pow(ModuleBase::BOHR_TO_A, 6) * ModuleBase::Ry_to_eV;
    EXPECT_NEAR(vdwd2_test->parameter().C6_["Si"], Si_C6,1e-13);
}

TEST_F(vdwd2Test, D2r0UnitBohr)
{
    input.vdw_R0_unit = "Bohr";
    std::unique_ptr<vdw::Vdwd2> vdwd2_test = std::make_unique<vdw::Vdwd2>(ucell);

    vdwd2_test->parameter().initial_parameters(input);
    EXPECT_EQ(vdwd2_test->parameter().R0_["Si"], 1.716);
}

TEST_F(vdwd2Test, D2WrongUnit)
{
    input.vdw_R0_unit = "B";
    input.vdw_C6_unit = "eV";
    std::unique_ptr<vdw::Vdwd2> vdwd2_test = std::make_unique<vdw::Vdwd2>(ucell);

    testing::internal::CaptureStdout();
    EXPECT_EXIT(vdwd2_test->parameter().C6_input(input.vdw_C6_file, input.vdw_C6_unit), ::testing::ExitedWithCode(0), "");
    EXPECT_EXIT(vdwd2_test->parameter().R0_input(input.vdw_R0_file, input.vdw_R0_unit), ::testing::ExitedWithCode(0), "");
    std::string output = testing::internal::GetCapturedStdout();
}

TEST_F(vdwd2Test, D2RadiusUnitAngstrom)
{
    input.vdw_cutoff_radius = "56.6918";
    input.vdw_radius_unit = "Angstrom";
    
    std::unique_ptr<vdw::Vdwd2> vdwd2_test = std::make_unique<vdw::Vdwd2>(ucell);
    vdwd2_test->parameter().initial_parameters(input);
    EXPECT_EQ(vdwd2_test->parameter().radius_, 56.6918/ModuleBase::BOHR_TO_A);
}

TEST_F(vdwd2Test, D2CutoffTypePeriod)
{
    input.vdw_cutoff_type = "period";
    input.vdw_cutoff_period = {3,3,3};
    
    std::unique_ptr<vdw::Vdwd2> vdwd2_test = std::make_unique<vdw::Vdwd2>(ucell);
    vdwd2_test->parameter().initial_parameters(input);
    EXPECT_EQ(vdwd2_test->parameter().period(), input.vdw_cutoff_period);
}

TEST_F(vdwd2Test, D2R0ZeroQuit)
{   
    std::unique_ptr<vdw::Vdwd2> vdwd2_test = std::make_unique<vdw::Vdwd2>(ucell);
    vdwd2_test->parameter().initial_parameters(input);
    vdwd2_test->parameter().R0_["Si"] = 0.0;
    
    testing::internal::CaptureStdout();
    EXPECT_EXIT(vdwd2_test->get_energy(), ::testing::ExitedWithCode(0), "");
    std::string output = testing::internal::GetCapturedStdout();
}

TEST_F(vdwd2Test, D2GetEnergy)
{   
    auto vdw_solver = vdw::make_vdw(ucell, input);
    double ene = vdw_solver->get_energy();
    EXPECT_NEAR(ene,-0.034526673470525196,1E-10);
}

TEST_F(vdwd2Test, D2GetForce)
{   
    auto vdw_solver = vdw::make_vdw(ucell, input);
    std::vector<ModuleBase::Vector3<double>> force = vdw_solver->get_force();
    EXPECT_NEAR(force[0].x, -0.00078824525563651242,1e-12);
    EXPECT_NEAR(force[0].y,  2.6299822052061785e-08,1e-12);
    EXPECT_NEAR(force[0].z,  2.6299822050796364e-08,1e-12);
    EXPECT_NEAR(force[1].x,  0.00078824525563651025,1e-12);
    EXPECT_NEAR(force[1].y, -2.6299822051961511e-08,1e-12);
    EXPECT_NEAR(force[1].z, -2.6299822050777911e-08,1e-12);
}

TEST_F(vdwd2Test, D2GetStress)
{   
    auto vdw_solver = vdw::make_vdw(ucell, input);
    ModuleBase::Matrix3 stress = vdw_solver->get_stress();
    EXPECT_NEAR(stress.e11, -0.00020532319044269705,1e-12);
    EXPECT_NEAR(stress.e12, -3.5642821939401251e-08,1e-12);
    EXPECT_NEAR(stress.e13, -3.5642821939437223e-08,1e-12);
    EXPECT_NEAR(stress.e21, -3.5642821939401251e-08,1e-12);
    EXPECT_NEAR(stress.e22, -0.00020491799872060745,1e-12);
    EXPECT_NEAR(stress.e23, -4.4369449138539291e-06,1e-12);
    EXPECT_NEAR(stress.e31, -3.5642821939437223e-08,1e-12);
    EXPECT_NEAR(stress.e32, -4.4369449138539291e-06,1e-12);
    EXPECT_NEAR(stress.e33, -0.00020491799872060734,1e-12);
}



class vdwd3Test: public testing::Test
{
    protected:
    UnitCell ucell;
    Input input;

    void SetUp(){
        stru_ structure{std::vector<double>{0.5,0.5,0.0,0.5,0.0,0.5,0.0,0.5,0.5},
                        std::vector<atomtype_>{atomtype_{"Si",
                                                         std::vector<std::vector<double>>{
                                                             {0., 0., 0.},
                                                             {0.3, 0.25, 0.25}
                                                         }}}};
        construct_ucell(structure,ucell);

        input.vdw_method = "d3_0";
        input.vdw_s6 = "1.0";
        input.vdw_s8 = "0.7875";
        input.vdw_a1 = "0.4289";
        input.vdw_a2 = "4.4407";
        input.vdw_abc = false;
        input.vdw_cutoff_type = "radius";
        input.vdw_radius_unit = "Bohr";
        input.vdw_cutoff_radius = "95";
        input.vdw_cn_thr_unit = "Bohr";
        input.vdw_cn_thr = 40;
        input.vdw_cutoff_period = {3,3,3};
    }

    void TearDown(){
        ClearUcell(ucell);
    }
};

TEST_F(vdwd3Test, D30Default)
{
    std::unique_ptr<vdw::Vdwd3> vdwd3_test = std::make_unique<vdw::Vdwd3>(ucell);

    vdwd3_test->parameter().initial_parameters(input);

    EXPECT_EQ(vdwd3_test->parameter().s6(), 1.0);
    EXPECT_EQ(vdwd3_test->parameter().s18(), 0.7875);
    EXPECT_EQ(vdwd3_test->parameter().rs6(), 0.4289);
    EXPECT_EQ(vdwd3_test->parameter().rs18(), 4.4407);
    EXPECT_EQ(vdwd3_test->parameter().abc(), false);
    EXPECT_EQ(vdwd3_test->parameter().version(), "d3_0");
    EXPECT_EQ(vdwd3_test->parameter().model(), "radius");
    EXPECT_EQ(vdwd3_test->parameter().rthr2(), std::pow(95, 2));
    EXPECT_EQ(vdwd3_test->parameter().cn_thr2(), std::pow(40, 2));   
}

TEST_F(vdwd3Test, D30UnitA)
{
    input.vdw_radius_unit = "A";
    input.vdw_cn_thr_unit = "A";
    std::unique_ptr<vdw::Vdwd3> vdwd3_test = std::make_unique<vdw::Vdwd3>(ucell);

    vdwd3_test->parameter().initial_parameters(input);

    EXPECT_EQ(vdwd3_test->parameter().rthr2(), std::pow(95/ModuleBase::BOHR_TO_A, 2));
    EXPECT_EQ(vdwd3_test->parameter().cn_thr2(), std::pow(40/ModuleBase::BOHR_TO_A, 2));   
}

TEST_F(vdwd3Test, D30Period)
{
    input.vdw_cutoff_type = "period";
    std::unique_ptr<vdw::Vdwd3> vdwd3_test = std::make_unique<vdw::Vdwd3>(ucell);

    vdwd3_test->parameter().initial_parameters(input);
    vdwd3_test->init();
    std::vector<int> rep_vdw_ref = {input.vdw_cutoff_period.x, input.vdw_cutoff_period.y, input.vdw_cutoff_period.z};

    EXPECT_EQ(vdwd3_test->parameter().period(), input.vdw_cutoff_period);
    EXPECT_EQ(vdwd3_test->rep_vdw_, rep_vdw_ref);  
}

TEST_F(vdwd3Test, D30GetEnergy)
{   
    auto vdw_solver = vdw::make_vdw(ucell, input);
    double ene = vdw_solver->get_energy();
    EXPECT_NEAR(ene,-0.20932367230529664,1E-10);
}

TEST_F(vdwd3Test, D30GetForce)
{   
    auto vdw_solver = vdw::make_vdw(ucell, input);
    std::vector<ModuleBase::Vector3<double>> force = vdw_solver->get_force();
    EXPECT_NEAR(force[0].x, -0.032450975169023302,1e-12);
    EXPECT_NEAR(force[0].y, 0.0,1e-12);
    EXPECT_NEAR(force[0].z, 0.0,1e-12);
    EXPECT_NEAR(force[1].x, 0.032450975169023302,1e-12);
    EXPECT_NEAR(force[1].y, 0.0,1e-12);
    EXPECT_NEAR(force[1].z, 0.0,1e-12);
}

TEST_F(vdwd3Test, D30GetStress)
{   
    auto vdw_solver = vdw::make_vdw(ucell, input);
    ModuleBase::Matrix3 stress = vdw_solver->get_stress();
    EXPECT_NEAR(stress.e11, -0.0011141545452036336,1e-12);
    EXPECT_NEAR(stress.e12, 0.0,1e-12);
    EXPECT_NEAR(stress.e13, 0.0,1e-12);
    EXPECT_NEAR(stress.e21, 0.0,1e-12);
    EXPECT_NEAR(stress.e22, -0.0012740017248971929,1e-12);
    EXPECT_NEAR(stress.e23, 0.00049503596239307496,1e-12);
    EXPECT_NEAR(stress.e31, 0.0,1e-12);
    EXPECT_NEAR(stress.e32, 0.00049503596239307496,1e-12);
    EXPECT_NEAR(stress.e33, -0.0012740017248971936,1e-12);
}

TEST_F(vdwd3Test, D3bjGetEnergy)
{  
    input.vdw_method = "d3_bj"; 
    auto vdw_solver = vdw::make_vdw(ucell, input);
    double ene = vdw_solver->get_energy();
    EXPECT_NEAR(ene,-0.047458675421836918,1E-10);
}

TEST_F(vdwd3Test, D3bjGetForce)
{   
    input.vdw_method = "d3_bj";
    auto vdw_solver = vdw::make_vdw(ucell, input);
    std::vector<ModuleBase::Vector3<double>> force = vdw_solver->get_force();
    EXPECT_NEAR(force[0].x, -0.0026006968781200602,1e-12);
    EXPECT_NEAR(force[0].y, 0.0,1e-12);
    EXPECT_NEAR(force[0].z, 0.0,1e-12);
    EXPECT_NEAR(force[1].x, 0.0026006968781200602,1e-12);
    EXPECT_NEAR(force[1].y, 0.0,1e-12);
    EXPECT_NEAR(force[1].z, 0.0,1e-12);
}

TEST_F(vdwd3Test, D3bjGetStress)
{   
    input.vdw_method = "d3_bj";
    auto vdw_solver = vdw::make_vdw(ucell, input);
    ModuleBase::Matrix3 stress = vdw_solver->get_stress();
    EXPECT_NEAR(stress.e11, -0.00014376286737216365,1e-12);
    EXPECT_NEAR(stress.e12, 0.0,1e-12);
    EXPECT_NEAR(stress.e13, 0.0,1e-12);
    EXPECT_NEAR(stress.e21, 0.0,1e-12);
    EXPECT_NEAR(stress.e22, -0.00015350088004991452,1e-12);
    EXPECT_NEAR(stress.e23, 1.8204947825641812e-05,1e-12);
    EXPECT_NEAR(stress.e31, 0.0,1e-12);
    EXPECT_NEAR(stress.e32, 1.8204947825641816e-05,1e-12);
    EXPECT_NEAR(stress.e33, -0.0001535008800499145,1e-12);
}


class vdwd3abcTest: public testing::Test
{
    protected:
    UnitCell ucell;
    Input input;

    void SetUp(){
        stru_ structure{
            std::vector<double>{0.75, 0.75, 0.0, 0.75, 0.0, 0.75, 0.0, 0.75, 0.75},
            std::vector<atomtype_>{atomtype_{"Si", std::vector<std::vector<double>>{{0., 0., 0.}, {0.3, 0.25, 0.25}}},
                                   atomtype_{"C", std::vector<std::vector<double>>{{0.5, 0.5, 0.5}}}}};
        construct_ucell(structure,ucell);

        input.vdw_method = "d3_0";
        input.vdw_s6 = "1.0";
        input.vdw_s8 = "0.7875";
        input.vdw_a1 = "0.4289";
        input.vdw_a2 = "4.4407";
        input.vdw_abc = true;
        input.vdw_cutoff_type = "radius";
        input.vdw_radius_unit = "Bohr";
        input.vdw_cutoff_radius = "95";
        input.vdw_cn_thr_unit = "Bohr";
        input.vdw_cn_thr = 40;
        input.vdw_cutoff_period = {3,3,3};
    }

    void TearDown(){
        ClearUcell(ucell);
    }
};


TEST_F(vdwd3abcTest, D30GetEnergy)
{  
    auto vdw_solver = vdw::make_vdw(ucell, input);
    double ene = vdw_solver->get_energy();
    EXPECT_NEAR(ene,-0.11487062308916372,1E-10);
}

TEST_F(vdwd3abcTest, D30GetForce)
{   
    auto vdw_solver = vdw::make_vdw(ucell, input);
    std::vector<ModuleBase::Vector3<double>> force = vdw_solver->get_force();
    EXPECT_NEAR(force[0].x, 0.030320738678429094,1e-12);
    EXPECT_NEAR(force[0].y, 0.025570534655235538,1e-12);
    EXPECT_NEAR(force[0].z, 0.025570534655235538,1e-12);
    EXPECT_NEAR(force[1].x, -0.0067036811361536061,1e-12);
    EXPECT_NEAR(force[1].y, 0.0037813111009633712,1e-12);
    EXPECT_NEAR(force[1].z, 0.0037813111009634614,1e-12);
}

TEST_F(vdwd3abcTest, D30GetStress)
{   
    auto vdw_solver = vdw::make_vdw(ucell, input);
    ModuleBase::Matrix3 stress = vdw_solver->get_stress();
    EXPECT_NEAR(stress.e11, -0.00023421562840819491,1e-12);
    EXPECT_NEAR(stress.e12, -0.00015112406243413323,1e-12);
    EXPECT_NEAR(stress.e13, -0.00015112406243413302,1e-12);
    EXPECT_NEAR(stress.e21, -0.00015112406243413323,1e-12);
    EXPECT_NEAR(stress.e22, -0.00023139547090668657,1e-12);
    EXPECT_NEAR(stress.e23, -0.00014931418741042754,1e-12);
    EXPECT_NEAR(stress.e31, -0.00015112406243413302,1e-12);
    EXPECT_NEAR(stress.e32, -0.00014931418741042754,1e-12);
    EXPECT_NEAR(stress.e33, -0.00023139547090668714,1e-12);
}

TEST_F(vdwd3abcTest, D3bjGetEnergy)
{  
    input.vdw_method = "d3_bj";
    auto vdw_solver = vdw::make_vdw(ucell, input);
    double ene = vdw_solver->get_energy();
    EXPECT_NEAR(ene,-0.030667806197006021,1E-10);
}

TEST_F(vdwd3abcTest, D3bjGetForce)
{   
    input.vdw_method = "d3_bj";
    auto vdw_solver = vdw::make_vdw(ucell, input);
    std::vector<ModuleBase::Vector3<double>> force = vdw_solver->get_force();
    EXPECT_NEAR(force[0].x, -0.0010630099217696475,1e-12);
    EXPECT_NEAR(force[0].y, -0.0010031953309458587,1e-12);
    EXPECT_NEAR(force[0].z, -0.0010031953309458642,1e-12);
    EXPECT_NEAR(force[1].x, 0.00015471729604904047,1e-12);
    EXPECT_NEAR(force[1].y,-0.00010902508913277635,1e-12);
    EXPECT_NEAR(force[1].z, -0.00010902508913277528,1e-12);
}

TEST_F(vdwd3abcTest, D3bjGetStress)
{   
    input.vdw_method = "d3_bj";
    auto vdw_solver = vdw::make_vdw(ucell, input);
    ModuleBase::Matrix3 stress = vdw_solver->get_stress();
    EXPECT_NEAR(stress.e11, -3.3803329202372578e-05,1e-12);
    EXPECT_NEAR(stress.e12, 5.1291622417145846e-06,1e-12);
    EXPECT_NEAR(stress.e13, 5.1291622417145889e-06,1e-12);
    EXPECT_NEAR(stress.e21, 5.1291622417145863e-06,1e-12);
    EXPECT_NEAR(stress.e22, -3.427844212559098e-05,1e-12);
    EXPECT_NEAR(stress.e23, 4.3904235877576825e-06,1e-12);
    EXPECT_NEAR(stress.e31, 5.1291622417145914e-06,1e-12);
    EXPECT_NEAR(stress.e32, 4.3904235877576833e-06,1e-12);
    EXPECT_NEAR(stress.e33, -3.4278442125590892e-05,1e-12);
}

#undef private


int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    testing::InitGoogleTest(&argc, argv);
    int result = RUN_ALL_TESTS();
    MPI_Finalize();
    return 0;
}