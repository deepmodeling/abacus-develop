#include "gtest/gtest.h"
#include "../relax.h"
#include "module_cell/unitcell.h"
#include "module_io/input.h"
#include "relax_test.h"
#include <fstream>

class Test_SETGRAD : public testing::Test
{
    protected:
        Relax rl;
        std::vector<double> result;

        void SetUp()
        {
            GlobalV::FORCE_THR = 0.001;
            GlobalV::CALCULATION = "cell-relax";

            ModuleBase::matrix force_in, stress_in;
            int nat = 3;

            force_in.create(nat,3);
            stress_in.create(3,3);

            force_in(0,0) = 1; force_in(0,1) = 2; force_in(0,2)= 3;
            force_in(1,0) = 4; force_in(1,1) = 5; force_in(1,2)= 6;
            force_in(2,0) = 7; force_in(2,1) = 8; force_in(2,2)= 9;

            stress_in(0,0) = 1; stress_in(0,1) = 2; stress_in(0,2)= 3;
            stress_in(1,0) = 4; stress_in(1,1) = 5; stress_in(1,2)= 6;
            stress_in(2,0) = 7; stress_in(2,1) = 8; stress_in(2,2)= 9;

            GlobalC::ucell.ntype = 1;
            GlobalC::ucell.nat = nat;
            GlobalC::ucell.atoms = new Atom[1];
            GlobalC::ucell.atoms[0].na = nat;
            GlobalC::ucell.omega = 1.0;
            GlobalC::ucell.lat0 = 1.0;
            
            GlobalC::ucell.iat2it = new int[nat];
            GlobalC::ucell.iat2ia = new int[nat];
            GlobalC::ucell.atoms[0].mbl = new ModuleBase::Vector3<int>[nat];
            GlobalC::ucell.atoms[0].taud = new ModuleBase::Vector3<double>[nat];
            GlobalC::ucell.lc = new int[3];

            GlobalC::ucell.iat2it[0] = 0;
            GlobalC::ucell.iat2it[1] = 0;
            GlobalC::ucell.iat2it[2] = 0;

            GlobalC::ucell.iat2ia[0] = 0;
            GlobalC::ucell.iat2ia[1] = 1;
            GlobalC::ucell.iat2ia[2] = 2;

            GlobalC::ucell.atoms[0].mbl[0].x = 0;
            GlobalC::ucell.atoms[0].mbl[0].y = 0;
            GlobalC::ucell.atoms[0].mbl[0].z = 1;

            GlobalC::ucell.atoms[0].mbl[1].x = 0;
            GlobalC::ucell.atoms[0].mbl[1].y = 1;
            GlobalC::ucell.atoms[0].mbl[1].z = 0;

            GlobalC::ucell.atoms[0].mbl[2].x = 1;
            GlobalC::ucell.atoms[0].mbl[2].y = 0;
            GlobalC::ucell.atoms[0].mbl[2].z = 0;

            GlobalC::ucell.atoms[0].taud[0] = 0.0;
            GlobalC::ucell.atoms[0].taud[1] = 0.0;
            GlobalC::ucell.atoms[0].taud[2] = 0.0;

            GlobalC::ucell.lc[0] = 1;
            GlobalC::ucell.lc[1] = 1;
            GlobalC::ucell.lc[2] = 1;

            rl.init_relax(nat);
            rl.relax_step(force_in,stress_in,0.0);

            for(int i=0;i<3;i++)
            {
                result.push_back(GlobalC::ucell.atoms[0].taud[i].x);
                result.push_back(GlobalC::ucell.atoms[0].taud[i].y);
                result.push_back(GlobalC::ucell.atoms[0].taud[i].z);
            }
            push_result();

            //reset lattice vector
            GlobalC::ucell.latvec.Identity();
            INPUT.fixed_axes = "shape";
            rl.init_relax(nat);
            rl.relax_step(force_in,stress_in,0.0);
            push_result();

            //reset lattice vector
            GlobalC::ucell.latvec.Identity();
            INPUT.fixed_axes = "volume";
            rl.init_relax(nat);
            rl.relax_step(force_in,stress_in,0.0);
            push_result();

            //reset lattice vector
            GlobalC::ucell.latvec.Identity();
            INPUT.fixed_axes = "a"; //anything other than "None"
            INPUT.fixed_ibrav = 1;
            GlobalC::ucell.lc[0] = 0;
            GlobalC::ucell.lc[1] = 0;
            GlobalC::ucell.lc[2] = 0;
            rl.init_relax(nat);
            rl.relax_step(force_in,stress_in,0.0);
            push_result();
        }

        void push_result()
        {
            result.push_back(GlobalC::ucell.latvec.e11);
            result.push_back(GlobalC::ucell.latvec.e12);
            result.push_back(GlobalC::ucell.latvec.e13);
            result.push_back(GlobalC::ucell.latvec.e21);
            result.push_back(GlobalC::ucell.latvec.e22);
            result.push_back(GlobalC::ucell.latvec.e23);
            result.push_back(GlobalC::ucell.latvec.e31);
            result.push_back(GlobalC::ucell.latvec.e32);
            result.push_back(GlobalC::ucell.latvec.e33);             
        }

};

TEST_F(Test_SETGRAD, relax_new)
{
    std::vector<double> result_ref = 
    {
        0,0,0.1709672056,0,0.2849453427,0,0.3989234797,0,0,1.005319517,
        0.01063903455,0.01595855183,0.0212780691,1.026597586,0.03191710366,
        0.03723662093,0.04255613821,1.047875655,1.059181731,0,0,0,1.059181731,
        0,0,0,1.059181731,1.034363264,0.01301504537,0.01952256806,
        0.02603009074,1.060393355,0.03904513611,0.0455526588,0.05206018148,
        1.086423445,1,0,0,0,1,0,0,0,1
    };

    for(int i=0;i<result.size();i++)
    {
        EXPECT_NEAR(result[i],result_ref[i],1e-8);
    }
}

class Test_RELAX : public testing::Test
{
    protected:
        Relax rl;
        std::vector<double> result;

        void SetUp()
        {
            std::ifstream force_file("./support/force.txt");
            std::ifstream stress_file("./support/stress.txt");
            std::ifstream energy_file("./support/stress.txt");

            int nstep = 66;
            int nat = 5;
            double energy;
            GlobalV::STRESS_THR = 0.01;

            this->setup_cell();

            ModuleBase::matrix force_in, stress_in;

            force_in.create(nat,3);
            stress_in.create(3,3);

            rl.init_relax(nat);

            for(int istep=0;istep<nstep;istep++)
            {
                for(int i=0;i<nat;i++)
                {
                    for(int j=0;j<3;j++)
                    {
                        force_file >> force_in(i,j);
                    }
                }
                for(int i=0;i<3;i++)
                {
                    for(int j=0;j<3;j++)
                    {
                        stress_file >> stress_in(i,j);
                    }
                }

                energy_file >> energy;

                rl.relax_step(force_in,stress_in,energy);

                result.push_back(GlobalC::ucell.atoms[0].taud[0].x);
                result.push_back(GlobalC::ucell.atoms[0].taud[0].y);
                result.push_back(GlobalC::ucell.atoms[0].taud[0].z);
                result.push_back(GlobalC::ucell.atoms[1].taud[0].x);
                result.push_back(GlobalC::ucell.atoms[1].taud[0].y);
                result.push_back(GlobalC::ucell.atoms[1].taud[0].z);
                result.push_back(GlobalC::ucell.atoms[2].taud[0].x);
                result.push_back(GlobalC::ucell.atoms[2].taud[0].y);
                result.push_back(GlobalC::ucell.atoms[2].taud[0].z);
                result.push_back(GlobalC::ucell.atoms[2].taud[1].x);
                result.push_back(GlobalC::ucell.atoms[2].taud[1].y);
                result.push_back(GlobalC::ucell.atoms[2].taud[1].z);
                result.push_back(GlobalC::ucell.atoms[2].taud[2].x);
                result.push_back(GlobalC::ucell.atoms[2].taud[2].y);
                result.push_back(GlobalC::ucell.atoms[2].taud[2].z);
                result.push_back(GlobalC::ucell.latvec.e11);
                result.push_back(GlobalC::ucell.latvec.e12);
                result.push_back(GlobalC::ucell.latvec.e13);
                result.push_back(GlobalC::ucell.latvec.e21);
                result.push_back(GlobalC::ucell.latvec.e22);
                result.push_back(GlobalC::ucell.latvec.e23);
                result.push_back(GlobalC::ucell.latvec.e31);
                result.push_back(GlobalC::ucell.latvec.e32);
                result.push_back(GlobalC::ucell.latvec.e33);
            }
        }

        void setup_cell()
        {
            int ntype = 3, nat = 5;
            GlobalC::ucell.ntype = ntype;
            GlobalC::ucell.nat = nat;

            GlobalC::ucell.omega = 452.590903143121;
            GlobalC::ucell.lat0 = 1.8897259886;
            GlobalC::ucell.iat2it = new int[nat];
            GlobalC::ucell.iat2ia = new int[nat];
            GlobalC::ucell.iat2it[0] = 0;
            GlobalC::ucell.iat2it[1] = 1;
            GlobalC::ucell.iat2it[2] = 2;
            GlobalC::ucell.iat2it[3] = 2;
            GlobalC::ucell.iat2it[4] = 2;

            GlobalC::ucell.iat2ia[0] = 0;
            GlobalC::ucell.iat2ia[1] = 0;
            GlobalC::ucell.iat2ia[2] = 0;
            GlobalC::ucell.iat2ia[3] = 1;
            GlobalC::ucell.iat2ia[4] = 2;

            GlobalC::ucell.atoms = new Atom[ntype];
            GlobalC::ucell.atoms[0].na = 1;
            GlobalC::ucell.atoms[1].na = 1;
            GlobalC::ucell.atoms[2].na = 3;
            
            for(int i=0;i<ntype;i++)
            {
                int na = GlobalC::ucell.atoms[i].na;
                GlobalC::ucell.atoms[i].mbl = new ModuleBase::Vector3<int>[na];
                GlobalC::ucell.atoms[i].taud = new ModuleBase::Vector3<double>[na];
                for (int j=0;j<na;j++)
                {
                    GlobalC::ucell.atoms[i].mbl[j] = {1,1,1};
                }
            }
            GlobalC::ucell.atoms[0].taud[0] = {0.5,0.5,0.00413599999956205};
            GlobalC::ucell.atoms[1].taud[0] = {0  ,0  ,0.524312999999893  };
            GlobalC::ucell.atoms[2].taud[0] = {0  ,0.5,0.479348999999274  };
            GlobalC::ucell.atoms[2].taud[1] = {0.5,0  ,0.479348999999274  };
            GlobalC::ucell.atoms[2].taud[2] = {0  ,0  ,0.958854000000429  };
            
            GlobalC::ucell.lc = new int[3];
            GlobalC::ucell.lc[0] = 1;
            GlobalC::ucell.lc[1] = 1;
            GlobalC::ucell.lc[2] = 1;

            GlobalC::ucell.latvec.e11 = 3.96;
            GlobalC::ucell.latvec.e12 = 0;
            GlobalC::ucell.latvec.e13 = 0;
            GlobalC::ucell.latvec.e21 = 0;
            GlobalC::ucell.latvec.e22 = 3.96;
            GlobalC::ucell.latvec.e23 = 0;
            GlobalC::ucell.latvec.e31 = 0;
            GlobalC::ucell.latvec.e32 = 0;
            GlobalC::ucell.latvec.e33 = 4.2768;
        }
};

TEST_F(Test_RELAX, relax_new)
{
    int size = 1584;
    double tmp;
    std::ifstream result_ref("./support/result_ref.txt");

    for(int i=0;i<size;i++)
    {
        result_ref >> tmp;
        EXPECT_NEAR(tmp,result[i],1e-8);
    }
}