#include "gtest/gtest.h"
#include <iomanip>
#define private public
#include "module_parameter/parameter.h"
#undef private
#include "../relax.h"
#include "module_cell/unitcell.h"
#include "relax_test.h"
#include <fstream>


class Test_SETGRAD : public testing::Test
{
    protected:
        Relax rl;
        std::vector<double> result;
        Input_para& input = PARAM.input;
        UnitCell ucell;

        void SetUp()
        {
            PARAM.input.force_thr  = 0.001;
            PARAM.input.calculation = "cell-relax";

            ModuleBase::matrix force_in, stress_in;
            int nat = 3;

            force_in.create(nat,3);
            stress_in.create(3,3);

            force_in(0,0) = 0.1; force_in(0,1) = 0.1; force_in(0,2)= 0.1;
            force_in(1,0) = 0; force_in(1,1) = 0.1; force_in(1,2)= 0.1;
            force_in(2,0) = 0; force_in(2,1) = 0; force_in(2,2)= 0.1;

            stress_in(0,0) = 1; stress_in(0,1) = 1; stress_in(0,2)= 1;
            stress_in(1,0) = 0; stress_in(1,1) = 1; stress_in(1,2)= 1;
            stress_in(2,0) = 0; stress_in(2,1) = 0; stress_in(2,2)= 1;

            ucell.ntype = 1;
            ucell.nat = nat;
            ucell.atoms = new Atom[1];
            ucell.atoms[0].na = nat;
            ucell.omega = 1.0;
            ucell.lat0 = 1.0;
            
            ucell.iat2it = new int[nat];
            ucell.iat2ia = new int[nat];
            ucell.atoms[0].mbl.resize(nat);
            ucell.atoms[0].taud.resize(nat);
            ucell.atoms[0].tau.resize(nat);
            ucell.atoms[0].dis.resize(nat);
            ucell.lc = new int[3];

            ucell.iat2it[0] = 0;
            ucell.iat2it[1] = 0;
            ucell.iat2it[2] = 0;

            ucell.iat2ia[0] = 0;
            ucell.iat2ia[1] = 1;
            ucell.iat2ia[2] = 2;

            ucell.atoms[0].mbl[0].x = 0;
            ucell.atoms[0].mbl[0].y = 0;
            ucell.atoms[0].mbl[0].z = 1;

            ucell.atoms[0].mbl[1].x = 0;
            ucell.atoms[0].mbl[1].y = 1;
            ucell.atoms[0].mbl[1].z = 0;

            ucell.atoms[0].mbl[2].x = 1;
            ucell.atoms[0].mbl[2].y = 0;
            ucell.atoms[0].mbl[2].z = 0;

            ucell.atoms[0].taud[0] = 0.0;
            ucell.atoms[0].taud[1] = 0.0;
            ucell.atoms[0].taud[2] = 0.0;

            ucell.atoms[0].tau.resize(nat);

            ucell.lc[0] = 1;
            ucell.lc[1] = 1;
            ucell.lc[2] = 1;

            rl.init_relax(nat);
            rl.relax_step(ucell,force_in,stress_in,0.0);

            for(int i=0;i<3;i++)
            {
                result.push_back(ucell.atoms[0].taud[i].x);
                result.push_back(ucell.atoms[0].taud[i].y);
                result.push_back(ucell.atoms[0].taud[i].z);
            }
            push_result();

            //reset lattice vector
            ucell.latvec.Identity();
            input.fixed_axes = "shape";
            rl.init_relax(nat);
            rl.relax_step(ucell,force_in,stress_in,0.0);
            push_result();

            //reset lattice vector
            ucell.latvec.Identity();
            input.fixed_axes = "volume";
            rl.init_relax(nat);
            rl.relax_step(ucell,force_in,stress_in,0.0);
            push_result();

            //reset lattice vector
            ucell.latvec.Identity();
            input.fixed_axes = "a"; //anything other than "None"
            input.fixed_ibrav = true;
            ucell.latName = "sc";
            ucell.lc[0] = 0;
            ucell.lc[1] = 0;
            ucell.lc[2] = 0;
            rl.init_relax(nat);
            rl.relax_step(ucell,force_in,stress_in,0.0);
            
            push_result();
        }

        void push_result()
        {
            result.push_back(ucell.latvec.e11);
            result.push_back(ucell.latvec.e12);
            result.push_back(ucell.latvec.e13);
            result.push_back(ucell.latvec.e21);
            result.push_back(ucell.latvec.e22);
            result.push_back(ucell.latvec.e23);
            result.push_back(ucell.latvec.e31);
            result.push_back(ucell.latvec.e32);
            result.push_back(ucell.latvec.e33);             
        }

};

TEST_F(Test_SETGRAD, relax_new)
{
    std::vector<double> result_ref = 
    {
        0,           0,            0.24293434145, 
        0,           0.242934341453,           0,
        0,           0,           0,
        //paramter for taud
        1.2267616333,0.2267616333,0.22676163333, 
        0,            1.2267616333 ,0.2267616333,
        0,            0,           1.22676163333,
        // paramter for fisrt time after relaxation
        1.3677603495, 0,            0,
        0,            1.36776034956, 0,
        0,            0,            1.36776034956,
        // paramter for second time after relaxation
        1.3677603495  ,0.3633367476,0.36333674766,
        0,            1.3677603495 ,0.36333674766,
        0,            0,            1.3677603495 ,
        // paramter for third time after relaxation
        1,0,0,0,1,0,0,0,1
        // paramter for fourth time after relaxation
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
        UnitCell ucell;

        void SetUp()
        {
            std::ifstream force_file("./support/force.txt");
            std::ifstream stress_file("./support/stress.txt");
            std::ifstream energy_file("./support/stress.txt");

            int nstep = 3;
            int nat = 5;
            double energy;
            PARAM.input.stress_thr = 0.01;

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

                PARAM.input.fixed_ibrav = false;
                rl.relax_step(ucell,force_in,stress_in,energy);
                result.push_back(ucell.atoms[0].taud[0].x);
                result.push_back(ucell.atoms[0].taud[0].y);
                result.push_back(ucell.atoms[0].taud[0].z);
                result.push_back(ucell.atoms[1].taud[0].x);
                result.push_back(ucell.atoms[1].taud[0].y);
                result.push_back(ucell.atoms[1].taud[0].z);
                result.push_back(ucell.atoms[2].taud[0].x);
                result.push_back(ucell.atoms[2].taud[0].y);
                result.push_back(ucell.atoms[2].taud[0].z);
                result.push_back(ucell.atoms[2].taud[1].x);
                result.push_back(ucell.atoms[2].taud[1].y);
                result.push_back(ucell.atoms[2].taud[1].z);
                result.push_back(ucell.atoms[2].taud[2].x);
                result.push_back(ucell.atoms[2].taud[2].y);
                result.push_back(ucell.atoms[2].taud[2].z);
                result.push_back(ucell.latvec.e11);
                result.push_back(ucell.latvec.e12);
                result.push_back(ucell.latvec.e13);
                result.push_back(ucell.latvec.e21);
                result.push_back(ucell.latvec.e22);
                result.push_back(ucell.latvec.e23);
                result.push_back(ucell.latvec.e31);
                result.push_back(ucell.latvec.e32);
                result.push_back(ucell.latvec.e33);
            }
        }

        void setup_cell()
        {
            int ntype = 3, nat = 5;
            ucell.ntype = ntype;
            ucell.nat = nat;

            ucell.omega = 452.590903143121;
            ucell.lat0 = 1.8897259886;
            ucell.iat2it = new int[nat];
            ucell.iat2ia = new int[nat];
            ucell.iat2it[0] = 0;
            ucell.iat2it[1] = 1;
            ucell.iat2it[2] = 2;
            ucell.iat2it[3] = 2;
            ucell.iat2it[4] = 2;

            ucell.iat2ia[0] = 0;
            ucell.iat2ia[1] = 0;
            ucell.iat2ia[2] = 0;
            ucell.iat2ia[3] = 1;
            ucell.iat2ia[4] = 2;

            ucell.atoms = new Atom[ntype];
            ucell.atoms[0].na = 1;
            ucell.atoms[1].na = 1;
            ucell.atoms[2].na = 3;
            
            for(int i=0;i<ntype;i++)
            {
                int na = ucell.atoms[i].na;
                ucell.atoms[i].mbl.resize(na);
                ucell.atoms[i].taud.resize(na);
                ucell.atoms[i].tau.resize(na);
                ucell.atoms[i].dis.resize(na);
                for (int j=0;j<na;j++)
                {
                    ucell.atoms[i].mbl[j] = {1,1,1};
                }
            }
            ucell.atoms[0].taud[0] = {0.5,0.5,0.00413599999956205};
            ucell.atoms[1].taud[0] = {0  ,0  ,0.524312999999893  };
            ucell.atoms[2].taud[0] = {0  ,0.5,0.479348999999274  };
            ucell.atoms[2].taud[1] = {0.5,0  ,0.479348999999274  };
            ucell.atoms[2].taud[2] = {0  ,0  ,0.958854000000429  };
            
            ucell.lc = new int[3];
            ucell.lc[0] = 1;
            ucell.lc[1] = 1;
            ucell.lc[2] = 1;

            ucell.latvec.e11 = 3.96;
            ucell.latvec.e12 = 0;
            ucell.latvec.e13 = 0;
            ucell.latvec.e21 = 0;
            ucell.latvec.e22 = 3.96;
            ucell.latvec.e23 = 0;
            ucell.latvec.e31 = 0;
            ucell.latvec.e32 = 0;
            ucell.latvec.e33 = 4.2768;
        }
};

TEST_F(Test_RELAX, relax_new)
{
    int size = 72;
    double tmp;
    std::vector<double> result_ref=
    {
        0.5000000586,0.4999998876,0.009364595811,
        0.9999999281,1.901333279e-07,0.5476035454,
        3.782097706e-07,0.4999999285,0.4770375874,
        0.4999997826,2.072311863e-07,0.477037871,
        0.9999998523,0.9999997866,0.9349574003,
        // paramter for taud after first relaxation
        4.006349654,-1.93128788e-07,5.793863639e-07,
        -1.93128788e-07,4.006354579,-3.86257576e-07,
        6.757962549e-07,-4.505308366e-07,3.966870038,
        // paramter for latvec after first relaxation
        0.5000000566,0.4999998916,0.009177239183,
        0.9999999308,1.832935626e-07,0.5467689737,
        3.647769323e-07,0.4999999311,0.4771204124,
        0.4999997903,1.998879545e-07,0.4771206859,
        0.9999998574,0.9999997943,0.9358136888,
        // paramter for taud after second relaxation
        3.999761277,-1.656764727e-07,4.97029418e-07,
        -1.656764727e-07,3.999765501,-3.313529453e-07,
        5.797351131e-07,-3.864900754e-07,4.010925071,
        // paramter for latvec after second relaxation
        0.500000082,0.4999999574,0.01057784352,
        0.9999999149,1.939640249e-07,0.5455830599,
        3.795967781e-07,0.4999998795,0.4765373919,
        0.4999998037,2.756298268e-07,0.4765374602,
        0.9999998196,0.9999996936,0.9367652445,
        // paramter for taud after third relaxation
        4.017733155,-1.420363309e-07,2.637046077e-07,
        -1.420364243e-07,4.017735987,3.126225134e-07,
        3.479123171e-07,2.578467568e-07,4.011674933
        // paramter for latvec after third relaxation
    };
    for(int i=0;i<size;i++)
    {
        EXPECT_NEAR(result_ref[i],result[i],1e-8);
    }
}