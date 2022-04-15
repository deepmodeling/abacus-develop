#include"../math_integral.h"
#include"gtest/gtest.h"

#include<math.h>

#define doublethreshold 1e-12



/************************************************
*  unit test of class Integral
***********************************************/

/**
 * Tested functions:
 *  - function Simpson_Integral with rab as input
 *  - function Simpson_Integral with dr as input
 *  - function Simpson_Integral_0toall
 *  - function Simpson_Integral_alltoinf
 * 
 */

class SimpsonIntegralSinx : public testing::Test
{
    /**
     * test the integral of sinx between [0,PI],
     * devide to mesh-1 parts
     */
    
    protected:

    double*     func;
    double*     rab;
    int         mesh = 10001;
    double      dr   = M_PI/(mesh-1);
    double      asum;
    double*     asumlist;
    double      expectvalue = 2.0;

    void SetUp()
    {
        func        = new double[mesh];
        rab         = new double[mesh];
        asumlist    = new double[mesh];

        for (int i=0;i<mesh;i++)
        {
            func[i] = sin(M_PI*i/(mesh-1));
            rab[i] = dr;
        }
    }

    void TearDown()
    {
        delete func;
        delete rab;
        delete asumlist;
    }
};



TEST_F(SimpsonIntegralSinx,SimpsonIntegralRab)
{
    ModuleBase::Integral::Simpson_Integral(mesh,func,rab,asum);
    EXPECT_NEAR(asum,expectvalue,doublethreshold);
    EXPECT_DEATH(ModuleBase::Integral::Simpson_Integral(mesh-1,func,rab,asum),"");
}

TEST_F(SimpsonIntegralSinx,SimpsonIntegralDr)
{
    ModuleBase::Integral::Simpson_Integral(mesh,func,dr,asum);
    EXPECT_NEAR(asum,expectvalue,doublethreshold);
    EXPECT_DEATH(ModuleBase::Integral::Simpson_Integral(mesh-1,func,dr,asum),"");
}

TEST_F(SimpsonIntegralSinx,SimpsonIntegral0toall)
{
    int halfmesh = round(mesh/2);
    ModuleBase::Integral::Simpson_Integral_0toall(mesh,func,rab,asumlist);
    EXPECT_NEAR(asumlist[mesh-1],expectvalue,doublethreshold);
    EXPECT_NEAR(asumlist[mesh-2],1.9999999506519738901,doublethreshold);
    EXPECT_NEAR(asumlist[halfmesh],expectvalue/2.0,doublethreshold);
    EXPECT_NEAR(asumlist[halfmesh-1],0.99968584073722277505,doublethreshold);
}


TEST_F(SimpsonIntegralSinx,SimpsonIntegralAlltoinf)
{
    int halfmesh = round(mesh/2);
    ModuleBase::Integral::Simpson_Integral_alltoinf(mesh,func,rab,asumlist);
    EXPECT_NEAR(asumlist[0],expectvalue,doublethreshold);
    EXPECT_NEAR(asumlist[1],1.9999999506519754444,doublethreshold);
    EXPECT_NEAR(asumlist[halfmesh],expectvalue/2.0,doublethreshold);
    EXPECT_NEAR(asumlist[halfmesh-1],1.0003141592627740053,doublethreshold);
}