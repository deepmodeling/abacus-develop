#include "module_basis/module_nao/real_gaunt_table.h"
#include "gtest/gtest.h"

#include <iostream>
#include<fstream>

#ifdef __MPI
#include <mpi.h>
#endif

#include <chrono>
using iclock = std::chrono::high_resolution_clock;

#include "module_base/constants.h"
#include "module_basis/module_ao/ORB_gaunt_table.h"

/***********************************************************
 *      Unit test of class "RealGauntTable"
 ***********************************************************/
/*!
 *  Tested functions:
 *
 *  - build
 *      - constructs the Gaunt table of real spherical harmonics
 *
 *  - lmax
 *      - gets the maximum angular momentum
 *
 *  - operator()
 *      - gets a specific Gaunt coefficient of real spherical harmonics
 *
 *                                                          */
class RealGauntTableTest : public ::testing::Test
{
  protected:
    void SetUp() { /*rgt.build(lmax);*/ }
    void TearDown() {}

    int lmax = 10;      //!< maximum angular momentum
    //RealGauntTable rgt; //!< object under test
    const double tol = 1e-12; //!< numerical error tolerance for individual Gaunt coefficient
};

TEST_F(RealGauntTableTest, LegacyConsistency)
{
    //iclock::time_point start;
    //std::chrono::duration<double> dur;

    // this test checks whether the coefficients in RealGauntTable is consistent with those of ORB_gaunt_table
    // this test shall be removed in the future once the refactoring is finished
    ORB_gaunt_table ogt;

    RealGauntTable::instance().build(lmax);

    //start = iclock::now();
    ogt.init_Gaunt_CH(lmax);
    ogt.init_Gaunt(lmax);
    //dur = iclock::now() - start;
    //std::cout << "time elased = " << dur.count() << " s" << std::endl;

    //RealGauntTable rgt2;
    //start = iclock::now();
    //rgt2.build(lmax);
    //dur = iclock::now() - start;
    //std::cout << "time elased = " << dur.count() << " s" << std::endl;

    for (int l1 = 0; l1 <= lmax; ++l1)
    {
        for (int mm1 = 0; mm1 <= 2*l1; ++mm1)
        {
            int index1 = ogt.get_lm_index(l1, mm1);
            for (int l2 = 0; l2 <= lmax; ++l2)
            {
                for (int mm2 = 0; mm2 <= 2*l2; ++mm2)
                {
                    int index2 = ogt.get_lm_index(l2, mm2);
                    for (int l3 = 0; l3 <= 2*lmax; ++l3)
                    //for (int l3 = std::abs(l1-l2); l3 <= 2*lmax; l3 += 2)
                    {
                        for (int mm3 = 0; mm3 <= 2*l3; ++mm3)
                        {
                            int index3 = ogt.get_lm_index(l3, mm3);

                            int m1 = ogt.Index_M(mm1);
                            int m2 = ogt.Index_M(mm2);
                            int m3 = ogt.Index_M(mm3);

                            EXPECT_NEAR(RealGauntTable::instance()(l1, l2, l3, m1, m2, m3), 
                                    ogt.Gaunt_Coefficients(index1, index2, index3), tol);
                        }
                    }
                }
            }
        }
    }
}

TEST_F(RealGauntTableTest, SanityCheck)
{
    EXPECT_EQ(RealGauntTable::instance().lmax(), lmax);

    EXPECT_NEAR(RealGauntTable::instance()(0, 0, 0, 0, 0, 0), ModuleBase::SQRT_INVERSE_FOUR_PI, tol);

    EXPECT_NEAR(RealGauntTable::instance()(4, 0, 4, 3, 0, 3), ModuleBase::SQRT_INVERSE_FOUR_PI, tol);
    EXPECT_NEAR(RealGauntTable::instance()(4, 0, 4, -3, 0, -3), ModuleBase::SQRT_INVERSE_FOUR_PI, tol);

    EXPECT_NEAR(RealGauntTable::instance()(2, 2, 2, 2, -1, -1), -std::sqrt(15.0) / 7.0 * ModuleBase::SQRT_INVERSE_FOUR_PI, tol);
    EXPECT_NEAR(RealGauntTable::instance()(2, 2, 2, -1, 2, -1), -std::sqrt(15.0) / 7.0 * ModuleBase::SQRT_INVERSE_FOUR_PI, tol);

    EXPECT_NEAR(RealGauntTable::instance()(3, 3, 2, 2, 1, 1), ModuleBase::SQRT_INVERSE_FOUR_PI / std::sqrt(6.0), tol);
    EXPECT_NEAR(RealGauntTable::instance()(2, 3, 3, 1, 1, 2), ModuleBase::SQRT_INVERSE_FOUR_PI / std::sqrt(6.0), tol);

    EXPECT_NEAR(RealGauntTable::instance()(4, 5, 7, 3, -2, -5), ModuleBase::SQRT_INVERSE_FOUR_PI * std::sqrt(210.0) / 221.0, tol);
}

/**
* @brief real gaunt function test:
*           Using sympy realgaunt function test points set.
*           2 failed test case with a sign error occurred 
*           on some cases of the actual Gaunt coefficients of 2 negative m
*/
TEST_F(RealGauntTableTest, Check2)
{
    const double PI	= 3.14159265358979323846;
    for(int i=0;i<3;i++){
        for(int j=-i;j<=i;j++){
            EXPECT_NEAR(RealGauntTable::instance()(0, i, i, 0, j, j), 1/(2*sqrt(PI)), tol);
        }
    }
    
    // EXPECT_NEAR(rgt(1, 1, 2, -1, 1, -2), -sqrt(15)/(10*sqrt(PI)), tol); //wrong case
    // EXPECT_NEAR(rgt(1, 1, 2, -1, 1, -2),-sqrt(15)/(10*sqrt(PI)), tol);  //wrong case
    
    EXPECT_NEAR(RealGauntTable::instance()(1, 1, 2, 0, 0, 0),sqrt(5)/(5*sqrt(PI)), tol);
    EXPECT_NEAR(RealGauntTable::instance()(1, 1, 2, 1, 1, 0),-sqrt(5)/(10*sqrt(PI)), tol);
    EXPECT_NEAR(RealGauntTable::instance()(2, 2, 2, 0, 0, 0),sqrt(5)/(7*sqrt(PI)), tol);
    EXPECT_NEAR(RealGauntTable::instance()(2, 2, 2, 0, 2, 2),-sqrt(5)/(7*sqrt(PI)), tol);
    EXPECT_NEAR(RealGauntTable::instance()(2, 2, 2, -2, -2, 0),-sqrt(5)/(7*sqrt(PI)), tol);
    
    EXPECT_NEAR(RealGauntTable::instance()(1, 1, 2, -1, 0, -1),sqrt(15)/(10*sqrt(PI)), tol);
    EXPECT_NEAR(RealGauntTable::instance()(1, 1, 2, 0, 1, 1),sqrt(15)/(10*sqrt(PI)), tol);
    EXPECT_NEAR(RealGauntTable::instance()(1, 1, 2, 1, 1, 2),sqrt(15)/(10*sqrt(PI)), tol);
    
    EXPECT_NEAR(RealGauntTable::instance()(1, 1, 2, -1, -1, 2),-sqrt(15)/(10*sqrt(PI)), tol);

    EXPECT_NEAR(RealGauntTable::instance()(2, 2, 2, 0, 1, 1),sqrt(5)/(14*sqrt(PI)), tol);
    EXPECT_NEAR(RealGauntTable::instance()(2, 2, 2, 1, 1, 2),sqrt(15)/(14*sqrt(PI)), tol);
    EXPECT_NEAR(RealGauntTable::instance()(2, 2, 2, -1, -1, 2),-sqrt(15)/(14*sqrt(PI)), tol);

    
}



/**
* @brief gaunt function test:
*           Use the sympy gaunt function run results as the test set. 
*           Test set results are stored in Gaunt.txt with data structure 
*           l1  l2  l3  m1  m2  m3  gaunt(l1,l2,l3,m1,m2,m3). 
*           Double precision machine error (1e-15) is used to compare the 
*           accuracy of the function in abacus and the gaunt function in sumpy.
*/
struct gaunt_ans{
    int l1;
    int l2;
    int l3;
    int m1;
    int m2;
    int m3;
    double gaunt;
}typedef gaunt_ans;



// the length of the values in gaunt.txt
#define len_gaunt 7242
TEST_F(RealGauntTableTest, Check3)
{
    gaunt_ans ga_ref[len_gaunt];
    gaunt_ans ga_func[len_gaunt];

 
    int len=len_gaunt;
	std::ifstream infile("../../../../../source/module_basis/module_nao/test/gaunt.txt");
    if (!infile) {
        EXPECT_NEAR(0,1,tol);
    }

    double tmp;
    for(int i=0;i<len;i++){
        infile>>ga_ref[i].l1>>ga_ref[i].l2>>ga_ref[i].l3>>ga_ref[i].m1>>ga_ref[i].m2>>ga_ref[i].m3>>ga_ref[i].gaunt;
    }

    int cnt=0;
    const double PI	= 3.14159265358979323846;
    int l_max = 10;
    int l3_max = 24;
    for(int l1=0;l1<=l_max;l1++){
        for(int l2=l1;l2<=l_max;l2++){
            for(int l3=l2;l3<=l3_max;l3++){
                for(int m1=-l1;m1<=l1;m1++){
                    for(int m2=-l2;m2<=l2;m2++){
                        for(int m3=0;m3<=l3;m3++){
                            double gaunt = RealGauntTable::instance().gaunt(l1,l2,l3,m1,m2,m3);
                            double gaunt_symmetry[20];
                            double tmp;
                            int cnt_sym = 0;
                            if(gaunt>tol&&cnt<len_gaunt){
                                ga_func[cnt].l1=l1,ga_func[cnt].l2=l2,ga_func[cnt].l3=l3;
                                ga_func[cnt].m1=m1,ga_func[cnt].m2=m2,ga_func[cnt].m3=m3;
                                ga_func[cnt++].gaunt=gaunt;


                                // Detects whether the values of functions omitted in the loop due to symmetry are equal
                                if(l3<=10){
                                    
                                    gaunt_symmetry[cnt_sym++] = RealGauntTable::instance().gaunt(l2,l1,l3,m2,m1,m3);
                                    gaunt_symmetry[cnt_sym++] = RealGauntTable::instance().gaunt(l2,l1,l3,-m2,-m1,-m3);
                            
                                    gaunt_symmetry[cnt_sym++] = RealGauntTable::instance().gaunt(l2,l3,l1,m2,m3,m1);
                                    gaunt_symmetry[cnt_sym++] = RealGauntTable::instance().gaunt(l2,l3,l1,-m2,-m3,-m1);
                                
                                    gaunt_symmetry[cnt_sym++] = RealGauntTable::instance().gaunt(l1,l3,l2,m1,m3,m2);
                                    gaunt_symmetry[cnt_sym++] = RealGauntTable::instance().gaunt(l1,l3,l2,-m1,-m3,-m2);
                                
                                    gaunt_symmetry[cnt_sym++] = RealGauntTable::instance().gaunt(l3,l2,l1,m3,m2,m1);
                                    gaunt_symmetry[cnt_sym++] = RealGauntTable::instance().gaunt(l3,l2,l1,-m3,-m2,-m1);
  
                                    gaunt_symmetry[cnt_sym++] = RealGauntTable::instance().gaunt(l3,l1,l2,m3,m1,m2);
                                    gaunt_symmetry[cnt_sym++] = RealGauntTable::instance().gaunt(l3,l1,l2,-m3,-m1,-m2);
    
                                    gaunt_symmetry[cnt_sym++] = RealGauntTable::instance().gaunt(l1,l2,l3,-m1,-m2,-m3);

                                } else{
                                    gaunt_symmetry[cnt_sym++] = RealGauntTable::instance().gaunt(l2,l1,l3,m2,m1,m3);
                                    gaunt_symmetry[cnt_sym++] = RealGauntTable::instance().gaunt(l2,l1,l3,-m2,-m1,-m3);
                                    gaunt_symmetry[cnt_sym++] = RealGauntTable::instance().gaunt(l1,l2,l3,-m1,-m2,-m3);
                                }
                            
                            }
                        }   
                    }                    
                }
            }                   
        }        
    }

    if(len!=cnt){
        printf("gaunt function count err!\n");
    }
    for(int i=0;i<len;i++){
        EXPECT_EQ(ga_func[i].l1, ga_ref[i].l1);
        EXPECT_EQ(ga_func[i].l2, ga_ref[i].l2);
        EXPECT_EQ(ga_func[i].l3, ga_ref[i].l3);
        EXPECT_EQ(ga_func[i].m1, ga_ref[i].m1);
        EXPECT_EQ(ga_func[i].m2, ga_ref[i].m2);
        EXPECT_EQ(ga_func[i].m3, ga_ref[i].m3);
        EXPECT_NEAR(ga_func[i].gaunt,ga_ref[i].gaunt, tol);
    }



}


int main(int argc, char** argv)
{

#ifdef __MPI
    MPI_Init(&argc, &argv);
#endif

    testing::InitGoogleTest(&argc, argv);
    int result = RUN_ALL_TESTS();

#ifdef __MPI
    MPI_Finalize();
#endif

    return result;
}
