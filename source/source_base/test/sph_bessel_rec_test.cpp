#include"../sph_bessel_recursive.h"
#include"gtest/gtest.h"

#define threshold 1e-12

/************************************************
*  unit test of class Sph_Bessel_Recursive
***********************************************/

/**
 * Note: this unit test try to ensure the invariance
 * of the spherical Bessel produced by class Sph_Bessel_Recursive,
 * and the reference results are produced by ModuleBase::Sph_Bessel_Recursive
 * at 2022-1-25.
 * 
 */

double mean(std::vector<double> &vect)
{
    double meanv = 0.0;

    int totN = vect.size();
    for (int i=0; i< totN; ++i) {meanv += vect[i]/totN;}

    return meanv;
}

TEST(SphBessel,D1)
{
    int lmax = 7;
    int rmesh = 700;
    double dx = 0.01;

    ModuleBase::Sph_Bessel_Recursive::D1 sphbesseld1;
    sphbesseld1.set_dx(dx);
    sphbesseld1.cal_jlx(lmax,rmesh);   
    std::vector<std::vector<double>> jlx = sphbesseld1.get_jlx();

    ASSERT_EQ(jlx.size(),static_cast<size_t>(lmax + 1));
    EXPECT_NEAR( mean(jlx[0])/0.2084468748396,    1.0, threshold);
    EXPECT_NEAR( mean(jlx[1])/0.12951635180384,   1.0, threshold);
    EXPECT_NEAR( mean(jlx[2])/0.124201140093879,  1.0, threshold);
    EXPECT_NEAR( mean(jlx[3])/0.118268654505568,  1.0, threshold);
    EXPECT_NEAR( mean(jlx[4])/0.0933871035384385, 1.0, threshold);
    EXPECT_NEAR( mean(jlx[5])/0.0603800487910689, 1.0, threshold);
    EXPECT_NEAR( mean(jlx[6])/0.0327117051555907, 1.0, threshold);
    EXPECT_NEAR( mean(jlx[7])/0.0152155566653926, 1.0, threshold);
}


TEST(SphBessel,D2)
{
    int lmax = 7;
    int rmesh = 700;
    int kmesh = 800;
    double dx = 0.0001;

    ModuleBase::Sph_Bessel_Recursive::D2 sphbesseld2;
    sphbesseld2.set_dx(dx);
    sphbesseld2.cal_jlx(lmax,rmesh,kmesh);   
    std::vector<std::vector<std::vector<double>>> jlxd2 = sphbesseld2.get_jlx();
    std::vector<std::vector<double>> jlx(lmax+1);

    ASSERT_EQ(jlxd2.size(),static_cast<size_t>(lmax + 1));

    //calculate the mean of jlxd2[i][j] and assign to jlx[i][j]
    for(int i=0; i<jlxd2.size(); ++i)
    {
        jlx[i].resize(jlxd2[i].size());
        for(int j=0; j< jlxd2[i].size(); ++j)
        {
            jlx[i][j] = mean(jlxd2[i][j]);
        }
    }  
    
    EXPECT_NEAR( mean(jlx[0])/0.130406547960426,  1.0, threshold);
    EXPECT_NEAR( mean(jlx[1])/0.0643093491554227, 1.0, threshold);
    EXPECT_NEAR( mean(jlx[2])/0.0434912807857165, 1.0, threshold);
    EXPECT_NEAR( mean(jlx[3])/0.0329463027246214, 1.0, threshold);
    EXPECT_NEAR( mean(jlx[4])/0.0264877891341284, 1.0, threshold);
    EXPECT_NEAR( mean(jlx[5])/0.0220804766801247, 1.0, threshold);
    EXPECT_NEAR( mean(jlx[6])/0.0188550846449362, 1.0, threshold);
    EXPECT_NEAR( mean(jlx[7])/0.0163891245775448, 1.0, threshold);
}
