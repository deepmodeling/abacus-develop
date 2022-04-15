#include"../math_sphbes.h"
#include<iostream>

#ifdef __MPI
#include"mpi.h"
#endif

#include"gtest/gtest.h"

#define doublethreshold 1e-7


/************************************************
*  unit test of class Integral
***********************************************/

/**
 * Note: this unit test try to ensure the invariance
 * of the spherical Bessel produced by class Sphbes,
 * and the reference results are produced by ABACUS 
 * at 2022-1-27.
 * 
 * Tested function: 
 *      - Spherical_Bessel.
 *      - Spherical_Bessel_Roots
 * 
 */

double mean(const double* vect, const int totN)
{
    double meanv = 0.0;
    for (int i=0; i< totN; ++i) {meanv += vect[i]/totN;}
    return meanv;
}

class Sphbes : public testing::Test
{
    protected:

    int     msh =   700;
    int     l0  =   0;
    int     l1  =   1;
    int     l2  =   2;
    int     l3  =   3;
    int     l4  =   4;
    int     l5  =   5;
    int     l6  =   6;
    int     l7  =   7;
    double  q   =   1.0;
    double  *r  =   new double[msh];       
    double  *jl =   new double[msh];

    void SetUp()
    {
        for(int i=0; i<msh; ++i) {r[i] = 0.01*(i);}
    }

    void TearDown()
    {
        delete [] r;
        delete [] jl;
    }       
};


TEST_F(Sphbes,SphericalBessel)
{
    //int l = 0;
    ModuleBase::Sphbes::Spherical_Bessel(msh,r,q,l0,jl);
    //reference result is from bessel_test.cpp which is calculated by
    //ModuleBase::Sph_Bessel_Recursive::D1
    EXPECT_NEAR(mean(jl,msh)/0.2084468748396,   1.0,doublethreshold); 


    //int l = 1;
    ModuleBase::Sphbes::Spherical_Bessel(msh,r,q,l1,jl);
    //reference result is from bessel_test.cpp which is calculated by
    //ModuleBase::Sph_Bessel_Recursive::D1
    EXPECT_NEAR(mean(jl,msh)/0.12951635180384,   1.0,doublethreshold);


    //int l = 2;
    ModuleBase::Sphbes::Spherical_Bessel(msh,r,q,l2,jl);
    //the result from bessel_test.cpp calculated by
    //ModuleBase::Sph_Bessel_Recursive::D1 is 0.124201140093879
    //reference result is calculated by Sphbes::Spherical_Bessel(msh,r,q,l,jl)
    EXPECT_NEAR(mean(jl,msh)/0.12420114009271901456,   1.0,doublethreshold);

    //int l = 3;
    ModuleBase::Sphbes::Spherical_Bessel(msh,r,q,l3,jl);
    //the result from bessel_test.cpp calculated by
    //ModuleBase::Sph_Bessel_Recursive::D1 is 0.118268654505568
    //reference result is calculated by Sphbes::Spherical_Bessel(msh,r,q,l,jl)
    EXPECT_NEAR(mean(jl,msh)/0.11826865448477598408,   1.0,doublethreshold);


    //int l = 4;
    ModuleBase::Sphbes::Spherical_Bessel(msh,r,q,l4,jl);
    //the result from bessel_test.cpp calculated by
    //ModuleBase::Sph_Bessel_Recursive::D1 is 0.0933871035384385
    //reference result is calculated by Sphbes::Spherical_Bessel(msh,r,q,l,jl)
    EXPECT_NEAR(mean(jl,msh)/0.093387100084701621383,   1.0,doublethreshold);


    //int l = 5;
    ModuleBase::Sphbes::Spherical_Bessel(msh,r,q,l5,jl);
    //the result from bessel_test.cpp calculated by
    //ModuleBase::Sph_Bessel_Recursive::D1 is 0.0603800487910689
    //reference result is calculated by Sphbes::Spherical_Bessel(msh,r,q,l,jl)
    EXPECT_NEAR(mean(jl,msh)/0.060380048719821471925,   1.0,doublethreshold);


    //int l = 6;
    ModuleBase::Sphbes::Spherical_Bessel(msh,r,q,l6,jl);
    //the result from bessel_test.cpp calculated by
    //ModuleBase::Sph_Bessel_Recursive::D1 is 0.0327117051555907
    //reference result is calculated by Sphbes::Spherical_Bessel(msh,r,q,l,jl)
    EXPECT_NEAR(mean(jl,msh)/0.032711705053977857549,   1.0,doublethreshold);


    //int l = 7;
    ModuleBase::Sphbes::Spherical_Bessel(msh,r,q,l7,jl);
    //the result from bessel_test.cpp calculated by
    //ModuleBase::Sph_Bessel_Recursive::D1 is 0.0152155566653926
    //reference result is calculated by Sphbes::Spherical_Bessel(msh,r,q,l,jl)
    EXPECT_NEAR(mean(jl,msh)/0.015215556095798710851,   1.0,doublethreshold);
}

TEST_F(Sphbes,SphericalBesselRoots)
{
    int neign = 100;
    double **eign = new double*[8];
    for(int i=0;i<8;++i)
    {
        eign[i] = new double[neign];
        ModuleBase::Sphbes::Spherical_Bessel_Roots(neign,i,1.0e-12,eign[i],10.0);
    }

    EXPECT_NEAR(eign[0][0]/0.31415926535899563188,       1.0,doublethreshold);
    EXPECT_NEAR(eign[0][99]/31.415926535896932847,       1.0,doublethreshold);
    EXPECT_NEAR(mean(eign[0],100)/15.865042900628463229, 1.0,doublethreshold);
    EXPECT_NEAR(eign[1][0]/0.44934094579091843347,       1.0,doublethreshold);
    EXPECT_NEAR(eign[1][99]/31.572689440204385392,       1.0,doublethreshold);
    EXPECT_NEAR(mean(eign[1],100)/16.020655759558295017, 1.0,doublethreshold);
    EXPECT_NEAR(eign[2][0]/0.57634591968946913276,       1.0,doublethreshold);
    EXPECT_NEAR(eign[2][99]/31.729140298172534784,       1.0,doublethreshold);
    EXPECT_NEAR(mean(eign[2],100)/16.175128483074864505, 1.0,doublethreshold);
    EXPECT_NEAR(eign[3][0]/0.69879320005004752492,       1.0,doublethreshold);
    EXPECT_NEAR(eign[3][99]/31.885283678838447941,       1.0,doublethreshold);
    EXPECT_NEAR(mean(eign[3],100)/16.328616567969248763, 1.0,doublethreshold);
    EXPECT_NEAR(eign[4][0]/0.81825614525711076741,       1.0,doublethreshold);
    EXPECT_NEAR(eign[4][99]/32.041124042016576823,       1.0,doublethreshold);
    EXPECT_NEAR(mean(eign[4],100)/16.481221742387987206, 1.0,doublethreshold);
    EXPECT_NEAR(eign[5][0]/0.93558121110426506473,       1.0,doublethreshold);
    EXPECT_NEAR(eign[5][99]/32.196665741899131774,       1.0,doublethreshold);
    EXPECT_NEAR(mean(eign[5],100)/16.633019118735202113, 1.0,doublethreshold);
    EXPECT_NEAR(eign[6][0]/1.051283540809391015,         1.0,doublethreshold);
    EXPECT_NEAR(eign[6][99]/32.351913030537232885,       1.0,doublethreshold);
    EXPECT_NEAR(mean(eign[6],100)/16.784067905062840964, 1.0,doublethreshold);
    EXPECT_NEAR(eign[7][0]/1.1657032192516516567,        1.0,doublethreshold);
    EXPECT_NEAR(eign[7][99]/32.506870061157627561,       1.0,doublethreshold);
    EXPECT_NEAR(mean(eign[7],100)/16.934416735327332049, 1.0,doublethreshold);
    
    for(int i=0;i<8;++i) delete [] eign[i];
    delete [] eign;
}


int main(int argc, char **argv)
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